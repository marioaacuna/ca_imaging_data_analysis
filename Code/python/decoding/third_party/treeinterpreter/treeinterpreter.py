 #String manipulation
import re
# Numerical
import numpy as np
import pandas as pd
# joblib is optional
try:
    from joblib import Parallel, delayed
    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False
    Parallel = None
    delayed = None


class TreeInterpreter(object):
    def __init__(self, model):
        """Checks that model is supported and has been trained.

        :param model: The input model to interpret.
        """

        # Check that model type is supported
        model_class, model_type, implementation, estimator_type = _validate_model_type(model)
        # Check that model has been fit
        _validate_model_is_trained(model, model_class)

        # Store model and info on it
        self.model = model
        self.__internals = dict(model_class=model_class,
                                model_type=model_type,
                                implementation=implementation,
                                estimator_type=estimator_type)

        # Initialize all attributes
        self.target_names = None
        self.feature_names = None
        self.predictions = None
        self.biases = None
        self.contributions = None
        self.conditional_contributions = dict()


    def interpret(self, X_test, target_names=None, feature_names=None,
                  compute_conditional_contribution=False, n_jobs=None,
                  verbose=False):
        """Main method to interpret the provided model. It computes the feature
        contributions , such that predictions ≈ bias + feature_contributions.

        :param X_test: [ndarray or DataFrame] Data on which to test feature
            contributions. It must have the same number of features of the
            dataset used to train the model.
        :param target_names: [list] Name of each target. Used for plots, output
            DataFrames and printing. If model is a classifier, these are the
            names of the classes, if a regressor these are the names of the
            target variables.
        :param feature_names: [list] Name of each feature. Used for plots,
            output DataFrames and printing. If X_test is a DataFrame, names are
            silently inferred from it, or they can be overridden with this
            argument.
        :param compute_conditional_contribution: [bool] Whether to also compute
            all conditional contributions along the path.
        :param n_jobs: [int or None] The number of parallel processes to use if
            joblib is installed. If None, trees are processed sequentially.
        :param verbose: [bool] Whether to print messages to the console regarding
            progress and outcomes.

        :return self: This allows the user to call this method together with
            initialization, and return the object in a variable, that is,
            TI = TreeInterpreter(model).interpret(X_test)

        The following attributes are stored in self:
        predictions: [ndarray] Contains the prediction of each feature to
            each observation and class, averaged across trees.
        biases: [ndarray] Contains the baseline prediction of each feature
            to each observation and class, averaged across trees.
        contributions: [ndarray] Contains the contribution of each feature
            to each observation and class, averaged across trees.
        conditional_contributions: [dictionary] (optional if
            `compute_conditional_contribution` == True) A dictionary containing
            the values of contribution of each feature conditioned on previous
            features. Each key contains the list of features along the decision
            path, from the root down to the leaf. Each dictionary value contains
            a numpy array containing the conditional contribution of each
            feature to each observation, averaged across trees.
        """

        # Set class names
        if target_names is None:
            target_names = ['target_%i' % (i + 1) for i in range(self.model.n_classes_)]
        if len(target_names) != self.model.n_classes_:
            raise ValueError('The number of target names provided does not match the number of targets on which the model was trained.')
        # Store class names
        self.target_names = target_names

        # Set feature names
        if isinstance(X_test, pd.DataFrame):
            if feature_names is None:
                feature_names = X_test.columns
            # Convert X_test to numpy array
            X_test = X_test.values.copy()

        elif isinstance(X_test, np.ndarray):
            if feature_names is None:
                feature_names = ['feature_%i' % (i + 1) for i in range(self.model.n_features_)]

        else:
            raise TypeError('X_test should be either a numpy array or a pandas DataFrame')

        # Check that number of feature names is correct
        if len(feature_names) != self.model.n_features_:
            raise ValueError('The number of feature names provided does not match the number of features on which the model was trained.')

        # Store feature names
        self.feature_names = feature_names

        # Initialize output variables
        predictions = None
        biases = None
        contributions = None
        n_evaluations = None
        conditional_contributions = dict()
        conditional_contributions_samples = dict()
        # Gather some info
        n_trees = len(self.model.estimators_)

        ########################################################################
        # Process trees in parallel
        ########################################################################
        if HAS_JOBLIB and n_jobs is not None:
            results = Parallel(n_jobs=n_jobs)(
                    delayed(self.__compute_feature_contributions_from_tree)(
                        estimator=tree, X_test=X_test,
                        compute_conditional_contribution=compute_conditional_contribution)
                    for tree in self.model.estimators_)

            # Unpack outputs and average across trees
            predictions = np.nanmean([i[0] for i in results], axis=0)
            biases = np.nanmean([i[1] for i in results], axis=0)
            contributions = np.nanmean([i[2] for i in results], axis=0)

            if compute_conditional_contribution:
                # Unpack outputs
                all_conditional_contributions = [i[3] for i in results]

                # Collect conditional contributions from all trees
                for i_tree in range(len(self.model.estimators_)):
                    # Get conditional feature contributions for this tree
                    this_conditional_contributions = all_conditional_contributions[i_tree]
                    if this_conditional_contributions is None:
                        continue

                    # Get list of feature sets used by this tree
                    feature_sets = list(this_conditional_contributions.keys())
                    for features in feature_sets:
                        if features in conditional_contributions.keys():
                            conditional_contributions[features] = np.dstack((conditional_contributions[features], this_conditional_contributions[features]))
                        else:
                            conditional_contributions[features] = this_conditional_contributions[features]

                # Average across all trees
                feature_sets = list(conditional_contributions.keys())
                for features in feature_sets:
                    if conditional_contributions[features].ndim == 3:
                        conditional_contributions[features] = np.mean(
                                conditional_contributions[features], axis=2)


        ########################################################################
        # Process trees sequentially
        ########################################################################
        else:
            for i_tree, tree in enumerate(self.model.estimators_):
                if verbose:
                    print('Processing tree %i/%i' % (i_tree + 1, n_trees))

                # Compute predictions
                results = self.__compute_feature_contributions_from_tree(
                    estimator=tree, X_test=X_test,
                    compute_conditional_contribution=compute_conditional_contribution)

                # Extract values of contributions, replace NaN with 0 and count number of non-NaN
                n = np.logical_not(np.isnan(results[2])).astype(int)
                c = results[2]
                c[np.logical_not(np.isfinite(c))] = 0

                # Unpack outputs and average trees processed so far
                if i_tree == 0:
                    predictions = results[0]
                    biases = results[1]
                    contributions = c
                    n_evaluations = n

                else:
                    # Simply add values (and divide at the end)
                    predictions = predictions + results[0]
                    biases = biases + results[1]
                    # Simply add the contribution values, and we'll divide at the end
                    n_evaluations = n_evaluations + n
                    contributions = contributions + c

                if compute_conditional_contribution:
                    # Get conditional feature contributions for this tree
                    this_conditional_contributions = results[3]
                    if this_conditional_contributions is None:
                        continue

                    # Get list of feature sets used by this tree
                    feature_sets = list(this_conditional_contributions.keys())
                    for features in feature_sets:
                        if features in conditional_contributions.keys():
                            i_iter = conditional_contributions_samples[features]
                            if this_conditional_contributions[
                                features].ndim == 3:
                                for ii in range(this_conditional_contributions[features].shape[2]):
                                    conditional_contributions[features] = _iterative_mean(i_iter + ii, conditional_contributions[features], this_conditional_contributions[features][:, :, ii])
                                    conditional_contributions_samples[features] += 1
                            else:
                                conditional_contributions[
                                    features] = _iterative_mean(i_iter, conditional_contributions[features], this_conditional_contributions[features])
                                conditional_contributions_samples[features] += 1

                        else:
                            if this_conditional_contributions[features].ndim == 3:
                                conditional_contributions[features] = np.mean(this_conditional_contributions[features], axis=2)
                            else:
                                conditional_contributions[features] = this_conditional_contributions[features]
                            conditional_contributions_samples[features] = 1

            # Compute mean values
            predictions = predictions / n_trees
            biases = biases / n_trees
            # Divide contributions only by the number of times the feature was evaluated.
            # Features that were never evaluated will return NaN
            contributions = _divide0(contributions, n_evaluations, replace_with=np.nan)

            if verbose:
                if np.any(n_evaluations == 0):
                    n_not_evaluated_features = np.unique(np.where(n_evaluations == 0)[1]).shape[0]
                    raise UserWarning('%i out %i (%.1f%%) features were never evaluated by the classifier.\nConsider increasing the number of estimators' % (
                            n_not_evaluated_features, n_evaluations.shape[1], n_not_evaluated_features / n_evaluations.shape[1] * 100))

        # Store results
        self.contributions = contributions
        self.predictions = predictions
        self.biases = biases

        # Compute mean feature contribution of each feature to each target
        print()
        self.mean_contribution_per_target = pd.DataFrame(np.mean(self.contributions, axis=0),
                                                         index=self.feature_names, columns=self.target_names)









        if compute_conditional_contribution:
            self.conditional_contributions = conditional_contributions

        return self


    def __compute_feature_contributions_from_tree(self, estimator, X_test,
                                                 compute_conditional_contribution=False):
        """For a given DecisionTreeRegressor, DecisionTreeClassifier,
        ExtraTreeRegressor, or ExtraTreeClassifier, returns [prediction, bias,
        feature_contributions], such that prediction ≈ bias + feature_contributions.

        :param estimator: The tree from which to calculate feature contributions.
        :param X_test: [ndarray] Data on which to test feature contributions. It must
            have the same number of features of the dataset used to train the model.
        :param compute_conditional_contribution [bool]: Whether to also return all
            the conditional contributions along the path.

        :return predictions: [ndarray] Prediction of each feature to each
            observation and target.
        :return biases: [ndarray] Baseline prediction of each feature to each
            observation and target.
        :return contributions: [ndarray] Contribution of each feature to each
            observation and target.
        :return feature_depth: [list of lists] Each item contains a feature and
            the depth at which it was used for a split. Unused features are not
            listed.

        :return conditional_contributions: [dictionary] (optional if `compute_conditional_contribution`
            == True) A dictionary containing the values of contribution of each
            feature conditioned on previous features. Each key contains the list of
            features along the decision path, from the root down to the leaf. Each
            dictionary value contains a numpy array containing the conditional
            contribution of each feature to each observation.
        :return conditional_contributions_samples: [dictionary] (optional if
            `compute_conditional_contribution` == True) A dictionary with the same
            keys of conditional_contributions. Instead of conditional contribution
            values, it contains the index of the observations where the values of
            conditional distributions were calculated. This dictionary could be used
            to select only observations from the true class of interest.
        """

        # Get number of test samples
        n_samples = X_test.shape[0]
        # Get list of the feature used at each node
        features = estimator.tree_.feature

        # Initialize output variables
        contributions = np.zeros((n_samples, estimator.n_features_,
                                  estimator.n_classes_)) * np.nan
        biases = np.zeros((n_samples, estimator.n_classes_)) * np.nan
        conditional_contributions = None
        conditional_contributions_samples = None
        # Initialize intermediate variables
        path_to_features = None

        # Retrieve leaves and paths
        leaves_X = estimator.apply(X_test)
        paths = _get_tree_paths(estimator.tree_, node_id=0)
        # Reverse direction of path, and convert to tuple
        paths = tuple([np.array(path)[::-1] for path in paths])
        # Map leaves to paths
        leaf_to_path = {path[-1]: path for path in paths}

        # Remove single-dimensional inner arrays
        values = estimator.tree_.value.squeeze(axis=1)
        # Reshape if squeezed into a single float
        if len(values.shape) == 0:
            values = np.array([values])

        # Compute bias per class
        if self.__internals['estimator_type'] == 'DecisionTreeRegressor':
            biases = np.full(X_test.shape[0], values[paths[0][0]])

        elif self.__internals['estimator_type'] == 'DecisionTreeClassifier':
            # scikit stores category counts, we turn them into probabilities
            normalizer = values.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            values /= normalizer
            biases = np.tile(values[paths[0][0]], (X_test.shape[0], 1))

        # Get the predictions of the tree for each observation
        predictions = values[leaves_X]

        # If the tree did not perform any split, return immediately
        if features[features >= 0].shape[0] == 0:
            if compute_conditional_contribution:
                return predictions, biases, contributions, None, None
            else:
                return predictions, biases, contributions

        if compute_conditional_contribution:
            # Map each path to the whole conditional set of features up to that node
            path_to_features = {tuple(path): tuple(features[path[:-1]]) for path
                                in paths}
            # Create the dictionaries containing the actual values of conditional
            # contribution, and the index of the samples that ended up in each array
            conditional_contributions = dict()
            conditional_contributions_samples = dict()

        # Compute contributions
        for i_obs, leaf in enumerate(leaves_X):
            path = leaf_to_path[leaf]
            # Compute absolute contribution of each feature at each node (that is, at each split)
            contrib_features = values[path[1:], :] - values[path[:-1], :]
            contrib_features_index = features[path[:-1]]
            # Store data
            contributions[i_obs, contrib_features_index, :] = np.sum(
                contrib_features, axis=0)

            # Compute conditional contributions
            if compute_conditional_contribution:
                # Compute incremental contributions down the path, due by conditioning
                # a feature preceding ones in the path
                contrib_features_joint = np.cumsum(contrib_features, axis=0)

                # Store values
                features_at_this_node = path_to_features[tuple(path)]
                if features_at_this_node in conditional_contributions.keys():
                    conditional_contributions[
                        features_at_this_node] = np.dstack((conditional_contributions[features_at_this_node], contrib_features_joint))
                    conditional_contributions_samples[
                        features_at_this_node].append(i_obs)
                else:
                    conditional_contributions[features_at_this_node] = contrib_features_joint
                    conditional_contributions_samples[features_at_this_node] = [i_obs]

        # Return results
        if compute_conditional_contribution:
            # Add a tuple of feature_map, contribution values, and sample index of
            # each contribution
            return predictions, biases, contributions, conditional_contributions, conditional_contributions_samples
        else:
            return predictions, biases, contributions


################################################################################
# Navigate scikit-learn tree model
################################################################################
def _get_tree_paths(tree, node_id, depth=0):
    """Recursively navigate a tree model to gather the node_ids of each decision
    path.
    """
    # The following variables mark leaves and undefined nodes. They can be obtained with:
    # from sklearn.tree._tree import TREE_LEAF, TREE_UNDEFINED
    TREE_LEAF = -1
    # For the records, TREE_UNDEFINED = -2

    if node_id == TREE_LEAF:
        raise ValueError("Invalid node_id %s" % TREE_LEAF)

    left_child = tree.children_left[node_id]
    right_child = tree.children_right[node_id]

    if left_child != TREE_LEAF:
        left_paths = _get_tree_paths(tree, left_child, depth=depth + 1)
        right_paths = _get_tree_paths(tree, right_child, depth=depth + 1)

        for path in left_paths:
            path.append(node_id)
        for path in right_paths:
            path.append(node_id)
        paths = left_paths + right_paths
    else:
        paths = [[node_id]]

    return paths


################################################################################
# Utilities
################################################################################
def _validate_model_type(model):
    """Check whether TreeInterpreter can work on this model. If not recognized
    this function returns a NotImplementedError error.

    :param model: The input model.
    :return model_class: [str] The model class (e.g., RandomForestClassifier).
    :return model_type: [str] The type of model (i.e., 'classifier' or
        'regressor').
    :return implementation: [str] The package implementing the model (e.g.,
        sklearn).
    :return estimator_type: [str] The class of the estimator (e.g.,
        DecisionTreeClassifier).
    """

    # List of known model types
    KNOWN_MODEL_TYPES = [
        # RandomForests
        'sklearn.ensemble.forest.RandomForestClassifier',
        'sklearn.ensemble.forest.RandomForestRegressor',
        # Ensemble of extremely randomized trees
        'sklearn.ensemble.forest.ExtraTreesClassifier',
        'sklearn.ensemble.forest.ExtraTreesRegressor',

        # Decision trees
        'sklearn.tree.DecisionTreeClassifier',
        'sklearn.tree.DecisionTreeRegressor',
        # Extremely randomized trees
        'sklearn.tree.ExtraTreeClassifier',
        'sklearn.tree.ExtraTreeRegressor'
        ]

    # The class of the current model
    model_class = str(type(model))
    # Check class against known types
    result = re.search('\'(.*)\'', model_class)
    model_class_str = result.group(1)
    if model_class_str in KNOWN_MODEL_TYPES:
        # Infer class and package of origin of the model
        model_class_parts = model_class_str.split('.')
        model_class = model_class_parts[-1]
        implementation = model_class_parts[0]

        # Get type of estimator
        if model_class == 'RandomForestClassifier':
            estimator_type = 'DecisionTreeClassifier'
        elif model_class == 'RandomForestRegressor':
            estimator_type = 'DecisionTreeRegressor'
        elif model_class == 'ExtraTreesClassifier':
            estimator_type = 'ExtraTreeClassifier'
        elif model_class == 'ExtraTreesRegressor':
            estimator_type = 'ExtraTreeRegressor'
        else:
            estimator_type = model_class

        # Get type of model
        if 'classifier' in estimator_type.lower():
            model_type = 'classifier'
        elif 'regressor' in estimator_type.lower():
            model_type = 'regressor'
        else:
            raise NotImplementedError('Not clear whether \'%s\' is a classifier or a regressor')

        return model_class, model_type, implementation, estimator_type

    else:
        raise NotImplementedError('Class \'%s\' is not supported by TreeInterpreter' % model_class_str)


def _validate_model_is_trained(model, model_type):
    """Check whether the model has been already trained.

    :param model: The input model.
    :param model_type: [str] The model type (e.g., RandomForestClassifier).
    :return Either nothing or an AttributeError.
    """
    if model_type in ['RandomForestClassifier', 'RandomForestRegressor',
                      'ExtraTreesClassifier', 'ExtraTreesRegressor']:
        has_this_attribute = 'n_estimators'
        has_not_this_attribute = 'estimators_'

    else:
        raise NotImplementedError('Don\'t know what to do with \'%s\'' % model_type)

    if hasattr(model, has_this_attribute) and not hasattr(model, has_not_this_attribute):
        raise AttributeError('The model has not been trained yet, and thus cannot be interpreted.')


def _iterative_mean(iter, current_mean, x):
    """
    Iteratively calculates mean using http://www.heikohoffmann.de/htmlthesis/node134.html
    Originally implemented in https://github.com/andosa/treeinterpreter/pull/24

    :param iter: [int > 0] Current iteration.
    :param current_mean: [ndarray] Current value of mean.
    :param x: [ndarray] New value to be added to mean.

    :return: [ndarray] Updated mean.
    """

    return current_mean + ((x - current_mean) / (iter + 1))


def _divide0(a, b, replace_with):
    """Divide two numbers but replace its result if division is not possible,
    e.g., when dividing a number by 0. No type-checking or agreement between
    dimensions is performed. Be careful!

    :param a: [ndarray or int or float] Numerator.
    :param b: [ndarray or int or float] Denominator.
    :param replace_with: [int or float] Return this number if a/b is not defined.

    :return: [ndarray or int or float] Result of division, cast by numpy to the
        best data type to hold it.
    """

    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        if isinstance(c, np.ndarray):
            c[np.logical_not(np.isfinite(c))] = replace_with
        else:
            if not np.isfinite(c):
                c = replace_with

    return c






if __name__ == '__main__':
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.datasets import make_classification
    from sklearn.model_selection import train_test_split

    X, y = make_classification(n_samples=1000, n_features=10, n_informative=2, n_redundant=2, n_repeated=0)
    RF = RandomForestClassifier(n_estimators=100)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=.20)
    RF.fit(X_train, y_train)

    TI = TreeInterpreter(RF).interpret(X_test)

    print()
