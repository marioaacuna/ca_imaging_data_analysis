# Imports
import numpy as np
from sklearn.tree import DecisionTreeRegressor, DecisionTreeClassifier
# joblib is optional
try:
    from joblib import Parallel, delayed
    HAS_JOBLIB = True
except ImportError:
    HAS_JOBLIB = False


# The following variables mark leaves and undefined nodes. They can be obtained with:
# from sklearn.tree._tree import TREE_LEAF, TREE_UNDEFINED
TREE_LEAF = -1
TREE_UNDEFINED = -2


################################################################################
# Compute feature contributions
################################################################################
def compute_feature_contributions_ensemble(clf, X_test, compute_conditional_contribution=False, n_jobs=None, verbose=False):
    """For a given RandomForestRegressor, RandomForestClassifier, ExtraTreesRegressor,
    or ExtraTreesClassifier returns the mean [prediction, bias, feature_contributions],
    such that prediction ≈ bias + feature_contributions.

    :param clf: [RandomForestRegressor or RandomForestClassifier or
        ExtraTreesRegressor or ExtraTreesClassifier] A trained ensemble classifier.
    :param X_test: [ndarray] Data on which to test feature contributions. It must
        have the same number of features of the dataset used to train the model.
    :param compute_conditional_contribution [bool]: Whether to also return all
        the conditional contributions along the path.
    :param n_jobs: [int or None] The number of parallel processes to use if joblib
        is installed. If None, trees in the ensemble are processed sequentially.
    :param verbose [bool] Whether to print messages to the console regarding
        progress.

    :return predictions: [ndarray] Contains the prediction of each feature to
        each observation and class, averaged across trees.
    :return biases: [ndarray] Contains the baseline prediction of each feature
        to each observation and class, averaged across trees.
    :return contributions: [ndarray] Contains the contribution of each feature
        to each observation and class, averaged across trees.
    :return conditional_contributions: [dictionary] (optional if `compute_conditional_contribution`
        == True) A dictionary containing the values of contribution of each
        feature conditioned on previous features. Each key contains the list of
        features along the decision path, from the root down to the leaf. Each
        dictionary value contains a numpy array containing the conditional
        contribution of each feature to each observation, averaged across trees.
    """

    # Initialize output variables
    predictions = None
    biases = None
    contributions = None
    n_evaluations = None
    conditional_contributions = dict()
    conditional_contributions_samples = dict()
    # Gather some info
    n_trees = len(clf.estimators_)


    ############################################################################
    # Process trees in parallel
    ############################################################################
    if HAS_JOBLIB and n_jobs is not None:
        results = Parallel(n_jobs=n_jobs)(
                delayed(_compute_feature_contributions_from_tree)(estimator=tree, X_test=X_test,
                                                                  compute_conditional_contribution=compute_conditional_contribution)
                for tree in clf.estimators_)

        # Unpack outputs and average across trees
        predictions = np.nanmean([i[0] for i in results], axis=0)
        biases = np.nanmean([i[1] for i in results], axis=0)
        contributions = np.nanmean([i[2] for i in results], axis=0)

        if compute_conditional_contribution:
            # Unpack outputs
            all_conditional_contributions = [i[3] for i in results]

            # Collect conditional contributions from all trees
            for i_tree in range(len(clf.estimators_)):
                # Get conditional feature contributions for this tree
                this_conditional_contributions = all_conditional_contributions[i_tree]
                if this_conditional_contributions is None:
                    continue

                # Get list of feature sets used by this tree
                feature_sets = list(this_conditional_contributions.keys())
                for features in feature_sets:
                    if features in conditional_contributions.keys():
                        conditional_contributions[features] = np.dstack((conditional_contributions[features],
                                                                        this_conditional_contributions[features]))
                    else:
                        conditional_contributions[features] = this_conditional_contributions[features]

            # Average across all trees
            feature_sets = list(conditional_contributions.keys())
            for features in feature_sets:
                if conditional_contributions[features].ndim == 3:
                    conditional_contributions[features] = np.mean(conditional_contributions[features], axis=2)


    ############################################################################
    # Process trees sequentially
    ############################################################################
    else:
        for i_tree, tree in enumerate(clf.estimators_):
            if verbose:
                print('Processing tree %i/%i' % (i_tree + 1, n_trees))

            # Compute predictions
            results = _compute_feature_contributions_from_tree(estimator=tree, X_test=X_test,
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
                        if this_conditional_contributions[features].ndim == 3:
                            for ii in range(this_conditional_contributions[features].shape[2]):
                                conditional_contributions[features] = _iterative_mean(i_iter + ii, conditional_contributions[features],
                                                                                      this_conditional_contributions[features][:, :, ii])
                                conditional_contributions_samples[features] += 1
                        else:
                            conditional_contributions[features] = _iterative_mean(i_iter, conditional_contributions[features],
                                                                                  this_conditional_contributions[features])
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


    # Return outputs
    if compute_conditional_contribution:
        return predictions, biases, contributions, conditional_contributions
    else:
        return predictions, biases, contributions


def _compute_feature_contributions_from_tree(estimator, X_test, compute_conditional_contribution=False):
    """For a given DecisionTreeRegressor, DecisionTreeClassifier,
    ExtraTreeRegressor, or ExtraTreeClassifier, returns [prediction, bias,
    feature_contributions], such that prediction ≈ bias + feature_contributions.

    :param estimator: [DecisionTreeRegressor or DecisionTreeClassifier or
        ExtraTreeRegressor or ExtraTreeClassifier] The tree from which to calculate
        feature contributions.
    :param X_test: [ndarray] Data on which to test feature contributions. It must
        have the same number of features of the dataset used to train the model.
    :param compute_conditional_contribution [bool]: Whether to also return all
        the conditional contributions along the path.

    :return predictions: [ndarray] Contains the prediction of each feature to
        each observation and class.
    :return biases: [ndarray] Contains the baseline prediction of each feature
        to each observation and class.
    :return contributions: [ndarray] Contains the contribution of each feature
        to each observation and class.
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
    contributions = np.zeros((n_samples, estimator.n_features_, estimator.n_classes_)) * np.nan
    biases = np.zeros((n_samples, estimator.n_classes_)) * np.nan
    conditional_contributions = None
    conditional_contributions_samples = None
    # Initialize intermediate variables
    path_to_features = None

    # Retrieve leaves and paths
    leaves_X = estimator.apply(X_test)
    paths = _get_tree_paths(estimator.tree_, 0)
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
    if isinstance(estimator, DecisionTreeRegressor):
        biases = np.full(X_test.shape[0], values[paths[0][0]])

    elif isinstance(estimator, DecisionTreeClassifier):
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
        path_to_features = {tuple(path): tuple(features[path[:-1]]) for path in paths}
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
        contributions[i_obs, contrib_features_index, :] = np.sum(contrib_features, axis=0)

        # Compute conditional contributions
        if compute_conditional_contribution:
            # Compute incremental contributions down the path, due by conditioning
            # a feature preceding ones in the path
            contrib_features_joint = np.cumsum(contrib_features, axis=0)

            # Store values
            features_at_this_node = path_to_features[tuple(path)]
            if features_at_this_node in conditional_contributions.keys():
                conditional_contributions[features_at_this_node] = \
                    np.dstack((conditional_contributions[features_at_this_node],
                               contrib_features_joint))
                conditional_contributions_samples[features_at_this_node].append(i_obs)
            else:
                conditional_contributions[features_at_this_node] = contrib_features_joint
                conditional_contributions_samples[features_at_this_node] = [i_obs]

    # Return results
    if compute_conditional_contribution:
        # Add a tuple of feature_map, contribution values, and sample index of
        # each contribution
        return predictions, biases, contributions, conditional_contributions, \
               conditional_contributions_samples
    else:
        return predictions, biases, contributions


################################################################################
# Navigate tree model
################################################################################
def _get_tree_paths(tree, node_id, depth=0):
    """Recursively navigate a tree model to gather the node_ids of each decision
    path.
    """
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
