# Numerical packages
import numpy as np
# import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler, label_binarize
from sklearn.model_selection import GridSearchCV, StratifiedKFold, RepeatedStratifiedKFold, StratifiedShuffleSplit, cross_val_predict
from sklearn.metrics import make_scorer, average_precision_score
from sklearn.calibration import CalibratedClassifierCV, _SigmoidCalibration, clone
from sklearn.isotonic import IsotonicRegression
from sklearn.utils import resample
# from mlxtend.evaluate import mcnemar, mcnemar_table
from scipy.stats import binom_test, ranksums

# Parallel computing
from joblib import Parallel, delayed

# Local repository
from decoding.third_party.pycm import ConfusionMatrix
# from decoding.third_party.treeinterpreter import treeinterpreter as ti
from Utilities.IO_operations import log


################################################################################
### Default parameters
################################################################################
N_JOBS = -1  # number of cores used for parallel computations


################################################################################
### Wrapper of RandomForestClassifier class that supports calibration and significance testing
################################################################################
class MultiClassClassifier(object):
    def __init__(self, optimized_parameters=None, verbose=False):
        """The object can be initialized providing optimized hyper-parameters.

        :param optimized_parameters: [dict or None] Contains the optimized
            parameters to pass the RandomForestClassifier constructor. It has to
            contain key-value combinations that are supported by that class.
            If None, hyperparameters will be optimized via cross-validation.
        :param verbose: [bool] Flag indicating whether to print messages to console.
        """
        # Set default parameters
        self.default_clf_parameters = dict(criterion='entropy', n_jobs=N_JOBS, bootstrap=True, oob_score=False, class_weight='balanced_subsample')
        self.optimization_CV_k_folds = 10

        # Set additional user parameters (these can overwrite default parameters)
        self.optimized_parameters = optimized_parameters
        self._clf_params_need_tuning = self.optimized_parameters is None

        # Store verbose state
        self._verbose = verbose

        # Initialize all other attributes
        self.calibration_method = None
        self._is_calibrated = False
        self._computed_p_values = False
        self._is_fit = False
        self.clf = None
        self.base_estimator = None
        self.X = None
        self.y = None
        self.n_samples = None
        self.n_features = None
        self.classes = None
        self.n_classes = None
        self.posterior_probabilities = None
        self.predicted_labels = None
        self.p_value = None


    def fit_predict(self, data, labels, calibrate=True, calibration_method='sigmoid', CV_k_folds=10):
        """Public method to fit and predict data via cross-validation splitting.

        :param data: [2D numpy array] Contains the observations in Tidy format,
            which is observations on rows and features on columns.
        :param labels: [1D numpy array] The true label of each observation.
        :param calibrate: [bool] Whether to perform calibration of the estimated
            probabilities, which is used to estimate significance of the prediction
            on each observation sample.
        :param calibration_method: [str] The default method to calibrate posterior
            probabilities output by RandomForestClassifier is 'sigmoid', which
            corresponds to Platt's method. The alternative is 'isotonic', which
            corresponds to a non-parametric method. This variable is used only if
            `calibrate` is True.
        :param CV_k_folds: [int] Number of k-fold cross-validation splits.

        :return self
        """
        # Copy data and labels in object attributes
        self.X = data
        self.y = labels
        # Get the structure of data and labels
        self.classes = np.unique(self.y)
        self.n_classes = len(self.classes)
        self.n_samples = self.X.shape[0]
        self.n_features = self.X.shape[1]
        # Get number of splits allowed to do
        n_splits = self._estimate_n_splits_CV(CV_k_folds)

        # Tune hyper-parameters, if not provided by user
        if self._clf_params_need_tuning:
            self._tune_hyperparameters()

        # Initialize RandomForestClassifier object with default parameters
        base_estimator = RandomForestClassifier(**self.default_clf_parameters)
        # Override default parameters with optimized parameters
        base_estimator.set_params(**self.optimized_parameters)
        # Insert classifier in a pipeline together with a StandardScaler object.
        # Store the entire pipeline as the `_base_estimator`.
        self.base_estimator = Pipeline([('scaler', StandardScaler(with_mean=True, with_std=True, copy=True)),
                                        ('classifier', base_estimator)])
        # Create partition for cross-validation
        CV_partition = StratifiedKFold(n_splits=n_splits, shuffle=False)

        # Code block if user chooses to calibrate posterior probabilities
        if calibrate:
            if self._verbose:
                log('Training, calibrating and testing Random Forest classifier')

            # Store calibration method
            self.calibration_method = calibration_method
            # Perform calibration based on the CalibratedClassifierCV class in scikit-learn
            self.clf = CalibratedClassifierCV(base_estimator=self.base_estimator, method=self.calibration_method, cv=CV_partition)
            # Fit and predict posterior probabilities
            self.clf.fit(self.X, self.y)
            self.posterior_probabilities = self.clf.predict_proba(self.X)
            # Toggle internal state
            self._is_calibrated = True

        # Code block if user does not choose to calibrate posterior probabilities
        else:
            if self._verbose:
                log('Training and testing Random Forest classifier')

            # Estimate labels by cross validation. Posterior probabilities are not calibrated
            self.clf = clone(self.base_estimator)
            self.posterior_probabilities = cross_val_predict(self.base_estimator, self.X, self.y, cv=CV_partition, n_jobs=N_JOBS, method='predict_proba')

        # Assign label to each observations according to highest posterior probability
        self.predicted_labels = self.classes[np.argmax(self.posterior_probabilities, axis=1)]
        # Toggle internal state
        self._is_fit = True


    def permute_labels(self, n_significance_permutation=100, test_size=0.1):
        """Method to assess significance of each observation by permuting labels.
        The null distribution is built via stratified shuffle split and not
        cross-validation, which creates a process more similar to bootstrap.

        :param n_significance_permutation: [int] Number of label permutations to
            build a consistent null distribution.
        :param test_size: [float] Fraction of data that will end in the test data.

        :return self
        """
        # Can continue only if performed calibration on empirical data
        if not self._is_calibrated:
            raise Exception('Can compute p-values only on calibrated classifiers')

        # Compute null-distribution of posterior probabilities
        CV_partition = StratifiedShuffleSplit(n_splits=n_significance_permutation, test_size=test_size)
        null_distribution = Parallel(n_jobs=N_JOBS)(
                delayed(self._calibrate_train_test_split)(train_indices, test_indices)
                for train_indices, test_indices in CV_partition.split(self.X, self.y))
        # Stack values along third dimension
        null_distribution = np.dstack(null_distribution)

        # Check the probability of each prediction against the null distribution
        difference = np.array(null_distribution >= self.posterior_probabilities[:, :, None], dtype=int)
        # Count the number of times the relationship is true
        count_per_resample = np.sum(difference, axis=2)
        # Get the count only in the target class
        count_per_resample = count_per_resample[np.arange(self.predicted_labels.shape[0]), self.predicted_labels]
        # Estimate a one-sided p-value
        self.p_value = (count_per_resample + 1) / (n_significance_permutation + 1)
        # Toggle internal state
        self._computed_p_values = True


    def _calibrate_train_test_split(self, train_indices, test_indices):
        """Private method called iteratively by the method `permute_labels`. This
        method represents the true novelty of this class compared to sklearn's
        CalibratedClassifierCV. Here, in fact, we calibrate the classifier via
        cross-validation on the train set and then predict the labels on the
        independent test set.

        :param train_indices: [numpy array] Indices of samples of the train set.
        :param test_indices: [numpy array] Indices of samples of the test set.

        :return calibrated_probabilities: [numpy array] Posterior probabilities
            of each sample to belong to each class. It has shape [n_samples x
            n_classes].
        """
        # Initialize variables used in this function
        calibrator = None

        # Shuffle labels without replacement (i.e., permute labels)
        y = resample(self.y, n_samples=self.n_samples, replace=False)

        # Get train and test subsets
        X_train = self.X[train_indices, :]
        y_train = y[train_indices]
        y_test = y[test_indices]

        # Fit a Random Forest on the train set
        RF = clone(self.base_estimator)
        RF.fit(X_train, y_train)
        # Predict uncalibrated probabilities of all samples
        uncalibrated_probabilities_all = RF.predict_proba(self.X)
        # Slice probabilities of test set
        uncalibrated_probabilities_test = uncalibrated_probabilities_all[test_indices, :].copy()
        if self.n_classes == 2:
            uncalibrated_probabilities_all = uncalibrated_probabilities_all[:, 1:]
            uncalibrated_probabilities_test = uncalibrated_probabilities_test[:, 1:]

        # Convert class labels to binary labels
        binary_y_test_labels = label_binarize(y_test, self.classes)
        calibrated_probabilities = np.zeros((self.n_samples, self.n_classes))

        # Calibrate posterior probabilities of test set. Note that these samples
        # have not been used to fit the RandomForest that predicted them.
        for k, p_test, p_test_all in zip(self.classes, uncalibrated_probabilities_test.T, uncalibrated_probabilities_all.T):
            # Initialize calibrator object according to the calibration method
            # chosen by the user
            if self.calibration_method == 'isotonic':
                calibrator = IsotonicRegression(out_of_bounds='clip')
            elif self.calibration_method == 'sigmoid':
                calibrator = _SigmoidCalibration()
            # Fit calibrator on uncalibrated probabilities of "test" samples
            calibrator.fit(p_test, binary_y_test_labels[:, k])
            if self.n_classes == 2:
                k += 1  # This makes the probabilities be assigned to the second column and terminate the for-loop
            # Adjust uncalibrated probabilities of all samples
            calibrated_probabilities[:, k] = calibrator.predict(p_test_all)

        # Normalize probabilities so that the sum for each observation is 1
        if self.n_classes == 2:
            calibrated_probabilities[:, 0] = 1. - calibrated_probabilities[:, 1]
        else:
            calibrated_probabilities /= np.sum(calibrated_probabilities, axis=1)[:, np.newaxis]
        # Fix when all probabilities were 0 and we ended up having a nan. In that
        # case, assign a uniform probability to each class
        calibrated_probabilities[np.isnan(calibrated_probabilities)] = 1. / self.n_classes
        # Deal with cases where the predicted probability minimally exceeds 1.0
        calibrated_probabilities[(1.0 < calibrated_probabilities) & (calibrated_probabilities <= 1.0 + 1e-5)] = 1.0

        return calibrated_probabilities


    def _estimate_n_splits_CV(self, CV_k_folds=10):
        """Utility method which makes sure that there are enough observations per
        fold in a cross-validation scheme.

        :param CV_k_folds: [int] Number of desired k-fold cross-validation splits.

        :return [int] number of allowed k-fold cross-validation splits given the
            composition of the data
        """
        # Set the number of folds for cross-validation
        _, n_observations_per_class = np.unique(self.y, return_counts=True)
        return min(n_observations_per_class.min(), CV_k_folds)


    def _tune_hyperparameters(self):
        """Private method used to optimize the hyper-parameters of the base
        estimator.
        """
        if self._verbose:
            log('Random Forest optimization started')

        # Set how to score performance
        scoring = make_scorer(self.compute_AUC, greater_is_better=True, needs_proba=True, needs_threshold=False, metric='PR', return_average=True)
        # Get number of stratified splits that one can do with these labels
        n_splits = self._estimate_n_splits_CV(self.optimization_CV_k_folds)
        # Initialize output variable
        optimized_parameters = dict()

        # 1. Tune the number of candidate features for splitting
        if self.n_features == 1:
            optimized_parameters['max_features'] = None
        else:
            # Initialize partition object and a RandomForestClassifier
            CV_partition = StratifiedKFold(n_splits=n_splits, shuffle=False)
            RF = RandomForestClassifier(n_estimators=300)  # use many trees to start with
            RF.set_params(**self.default_clf_parameters)
            RF.set_params(min_samples_leaf=1)  # Maximal depth
            params_grid = dict(max_features=[0.33, 'sqrt', 'log2'])  # Set rules-of-thumb
            # Perform an exhaustive grid search over the parameter grid
            random_search = GridSearchCV(RF, param_grid=params_grid, scoring=scoring, cv=CV_partition, refit=False, iid=False, n_jobs=N_JOBS).fit(self.X, self.y)
            # Pick the best number of features
            optimized_parameters['max_features'] = random_search.best_params_['max_features']
        if self._verbose:
            log('Random Forest optimization: max_features \'%s\'' % optimized_parameters['max_features'])

        # 2. Tune the depth of trees
        # Initialize partition object and a RandomForestClassifier
        CV_partition = StratifiedKFold(n_splits=n_splits, shuffle=False)
        RF = RandomForestClassifier(n_estimators=300)
        RF.set_params(**self.default_clf_parameters)
        params_grid = dict(min_samples_leaf=[1, 2, 4], min_samples_split=[2, 4, 10])  # Set values that are common in the literature
        # Perform an exhaustive grid search over the parameter grid
        random_search = GridSearchCV(RF, param_grid=params_grid, scoring=scoring, cv=CV_partition, refit=False, iid=False, n_jobs=N_JOBS).fit(self.X, self.y)
        # Pick the best number of samples per leaf and splits to perform by each
        # decision tree
        optimized_parameters['min_samples_leaf'] = random_search.best_params_['min_samples_leaf']
        optimized_parameters['min_samples_split'] = random_search.best_params_['min_samples_split']
        if self._verbose:
            log('Random Forest optimization: min_samples_leaf \'%s\'' % optimized_parameters['min_samples_leaf'])
            log('Random Forest optimization: min_samples_split \'%s\'' % optimized_parameters['min_samples_split'])

        # 3. Tune the number of trees
        # Initialize partition object and a RandomForestClassifier
        CV_partition = StratifiedKFold(n_splits=n_splits, shuffle=False)
        RF = RandomForestClassifier()
        RF.set_params(**self.default_clf_parameters)
        # Set a range of number of trees that starts from 300 (what we have used
        # in step 1) and goes in increments of 100 until reaching a maximum value
        # of 200 trees per class.
        min_n_trees = 300
        n_classes = np.unique(self.y).shape[0]
        max_n_trees = (n_classes + 1) * 200
        params_grid = dict(n_estimators=np.arange(min_n_trees, max_n_trees, 100, dtype=int))
        # Perform an exhaustive grid search over the parameter grid
        random_search = GridSearchCV(RF, param_grid=params_grid, scoring=scoring, cv=CV_partition, refit=False, iid=False, n_jobs=N_JOBS).fit(self.X, self.y)
        # Pick the best number of trees
        optimized_parameters['n_estimators'] = random_search.best_params_['n_estimators']
        if self._verbose:
            log('Random Forest optimization: n_estimators \'%i\'' % optimized_parameters['n_estimators'])

        # Store optimized parameters
        self.optimized_parameters = optimized_parameters

        # Toggle off internal flag
        self._clf_params_need_tuning = False


    def compute_performance(self, class_names=None, performance_estimate='AUPRC', only_significant=True):
        """Compute the performance of the classifier.

        :param class_names: [list or None] Names of classes so that labels used
            internally by algorithm are now translated to meaningful names.
        :param performance_estimate: [list or str] A list or a name of a performance
            estimate, such as 'accuracy', 'AUROC' (the area under the ROC curve),
            or 'AUPRC' (the area under the precision-recall curve).
        :param only_significant: [bool] Whether to consider only the performance
            of the significantly classified samples. Note that these also include
            wrongly but confidently classified samples.

        :return performance: [dict] Dictionary of performance estimates, stored
            both per class and the average across classes.
        """
        # Initialize output variable
        performance = dict()
        # Make a dictionary to convert labels to strings
        class_label_to_name = dict({key: val for key, val in zip(self.classes, class_names)})

        # Make sure that performance_estimate is iterable
        if not isinstance(performance_estimate, (list, np.ndarray)):
            performance_estimate = [performance_estimate]

        # Translate labels to class names
        y_true = self.y
        y_true_label = np.array([class_label_to_name[i] for i in y_true])
        y_pred = np.array([class_label_to_name[i] for i in self.predicted_labels])
        y_prob = self.posterior_probabilities
        # Keep only observations that have been confidently classified
        if self._computed_p_values and only_significant:
            indices = np.where(self.p_value <= 0.05)[0]  # Note that p-values are not corrected, but here I do not want to be too conservative
            y_true = y_true[indices]
            y_true_label = y_true_label[indices]
            y_pred = y_pred[indices]
            y_prob = y_prob[indices, :]

        # Compute confusion matrix
        cm = ConfusionMatrix(actual_vector=y_true_label, predict_vector=y_pred)
        # Format confusion matrix as a list of lists that can reconstruct a table
        conf_mat = [[''] + list(cm.table.keys())]
        for key, value in cm.table.items():
            row = [key] + list(value.values())
            conf_mat.append(row)
        performance['confusion_matrix'] = conf_mat

        # Iterate through the performance estimates chosen by the user
        if 'accuracy' in performance_estimate:
            # Accuracy
            performance['accuracy'] = {key: value for key, value in cm.ACC.items()}
            performance['accuracy_mean'] = cm.Overall_ACC
            # Compute p-value of accuracy as in caret's ConfusionMatrix.R
            y = np.vstack(conf_mat)[1:, 1:].astype(int)
            no_information_rate = (y.sum(axis=0) / y.sum()).max()
            performance['accuracy_p_value'] = binom_test(np.diag(y).sum(), y.sum(), no_information_rate, alternative='greater')

        if 'AUROC' in performance_estimate:
            # Area under the ROC curve
            performance['AUROC'] = {key: value for key, value in cm.AUC.items()}
            performance['AUROC_mean'] = cm.AUNP

        if 'AUPRC' in performance_estimate:
            # Area under the precision-recall curve
            AUPRC = self.compute_AUC(y_true, y_prob, metric='PR', return_average=False)
            performance['AUPRC'] = {key: value for key, value in zip(class_names, AUPRC)}
            performance['AUPRC_mean'] = AUPRC.mean()

        return performance


    # def compute_feature_importance(self, only_correctly_classified=True, only_significant=True):
    #     if not self._is_fit:
    #         raise Exception('Can only compute feature importance on fit classifier')
    #
    #     if not self._computed_p_values:
    #         only_significant = False
    #     # Get indices of correctly classified bins
    #     if only_correctly_classified:
    #         indices_cor = np.where(self.predicted_labels == self.y)[0]
    #     else:
    #         indices_cor = np.arange(self.n_samples)
    #     # Get indices of significant bins
    #     if only_significant:
    #         indices_sig = np.where(self.p_value <= 0.05)[0]
    #     else:
    #         indices_sig = np.arange(self.n_samples)
    #
    #     # Select data for analysis
    #     indices = np.intersect1d(indices_cor, indices_sig)
    #     data = self.X[indices, :]
    #     labels = self.y[indices]
    #     # Scale data
    #     data = StandardScaler(with_mean=True, with_std=True, copy=True).fit_transform(data)
    #
    #     # Initialize MultiSurf to select all features from continuous predictors
    #     ms = MultiSURF(self.n_features, self.n_samples + 1, n_jobs=N_JOBS)
    #     # Run algorithm
    #     ms.fit(data, labels)
    #     # Make DataFrame and sort features by importance
    #     feature_importance = pd.DataFrame(ms.feature_importances_, columns=['importance'])
    #     feature_importance['feature'] = np.arange(feature_importance.shape[0])
    #     feature_importance.sort_values(by='importance', ascending=False, inplace=True)
    #     feature_importance.reset_index(drop=True, inplace=True)
    #
    #     return feature_importance

    # def perform_recursive_feature_selection(self):
    #     pass

    # def compute_feature_contribution(self, n_repeats=5, n_splits=2, only_significant=False):
    #     if not self._is_fit:
    #         raise Exception('Can only compute feature contribution on fit classifier')
    #
    #     if not self._computed_p_values:
    #         only_significant = False
    #     # Get data
    #     if only_significant:
    #         indices = np.where(self.p_value <= 0.05)[0]
    #     else:
    #         indices = np.arange(self.n_samples)
    #     n_indices = indices.shape[0]
    #     data = self.X[indices, :]
    #     labels = self.y[indices]
    #
    #     # Allocate variables
    #     contributions_training = np.zeros((n_indices, self.n_features, self.n_classes))
    #     contributions_testing = np.zeros((n_indices, self.n_features, self.n_classes))
    #     n_trained = np.zeros((n_indices, ), dtype=int)
    #     n_tested = np.zeros((n_indices,), dtype=int)
    #     # Extract classifier from pipeline and initialize cross-validation scheme
    #     classifier = clone(self._base_estimator.named_steps['classifier'])
    #     CV_partition = RepeatedStratifiedKFold(n_repeats=n_repeats, n_splits=n_splits)
    #     for train_indices, test_indices in CV_partition.split(data, labels):
    #         # Normalize data to unit variance
    #         X_train = data[train_indices, :]
    #         scl = StandardScaler(with_mean=True, with_std=True, copy=True).fit(X_train)
    #         X_train = scl.transform(X_train)
    #         X_test = scl.transform(data[test_indices, :])
    #         # Compute contribution of each feature to each observation in the training set
    #         classifier.fit(X_train, labels[train_indices])
    #         _, _, c = ti.compute_feature_contributions_ensemble(classifier, X_train, joint_contribution=False, use_original_implementation=False)
    #         contributions_training[train_indices, :, :] += c
    #         n_trained[train_indices] += 1
    #
    #         # Compute contribution of each feature to each observation in the training set
    #         _, _, c = ti.compute_feature_contributions_ensemble(classifier, X_test, joint_contribution=False, use_original_implementation=False)
    #         contributions_testing[test_indices, :, :] += c
    #         n_tested[test_indices] += 1
    #
    #     # Divide each value by the number of times it has been used for estimation
    #     contributions_training /= n_trained.reshape(-1, 1, 1)
    #     contributions_testing /= n_tested.reshape(-1, 1, 1)
    #
    #     return contributions_training, contributions_testing


    # def cluster_significant_observations(self, class_names):
    #     if not self._computed_p_values:
    #         raise Exception('Cannot cluster observations before computing p-values')
    #
    #     raise NotImplementedError


    # def estimate_feature_selectivity(self, contributions, class_names, reference_class='SP'):
    #     raise NotImplementedError
    #
    #     indices = np.where(self.p_value <= 0.05)[0]
    #     data = self.X[indices, :]
    #     labels = self.y[indices]
    #     # Get size of data
    #     n_features = data.shape[1]
    #     n_classes = np.unique(labels).shape[0]
    #     # Split data in positive and negative contributions, and look at the activity
    #     # during these observations
    #     avg_activity = np.zeros((n_features, n_classes, 2))  # 2 layers for positive and negative contributions
    #     avg_activity_p = np.zeros((n_features, n_classes))  # p-values
    #     for i_cell in range(n_features):
    #         for ci in range(n_classes):
    #             # Values for positive contributions
    #             idx = np.where(contributions[:, i_cell, ci] > 0)[0]
    #             act_pos = data[idx, i_cell]
    #             # Values for negative contributions
    #             idx = np.where(contributions[:, i_cell, ci] <= 0)[0]
    #             act_neg = data[idx, i_cell]
    #             # Concatenate data and perform Wilcoxon rank-sum test against
    #             # null-hypothesis of median = 0
    #             Wrs = ranksums(act_pos, act_neg)
    #             avg_activity_p[i_cell, ci] = getattr(Wrs, 'pvalue')
    #             avg_activity[i_cell, ci, 0] = np.mean(act_pos)
    #             avg_activity[i_cell, ci, 1] = np.mean(act_neg)
    #     # Apply Bonferroni correction on p-values
    #     avg_activity_p *= n_classes
    #
    #     # Stack results in an array
    #     selectivity_strength = np.hstack((avg_activity[:, :, 0] - avg_activity[:, :, 1], avg_activity_p, np.zeros(n_features).reshape(-1, 1))).astype(object)
    #     selectivity_column = n_classes * 2
    #     # Get indices of reference and studied classes
    #     ref_class_idx = class_names.index(reference_class)
    #     stim_class_idx = np.setdiff1d(range(n_classes), ref_class_idx)
    #     # Loop through features
    #     for i_cell in range(n_features):
    #         # Get which classes were significantly contributed to
    #         signif = np.zeros((n_classes, 2), dtype=bool)
    #         for ci in range(n_classes):
    #             signif[ci, 0] = selectivity_strength[i_cell, ci] > 0
    #             signif[ci, 1] = selectivity_strength[i_cell, ci + n_classes] <= 0.05
    #
    #         # Preference toward reference class
    #         ref_strength = ''
    #         if selectivity_strength[i_cell, ref_class_idx] > 0 and selectivity_strength[i_cell, ref_class_idx + n_classes] <= 0.05:
    #             ref_strength = 'prefers %s' % reference_class
    #         elif selectivity_strength[i_cell, ref_class_idx] <= 0 and selectivity_strength[i_cell, ref_class_idx + n_classes] <= 0.05:
    #             ref_strength = 'prefers stim'
    #         elif selectivity_strength[i_cell, ref_class_idx + n_classes] > 0.05:
    #             ref_strength = 'no %s pref' % reference_class
    #
    #         # Preference toward stimuli
    #         stim_strength = np.empty((len(stim_class_idx), ), dtype=object)
    #         for idx, i in enumerate(stim_class_idx):
    #             if selectivity_strength[i_cell, i] > 0 and selectivity_strength[i_cell, i + n_classes] <= 0.05:
    #                 stim_strength[idx] = 'exc'
    #             elif selectivity_strength[i_cell, i] <= 0 and selectivity_strength[i_cell, i + n_classes] <= 0.05:
    #                 stim_strength[idx] = 'inh'
    #             elif selectivity_strength[i_cell, i + n_classes] > 0.05:
    #                 stim_strength[idx] = 'non-sel'
    #
    #         # Combine information from reference and other classes to assess selectivity of this feature
    #         if np.all(stim_strength == 'non-sel'):
    #             selectivity_strength[i_cell, selectivity_column] = 'non-selective'
    #         elif np.all(stim_strength == 'exc'):
    #             if ref_strength == 'prefers %s' % reference_class:
    #                 selectivity_strength[i_cell, selectivity_column] = 'weak exc, mixed'
    #             else:
    #                 selectivity_strength[i_cell, 6] = 'strong exc, mixed'
    #         elif np.any(stim_strength == 'exc'):
    #             stims = stim_class_idx[np.where(stim_strength == 'exc')[0]]
    #             if stims.shape[0] == 1:
    #                 stim_name = class_names[stims[0]]
    #             else:
    #                 stim_name = 'mixed'
    #             if ref_strength == 'prefers stim' or ref_strength == 'non-sel':
    #                 selectivity_strength[i_cell, selectivity_column] = 'strong exc, %s' % stim_name
    #             else:
    #                 selectivity_strength[i_cell, selectivity_column] = 'weak exc, %s' % stim_name
    #
    #         elif np.all(stim_strength == 'inh'):
    #             if ref_strength == 'prefers %s' % reference_class:
    #                 selectivity_strength[i_cell, selectivity_column] = 'non-selective'
    #             else:
    #                 selectivity_strength[i_cell, 6] = 'strong inh, mixed'
    #
    #         elif np.any(stim_strength == 'inh'):
    #             stims = stim_class_idx[np.where(stim_strength == 'inh')[0]]
    #             if stims.shape[0] == 1:
    #                 stim_name = class_names[stims[0]]
    #             else:
    #                 stim_name = 'mixed'
    #             if ref_strength == 'prefers %s' % reference_class:
    #                 selectivity_strength[i_cell, selectivity_column] = 'strong inh, %s' % stim_name
    #             else:
    #                 selectivity_strength[i_cell, selectivity_column] = 'weak inh, %s' % stim_name
    #
    #     mean_contribution = np.median(contributions, axis=0)
    #     results = pd.DataFrame(np.hstack((mean_contribution, selectivity_strength[:, :n_classes], selectivity_strength[:, selectivity_column].reshape(-1, 1))), columns=['c %s' % i for i in class_names] + ['a %s' % i for i in class_names] + ['selectivity']).sort_values(by='selectivity').reset_index()
    #
    #     return results

    # Own implementation of ROC analysis
    def compute_AUC(self, y, y_proba, **kwargs):
        """Compute the ara under the curve (AUC) of the receiving-operator
        characteristics (ROC) curve or of the precision-call curve, while comparing
        the data in two distributions. Distributions can have unequal size, as this
        is taken into account when computing the performance of the separation
        between the two.

        :param y: [numpy array] True labels.
        :param y_proba: [numpy array] Predicted probabilities of each sample to
            belong to each class.
        :param kwargs: [dict] Additional arguments. Supported inputs are:
            `return_average` [bool]: Whether to return the average precision score
                of all classes (one value). Default is True to allow scoring
                functions to use this method for optimization. If user wants to
                return the value for each class, this parameter should be False.
            `metric` [str]: Either 'ROC' or 'PR' (for precision-recall). Default
                is 'PR'.

        :return: The value of AUC minus the performance that a random classifier would
            achieve. This means that AUC can range [-1, 1], with negative values
            corresponding to the positive class having a mean lower than the negative
            class, and vice versa when AUC is positive.
        """
        # Unpack inputs
        metric = kwargs.pop('metric', 'PR')
        return_average = kwargs.pop('return_average', True)

        # Make sure inputs are 2D
        if self.n_classes == 2 and y_proba.ndim != 2:
            y_proba = y_proba.reshape(-1, 1)  # Reshape to column vector
            # Concatenate with probabilities for other class. According to
            # sklearn.metrics.scorer._ProbaScorer.call():
            # if y_type == "binary":
            #     y_pred = y_pred[:, 1]
            # This means that we keep only the probability for the second class.
            # Restore probabilities for the first class as 1-p.
            y_proba = np.hstack((1 - y_proba, y_proba))

        # Set performance of random classifier for the binary classification task
        r = 1 / 2

        # Compute AUC score per class
        AUC = np.zeros((self.n_classes,), dtype=float)
        for idx, label in enumerate(self.classes):
            # Separate probabilities for this and other classes
            p_this_class = y_proba[y == label, idx]
            p_not_this_class = y_proba[y != label, idx]

            # Get number of elements in each group
            n_pos_observations = p_this_class.size
            n_neg_observations = p_not_this_class.size
            # Determine accuracy of a random classifier
            random_classifier_AUC = divide0(n_pos_observations, n_pos_observations + n_neg_observations, replace_with=0)

            # Combine all observations
            data = np.hstack((p_this_class.ravel(), p_not_this_class.ravel()))
            # Calculate the threshold values between data points
            s_data = np.sort(np.unique(data))
            d_data = np.diff(s_data)
            # If there are insufficient data points, return the AUC of a random
            # classifier
            if d_data.size == 0:
                AUC[idx] = random_classifier_AUC
                continue

            # Compute range of thresholds
            d_data = np.hstack((d_data, d_data[-1]))
            thresholds = np.hstack((s_data[0] - d_data[0], s_data + d_data / 2))

            # Calculate hits and misses
            TP = np.sum((p_this_class.ravel()[:, None] >= thresholds).astype(int), axis=0)
            FP = np.sum((p_not_this_class.ravel()[:, None] >= thresholds).astype(int), axis=0)
            # Compute the rest of the confusion matrix
            FN = n_pos_observations - TP  # False negatives

            if metric == 'ROC':
                TN = n_neg_observations - FP  # True negatives
                # Compute the area under the ROC curve in the ROC space
                # https://en.wikipedia.org/wiki/Receiver_operating_characteristic
                TPR = divide0(TP, TP + FN, 0)  # true positive rate
                FPR = divide0(FP, FP + TN, 0)  # false positive rate
                AUC[idx] = np.abs(np.trapz(x=FPR, y=TPR))

            elif metric == 'PR':
                precision = divide0(TP, TP + FP, 0)
                recall = divide0(TP, TP + FN, 0)
                AUC[idx] = np.abs(np.trapz(x=recall, y=precision))

            # Normalize AUC to [0, 1] interval
            AUC[idx] = r + (AUC[idx] - random_classifier_AUC) * (1 - r) / (np.abs(random_classifier_AUC - r) + r)

        # Compute average, if user requested it
        if return_average:
            AUC = AUC.mean()

        return AUC


def divide0(a, b, replace_with):
    """Divide two numbers but replace its result if division is not possible,
    e.g., when dividing a number by 0. No type-checking or agreement between
    dimensions is performed. Be careful!

    :param a: Numerator.
    :param b: Denominator.
    :param replace_with: Return this number if a/b is not defined.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        if isinstance(c, np.ndarray):
            c[np.logical_not(np.isfinite(c))] = replace_with
        else:
            if not np.isfinite(c):
                c = replace_with

    return c


################################################################################
### Statistics
################################################################################
# def compare_classifiers_performance(y_true, y_mod1, y_mod2):
#     tb = mcnemar_table(y_target=y_true, y_model1=y_mod1, y_model2=y_mod2)
#     _, p = mcnemar(ary=tb, corrected=True)
#     return p
