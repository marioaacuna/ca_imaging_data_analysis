# System packages
import sys
import warnings

# Numerical packages
import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.calibration import CalibratedClassifierCV, clone
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold, StratifiedShuffleSplit, train_test_split
from sklearn.metrics import make_scorer
from scipy.stats import randint
import pycm

# Parallel computing
from joblib import Parallel, delayed
# Suppress all warnings
warnings.filterwarnings('ignore')


N_JOBS = -1  # number of cores used for parallel computations


class EnsembleTreeClassifier(object):
    def __init__(self, ensemble_type='RandomForest', tree_split_criterion='gini', verbose=False):
        """
        :param tree_split_criterion: [str] Refers to 'criterion' parameter of the
            RandomForestClassifier.
        :param verbose: [bool] Flag indicating whether to print messages to console.
        """
        # Unpack inputs
        self.ensemble_type = ensemble_type.lower()
        if self.ensemble_type == 'randomforest':
            self._base_estimator_class = RandomForestClassifier
        elif self.ensemble_type == 'extratrees':
            self._base_estimator_class = ExtraTreesClassifier
        else:
            raise AttributeError('\'ensemble_type\' can only be one of \'RandomForest\' or \'ExtraTrees\'')
        self._default_clf_parameters = dict(criterion=tree_split_criterion,
                                            bootstrap=True, oob_score=False,
                                            class_weight='balanced_subsample',
                                            n_jobs=N_JOBS)
        self._fallback_optimized_parameters = dict()
        self._verbose = verbose

        # Initialize internal states
        self.is_fit = False
        self._computed_feature_contributions = False
        self._computed_feature_importances = False

        # Initialize all other attributes
        self.trained_classifiers = None
        self._base_estimator = None
        self.X = None
        self.y = None
        self.n_samples = None
        self.n_features = None
        self.classes = None
        self.n_classes = None
        self.class_names = None
        self.posterior_probabilities = None
        self.crossval_posterior_probabilities = None
        self.posterior_probabilities = None
        self.predicted_labels = None
        self.sig_predicted_labels = None
        self.crossval_performance = None
        self.performance = None
        self.performance_mean = None
        self.all_confusion_matrix = None
        self.confusion_matrix = None
        self.sig_crossval_performance = None
        self.sig_performance = None
        self.sig_performance_mean = None
        self.sig_all_confusion_matrix = None
        self.sig_confusion_matrix = None
        self.p_values = None
        self.correct_classifications = None
        self.correct_classifications_perc = None
        self.train_indices = None
        self.test_indices = None


    def fit_predict(self, data, labels, class_names=None,
                    CV_scheme='k-fold',
                    fraction_training=0.8,
                    max_n_splits_training=20, v_fraction=0.1, shuffle_split_repeats=5,
                    fraction_training_for_tuning=0,
                    max_n_splits_tuning=10,
                    fraction_training_for_calibration=0,
                    max_n_splits_calibration=5,
                    calibration_method='sigmoid',
                    data_shape=None, data_scaling=None,
                    n_permutations=0):
        """Public method to fit and predict data via cross-validation splitting.

        :param data: [2D numpy array] Contains the observations in Tidy format,
            which is observations on rows and features on columns.
        :param labels: [1D numpy array] The true label of each observation.
        :param class_names: [list or None] Names of classes.

        :param CV_scheme: [str] Either 'k-fold', 'hv-block' or 'balanced k-fold'.
        :param v_fraction: [float] Fraction of data to discard in the 'hv-block'
            scheme.
        :param shuffle_split_repeats: [int] Number of splits to perform with the
            'shuffle split' CV scheme.
        :param fraction_training: [float] Fraction of data to use
            for training in the 'shuffle split' scheme.
        :param max_n_splits_training:
        :param fraction_training_for_tuning: [float] Out of the training set, how
            much should be allocated for hyperparameter tuning. The rest will be
            left for fitting (and possibly calibrating) the model.
        :param max_n_splits_tuning:
        :param fraction_training_for_calibration: [float] Out of the training set,
            how much should be allocated for probability calibration. The rest will
            be left for fitting the model.
        :param max_n_splits_calibration:

        :param calibration_method: [str] The default method to calibrate posterior
            probabilities output by RandomForestClassifier is 'sigmoid', which
            corresponds to Platt's method. The alternative is 'isotonic', which
            corresponds to a non-parametric method. If None, no calibration is
            performed.

        :param data_shape: [dict] Original data shape of each class' data, indicated
            as [n_bins, n_trials].
        :param data_scaling: [None or str] Whether to standardize data for training.
            No scaling will be performed if None; 'minmax' will apply the Min-Max
            scaler, and 'standard' to standardize the data.

        :param n_permutations: [int] Number of permutations to compute null
            distribution and p-value of each posterior probability.
        """
        # Check inputs
        if CV_scheme not in ['k-fold', 'hv-block', 'balanced k-fold', 'shuffle split']:
            raise AttributeError('\'CV_scheme\' has to be one of [\'%s\']' % '\', \''.join(['k-fold', 'hv-block', 'balanced k-fold']))
        if fraction_training < 0 or fraction_training >= 1:
            raise ValueError('\'fraction_training\' should be between 0 and 1 (1 excluded)')
        if fraction_training_for_tuning < 0 or fraction_training_for_tuning >= 1:
            raise ValueError('\'fraction_training_for_tuning\' should be between 0 and 1 (1 excluded)')
        if fraction_training_for_calibration < 0 or fraction_training_for_calibration >= 1:
            raise ValueError('\'fraction_training_for_calibration\' should be between 0 and 1 (1 excluded)')

        # Set which splits will be performed
        do_hyperparameter_tuning = fraction_training_for_tuning > 0
        do_calibration = calibration_method is not None and fraction_training_for_calibration > 0

        # If no tuning will be performed, set rule-of-thumb parameters
        self._fallback_optimized_parameters = dict()
        if not do_hyperparameter_tuning:
            self._fallback_optimized_parameters['classifier__max_features'] = int(np.clip(np.floor(np.log2(data.shape[1] + 1)), 1, data.shape[1]))
            self._fallback_optimized_parameters['classifier__min_samples_leaf'] = 1
            self._fallback_optimized_parameters['classifier__min_samples_split'] = 2
            self._fallback_optimized_parameters['classifier__n_estimators'] = 128

        # Count number of classes
        classes = np.unique(labels)
        n_classes = len(classes)
        if class_names is not None:
            if not isinstance(class_names, (list, np.ndarray)):
                class_names = [class_names]
            if len(class_names) != n_classes:
                raise ValueError('\'class_names\' should be as long as the number of classes in \'labels\'')
            if not all([type(i) == str for i in class_names]):
                raise ValueError('\'class_names\' should contain only strings')
        self.class_names = class_names
        # Copy data and labels in object attributes
        self.X = data
        self.y = labels
        # Get the structure of data and labels
        self.classes = classes
        self.n_classes = n_classes
        self.n_samples = self.X.shape[0]
        self.n_features = self.X.shape[1]
        # Initialize classifier with default parameters
        base_estimator = self._base_estimator_class(**self._default_clf_parameters)

        # Initialize local variables
        CV_partition = None
        inner_CV_partition = None
        tuning_CV_partition = None

        # Get number of splits allowed to do
        fraction_training_outer_CV, n_splits = self._estimate_fraction_training_CV(labels=self.y, fraction_training=fraction_training, CV_scheme=CV_scheme)
        # Create partition for cross-validation
        if n_splits >= 2:
            if CV_scheme == 'k-fold':
                CV_partition = StratifiedKFold(n_splits=n_splits, shuffle=False)
            elif CV_scheme == 'hv-block':
                CV_partition = StratifiedHVBlock(v_fraction=v_fraction, n_splits=n_splits, shuffle=False)
            elif CV_scheme == 'balanced k-fold':
                CV_partition = BalancedStratifiedKFold(n_splits=n_splits)
            elif CV_scheme == 'shuffle split':
                CV_partition = StratifiedShuffleSplit(n_splits=np.clip(n_splits * shuffle_split_repeats, 2, max_n_splits_training), train_size=fraction_training_outer_CV)

        else:
            log('Not enough data in at least one class to perform cross-validation. Skipping dataset')
            return

        # Insert classifier in a pipeline together with a StandardScaler object.
        pipeline_named_steps = list()
        if data_scaling is not None:
            scaler = None
            if data_scaling == 'standard':
                scaler = StandardScaler(with_mean=True, with_std=True, copy=True)
            elif data_scaling == 'minmax':
                scaler = MinMaxScaler(feature_range=(0, 1), copy=True)
            pipeline_named_steps.append(('scaler', scaler))
        pipeline_named_steps.append(('classifier', base_estimator))
        self._base_estimator = Pipeline(pipeline_named_steps)

        # Allocate output variables
        classifiers_cv = list()
        AUPRC_cv = np.zeros((n_splits, self.n_classes)) * np.nan
        sig_AUPRC_cv = np.zeros((n_splits, self.n_classes)) * np.nan
        all_confusion_matrix = dict()
        conf_mat = np.zeros((self.n_classes, self.n_classes), dtype=float)
        sig_all_confusion_matrix = dict()
        sig_conf_mat = np.zeros((self.n_classes, self.n_classes), dtype=float)
        posterior_probabilities = np.zeros((self.n_samples, self.n_classes, n_splits)) * np.nan
        p_values = np.ones((self.n_samples, n_splits)) * np.nan
        train_indices = list()
        test_indices = list()

        for i_split, (train_index, test_index) in enumerate(CV_partition.split(self.X, self.y)):
            # Copy fallback optimized parameters
            optimized_parameters = self._fallback_optimized_parameters.copy()

            # Split data and keep indices
            X_train, X_test = self.X[train_index, :], self.X[test_index, :]
            y_train, y_test = self.y[train_index], self.y[test_index]
            train_indices.append(train_index)
            test_indices.append(test_index)

            # Tune hyperparameters
            clf = clone(self._base_estimator)
            if do_hyperparameter_tuning:
                # Check whether there is enough data for tuning
                _, _, y_tuning, _ = train_test_split(X_train, y_train, test_size=1 - fraction_training_for_tuning)
                inner_fraction_training_for_tuning, n_splits_tuning = self._estimate_fraction_training_CV(labels=y_tuning, fraction_training=fraction_training, CV_scheme=CV_scheme)

                if inner_fraction_training_for_tuning > 0 and n_splits_tuning >= 2:
                    # Create partition for cross-validation
                    if CV_scheme == 'k-fold':
                        tuning_CV_partition = StratifiedKFold(n_splits=n_splits_tuning, shuffle=False)
                    elif CV_scheme == 'hv-block':
                        tuning_CV_partition = StratifiedHVBlock(v_fraction=v_fraction, n_splits=n_splits_tuning, shuffle=False)
                    elif CV_scheme == 'balanced k-fold':
                        tuning_CV_partition = BalancedStratifiedKFold(n_splits=n_splits_tuning)
                    elif CV_scheme == 'shuffle split':
                        tuning_CV_partition = StratifiedShuffleSplit(n_splits=np.clip(n_splits_tuning * shuffle_split_repeats, 2, max_n_splits_tuning), train_size=inner_fraction_training_for_tuning)

                    # Cut out a validation set for hyperparameter tuning
                    X_tuning, X_train_split, y_tuning, y_train_split = train_test_split(X_train, y_train, test_size=1 - fraction_training_for_tuning)
                    try:
                        optimized_parameters = self.tune_hyperparameters(data=X_tuning, labels=y_tuning, CV_scheme=tuning_CV_partition)
                        X_train = X_train_split.copy()
                        y_train = y_train_split.copy()
                    except:
                        pass

            # Override default parameters with optimized parameters
            clf.set_params(**optimized_parameters)

            # Calibrate classifier
            if do_calibration:
                inner_fraction_training_for_calibration, n_splits_calibration = self._estimate_fraction_training_CV(labels=y_train, fraction_training=fraction_training_for_calibration, CV_scheme=CV_scheme)
                if inner_fraction_training_for_calibration > 0 and n_splits_calibration >= 2:
                    if self._verbose:
                        log('Calibrating classifier')

                    # Set new cross-validation scheme
                    if CV_scheme == 'k-fold':
                        inner_CV_partition = StratifiedKFold(n_splits=n_splits_calibration, shuffle=False)
                    elif CV_scheme == 'hv-block':
                        inner_CV_partition = StratifiedHVBlock(v_fraction=v_fraction, n_splits=n_splits_calibration, shuffle=False)
                    elif CV_scheme == 'balanced k-fold':
                        inner_CV_partition = BalancedStratifiedKFold(n_splits=n_splits_calibration)
                    elif CV_scheme == 'shuffle split':
                        inner_CV_partition = StratifiedShuffleSplit(n_splits=np.clip(n_splits_calibration * shuffle_split_repeats, 2, max_n_splits_calibration), train_size=inner_fraction_training_for_calibration)

                else:
                    inner_CV_partition = 'prefit'

                # Fit and predict posterior probabilities
                clf = CalibratedClassifierCV(base_estimator=clf, method=calibration_method, cv=inner_CV_partition)

            # Train classifier
            clf.fit(X_train, y_train)
            # Estimate posterior probabilities
            y_proba = clf.predict_proba(X_test)

            # Compute significance with a permutation test
            if n_permutations > 0:
                if self._verbose:
                    log('Permuting labels')

                null_distribution = Parallel(n_jobs=N_JOBS)(
                        delayed(self.__permute_labels)(clf, X_train, y_train, X_test)
                        for _ in range(n_permutations))
                # Stack values along third dimension
                null_distribution = np.dstack(null_distribution)
                # Check the probability of each prediction against the null distribution
                difference = np.array(null_distribution >= y_proba[:, :, None], dtype=int)
                # Count the number of times the relationship is true
                count_per_resample = np.sum(difference, axis=2)
                # Get the count only in the target class
                predicted_labels = self.classes[np.argmax(y_proba, axis=1)]
                count_per_resample = count_per_resample[np.arange(predicted_labels.shape[0]), predicted_labels]
                # Estimate a one-sided p-value
                p_value = (count_per_resample + 1) / (n_permutations + 1)
                p_values[test_index, i_split] = p_value
                # Keep only significant observations
                significant_samples = np.where(p_value <= 0.05)[0]
                if significant_samples.shape[0] > 0:
                    # Compute performance on significant samples
                    sig_AUPRC_cv[i_split, :], cm = self.compute_AUC(y_test[significant_samples], y_proba[significant_samples, :], return_average=False, return_confusion_matrix=True, missing_class_is_nan_AUPRC=True)
                    sig_all_confusion_matrix['split_%i' % (i_split + 1)] = pd.DataFrame(cm, columns=self.class_names, index=self.class_names)
                    cm = divide0(cm, np.sum(cm, axis=1).reshape(-1, 1), replace_with=0)
                    sig_conf_mat = sig_conf_mat + cm

            # Compute performance on all samples
            AUPRC_cv[i_split, :], cm = self.compute_AUC(y_test, y_proba, return_average=False, return_confusion_matrix=True, missing_class_is_nan_AUPRC=True)
            all_confusion_matrix['split_%i' % (i_split + 1)] = pd.DataFrame(cm, columns=self.class_names, index=self.class_names)
            cm = divide0(cm, np.sum(cm, axis=1).reshape(-1, 1), replace_with=0)
            conf_mat = conf_mat + cm

            # Store info on training
            classifiers_cv.append(clf)
            posterior_probabilities[test_index, :, i_split] = y_proba

            # Print progress to console
            if self._verbose:
                # if n_permutations > 0 and significant_samples.shape[0] > 0:
                #     print(sig_AUPRC_cv[i_split, :])
                # else:
                print('Split %i: %s' % (i_split + 1, AUPRC_cv[i_split, :]))

        # Average performance across all cross-validation iterations
        self.trained_classifiers = classifiers_cv
        self.crossval_performance = AUPRC_cv
        self.performance = np.nanmean(self.crossval_performance, axis=0)
        self.performance_mean = np.nanmean(self.performance)
        self.all_confusion_matrix = all_confusion_matrix
        self.confusion_matrix = conf_mat / n_splits
        # Performance on significant samples only
        self.sig_crossval_performance = sig_AUPRC_cv
        self.sig_performance = np.nanmean(self.sig_crossval_performance, axis=0)
        self.sig_performance_mean = np.nanmean(self.sig_performance)
        self.sig_all_confusion_matrix = sig_all_confusion_matrix
        self.sig_confusion_matrix = sig_conf_mat / n_splits

        self.crossval_posterior_probabilities = posterior_probabilities
        p = np.nanmean(self.crossval_posterior_probabilities, axis=2)
        p[np.isnan(p)] = 1 / self.n_classes
        self.posterior_probabilities = p / np.sum(p, axis=1).reshape(-1, 1)
        self.predicted_labels = self.posterior_probabilities.argmax(axis=1)
        self.p_values = p_values
        self.train_indices = train_indices
        self.test_indices = test_indices

        # Add predicted labels for significant observations
        self.sig_predicted_labels = self.predicted_labels * np.nan
        n_folds = len(self.test_indices)
        for i_fold in range(n_folds):
            idx = self.test_indices[i_fold]
            p_values = self.p_values[idx, i_fold]
            p = self.posterior_probabilities[idx, :]
            sig_idx = np.where(p_values < 0.5)[0]
            sig_p = p[sig_idx, :]
            self.sig_predicted_labels[idx[sig_idx]] = sig_p.argmax(axis=1)

        # Count correctly classified classes
        correct_classifications = np.zeros((self.n_classes, ), dtype=int)
        correct_classifications_perc = np.zeros((self.n_classes, ), dtype=float)
        for i_class in range(self.n_classes):
            # Count number of times the classifier correctly identifies class of interest in a trial
            idx_this_class = np.where(self.y == i_class)[0]
            label_per_observation = self.sig_predicted_labels[idx_this_class].reshape(data_shape[i_class, 1], data_shape[i_class, 0])
            this_correct_classifications = np.any(label_per_observation == i_class, axis=1).astype(int).sum()
            # Store results
            correct_classifications[i_class] = this_correct_classifications
            correct_classifications_perc[i_class] = this_correct_classifications / data_shape[i_class, 1] * 100
        # Store results
        self.correct_classifications = correct_classifications
        self.correct_classifications_perc = correct_classifications_perc

        # Toggle internal state
        self.is_fit = True


    def tune_hyperparameters(self, data, labels, CV_scheme):
        """Optimize hyper-parameters of the base estimator.

        :param data:
        :param labels:
        :param CV_scheme:
        """
        if self._verbose:
            log('Tuning hyperparameters')

        # Define a space for each parameter to optimize
        n_classes = np.unique(self.y).shape[0]
        min_n_trees = n_classes
        max_n_trees = 200#(n_classes + 1) * 200
        if self.ensemble_type == 'extratrees':
            min_n_trees *= 2
            max_n_trees *= 3
        n_features = data.shape[1]
        max_features = int(np.ceil(np.sqrt(n_features)))
        if max_features < 1:
            max_features = 1
        parameter_space = dict(n_estimators=randint(min_n_trees, max_n_trees),
                               min_samples_split=randint(2, 10),
                               min_samples_leaf=randint(1, 10))
        if max_features > 1:
            parameter_space['max_features'] = randint(1, max_features)

        # Set how to score performance
        scoring = make_scorer(self.compute_AUC, greater_is_better=True, needs_proba=True, needs_threshold=False, metric='PR', return_average=True)

        # Re-instantiate a copy of the base estimator
        clf = clone(self._base_estimator)
        parameter_space = {'classifier__' + k: v for k, v in parameter_space.items()}

        # Randomly search best hyperparameters
        random_search = RandomizedSearchCV(estimator=clf,
                                           param_distributions=parameter_space,
                                           scoring=scoring,
                                           n_iter=5,
                                           cv=CV_scheme,
                                           iid=False,
                                           n_jobs=N_JOBS,
                                           refit=False)
        random_search.fit(data, labels)

        return random_search.best_params_


    ############################################################################
    # Performance
    ############################################################################
    def compute_AUC(self, y_true, y_proba, **kwargs):
        """Compute the ara under the curve (AUC) of the receiving-operator
        characteristics (ROC) curve or of the precision-call curve, while comparing
        the data in two distributions. Distributions can have unequal size, as this
        is taken into account when computing the performance of the separation
        between the two.

        :param y_true: [numpy array] True labels.
        :param y_proba: [numpy array] Predicted probabilities of each sample to
            belong to each class.
        :param kwargs: [dictionary] Additional arguments. Supported inputs are:
            `return_average` [bool]: Whether to return the average precision score
                of all classes (one value). Default is True to allow scoring
                functions to use this method for optimization. If user wants to
                return the value for each class, this parameter should be False.
            `metric` [str]: Either 'ROC' or 'PR' (for precision-recall). Default
                is 'PR'.
            `return_confusion_matrix` [bool]: Whether to return the confusion
                matrix.
            `target_names` [list]: List of class names, in the same order as they
                appear in `y`.

        :return: The value of AUC minus the performance that a random classifier would
            achieve. This means that AUC can range [0, 1], with values < 0.5
            corresponding to the positive class having a mean lower than the negative
            class, and vice versa when AUC is > 0.5.
        """
        # Unpack inputs
        metric = kwargs.pop('metric', 'PR')
        missing_class_is_nan_AUPRC = kwargs.pop('missing_class_is_nan_AUPRC', False)
        return_average = kwargs.pop('return_average', True)
        return_confusion_matrix = kwargs.pop('return_confusion_matrix', False)

        # Initialize local variables
        conf_mat = np.zeros((self.n_classes, self.n_classes), dtype=int)
        AUC = (np.zeros((self.n_classes, )) - 1).astype(object)

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

        # Get predicted labels
        y_pred = np.argmax(y_proba, axis=1)
        # Get number of classes
        all_classes_analyzed = np.array(np.union1d(y_true, y_pred))
        n_classes = all_classes_analyzed.shape[0]
        # Cannot build a confusion matrix without at least 2 classes
        if n_classes >= 2:
            # Get index of present classes
            present_classes_idx = np.where(np.in1d(self.classes, all_classes_analyzed))[0]

            # Compute confusion matrix
            cm = pycm.ConfusionMatrix(actual_vector=y_true, predict_vector=y_pred)
            # Get metric
            if metric == 'PR':
                metric_values = list(getattr(cm, 'AUPR').values())
            elif metric == 'ROC':
                metric_values = list(getattr(cm, 'AUC').values())
            else:  # Accuracy
                metric_values = list(getattr(cm, 'ACC').values())
            AUC[present_classes_idx] = metric_values
            AUC = np.array([-1 if i in ['None', 0] else i for i in AUC], dtype=float)

            # Format confusion matrix for output
            if return_confusion_matrix:
                cm_matrix = np.vstack([list(i.values()) for i in cm.matrix.values()])
                # Format confusion matrix
                for row, idx in enumerate(present_classes_idx):
                    conf_mat[idx, present_classes_idx] = cm_matrix[row, :]

        # Set performance of random classifier for the binary classification task
        # r = 1 / 2
        #
        # # Compute AUC score per class
        # AUC = np.zeros((self.n_classes,), dtype=float)
        # for idx, label in enumerate(self.classes):
        #     # Get probabilities for this class
        #     y_score = y_proba[:, idx]
        #     # If this stimulus was not analyzed by the classifier, all probabilities will be equal to random
        #     if np.all(y_score == 0):
        #         AUC[idx] = r
        #         continue
        #
        #     # Area under the ROC curve
        #     if metric == 'ROC':
        #         fpr, tpr, _ = roc_curve(y_true, y_score, pos_label=label, drop_intermediate=True)
        #         AUC[idx] = auc(fpr, tpr)
        #
        #     # Area under the precision-recall curve
        #     elif metric == 'PR':
        #         precision, recall, _ = precision_recall_curve(y_true, y_score, pos_label=label)
        #         AUC[idx] = auc(recall, precision)
        #
        #         # Normalize AUC to [0, 1] interval
        #         if normalize_AUPRC_to_random:
        #             # Get number of elements in each group
        #             n_pos_observations = np.where(y_true == label)[0].shape[0]
        #             n_neg_observations = np.where(y_true != label)[0].shape[0]
        #             # Compute performance of random classifier
        #             random_classifier_AUC = divide0(n_pos_observations, n_pos_observations + n_neg_observations, replace_with=0)
        #             # Normalize value
        #             AUC[idx] = r + (AUC[idx] - random_classifier_AUC) * (1 - r) / (np.abs(random_classifier_AUC - r) + r)
        #
        #     # Accuracy
        #     elif metric == 'accuracy':
        #         n_observations_this_class = np.where(y_true == label)[0].shape[0]
        #         AUC[idx] = np.where(np.logical_and(y_true == label, y_pred == y_true))[0].shape[0] / n_observations_this_class

        # Replace missing class values with NaNs
        if missing_class_is_nan_AUPRC:
            AUC[AUC == -1] = np.nan

        # Compute average, if user requested it
        if return_average:
            AUC = np.nanmean(AUC)

        if return_confusion_matrix:
            return AUC, conf_mat
        else:
            return AUC


    ############################################################################
    # Utilities
    ############################################################################
    @staticmethod
    def _estimate_fraction_training_CV(labels, fraction_training, CV_scheme):
        """Utility method which makes sure that there are enough observations per
        fold in a cross-validation scheme.

        :param labels:
        :param fraction_training: [float] Fraction of samples to allocate to
            training.
        :param CV_scheme: [str] Name of CV scheme to apply.

        :return [float] Max fraction
        """
        # Get count in each class
        _, class_counts = np.unique(labels, return_counts=True)
        min_class_count = min(class_counts)

        # Compute max number of k-fold splits allowed
        fraction_test = 1 - fraction_training
        max_n_splits = int(np.floor(1 / fraction_test))

        # Check whether a quality criterion can be satisfied , according to the chosen CV scheme
        if CV_scheme == 'shuffle split':
            did_satisfy_criterion = min_class_count >= 2
        else:
            n_samples_per_fold = int(np.floor(min_class_count / max_n_splits))
            did_satisfy_criterion = n_samples_per_fold >= 1

        # Choose what to return
        if did_satisfy_criterion:
            return fraction_training, max_n_splits
        else:
            return 0, max_n_splits


    @ staticmethod
    def __permute_labels(clf, X_train, y_train, X_test):
        # Shuffle labels without replacement (i.e., permute labels)
        indices = np.arange(y_train.shape[0])
        np.random.shuffle(indices)
        indices = indices[:y_train.shape[0]]
        y_train_permuted = y_train[indices]
        # Train classifier
        this_clf = clone(clf)
        this_clf.fit(X_train, y_train_permuted)

        # Return posterior probabilities
        return this_clf.predict_proba(X_test)


################################################################################
# Implementation of hv-block cross-validation scheme
################################################################################
class StratifiedHVBlock(StratifiedKFold):
    def __init__(self, v_fraction, n_splits=5, shuffle=False, random_state=None):
        super(StratifiedHVBlock, self).__init__(n_splits=n_splits, shuffle=shuffle, random_state=random_state)
        self.v_fraction = v_fraction

    @staticmethod
    def _make_indices(start, end):
        lens = end - start
        np.cumsum(lens, out=lens)
        i = np.ones(lens[-1], dtype=int)
        i[0] = start[0]
        i[lens[:-1]] += start[1:]
        i[lens[:-1]] -= end[:-1]
        np.cumsum(i, out=i)

        return i

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, of length n_samples
            The target variable for supervised learning problems.

        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.

        Yields
        ------
        train : ndarray
            The training set indices for that split.

        test : ndarray
            The testing set indices for that split.
        """
        n_samples = X.shape[0]
        classes = np.unique(y)
        indices = np.arange(n_samples)
        for test_index in self._iter_test_masks(X, y, groups):
            train_index = indices[np.logical_not(test_index)]
            test_index = indices[test_index]
            # Initialize list with all train indices
            new_train_index = list()

            # For each class, take the samples that belong to each class
            for i_class in classes:
                samples_this_class = np.where(y == i_class)[0]
                train_index_this_class = np.array(np.intersect1d(train_index, samples_this_class))
                test_index_this_class = np.array(np.intersect1d(test_index, samples_this_class))
                # Find first and last element of this train block
                edges = idx2range(train_index_this_class)
                # Calculate number of samples to discard
                n_samples_in_this_test_set = test_index_this_class.shape[0]
                n_samples_to_discard = int(np.ceil(n_samples_in_this_test_set * self.v_fraction))
                # Shrink training blocks
                edges[:, 0] += n_samples_to_discard
                edges[:, 1] -= n_samples_to_discard
                # Make sure that there is at least one observation per group
                collapsed_blocks = edges[:, 1] <= edges[:, 0]
                if np.any(collapsed_blocks):
                    edges[collapsed_blocks, :2] = np.vstack((edges[collapsed_blocks, 0], edges[collapsed_blocks, 0])).T
                # Make new array of training indices
                new_train_index_this_class = self._make_indices(edges[:, 0], edges[:, 1] + 1)
                # Make sure that these indices only belong to this class
                new_train_index_this_class = np.intersect1d(train_index, new_train_index_this_class)
                # Store values
                new_train_index.append(new_train_index_this_class)

            yield np.hstack(new_train_index), test_index


################################################################################
# Implementation of a balanced k-fold cross-validation scheme
################################################################################
class BalancedStratifiedKFold(object):
    def __init__(self, n_splits=5, random_state=None):
        self.n_splits = n_splits
        self.random_state = random_state

    def get_n_splits(self):
        return self.n_splits

    def split(self, **kwargs):
        y = kwargs.pop('y')

        # Get count of samples in each class
        classes, class_counts = np.unique(y, return_counts=True)
        n_classes = len(classes)

        # Set number of samples in each fold per class
        min_class_count = min(class_counts)
        n_samples_per_fold = int(np.floor(min_class_count / self.n_splits))

        # Create list of indices in each class
        fold_indices = np.empty((n_classes, self.n_splits), dtype=object)
        for i_class in range(n_classes):
            samples_this_class = np.where(y == i_class)[0]
            n_samples_this_class = samples_this_class.shape[0]
            # Check whether all samples will be used in the cross-validation scheme
            all_samples_used = n_samples_per_fold * self.n_splits == n_samples_this_class
            # If not all samples will be used, select a random subsample
            if not all_samples_used:
                n_samples_to_pick = n_samples_per_fold * self.n_splits
                np.random.shuffle(samples_this_class)
                samples_this_class = np.sort(samples_this_class[:n_samples_to_pick])

            # Divide the samples this class in train and test indices, according
            # to a k-fold scheme
            fold_indices[i_class, :] = np.split(samples_this_class, self.n_splits)

        for i_fold in range(self.n_splits):
            # Set which folds will be used for either training or testing
            test_fold = i_fold
            train_folds = np.setdiff1d(np.arange(self.n_splits), test_fold)

            # Concatenate and yield indices
            train_indices = np.hstack(fold_indices[:, train_folds].ravel())
            test_indices = np.hstack(fold_indices[:, test_fold])

            yield train_indices, test_indices


################################################################################
# Utilities
################################################################################
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


def log(message):
    sys.stdout.write(str(message) + '\n')
    sys.stdout.flush()


def idx2range(idx):
    # Convert to numpy array
    if not type(idx) is np.ndarray:
        idx = np.array([idx], dtype=int).ravel()

    if idx.shape[0] > 1:
        # Find discontinuities in index
        dataIDX = np.atleast_2d(np.unique(np.hstack((0, np.where(np.diff(idx) > 1)[0]+1)))).transpose()
        dataIDX = np.hstack((dataIDX, np.atleast_2d(np.hstack((dataIDX[1:,0]-1, idx.shape[0]-1))).transpose()))
        # Get original values
        dataIDX = idx[dataIDX]

        # Add column for duration
        dataIDX = np.hstack((dataIDX, np.atleast_2d(dataIDX[:,1] - dataIDX[:,0] + 1).transpose()))

    else:
        dataIDX = np.empty((0, 3), dtype=int)

    return dataIDX


