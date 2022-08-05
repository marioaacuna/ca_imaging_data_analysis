# -*- coding: utf-8 -*-
# This module contains functions to decode activation periods in multi-trial
# imaging sessions. It takes inputs from MATLAB and returns a JSON file that can
# be read back in MATLAB.
# Threshold p value = 0.05

# Add Utilities folder to system path
import os, sys
script_path = os.path.realpath(__file__)
root_path = os.path.split(os.path.dirname(script_path))[0]
sys.path.insert(0, root_path)

# System packages
import json
from typing import Union, Tuple

# Numerical packages
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
# Parallel computing
from joblib import Parallel, delayed

# Local repository
from Utilities.IO_operations import log
from Utilities import matlab_file


N_JOBS = -1  # number of cores used for parallel computations


def run(input_filename: str, verbose: bool = True) -> None:
    """ Main function called when running this module.

    :param input_filename: The fullpath to the .mat filename containing the
            data and other parameters.
    :param verbose: Whether to print the outcome of some operations.
    :return: A JSON file on disk
    """

    # Set random number generator seed for replicability
    filename = os.path.splitext(os.path.basename(input_filename))[0]
    seed = int(np.clip(np.int64(np.sum([ord(char) for char in filename])), a_min=0, a_max=2**32 - 1))
    np.random.seed(seed)
    log('Random number generator seed set to %i' % seed)

    # Load parameter file
    if verbose:
        log('Analyzing \'%s\'' % input_filename)
    PARAMETERS = matlab_file.load(input_filename)
    # Unpack some parameters
    n_cells = PARAMETERS['data'].shape[1]
    CV_k_folds = PARAMETERS['CV_k_folds']

    # Make sure some parameters are lists
    if not isinstance(PARAMETERS['timestamps'], (list, np.ndarray)):
        PARAMETERS['timestamps'] = np.array([PARAMETERS['timestamps']])
    if not isinstance(PARAMETERS['timestamps_s'], (list, np.ndarray)):
        PARAMETERS['timestamps_s'] = np.array([PARAMETERS['timestamps_s']])

    # Make time axis for this stimulus, staring at the 1st timestamp (stimulus onset)
    n_bins = PARAMETERS['data_size'][0]
    timestamps_s = PARAMETERS['timestamps_s']
    times = np.linspace(PARAMETERS['trial_time'][0], PARAMETERS['trial_time'][1], n_bins + 1)
    times -= timestamps_s[0]

    # Get the observation markers
    bins_to_cumulate = PARAMETERS['bins_to_cumulate']
    bins_to_cumulate[np.isnan(bins_to_cumulate)] = -1
    bins_to_cumulate = bins_to_cumulate.astype(int)
    # Get the data and its original [trial x time] shape
    data = PARAMETERS['data']
    data_shape = PARAMETERS['data_size']

    # Compute selectivity of individual cells
    columns_to_copy = ['bin', 'response_probability', 'AUC']
    performance_all_cells = pd.DataFrame(columns=['cell'] + columns_to_copy)
    for i_cell in range(n_cells):
        X = data[:, i_cell].reshape(-1, 1).copy()
        performance_single_cell = _compute_accuracy_observation_subset(X, bins_to_cumulate, data_shape, subset_type='single',
                                                                       CV_k_folds=CV_k_folds, n_significance_permutation=PARAMETERS['significance_test_n_permutations'])

        # Store data
        performance_single_cell['bin'] = performance_single_cell.index + 1
        significant_activation = np.where(performance_single_cell['AUC'] != 0.5)[0]
        n_significant_activations = significant_activation.shape[0]
        if n_significant_activations > 0:
            for i in range(n_significant_activations):
                row = performance_all_cells.shape[0]
                performance_all_cells.at[row, 'cell'] = i_cell + 1  # 1-indexed
                performance_all_cells.loc[row, columns_to_copy] = performance_single_cell.loc[significant_activation[i], columns_to_copy]

            if verbose:
                log('Cell %i / %i: p(response): %.2f  |  AUC=%.2f' % (i_cell + 1, n_cells, float(np.mean(performance_single_cell.loc[significant_activation, 'response_probability'])), float(np.mean(performance_single_cell.loc[significant_activation, 'AUC']))))
        else:
            if verbose:
                log('Cell %i / %i: Non-selective' % (i_cell + 1, n_cells))

    # Store data in dictionary
    RESULTS = dict()
    performance = list([list(performance_all_cells.columns)])
    performance.append(performance_all_cells.values.tolist())
    RESULTS['performance'] = performance

    # Write results to a text file in JSON format
    with open(PARAMETERS['output_filename'], 'w') as f:
        json.dump(RESULTS, f)

    if verbose:
        log('Finished')


################################################################################
### Statistics
################################################################################
def _compute_accuracy_observation_subset(data: np.ndarray, bins_to_cumulate: np.ndarray,
                                         data_shape: np.ndarray, subset_type: str = 'single',
                                         CV_k_folds: int = 5, n_significance_permutation:int = 1000):
    """This function performs the main job of iterating through time bins and
    computing whether they differ from a baseline period.

    :param data: numpy-array obtained from [n_samples x n_features].ravel()
    :param bins_to_cumulate: It marks each time bin with either -1 (discard), 0
        (baseline) or an int > 0. Bins with the same marker will be grouped.
    :param data_shape: List of 2 items representing the number of bins and trials
        in the original data.
    :param subset_type: A string, either 'single' or 'cumulative'. Whether to
        consider bins individually or cumulate them. Only the 'single' mode is
        fully supported and tested.
    :param CV_k_folds: Number of folds for the stratified k-fold cross-validation.
    :param n_significance_permutation: Number of iterations of the permutation
        test to perform to assess significance.

    :return performance_single_cell [pandas DataFrame]: It contains info on each
        time bin: Whether it was significantly modulated (and on which trials),
        the corresponding response probability, and the activity difference with
        the baseline
    :return significant_activations [pandas DataFrame]: It contains info on each
        significant period of activation, including the AUC of the precision-recall
        curve, calculated between the baseline and the period of activation.
    """

    # Initialize variables
    positive_bin_markers = None
    n_features = data.shape[1]
    n_bins, n_trials = data_shape

    # Expand array of bins to match the number of observations in the data matrix
    bins_to_cumulate_array = np.tile(bins_to_cumulate, (n_trials, ))

    # Initialize classifier
    clf = LogisticRegression(penalty='l2', solver='lbfgs', class_weight='balanced')

    # Get the number of times observations have to be cumulated
    bin_markers = np.unique(bins_to_cumulate)
    bin_markers = bin_markers[bin_markers > 0]
    n_iterations = bin_markers.shape[0]

    # Allocate output variable where to store response probability and strength, and whether it is selective
    performance = np.tile([0, 0, .5], (n_iterations, 1))
    for i_iter in range(n_iterations):
        # Get positive data
        if subset_type == 'cumulative':  # Keep adding bins
            positive_bin_markers = bin_markers[:i_iter + 1]

        elif subset_type == 'single':  # Consider each bin separately
            positive_bin_markers = bin_markers[i_iter]

        # Compute posterior probabilities
        posterior_probabilities = _fit_predict(clf, data, (n_trials, n_features), bins_to_cumulate_array, positive_bin_markers, CV_k_folds, shuffle=False)
        # Assess statistical significance with a bootstrap test
        null_distribution = Parallel(n_jobs=N_JOBS)(
                delayed(_fit_predict)(clf, data, (n_trials, n_features), bins_to_cumulate_array, positive_bin_markers, CV_k_folds, shuffle=True)
                for _ in range(n_significance_permutation))
        null_distribution = np.vstack(null_distribution)
        # Check the probability of each prediction against the null distribution
        difference = np.array(null_distribution >= posterior_probabilities, dtype=int)
        # Count the number of times the relationship is true
        count_per_resample = np.sum(difference, axis=0)
        # Estimate one-sided p-values for each observation
        p_value = (count_per_resample + 1) / (n_significance_permutation + 1)
        # Keep only significant observations
        significant_trials = np.where(p_value <= 0.05)[0]
        if significant_trials.shape[0] > 0:
            is_selective = 1
            probability_response = significant_trials.shape[0] / n_trials
            # Compute AUC
            X_positive = _get_bins(data, bins_to_cumulate_array, positive_bin_markers,     n_trials, n_features).reshape(-1, n_trials).T
            X_negative = _get_bins(data, bins_to_cumulate_array, np.array([0], dtype=int), n_trials, n_features).reshape(-1, n_trials).T
            # AUC = compute_AUC(X_positive[significant_trials, :], X_negative[significant_trials, :], metric='PR')
            AUC = compute_AUC(X_positive, X_negative, metric='PR')
            # Store info
            performance[i_iter, :] = [is_selective, probability_response, AUC]
        else:
            continue

    # Convert output to DataFrame
    performance_single_cell = pd.DataFrame(performance, columns=['is_selective', 'response_probability', 'AUC'])
    # Store response probability
    return performance_single_cell


### Helper functions ###
def _fit_predict(clf, data: np.ndarray, data_shape: Tuple, bins_to_cumulate_array: np.ndarray, positive_bin_markers: np.ndarray, CV_k_folds: int, shuffle: bool = False):
    """This function is called to compute the predicted posterior probabilities
    while classifying observations in data.

    :param clf: sklearn classifier with a predict_proba() method.
    :param data: The data.
    :param data_shape: Size of data.
    :param bins_to_cumulate_array: Array that marks each observation in data with
        the corresponding marker in bins_to_cumulate (not passed to this function).
    :param positive_bin_markers: List of bins corresponding to the positive class.
    :param CV_k_folds: Number of folds for the stratified k-fold cross-validation.
    :param shuffle: Whether to shuffle the bin markers before computing the score.
        This makes for a permutation test in which the labels are randomized.
    """
    # Unpack input
    n_trials, n_features = data_shape

    # Shuffle bins to group
    if shuffle:
        data = np.random.choice(data.ravel(), size=data.shape, replace=False)
    # Get data
    X_positive = _get_bins(data, bins_to_cumulate_array, positive_bin_markers, n_trials, n_features)
    X_negative = _get_bins(data, bins_to_cumulate_array, np.array([0], dtype=int), n_trials, n_features)
    X = np.vstack((X_negative, X_positive))
    y = np.hstack((np.zeros(X_negative.shape[0], dtype=int), np.ones(X_positive.shape[0], dtype=int)))
    # Compute posterior probabilities
    CV_partition = StratifiedKFold(n_splits=CV_k_folds, shuffle=True)
    y_pred = cross_val_predict(clf, X, y, cv=CV_partition, method='predict_proba')[y == 1, 1]

    return y_pred


# Own implementation of ROC analysis
def compute_AUC(dist_pos: np.ndarray, dist_neg: np.ndarray, metric: str = 'PR') -> float:
    """Compute the ara under the curve (AUC) of the receiving-operator
    characteristics (ROC) curve or of the precision-call curve, while comparing
    the data in two distributions. Distributions can have unequal size, as this
    is taken into account when computing the performance of the separation
    between the two.

    :param dist_pos: Data of the positive class.
    :param dist_neg: Data of the negative class.
    :param metric: Either 'ROC' or 'PR' (for precision-recall)

    :return: The value of AUC minus the performance that a random classifier would
        achieve. This means that AUC can range [-1, 1], with negative values
        corresponding to the positive class having a mean lower than the negative
        class, and vice versa when AUC is positive.
    """
    # Initialize variables
    AUC = None

    # Get number of elements in each distribution
    n_pos_observations = dist_pos.size
    n_neg_observations = dist_neg.size
    # Determine accuracy of a random classifier
    random_classifier_AUC = divide0(n_pos_observations, n_pos_observations + n_neg_observations, replace_with=0)

    # Combine all observations
    data = np.hstack((dist_pos.ravel(), dist_neg.ravel()))
    # Calculate the threshold values between data points
    s_data = np.sort(np.unique(data))
    d_data = np.diff(s_data)
    # If there are insufficient data points, return the AUC of a random classifier
    if d_data.size == 0:
        return random_classifier_AUC

    # Compute range of thresholds
    d_data = np.hstack((d_data, d_data[-1]))
    thresholds = np.hstack((s_data[0] - d_data[0], s_data + d_data / 2))

    # Calculate hits and misses
    TP = np.sum((dist_pos.ravel()[:, None] >= thresholds).astype(int), axis=0)
    FP = np.sum((dist_neg.ravel()[:, None] >= thresholds).astype(int), axis=0)
    # Compute the rest of the confusion matrix
    FN = n_pos_observations - TP  # False negatives

    if metric == 'ROC':
        TN = n_neg_observations - FP  # True negatives
        # Compute the area under the ROC curve in the ROC space
        # https://en.wikipedia.org/wiki/Receiver_operating_characteristic
        TPR = divide0(TP, TP + FN, 0)  # true positive rate
        FPR = divide0(FP, FP + TN, 0)  # false positive rate
        AUC = np.abs(np.trapz(x=FPR, y=TPR))

    elif metric == 'PR':
        precision = divide0(TP, TP + FP, 0)
        recall = divide0(TP, TP + FN, 0)
        AUC = np.abs(np.trapz(x=recall, y=precision))

    # Normalize AUC to [0, 1] interval
    r = 0.5  # set performance of random classifier
    AUC = r + (AUC - random_classifier_AUC) * (1 - r) / (np.abs(random_classifier_AUC - r) + r)

    return AUC


def _get_bins(X: np.ndarray, bins_to_cumulate_array, bin_markers_to_select: np.ndarray, n_trials: int, n_features: int) -> np.ndarray:
    """Simple sub-function that selects and returns observations from a larger
    data array.

    :param X: The main data matrix.
    :param bins_to_cumulate_array [list or numpy array]: Array with markers for
        each observation to which bin they belong.
    :param bin_markers_to_select: List of bins to select.
    :param n_trials: Number of trials in the original data.
    :param n_features: Number of features in the original data. Data with more
        than one trial and feature are not supported yet.

    :return X_out: A 1d numpy array with the selected observations.
    """

    if n_trials > 1 and n_features > 1:
        raise Exception('Should not return raveled output. Return mean?')

    ### OLD ###
    # mean_X = np.zeros((n_trials, n_features))
    # # Average across bins of the same trial
    # for i_col in range(X_binned.shape[1]):
    #     this_col_X = X_binned[:, i_col].reshape(n_trials, -1)
    #     # Compute mean and standard deviation
    #     mean_X[:, i_col] = np.mean(this_col_X, axis=1)
    #
    # return mean_X

    ### NEW ###
    # Get the observations from a class
    X_out = X[np.where(np.in1d(bins_to_cumulate_array, bin_markers_to_select))[0], :].ravel().reshape(-1, 1)

    return X_out


def divide0(a: Union[float, np.ndarray], b: Union[float, np.ndarray], replace_with: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
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
### Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        # See more rows and columns of variables printed in console
        np.set_printoptions(suppress=True, linewidth=500)
        pd.set_option('display.max_columns', 500)
        pd.set_option('display.width', 1000)

        run(input_filename=r"D:\_MATLAB_2PI\FK_5\session160919_HPS_detection_data.mat", verbose=True)

