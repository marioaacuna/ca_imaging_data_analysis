# -*- coding: utf-8 -*-
# This module contains functions to decode activation periods in multi-trial
# imaging sessions. It takes inputs from MATLAB and returns a JSON file that can
# be read back in MATLAB.
# Threshold p value = 0.01

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
from Utilities.Arrays import idx2range


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

    # Get timestamps
    timestamps = PARAMETERS['timestamps']
    # Make time axis for this stimulus, staring at the 1st timestamp (stimulus onset)
    n_bins = PARAMETERS['data_size'][0]
    timestamps_s = PARAMETERS['timestamps_s']
    times = np.linspace(PARAMETERS['trial_time'][0], PARAMETERS['trial_time'][1], n_bins + 1)
    times -= timestamps_s[0]
    # Set maximum latency from last timestamp for an activation period to be
    # regarded as relevant. It is relative to the stimulus onset
    max_latency = PARAMETERS['evoked_activity_max_latency_bins'] + (timestamps[-1] - timestamps[0])

    # Get the observation markers
    bins_to_cumulate = PARAMETERS['bins_to_cumulate']
    bins_to_cumulate[np.isnan(bins_to_cumulate)] = -1
    bins_to_cumulate = bins_to_cumulate.astype(int)
    n_positive_bins = np.where(bins_to_cumulate > 0)[0].shape[0]
    # Get the data and its original [trial x time] shape
    data = PARAMETERS['data']
    data_shape = PARAMETERS['data_size']

    # Compute selectivity of individual cells
    columns_to_copy = ['start', 'finish', 'response_probability', 'peak_probability', 'activity_difference', 'AUC', 'p', 'effect_size']
    performance_all_cells = pd.DataFrame(columns=['cell'] + columns_to_copy)
    response_probability_traces = np.zeros((n_cells, n_positive_bins), dtype=float)
    significant_trials = list(np.empty((n_cells, ), dtype=object))
    for i_cell in range(n_cells):
        X = data[:, i_cell].reshape(-1, 1).copy()
        performance_single_cell, significant_activations = _compute_accuracy_observation_subset(X, bins_to_cumulate, data_shape, max_latency,
                                                                                                subset_type='single', CV_k_folds=CV_k_folds,
                                                                                                n_significance_permutation=PARAMETERS['significance_test_n_permutations'])

        # Store values
        response_probability_traces[i_cell, :] = performance_single_cell['response_probability']
        significant_trials[i_cell] = performance_single_cell['significant_trials'].values.tolist()
        if significant_activations is None:
            n_significant_activations = 0
        else:
            n_significant_activations = significant_activations.shape[0]
        if n_significant_activations > 0:
            for i in range(n_significant_activations):
                row = performance_all_cells.shape[0]
                performance_all_cells.at[row, 'cell'] = i_cell + 1  # 1-indexed
                performance_all_cells.loc[row, columns_to_copy] = significant_activations.loc[i, columns_to_copy]

            if verbose:
                log('Cell %i / %i: %i activations: Peak at %ss, AUC=%s  |  P_value %s' %
                    (i_cell + 1, n_cells, n_significant_activations, np.round(times[significant_activations['peak_probability'].values], 2), np.round(significant_activations['AUC'].values, 2), significant_activations['p'].values))
        else:
            if verbose:
                log('Cell %i / %i: Non-selective' % (i_cell + 1, n_cells))

    # Store data in dictionary
    RESULTS = dict()
    performance = list([list(performance_all_cells.columns)])
    performance.append(performance_all_cells.values.tolist())
    RESULTS['performance'] = performance
    RESULTS['response_probability_traces'] = response_probability_traces.tolist()
    RESULTS['significant_trials'] = significant_trials

    # Write results to a text file in JSON format
    with open(PARAMETERS['output_filename'], 'w') as f:
        json.dump(RESULTS, f)

    if verbose:
        log('Finished')


################################################################################
### Statistics
################################################################################
def _compute_accuracy_observation_subset(data: np.ndarray, bins_to_cumulate: np.ndarray,
                                         data_shape: np.ndarray, max_latency: int,
                                         subset_type: str = 'single', CV_k_folds: int = 5,
                                         n_significance_permutation:int = 1000):
    """This function performs the main job of iterating through time bins and
    computing whether they differ from a baseline period.

    :param data: numpy-array obtained from [n_samples x n_features].ravel()
    :param bins_to_cumulate: It marks each time bin with either -1 (discard), 0
        (baseline) or an int > 0. Bins with the same marker will be grouped.
    :param data_shape: List of 2 items representing the number of bins and trials
        in the original data.
    :param max_latency: Maximum number of bins after the last timestamp after
        which a new activation period will be discarded. Note that this does not
        exclude activations that started before this delay and might peak after.
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
    significant_activations = None
    n_significant_activations = 0
    stimulus_onset = np.where(bins_to_cumulate > 0)[0][0]
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
    performance = np.zeros((n_iterations, 3), dtype=float)
    P_values = np.zeros((n_trials, bin_markers.shape[0]), dtype=np.float32)
    data_diff = np.zeros((n_trials, bin_markers.shape[0]), dtype=np.float32)
    for i_iter in range(n_iterations):
        # Get positive data
        if subset_type == 'cumulative':  # Keep adding bins
            positive_bin_markers = bin_markers[:i_iter + 1]

        elif subset_type == 'single':  # Consider each bin separately
            positive_bin_markers = bin_markers[i_iter]

        # Compute posterior probabilities
        posterior_probabilities = _estimate_score(clf, data, (n_trials, n_features), bins_to_cumulate_array, positive_bin_markers, CV_k_folds, shuffle=False)
        # Assess statistical significance with a bootstrap test
        null_distribution = Parallel(n_jobs=-1, prefer='processes')(
                delayed(_estimate_score)(clf, data, (n_trials, n_features), bins_to_cumulate_array, positive_bin_markers, CV_k_folds, shuffle=True)
                for _ in range(n_significance_permutation))
        null_distribution = np.vstack(null_distribution)
        # Compute a one-sided p-value
        p_values = (1 + np.sum(null_distribution >= posterior_probabilities[None, :], axis=0)) / (n_significance_permutation + 1)
        P_values[:, i_iter] = p_values

        if not isinstance(posterior_probabilities, np.ndarray):
            p_values = float(p_values)
            is_selective = p_values <= 0.01
            probability_response = int(is_selective) * 100
            response_strength = posterior_probabilities

        else:
            # Take average difference between two epochs
            X_positive = _get_bins(data, bins_to_cumulate_array, positive_bin_markers, n_trials, n_features)
            X_negative = _get_bins(data, bins_to_cumulate_array, np.array([0], dtype=int), n_trials, n_features).reshape(n_trials, -1).mean(1).reshape(-1, 1)
            X_diff = (X_positive - X_negative).ravel()
            data_diff[:, i_iter] = X_diff
            # Take trials where there is a significant activation
            significant_trials = np.where(np.logical_and(p_values <= 0.01, X_diff != 0))[0]
            if significant_trials.size > 0:
                probability_response = significant_trials.shape[0] / p_values.shape[0]  * 100
                response_strength = np.mean(X_diff[significant_trials])
                is_selective = True
            else:
                probability_response = 0
                response_strength = 0
                is_selective = False

        # Store values
        performance[i_iter, :] = [is_selective, probability_response, response_strength]

    # Convert output to DataFrame
    performance_single_cell = pd.DataFrame(performance, columns=['is_selective', 'response_probability', 'response_strength'])
    # Store response probability
    this_response_probability_trace = performance_single_cell['response_probability'] * np.sign(performance_single_cell['response_strength'].values)
    # Get range is in which the response probability is statistically significant
    idx = np.empty((0, 4), dtype=int)
    idx_exc = idx2range(np.where(this_response_probability_trace > 0)[0])
    idx_inh = idx2range(np.where(this_response_probability_trace < 0)[0])
    if idx_exc.size > 0:
        idx = np.vstack((idx, np.hstack((idx_exc, np.ones((idx_exc.shape[0], 1))))))
    if idx_inh.size > 0:
        idx = np.vstack((idx, np.hstack((idx_inh, np.zeros((idx_inh.shape[0], 1)) - 1))))
    # Merge positive and negative modulations and sort them by latency
    if idx.shape[0] > 0:
        idx = idx.astype(int)
        idx = idx[idx[:, 0].argsort(), :]
        # Remove activation periods that started too late
        idx = idx[idx[:, 0] <= max_latency, :]
    if idx.shape[0] == 0:
        performance_single_cell['significant_trials'] = [list() for _ in range(performance_single_cell.shape[0])]

    else:
        # Make ranges of bins
        n_intervals = idx.shape[0]
        evoked_periods = [bin_markers[idx[i, 0]:idx[i, 1] + 1] for i in range(n_intervals)]
        # Compute AUC in each period compared to baseline period
        baseline_activity = data[np.where(bins_to_cumulate_array == 0)[0]].ravel()
        AUC = np.zeros((n_intervals, ), dtype=float)
        for i in range(n_intervals):
            evoked_activity = data[np.where(np.in1d(bins_to_cumulate_array, evoked_periods[i]))[0]].ravel()
            AUC[i] = compute_AUC(evoked_activity, baseline_activity, metric='PR')
        # Perform a permutation test
        null_distribution = list()
        for i_rep in range(n_significance_permutation):
            null_AUC = np.zeros((n_intervals,), dtype=float)
            for i in range(n_intervals):
                bins_to_cumulate_array_shuffled = np.random.choice(bins_to_cumulate_array, bins_to_cumulate_array.shape, replace=False)
                baseline_activity = data[np.where(bins_to_cumulate_array_shuffled == 0)[0]].ravel()
                evoked_activity = data[np.where(np.in1d(bins_to_cumulate_array_shuffled, evoked_periods[i]))[0]].ravel()
                null_AUC[i] = compute_AUC(evoked_activity, baseline_activity, metric='PR')
            null_distribution.append(null_AUC)
        null_distribution = np.vstack(null_distribution)
        # Compute two-tailed p-values
        p_left = (1 + np.sum(null_distribution <= AUC[None, :], axis=0)) / (n_significance_permutation + 1)
        p_right = (1 + np.sum(null_distribution >= AUC[None, :], axis=0)) / (n_significance_permutation + 1)
        p_values = np.clip(2.0 * np.min(np.vstack((p_left, p_right)), axis=0), a_min=0, a_max=1)
        # Find significant trials
        significant_trials_idx = [np.where(np.logical_and(P_values[:, i] <= 0.01, data_diff[:, i] != 0))[0].tolist() for i in range(bin_markers.shape[0])]
        performance_single_cell['significant_trials'] = significant_trials_idx
        # Find significant activations
        significant_activations_idx = np.where(np.logical_and(p_values <= 0.01, np.sign(AUC - np.mean(null_distribution, axis=0)) == idx[:, 3]))[0]
        if significant_activations_idx.shape[0] > 0:
            significant_activations = pd.DataFrame(np.atleast_2d(idx[significant_activations_idx, :2]), columns=['start', 'finish'])
            n_significant_activations = significant_activations.shape[0]
            # Get max response probability in each period
            significant_activations['response_probability'] = [performance_single_cell.loc[significant_activations.loc[i, 'start']:significant_activations.loc[i, 'finish'], 'response_probability'].max() for i in range(n_significant_activations)]
            significant_activations['peak_probability'] = [performance_single_cell.loc[significant_activations.loc[i, 'start']:significant_activations.loc[i, 'finish'], 'response_probability'].idxmax() for i in range(n_significant_activations)]
            # Add offset of stimulus onset
            significant_activations[['start', 'finish', 'peak_probability']] += stimulus_onset
            # Store AUC and its p-value
            significant_activations['AUC'] = AUC[significant_activations_idx]
            significant_activations['p'] = p_values[significant_activations_idx]

            # Get activity difference between evoked and baseline
            X = data.reshape(data_shape[::-1])
            baseline_activity = np.mean(X[:, np.where(bins_to_cumulate == 0)[0]], axis=1)
            for i in range(n_significant_activations):
                evoked_activity = np.mean(X[:, significant_activations.loc[i, 'start']:significant_activations.loc[i, 'finish'] + 1], axis=1)
                significant_activations.loc[i, 'activity_difference'] = np.mean(evoked_activity - baseline_activity)
                significant_activations.loc[i, 'effect_size'] = np.mean(divide0(evoked_activity - baseline_activity, evoked_activity + baseline_activity, replace_with=0))

            # Check whether decrease in response might correspond to a true suppression
            if np.any(significant_activations['AUC'] < 0):
                # Get other intervals that were not considered for analysis
                all_bin_idx = np.ones_like(bins_to_cumulate, dtype=bool)
                for i in range(n_significant_activations):
                    all_bin_idx[significant_activations.loc[i, 'start']:significant_activations.loc[i, 'finish'] + 1] = False
                all_bin_idx[bins_to_cumulate < 0] = False
                all_bin_idx = idx2range(np.where(all_bin_idx)[0])
                # Initialize flag
                significant_activations['confirmed_suppression'] = True

                # Loop through periods of putative suppression
                for i in range(n_significant_activations):
                    if significant_activations.loc[i, 'AUC'] < 0:
                        # Take the latest interval before the one under consideration
                        interval_pre = np.where(all_bin_idx[:, 0] < significant_activations.loc[i, 'start'])[0][-1]
                        bins_pre = np.unique(bins_to_cumulate[all_bin_idx[interval_pre, 0]:all_bin_idx[interval_pre, 1]])
                        activity_pre = data[np.where(np.in1d(bins_to_cumulate_array, bins_pre))[0]].reshape(n_trials, -1)
                        # Take the first interval following the one under consideration
                        interval_post = np.where(all_bin_idx[:, 1] > significant_activations.loc[i, 'finish'])[0]
                        if interval_post.shape[0] > 0:
                            interval_post = interval_post[0]
                            bins_post = np.unique(bins_to_cumulate[all_bin_idx[interval_post, 0]:all_bin_idx[interval_post, 1]])
                            activity_post = data[np.where(np.in1d(bins_to_cumulate_array, bins_post))[0]].reshape(n_trials, -1)
                            remove_edge_trials = False
                        else:
                            bins_post = [-1]
                            # Consider the period before the baseline
                            activity_post = data[np.where(np.in1d(bins_to_cumulate_array, bins_post))[0]].reshape(n_trials, -1).reshape(n_trials, -1)
                            # Remove first and last trial because there the pre/post relationship does not apply to these two intervals
                            activity_post = activity_post[1:-1, :]
                            activity_pre = activity_pre[1:-1, :]
                            remove_edge_trials = True
                        # Compute AUC
                        AUC = compute_AUC(activity_post.ravel(), activity_pre.ravel(), metric='PR')
                        # Perform a permutation test
                        null_distribution = list()
                        for i_rep in range(n_significance_permutation):
                            bins_to_cumulate_array_shuffled = np.random.choice(bins_to_cumulate_array, bins_to_cumulate_array.shape, replace=False)
                            activity_pre = data[np.where(np.in1d(bins_to_cumulate_array_shuffled, bins_pre))[0]]
                            activity_post = data[np.where(np.in1d(bins_to_cumulate_array_shuffled, bins_post))[0]]
                            if remove_edge_trials:
                                activity_pre = activity_pre.reshape(n_trials, -1)[1:-1, :]
                                activity_post = activity_post.reshape(n_trials, -1)[1:-1, :]
                            null_distribution.append(compute_AUC(activity_post.ravel(), activity_pre.ravel(), metric='PR'))
                        null_distribution = np.hstack(null_distribution)
                        # Compute one-tailed p-value
                        p_value = (1 + np.sum(null_distribution <= AUC)) / (n_significance_permutation + 1)
                        if p_value <= 0.01:
                            significant_activations.loc[i, 'confirmed_suppression'] = False

                # Remove activations where suppression could not be confirmed
                significant_activations.drop(np.where(np.logical_not(significant_activations['confirmed_suppression']))[0], inplace=True)
                significant_activations.reset_index(drop=True, inplace=True)
                significant_activations.drop(columns='confirmed_suppression', axis=1, inplace=True)
                n_significant_activations = significant_activations.shape[0]
                if n_significant_activations == 0:
                    significant_activations = None

        # Get indices of significant bins
        significant_bins_idx = np.zeros((performance_single_cell.shape[0], ), dtype=np.bool)
        for i in range(n_significant_activations):
            significant_bins_idx[significant_activations.loc[i, 'start'] - stimulus_onset:significant_activations.loc[i, 'finish'] + 1 - stimulus_onset] = True
        performance_single_cell['is_selective'] = significant_bins_idx
        # Remove identifications of significant trials that were not confirmed
        performance_single_cell['significant_trials'] = [list() if not performance_single_cell.loc[i, 'is_selective'] else performance_single_cell.loc[i, 'significant_trials'] for i in range(performance_single_cell.shape[0])]

    return performance_single_cell, significant_activations


### Helper functions ###
def _estimate_score(clf, data: np.ndarray, data_shape: Tuple, bins_to_cumulate_array: np.ndarray, positive_bin_markers: np.ndarray, CV_k_folds: int, shuffle: bool = False):
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
    # Compute performance score
    CV_partition = StratifiedKFold(n_splits=CV_k_folds, shuffle=True)
    y_pred = cross_val_predict(clf, X, y, cv=CV_partition, method='predict_proba')[:, 1]
    return y_pred[y == 1]


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

        run(input_filename=r"D:\_MATLAB_CaImaging\FK_4\session160729_HPS_detection_data", verbose=True)

