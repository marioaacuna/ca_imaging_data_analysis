# Add Utilities folder to system path
import os, sys
script_path = os.path.realpath(__file__)
root_path = os.path.split(os.path.dirname(script_path))[0]
sys.path.insert(0, root_path)

# System packages
import json
from collections import OrderedDict

# Numerical packages
import numpy as np
import pandas as pd
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.utils import resample
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.calibration import clone
from statsmodels.stats.multitest import multipletests
# Parallel computing
from joblib import Parallel, delayed

# Printing figures
from matplotlib.backends.backend_pdf import PdfPages

# Local repository
from decoding import visualization
from Utilities.IO_operations import log
from Utilities import matlab_file


overwrite_figures = False
reload_data_from_disk = False


def run(input_filename, debug=False):
    # Default settings
    CV_k_folds = 5

    # Set random number generator seed for replicability
    filename = os.path.splitext(os.path.basename(input_filename))[0]
    seed = int(np.clip(np.int64(np.sum([ord(char) for char in filename])), a_min=0, a_max=2**32 - 1))
    np.random.seed(seed)
    if debug:
        log('Random number generator seed set to %i' % seed)

    # Initialize output variable
    RESULTS = OrderedDict()
    RESULTS['bins_significant_activation'] = dict()

    # Load parameter file
    if debug:
        log('Analyzing \'%s\'' % input_filename)
    PARAMETERS = matlab_file.load(input_filename)
    # Unpack some parameters
    class_names = list(PARAMETERS['data'].keys())
    n_classes = len(class_names)
    n_cells = np.min([PARAMETERS['data'][c].shape[1] for c in class_names])

    # Make sure some parameters are lists
    if isinstance(PARAMETERS['output_figures']['activation_detection'], str):
        PARAMETERS['output_figures']['activation_detection'] = [PARAMETERS['output_figures']['activation_detection']]
    for cl in class_names:
        if not isinstance(PARAMETERS['timestamps'][cl], (list, np.ndarray)):
            PARAMETERS['timestamps'][cl] = np.array([PARAMETERS['timestamps'][cl]])

    # Use a linear classifier to separate activity bins before and after each timestamp
    if debug:
        log('Computing probability of response')
    # Make sure that output folder exists
    if not os.path.exists(PARAMETERS['output_figures']['folder']):
        os.mkdir(PARAMETERS['output_figures']['folder'])

    print()
    # Take the training stimulus
    training_stimulus = 'HPS'
    X_stim_raw = PARAMETERS['data'][training_stimulus]
    n_times, n_trials = PARAMETERS['data_size'][training_stimulus]
    X_stim = np.zeros((n_trials, n_cells, n_times))
    for i_trial in range(n_trials):
        X_stim[i_trial, :, :] = X_stim_raw[i_trial * n_times:i_trial * n_times + n_times, :].T
    # Take the spontaneous activity
    X_sp_raw = PARAMETERS['data']['SP']
    X_sp = np.zeros((n_trials, n_cells, n_times))
    for i_trial in range(n_trials):
        X_sp[i_trial, :, :] = X_sp_raw[i_trial * n_times:i_trial * n_times + n_times, :].T


    from mne.decoding import SlidingEstimator, LinearModel, cross_val_multiscore
    from sklearn.metrics import matthews_corrcoef, make_scorer
    scoring = 'roc_auc' # make_scorer(matthews_corrcoef)

    clf = LinearModel(LogisticRegression(penalty='l2', solver='lbfgs', class_weight='balanced'))
    time_decod = SlidingEstimator(clf, scoring=scoring, verbose=0)

    # Run cross-validated decoding analyses:
    X_all = np.vstack((X_sp, X_stim))
    y = np.hstack((np.zeros(X_sp.shape[0], dtype=int), np.ones(X_stim.shape[0], dtype=int)))
    neuron_n = 40
    X = X_all[:, neuron_n, :].reshape(n_trials * 2, 1, n_times).copy()
    # From empirical data
    CV_partition = StratifiedKFold(n_splits=10, shuffle=True)
    CV_scores = cross_val_multiscore(time_decod, X, y, cv=CV_partition, n_jobs=1, verbose=0)
    scores = np.mean(CV_scores, axis=0)
    # Permute labels
    n_significance_permutation = 100
    null_distribution = list()
    for i_rep in range(n_significance_permutation):
        print(i_rep + 1)
        X_shuffled = np.random.choice(X.ravel(), X.shape, replace=False)
        CV_partition = StratifiedKFold(n_splits=10, shuffle=True)
        null_CV_scores = cross_val_multiscore(time_decod, X_shuffled, y, cv=CV_partition, verbose=0)
        null_scores = np.mean(null_CV_scores, axis=0)
        null_distribution.append(null_scores)
    null_distribution = np.vstack(null_distribution)
    # Compute two-tailed p-value
    p_left = (1 + np.sum(null_distribution <= scores[None, :], axis=0)) / (n_significance_permutation + 1)
    p_right = (1 + np.sum(null_distribution >= scores[None, :], axis=0)) / (n_significance_permutation + 1)
    p_values = np.clip(2.0 * np.min(np.vstack((p_left, p_right)), axis=0), a_min=0, a_max=1)
    print(np.where(p_values <= 0.05)[0])

    plt.clf(); plt.plot(np.arange(n_times), scores, '-ok', label='score')
    plt.axhline(.5, color='k', linestyle='--', label='chance'); plt.axvline(20.5, color='r'); plt.legend()
    plt.show()

    plt.clf(); plt.imshow(X[:, 0, :])

    plt.clf(); plt.plot(X_stim.mean(0).mean(0)); plt.plot(X_sp.mean(0).mean(0))


    return

    for class_idx in range(n_classes):
        # Get timestamps
        timestamps = PARAMETERS['timestamps'][class_names[class_idx]]
        if timestamps.size == 0:
            continue

        # Create pdf file
        output_filename = os.path.join(PARAMETERS['output_figures']['folder'], PARAMETERS['output_figures']['base_filename'] + 'response_probability_single_cells_%s.pdf' % (class_names[class_idx]))
        try:
            if not os.path.exists(output_filename) or overwrite_figures:
                pdf_is_open = True
                PDF_file = PdfPages(output_filename)
            else:
                pdf_is_open = False
        except PermissionError as err:
            if err.strerror == 'Permission denied':
                raise PermissionError(
                    'Cannot continue while the PDF \'%s\'\nis open. Please close it and retry.' % output_filename)
            else:
                raise err

        # Get the observation markers
        bins_to_cumulate = PARAMETERS['bins_to_cumulate'][class_names[class_idx]]
        bins_to_cumulate[np.isnan(bins_to_cumulate)] = -1
        bins_to_cumulate = bins_to_cumulate.astype(int)
        n_positive_bins = np.where(bins_to_cumulate > 0)[0].shape[0]
        # Get the data and its original [trial x time] shape
        data = PARAMETERS['data'][class_names[class_idx]].copy()
        data_shape = PARAMETERS['data_size'][class_names[class_idx]]

        # Compute selectivity of individual cells
        performance = np.zeros((n_cells, 3), dtype=np.float64) * np.nan
        response_probability_traces = np.zeros((n_cells, n_positive_bins), dtype=float)
        for i_cell in range(n_cells):
            # Compute performance in cumulative bins
            X = data[:, i_cell].reshape(-1, 1).copy()
            performance_single_cell = _compute_accuracy_observation_subset(X, bins_to_cumulate, data_shape,
                                                                           subset_type='single', CV_k_folds=CV_k_folds,
                                                                           n_significance_permutation=500,
                                                                           debug=True)

            # Make figure
            if 'single_cells' in PARAMETERS['output_figures']['activation_detection']:
                # Reshape data
                X = X.reshape(data_shape[::-1])
                # Make figure
                fig = visualization.heatmap_response_probability_per_trial(X, performance_single_cell, bins_to_cumulate, timestamps,
                                                                     PARAMETERS['frame_rate'], PARAMETERS['n_frames_per_bin'], title='Cell %i' % (i_cell + 1))
                if pdf_is_open:
                    PDF_file.savefig(fig, bbox_inches='tight')
                plt.close(fig)

            # Store response probability
            response_probability_traces[i_cell, :] = performance_single_cell['response_probability'] * np.sign(performance_single_cell['response_strength'])
            # Find highest response probability
            peak_reliability = performance_single_cell['response_probability'].max()
            if peak_reliability > 0:
                peak_latency = np.where(performance_single_cell['response_probability'] == peak_reliability)[0]
                peak_strength = performance_single_cell.loc[peak_latency, 'response_strength'].values
                if peak_latency.shape[0] > 1:
                    highest_peak = peak_strength.argmax()
                    peak_latency = peak_latency[highest_peak] + 1
                    peak_strength = peak_strength[highest_peak]
                # Store values
                performance[i_cell, :] = [peak_reliability, peak_latency, peak_strength]
                if debug:
                    log('Cell %i / %i: Reliability %i%%  |  Latency %i bins  |  Difference %+.3f' % (i_cell + 1, n_cells, peak_reliability, peak_latency, peak_strength))

            else:
                if debug:
                    log('Cell %i / %i: Non-selective' % (i_cell + 1, n_cells))

        # Convert the array performance to a pandas DataFrame
        performance = pd.DataFrame(performance, columns=['peak_reliability', 'peak_latency', 'activity_difference'])
        performance['cell'] = np.arange(performance.shape[0]) + 1
        performance = performance.sort_values(by='peak_reliability', ascending=False, na_position='last').reset_index(drop=True)

        # Make figure with all cells
        fig = visualization.heatmap_response_probability_per_cell(response_probability_traces, timestamps,
                                                                  PARAMETERS['frame_rate'], PARAMETERS['n_frames_per_bin'], title=PARAMETERS['output_figures']['base_filename'][:-1] + ' %s' % class_names[class_idx])
        if pdf_is_open:
            PDF_file.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        # Flush images to disk
        if pdf_is_open:
            PDF_file.close()
            pdf_is_open = False

        # Temporarily store files on disk
        pickle.dump((performance, response_probability_traces), open(r'D:\_MATLAB_2PI\%s.p' % filename, 'wb'))

        print()

        performance, response_probability_traces = pickle.load(open(r'D:\_MATLAB_2PI\%s.p' % filename, 'rb'))















################################################################################
### Statistics
################################################################################
def _compute_accuracy_observation_subset(data, bins_to_cumulate, data_shape, subset_type='cumulative', CV_k_folds=10, n_significance_permutation=100, debug=False):
    # Initialize variables
    positive_bin_markers = None
    n_features = data.shape[1]
    n_bins, n_trials = data_shape

    # Expand array of bins to match the number of observations in the data matrix
    bins_to_cumulate_array = np.tile(bins_to_cumulate, (n_trials, ))

    from mne.decoding import LinearModel

    # Initialize classifier
    clf = Pipeline([('scaler', StandardScaler(with_mean=True, with_std=True, copy=True)),
                    ('classifier', LinearModel(LogisticRegression(penalty='l2', solver='lbfgs', class_weight='balanced')))])




    # Get the number of times observations have to be cumulated
    bin_markers = np.unique(bins_to_cumulate)
    bin_markers = bin_markers[bin_markers > 0]
    n_iterations = bin_markers.shape[0]

    # Allocate output variable where to store response probability and strength, and whether it is selective
    performance = np.zeros((n_iterations, 3), dtype=float)
    significant_trials_idx = list()
    for i_iter in range(n_iterations):
        # Get positive data
        if subset_type == 'cumulative':  # Keep adding bins
            positive_bin_markers = bin_markers[:i_iter + 1]

        elif subset_type == 'single':  # Consider each bin separately
            positive_bin_markers = bin_markers[i_iter]

        # Compute posterior probabilities
        posterior_probabilities = _estimate_score(clf, data, (n_trials, n_features), bins_to_cumulate_array, bins_to_cumulate, positive_bin_markers, CV_k_folds, shuffle=False)
        # Assess statistical significance with a bootstrap test
        null_distribution = Parallel(n_jobs=-1, prefer='processes')(
                delayed(_estimate_score)(clf, data, (n_trials, n_features), bins_to_cumulate_array, bins_to_cumulate, positive_bin_markers, CV_k_folds, shuffle=True)
                for _ in range(n_significance_permutation))
        if isinstance(posterior_probabilities, np.ndarray):
            null_distribution = np.vstack(null_distribution)
            # Compute two-tailed p-value
            p_left = (1 + np.sum(null_distribution <= posterior_probabilities[None, :], axis=0)) / (n_significance_permutation + 1)
            p_right = (1 + np.sum(null_distribution >= posterior_probabilities[None, :], axis=0)) / (n_significance_permutation + 1)

        else:
            null_distribution = np.hstack(null_distribution)
            p_left = (1 + np.sum(null_distribution <= posterior_probabilities, axis=0)) / (n_significance_permutation + 1)
            p_right = (1 + np.sum(null_distribution >= posterior_probabilities, axis=0)) / (n_significance_permutation + 1)

        p_values = np.clip(2.0 * np.min(np.vstack((p_left, p_right)), axis=0), a_min=0, a_max=1)
        # p-value correction with False Discovery Rate
        # _, p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

        if not isinstance(posterior_probabilities, np.ndarray):
            p_values = float(p_values)
            is_selective = p_values <= 0.05
            probability_response = int(is_selective) * 100
            response_strength = posterior_probabilities

        else:
            # Take average difference between two epochs
            X = data[np.where(bins_to_cumulate_array >= 0)[0], :].reshape(n_trials, -1).transpose()
            y = np.hstack((np.zeros(np.where(bins_to_cumulate == 0)[0].shape[0], dtype=int), np.ones(np.where(bins_to_cumulate > 0)[0].shape[0], dtype=int)))
            X_positive = np.mean(X[y == 1, :], axis=0)
            X_negative = np.mean(X[y == 0, :], axis=0)
            X_diff = X_positive - X_negative
            # Take trials where there is a significant activation
            # significant_trials = np.where(np.logical_and(p_values <= 0.05, X_diff != 0))[0]
            significant_trials = np.where(p_values <= 0.05)[0]
            significant_trials_idx.append(significant_trials)
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
        # Log outcome
        if debug:
            message = 'Iteration %i / %i: Selective? %s' % (i_iter + 1, n_iterations, is_selective)
            if is_selective:
                message += '  |  Response probability %i%%  |  Difference %+.3f' % (probability_response, response_strength)
            log(message)

    # Convert output to DataFrame
    performance = pd.DataFrame(performance, columns=['is_selective', 'response_probability', 'response_strength'])
    performance['significant_trials'] = significant_trials_idx

    return performance


def _estimate_score(clf, data, data_shape, bins_to_cumulate_array, bins_to_cumulate, positive_bin_markers, CV_k_folds, shuffle=False):
    # Unpack input
    n_trials, n_features = data_shape

    # Prepare data for classification
    X = data[np.where(bins_to_cumulate_array >= 0)[0], :].reshape(n_trials, -1).transpose()
    y = np.hstack((np.zeros(np.where(bins_to_cumulate == 0)[0].shape[0], dtype=int),
                   np.ones(np.where(bins_to_cumulate > 0)[0].shape[0], dtype=int)))
    # Shuffle observations
    if shuffle:
        X = np.random.choice(X.ravel(), size=X.shape, replace=True); plt.clf(); plt.imshow(X.T)

    # Compute performance score
    CV_partition = StratifiedKFold(n_splits=CV_k_folds, shuffle=True)
    # y_pred = cross_val_predict(classifier, X, y, cv=CV_partition, method='predict')
    y_pred = cross_val_predict(clf, X, y, cv=CV_partition, method='predict_proba')[:, 1]

    return y_pred[y == 1]

    # AUC = compute_AUC(X_positive, X_negative, metric='PR')
    # return AUC


# Own implementation of ROC analysis
def compute_AUROC(dist_pos, dist_neg, metric='PR'):
    # Get number of elements in each distribution
    n_pos_observations = dist_pos.size
    n_neg_observations = dist_neg.size

    # Combine all observations
    data = np.hstack((dist_pos.ravel(), dist_neg.ravel()))
    # Calculate the threshold values between data points
    s_data = np.sort(np.unique(data))
    d_data = np.diff(s_data)
    if d_data.size == 0:
        # Determine accuracy of random classifier
        random_classifier_AUC = n_pos_observations / (n_pos_observations + n_neg_observations)
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

    return AUC


def _get_bins(X, bins_to_cumulate_array, positive_bin_markers, n_trials, n_features):
    # Get the observations from a class
    X_binned = X[np.where(np.in1d(bins_to_cumulate_array, positive_bin_markers))[0], :]
    mean_X = np.zeros((n_trials, n_features))
    # Average across bins of the same trial
    for i_col in range(X_binned.shape[1]):
        this_col_X = X_binned[:, i_col].reshape(n_trials, -1)
        # Compute mean and standard deviation
        mean_X[:, i_col] = np.mean(this_col_X, axis=1)

    return mean_X


def divide0(a, b, rep):
    """Divide two numbers but replace its result if division is not
    possible, e.g., when dividing a number by 0.

    :param a: [numeric] A number.
    :param b: [numeric] A number.
    :param rep: [numeric] If a/b is not defined return this number instead.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = rep
    return c


################################################################################
### Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        import pickle
        import matplotlib
        matplotlib.use('Qt5Agg')
        from matplotlib import pyplot as plt

        np.set_printoptions(suppress=True)
        pd.set_option('display.max_columns', 500)
        pd.set_option('display.width', 1000)

        run(input_filename=r'D:\_MATLAB_2PI\FK_11_cond1_decoding_data.mat', debug=True)
