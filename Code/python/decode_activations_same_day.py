# System and IO packages
import os
import sys
import json
import pickle

# Numerical packages
import numpy as np

# Local repository
from decoding.EnsembleTreeClassifier import EnsembleTreeClassifier as ETC
from decoding.LogisticRegressionClassifier import LogisticRegressionClassifier as LRC
from Utilities.IO_operations import log
from Utilities import matlab_file


def run(input_filename, verbose=False):
    # Set random number generator seed for replicability
    filename = os.path.splitext(os.path.basename(input_filename))[0]
    seed = int(np.clip(np.int64(np.sum([ord(char) for char in filename])), a_min=0, a_max=2**32 - 1))
    np.random.seed(seed)
    if verbose:
        log('Random number generator seed set to %i' % seed)

    # Load parameter file
    if verbose:
        log('Analyzing \'%s\'' % input_filename)
    PARAMETERS = matlab_file.load(input_filename)
    # Unpack parameters
    class_names = list(PARAMETERS['data'].keys())
    n_classes = len(class_names)

    # Make sure that data is 2D
    for cl in class_names:
        if PARAMETERS['data'][cl].ndim == 1:
            PARAMETERS['data'][cl] = np.atleast_2d(PARAMETERS['data'][cl]).transpose()

    # Prepare data and labels for multi-class classifier
    LABELS = []  # array of labels
    for ci, cl in enumerate(class_names):
        n_samples = PARAMETERS['data'][cl].shape[0]
        these_labels = np.zeros((n_samples,), dtype=int) + ci
        LABELS.append(these_labels)
    # Concatenate the labels array, and make data array
    LABELS_response = np.hstack(LABELS)
    DATA_response = np.vstack(([PARAMETERS['data'][cl] for cl in class_names]))

    if verbose:
        log('Comparing classes: %s' % ', '.join(class_names))
    # Initialize classifier
    if PARAMETERS['decoder_implementation'] == 'random forest':
        mcc = ETC(ensemble_type='RandomForest',
                  tree_split_criterion='entropy',
                  verbose=verbose)
        # Train and test classifier
        mcc.fit_predict(DATA_response, LABELS_response, class_names=class_names,
                        CV_scheme=PARAMETERS['CV_scheme'],
                        v_fraction=PARAMETERS['v_fraction'],
                        shuffle_split_repeats=PARAMETERS['shuffle_split_repeats'],
                        fraction_training=PARAMETERS['fraction_training'],
                        fraction_training_for_tuning=PARAMETERS['fraction_training_for_tuning'],
                        fraction_training_for_calibration=PARAMETERS['fraction_training_for_calibration'],
                        max_n_splits_training=PARAMETERS['max_n_splits_training'],
                        max_n_splits_tuning=PARAMETERS['max_n_splits_tuning'],
                        max_n_splits_calibration=PARAMETERS['max_n_splits_calibration'],
                        calibration_method='sigmoid',
                        data_shape=PARAMETERS['data_shape'],
                        data_scaling=PARAMETERS['data_scaling'],
                        n_permutations=PARAMETERS['significance_test_n_permutations'])

    elif PARAMETERS['decoder_implementation'] == 'logistic regression':
        mcc = LRC(verbose=verbose)
        # Train and test classifier
        mcc.fit_predict(DATA_response, LABELS_response, class_names=class_names,
                        CV_scheme=PARAMETERS['CV_scheme'],
                        v_fraction=PARAMETERS['v_fraction'],
                        shuffle_split_repeats=PARAMETERS['shuffle_split_repeats'],
                        fraction_training=PARAMETERS['fraction_training'],
                        fraction_training_for_tuning=PARAMETERS['fraction_training_for_tuning'],
                        max_n_splits_training=PARAMETERS['max_n_splits_training'],
                        max_n_splits_tuning=PARAMETERS['max_n_splits_tuning'],
                        data_shape=PARAMETERS['data_shape'],
                        data_scaling=PARAMETERS['data_scaling'],
                        random_state=seed)

        print(mcc.performance)

    # Abort if classifier could not be fit
    if not mcc.is_fit:
        return

    # Count number of times the classifier correctly identified SP
    SP_idx = np.where(np.in1d(mcc.class_names, 'SP'))[0]
    if SP_idx.shape[0] > 0:
        SP_idx = SP_idx[0]

        # Get source of values to copy
        if PARAMETERS['significance_test_n_permutations'] > 0:
            predicted_labels = mcc.sig_predicted_labels
        else:
            predicted_labels = mcc.predicted_labels

        idx_this_class = np.where(mcc.y == SP_idx)[0]
        mcc.correct_classifications[SP_idx] = np.where(predicted_labels[idx_this_class] == SP_idx)[0].shape[0]
        mcc.correct_classifications_perc[SP_idx] = mcc.correct_classifications[SP_idx] / idx_this_class.shape[0]

    # Store model performance
    RESULTS = dict()
    RESULTS['performance'] = dict()

    # Get source of values to copy
    if PARAMETERS['significance_test_n_permutations'] > 0 and PARAMETERS['decoder_implementation'] == 'random forest':
        performance = mcc.sig_performance
        performance_mean = mcc.sig_performance_mean
        confusion_matrix = mcc.sig_confusion_matrix
    else:
        performance = mcc.performance
        performance_mean = mcc.performance_mean
        confusion_matrix = mcc.confusion_matrix

    # Make sure that all fields are converted to string
    fn = lambda x: x if isinstance(x, str) else str(x)  # Make lambda function to pass to dictionary- and list-comprehensions below
    RESULTS['performance']['AUPRC'] = {k: fn(v) for k, v in zip(class_names, performance)}
    RESULTS['performance']['AUPRC_mean'] = fn(performance_mean)
    conf_mat = [[''] + class_names]
    for i_class in range(n_classes):
        this_class = class_names[i_class]
        row = [this_class] + list(confusion_matrix[i_class, :])
        conf_mat.append(row)
    RESULTS['performance']['confusion_matrix'] = [[fn(sublist_item) for sublist_item in list_item] for list_item in conf_mat]
    RESULTS['n_correct'] = mcc.correct_classifications.tolist()
    RESULTS['percent_correct'] = mcc.correct_classifications_perc.tolist()
    RESULTS['class_names'] = mcc.class_names
    # Write results in JSON format
    with open(PARAMETERS['output_filenames']['results'], 'w+') as f:
        json.dump(RESULTS, f)

    # Store classifier to serialized object
    pickle.dump(mcc, open(PARAMETERS['output_filenames']['classifier'], 'wb'))

    # Log end of function
    if verbose:
        log('done')


################################################################################
# Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        run(input_filename=r'D:\_MATLAB_CaImaging\FK_5\FK_5_cond1_decoding_data.mat', verbose=True)
