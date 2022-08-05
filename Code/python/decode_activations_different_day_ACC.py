# System and IO packages
import os
import sys
import json
import pickle
from itertools import permutations

# Numerical packages
import numpy as np
import pandas as pd

# Local repository
from Utilities.IO_operations import log


def run(folder, output_filename_results, output_filename_confusion_matrix, metric, verbose=False):
    # Quit if folder does not exist
    if not os.path.exists(folder):
        log('\'%s\' does not exist' % folder)
        return
    # Look for and load stored classifiers
    files = os.listdir(folder)
    classifiers_to_load = sorted([i for i in files if i.endswith('_decoding_classifier.p')])
    n_classifiers = len(classifiers_to_load)
    if n_classifiers == 0:
        log('No classifiers in \'%s\'' % folder)
        exit(0)
    results_to_load = sorted([i for i in files if i.endswith('_decoding_results_ACC.json')])
    if verbose:
        log('Found %i classifiers in \'%s\'' % (n_classifiers, folder))
    classifiers = [pickle.load(open(os.path.join(folder, f), 'rb')) for f in classifiers_to_load]
    results = [json.load(open(os.path.join(folder, f), 'r')) for f in results_to_load]

    # Get names of all classes
    all_class_names = np.unique(np.hstack([i.class_names for i in classifiers]))
    # Initialize output variables
    PERFORMANCE = {key: list() for key in all_class_names}
    CONFUSION_MATRICES = list()

    # Iterate through pairs of classifiers
    if n_classifiers >= 2:
        pairs = list(permutations(np.arange(n_classifiers), 2))
        n_pairs = len(pairs)
        for i_pair in range(n_pairs):
            if verbose:
                log('Testing classifier #%i on data of condition #%i' % (pairs[i_pair][0] + 1, pairs[i_pair][1] + 1))
            # Get the classifier to test
            clf_day1 = classifiers[pairs[i_pair][0]]
            trained_classes_day1 = clf_day1.class_names
            # Get the data for testing from another classifier object
            clf_day2 = classifiers[pairs[i_pair][1]]
            data_day2 = clf_day2.X

            # Predict posterior probabilities
            posterior_probabilities_day2_with_clf_day1 = list()
            for clf in clf_day1.trained_classifiers:
                posterior_probabilities_day2_with_clf_day1.append(clf.predict_proba(data_day2))
            posterior_probabilities_day2_with_clf_day1 = np.dstack(posterior_probabilities_day2_with_clf_day1)
            posterior_probabilities_day2_with_clf_day1 = np.mean(posterior_probabilities_day2_with_clf_day1, axis=2)

            # Make sure there are probabilities for all classes, even those that the
            # classifier hasn't seen yet.
            posterior_probabilities_day2_with_clf_day1 = pd.DataFrame(posterior_probabilities_day2_with_clf_day1, columns=trained_classes_day1)
            # Add 0 probability to untrained stimuli
            tested_classes_day2 = clf_day2.class_names
            untrained_stimuli_day2 = list(set(trained_classes_day1).symmetric_difference(tested_classes_day2))
            if len(untrained_stimuli_day2) > 0:
                for st in untrained_stimuli_day2:
                    posterior_probabilities_day2_with_clf_day1[st] = 0
            # Re-order columns, as expected by other classifier
            posterior_probabilities_day2_with_clf_day1 = posterior_probabilities_day2_with_clf_day1[tested_classes_day2].values

            # Compute AUC of precision-recall curve
            labels_day2 = clf_day2.y
            AUPRC, cm = clf_day2.compute_AUC(labels_day2, posterior_probabilities_day2_with_clf_day1, metric=metric, return_average=False, return_confusion_matrix=True, normalize_AUPRC_to_random=True)
            AUPRC_day2_with_clf_day1 = {key: value for key, value in zip(tested_classes_day2, AUPRC)}
            AUPRC_day2_with_clf_day2 = results[pairs[i_pair][1]]#results[pairs[i_pair][1]]['performance']['AUPRC']

            # Format confusion matrix
            conf_mat = [[''] + tested_classes_day2]
            for i_class in range(len(tested_classes_day2)):
                this_class = tested_classes_day2[i_class]
                row = [this_class] + ['%i' % i for i in cm[i_class, :]]
                conf_mat.append(row)

            # Store data
            [PERFORMANCE[key].append([str(pairs[i_pair][0] + 1), str(pairs[i_pair][1] + 1), str(AUPRC_day2_with_clf_day2[key]), str(AUPRC_day2_with_clf_day1[key])]) for key in tested_classes_day2]
            CONFUSION_MATRICES.append([str(pairs[i_pair][0] + 1), str(pairs[i_pair][1] + 1), conf_mat])

    else:
        AUPRC_day2_with_clf_day2 = results[0]['performance']['AUPRC']
        [PERFORMANCE[key].append([None, '1', AUPRC_day2_with_clf_day2[key], '']) for key in all_class_names]

    # Write results in JSON format
    with open(output_filename_results, 'w+') as f:
        json.dump(PERFORMANCE, f)

    if n_classifiers >= 2:
        with open(output_filename_confusion_matrix, 'w+') as f:
            json.dump(CONFUSION_MATRICES, f)

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
        folder = r'V:\Ca_imaging_pain\6_data\ACC_CCI_anesth\\response_decoding\\FK_5\\'
        run(folder=folder, output_filename_results=r'V:\Ca_imaging_pain\6_data\ACC_CCI_anesth\\response_decoding\\FK_6\\FK_6_crossday_results_ACC.json',
            output_filename_confusion_matrix='',
            metric ='ACC',
            verbose=True)
