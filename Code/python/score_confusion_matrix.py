import sys
import os
import re
import pickle
import numpy as np
import json
from pycm import *

def run(root_folder, analysis_type, animal_list, metric,   verbose=False):
    # Split names of animals
    animal_list = [i.strip() for i in animal_list.split(',')]

    # Set default classes
    class_names =  ['HPS', 'FPS', 'temp_48', 'temp_43', 'pinprick', 'puff', 'touch', 'sound',]#epi : ['cold', 'heat', 'pinprick', 'touch'] # #'temp_38', 'odor_plus'
    n_classes = len(class_names)

    # Initialize local variables
    search_recursively = False
    confusion_matrix = None
    condition_names = None

    # Initialize confusion matrices and list of condition names
    if analysis_type == '':
        confusion_matrix = np.empty((n_classes, n_classes, 0), dtype=float)
        condition_names = list()
        search_recursively = False

    elif analysis_type == 'selective':
        confusion_matrix = {cls: np.empty((n_classes, n_classes, 0), dtype=float) for cls in class_names}
        condition_names = {cls: list() for cls in class_names}
        search_recursively = True

    elif analysis_type == 'stable':
        confusion_matrix = {cls: np.empty((n_classes, n_classes, 0), dtype=float) for cls in class_names}
        condition_names = {cls: list() for cls in class_names}
        search_recursively = True

    elif analysis_type == 'ensembles':
        confusion_matrix = np.empty((n_classes, n_classes, 0), dtype=float)
        condition_names = list()
        search_recursively = True

    elif str.startswith(analysis_type, 'H_'):
        confusion_matrix = {analysis_type: np.empty((n_classes, n_classes, 0), dtype=float)}
        condition_names = {analysis_type: list()}  # for cls in class_names}
        search_recursively = True

    # Visit each folder and subfolder to look for classifiers
    for animal_ID in animal_list:
        # Look for files in the appropriate folder
        folder = os.path.join(root_folder, animal_ID, analysis_type)
        if not os.path.exists(folder):
            continue

        # Get list of folders to analyze
        if search_recursively:
            list_of_folders = [os.path.join(folder, i) for i in os.listdir(folder)]
        else:
            list_of_folders = [folder]

        for folder in list_of_folders:
            this_confusion_matrix, this_condition_names = _read_confusion_matrices_in_folder(animal_ID, folder, class_names, verbose=verbose)

            # Concatenate data
            if analysis_type == '':
                confusion_matrix = np.dstack((confusion_matrix, this_confusion_matrix))
                condition_names += this_condition_names

            elif analysis_type == 'ensembles':
                confusion_matrix = np.dstack((confusion_matrix, this_confusion_matrix))
                condition_names += this_condition_names



            else:
                # Look for stimulus to which cells were selective
                folder_name_parts = os.path.basename(folder).split('__')
                selective_to = folder_name_parts[1]
                if selective_to == 'HPS_temp_48':
                    continue

                # Get when cells were selective
                if analysis_type == 'selective':
                    selective_when = int(folder_name_parts[0].split('session')[1])
                    [i.append(selective_when) for i in this_condition_names]

                elif analysis_type == 'stable':
                    selective_when = folder_name_parts[0].split('session')[1]
                    sessions = re.findall(r'\d+', selective_when)
                    sessions_sign = re.split(r'\d+', selective_when)
                    selective_when = [int('%s%s' % (i, j)) for i, j in zip(sessions_sign, sessions)]
                    this_condition_names = [i + selective_when for i in this_condition_names]

                elif str.startswith(analysis_type, 'H_'):
                    selective_when = folder_name_parts[0].split('session')[1]
                    [i.append(selective_when) for i in this_condition_names]
                    selective_to = analysis_type
                # elif analysis_type == 'ensembles':
                #     selective_when = int(folder_name_parts[0].split('session')[1])
                #     [i.append(selective_when) for i in this_condition_names]

                # Append data
                confusion_matrix[selective_to] = np.dstack((confusion_matrix[selective_to], this_confusion_matrix))
                condition_names[selective_to] += this_condition_names

        class_names_arr = np.array(class_names)
        if str.startswith(analysis_type, 'H_'):
            confusion_matrix = confusion_matrix[analysis_type]
            confusion_matrix.astype(int)

        # Perform score of the Confussion matrices
        for i_cond in range(np.size(confusion_matrix,2)):

            conf_mat = confusion_matrix[:,:,i_cond]
            mask_log = ~(np.all(np.isnan(conf_mat), axis = 0))
            mask = conf_mat[np.isnan(conf_mat) == False]
            conf_mat= mask.reshape(np.sqrt(len(mask)),np.sqrt(len(mask)))
            conf_mat_dict =  dict()
            these_class_names  = class_names_arr[mask_log]
            AUC = (np.zeros((len(these_class_names),)) - 1).astype(object)
            # Iterate through rows
            for i_row in range(conf_mat.shape[0]):
                row = conf_mat[i_row,].astype(np.int)
                row_dict = {class_name: value for class_name, value in
                            zip(these_class_names, row)}
                conf_mat_dict[these_class_names[i_row]] = row_dict

            matrix = conf_mat_dict
            cm = ConfusionMatrix(matrix= matrix)
            # Get metric values per class
            if metric == 'PR':
                metric_values = list(getattr(cm, 'AUPR').values())

            elif metric == 'ROC':
                metric_values = list(getattr(cm, 'AUC').values())

            elif metric == 'ACC':  # Accuracy
                metric_values = list(getattr(cm, 'ACC').values())

            elif metric == 'sInd':  # Similarity index
                metric_values = list(getattr(cm, 'sInd').values())

            elif metric == 'AGM':  # Similarity index
                metric_values = list(getattr(cm, 'AGM').values())

            elif metric == 'Kappa':  # Similarity index
                metric_values = getattr(cm, metric)

            else:
                metric_values = list(getattr(cm, metric).values())

            AUC = metric_values
            missing_value = np.nan
            AUC = np.array([missing_value if i == 'None' else i for i in AUC],
                           dtype=float)

            # Format output
            # Normalize AUC to [0, 1] interval
            r = 1 / 2
            if metric == 'PR':
                for idx, label in enumerate(cm.classes):
                    # Get number of elements in each group
                    n_pos_observations = cm.TP[label] #np.where(y_true == label)[0].shape[0]
                    n_neg_observations = cm.TN[label] #np.where(y_true != label)[0].shape[0]
                    # Compute performance of random classifier
                    random_classifier_AUC = divide0(n_pos_observations, n_pos_observations + n_neg_observations, replace_with=0)
                    # Normalize value
                    AUC[idx] = r + (AUC[idx] - random_classifier_AUC) * (1 - r) / (np.abs(random_classifier_AUC - r) + r)

            # # Set recursively back to false
            # if analysis_type == 'ensembles':
            #     search_recursively = False
            #
            # if search_recursively:
            #     # Remove stimuli not analyzed
            #     stimuli_analyzed = list(condition_names.keys())
            #     stimuli_to_discard = [i for i in stimuli_analyzed if len(condition_names[i]) == 0]
            #     [condition_names.pop(i) for i in stimuli_to_discard]
            #     [confusion_matrix.pop(i) for i in stimuli_to_discard]
            #     # Convert output
            #     stimuli_analyzed = list(confusion_matrix.keys())
            #     confusion_matrix = {i: confusion_matrix[i].tolist() for i in stimuli_analyzed}
            file_s =str(animal_ID + '_cond' + str(i_cond+1) +'_decoding_results_ACC.json')
            output_filename=os.path.join(root_folder,animal_ID, file_s)
            # Write results in JSON format
            if search_recursively:
                with open(output_filename, 'w+') as f:
                    json.dump((), f)
            else:
                with open(output_filename, 'w+') as f:
                    json.dump({key: value for key, value in
                            zip(cm.classes, AUC)}, f)
            print('kappa = ', str(cm.Kappa))

    if verbose:
        print('Completed')


def _read_confusion_matrices_in_folder(animal_ID, folder, class_names, verbose=False):
    # Initialize output variable
    n_classes = len(class_names)
    confusion_matrix = np.empty((n_classes, n_classes, 0), dtype=float)
    condition_names = list()

    files = os.listdir(folder)
    files = sorted([i for i in files if i.endswith('.p')])
    for f in files:
        if verbose:
            print(os.path.join(folder, f))

        # Load classifier
        mcc = pickle.load(open(os.path.join(folder, f), 'rb'))
        # Take confusion matrix
        if mcc.sig_all_confusion_matrix is not None:
            confusion_matrix_this_cond = mcc.sig_all_confusion_matrix
        else:
            confusion_matrix_this_cond = mcc.all_confusion_matrix
        n_samples = mcc.n_samples # np.array(np.round(confusion_matrix_this_cond * n_samples, 0), np.int)
        data = np.sum(np.dstack(confusion_matrix_this_cond.values()), 2)
        # Allocate confusion matrix
        this_confusion_matrix = np.empty((n_classes, n_classes), dtype=int) * np.nan

        # Get list of stimuli to consider (https://stackoverflow.com/a/23529016)
        stimuli_to_consider = [x for x in mcc.class_names if x in frozenset(class_names)]

        # Place stimuli in correct location of output variable
        for i_stim, stim_i in enumerate(stimuli_to_consider):
            row = class_names.index(stim_i)
            for j_stim, stim_j in enumerate(stimuli_to_consider):
                col = class_names.index(stim_j)
                # Copy value
                this_confusion_matrix[row, col] = data[i_stim, j_stim]

        # Concatenate results
        confusion_matrix = np.dstack((confusion_matrix, this_confusion_matrix))

        # Store name of condition
        name = f.split('_')[2]
        condition_names.append([animal_ID, name])

    return confusion_matrix, condition_names
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




################################################################################
# Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        run(root_folder='V:\\Ca_imaging_pain\\6_data\\ACC_CCI_anesth\\response_decoding',
            analysis_type='H_noci',
            animal_list= #'MA_15',#, FK_6, FK_7, FK_8, FK_9, FK_10, FK_11, FK_12, FK_14, FK_17, FK_18, FK_19, FK_20, MA_15, MA_16, MA_17, MA_23, MA_24, MA_25, MA_26, MA_27',
            'FK_10',
            #'FK_20',#, MA_30epi, MA_32epi, MA_33epi',
             #'FK_5, FK_6, FK_7, FK_8, FK_9, FK_10, FK_11, FK_12, FK_14, FK_17, FK_18, FK_19, FK_20, MA_15, MA_16, MA_17, MA_23, MA_24, MA_25, MA_26, MA_27',
            # animal_list='MA_23, MA_24, MA_25, MA_26, MA_27',
            # for TOUCH use these:       animal_list='FK_10, FK_11, FK_12, FK_17, FK_18, FK_19, FK_20, MA_15, MA_16, MA_17, MA_23, MA_24, MA_25, MA_26, MA_27', #
            metric = 'ACC',
            # class_names=['HPS', 'puff', 'temp_48'],
            verbose=True)

