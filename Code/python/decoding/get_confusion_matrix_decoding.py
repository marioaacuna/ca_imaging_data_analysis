import sys
import os
import re
import pickle
import numpy as np
import json


def run(root_folder, analysis_type, animal_list, output_filename,  verbose=False):
    # Split names of animals
    animal_list = [i.strip() for i in animal_list.split(',')]

    # Set default classes
    # check first if they are grouped in categories
    if str.endswith(analysis_type, 'categories'):
        class_names = ['noci', 'ave', 'inn' ]  # epi : ['cold', 'heat', 'pinprick', 'touch'] # #'temp_38', 'odor_plus'
    else:
        # class_names =  ['HPS', 'FPS', 'temp_48', 'temp_43', 'pinprick', 'puff', 'touch', 'sound']
        class_names =  ['cold', 'heat', 'pinprick', 'puff','touch'] # #'temp_38', 'odor_plus'

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

    if analysis_type =='all_categories':
        confusion_matrix = np.empty((n_classes, n_classes, 0), dtype=float)
        condition_names = list()
        search_recursively = False

    elif analysis_type == 'selective' :
        confusion_matrix = {cls: np.empty((n_classes, n_classes, 0), dtype=float) for cls in class_names}
        condition_names = {cls: list() for cls in class_names}
        search_recursively = True

    elif analysis_type ==  'non_selective':
        confusion_matrix = {cls: np.empty((n_classes, n_classes, 0), dtype=float) for cls in class_names}
        condition_names = {cls: list() for cls in class_names}
        search_recursively = True

    elif str.startswith(analysis_type, 'H_'):
        confusion_matrix = {analysis_type: np.empty((n_classes, n_classes, 0), dtype=float)}
        condition_names = {analysis_type: list()}# for cls in class_names}
        search_recursively = True

    elif analysis_type == 'stable':
        confusion_matrix = {cls: np.empty((n_classes, n_classes, 0), dtype=float) for cls in class_names}
        condition_names = {cls: list() for cls in class_names}
        search_recursively = True

    elif analysis_type == 'ensembles':
        confusion_matrix = {
        cls: np.empty((n_classes, n_classes, 0), dtype=float) for cls in
        class_names}
        condition_names = {cls: list() for cls in class_names}
        search_recursively = True
        #confusion_matrix = np.empty((n_classes, n_classes, 0), dtype=float)
        #condition_names = list()
        #search_recursively = True


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

            elif analysis_type == 'all_categories': # if
                confusion_matrix = np.dstack((confusion_matrix, this_confusion_matrix))
                condition_names += this_condition_names

            else:
                # Look for stimulus to which cells were selective
                folder_name_parts = os.path.basename(folder).split('__')

                # Get when cells were selective
                if analysis_type == 'selective' :
                    selective_to = folder_name_parts[1]
                    if selective_to == 'SP':  # 'HPS_temp_48':
                        continue
                    selective_when = int(folder_name_parts[0].split('session')[1])
                    [i.append(selective_when) for i in this_condition_names]

                elif analysis_type == 'non_selective' :
                    selective_to = folder_name_parts[1]
                    if selective_to == 'SP':  # 'HPS_temp_48':
                        continue
                    selective_when = int(
                            folder_name_parts[0].split('session')[1])
                    [i.append(selective_when) for i in this_condition_names]

                elif str.startswith(analysis_type, 'H_'):
                    selective_when = int(
                            folder_name_parts[0].split('session')[1])
                    [i.append(selective_when) for i in this_condition_names]

                elif analysis_type == 'stable':
                    selective_to = folder_name_parts[1]
                    if selective_to == 'SP':  # 'HPS_temp_48':
                        continue
                    selective_when = folder_name_parts[0].split('session')[1]
                    sessions = re.findall(r'\d+', selective_when)
                    sessions_sign = re.split(r'\d+', selective_when)
                    selective_when = [int('%s%s' % (i, j)) for i, j in zip(sessions_sign, sessions)]
                    this_condition_names = [i + selective_when for i in this_condition_names]

                elif analysis_type == 'ensembles':
                    selective_to = folder_name_parts[1]
                    if selective_to == 'SP':  # 'HPS_temp_48':
                        continue
                    selective_when = int(folder_name_parts[0].split('session')[1])
                    [i.append(selective_when) for i in this_condition_names]

                # Append data
                if  str.startswith(analysis_type, 'H_'):
                    selective_to   = analysis_type
                    confusion_matrix[selective_to] = np.dstack((confusion_matrix[selective_to], this_confusion_matrix))
                    condition_names[selective_to] += this_condition_names

                else:
                    confusion_matrix[selective_to] = np.dstack((
                                                               confusion_matrix[
                                                                   selective_to],
                                                               this_confusion_matrix))
                    condition_names[selective_to] += this_condition_names

    # Format output
    # Set recursively back to false
    if analysis_type == 'ensembles':
        stimuli_analyzed = list(confusion_matrix.keys())
        confusion_matrix = {i: confusion_matrix[i].tolist() for i in
                            stimuli_analyzed}
        search_recursively = False

    if search_recursively:
        # Remove stimuli not analyzed
        stimuli_analyzed = list(condition_names.keys())
        stimuli_to_discard = [i for i in stimuli_analyzed if len(condition_names[i]) == 0]
        [condition_names.pop(i) for i in stimuli_to_discard]
        [confusion_matrix.pop(i) for i in stimuli_to_discard]
        # Convert output
        stimuli_analyzed = list(confusion_matrix.keys())
        confusion_matrix = {i: confusion_matrix[i].tolist() for i in stimuli_analyzed}

    # Write results in JSON format
    if search_recursively or analysis_type == 'ensembles':
        with open(output_filename, 'w+') as f:
            json.dump((confusion_matrix, condition_names, class_names), f)

    else:
        with open(output_filename, 'w+') as f:
            json.dump((confusion_matrix.tolist(), condition_names, class_names), f)

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
        data = np.sum(np.dstack(confusion_matrix_this_cond.values()), 2)
        # Allocate confusion matrix
        this_confusion_matrix = np.empty((n_classes, n_classes), dtype=float) * np.nan

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
# Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        run(root_folder='V:\\Ca_imaging_pain\\6_data\\ACC_SNI_anxiety\\response_decoding',
            analysis_type='',
            animal_list= #'FK_5, FK_6, FK_7, FK_8, FK_9, FK_10, FK_11, FK_12, FK_14, FK_17,FK_18, FK_19, FK_20,MA_15, MA_16, MA_17, MA_23, MA_24, MA_25, MA_26, MA_27',
            #'FK_20',#, MA_30epi, MA_32epi, MA_33epi',
            # 'FK_5, FK_6, FK_7, FK_8, FK_9, FK_10, FK_11, FK_12, FK_14, FK_17, FK_18, FK_19, FK_20, MA_15, MA_16, MA_17, MA_23, MA_24, MA_25, MA_26, MA_27',
            # animal_list='FK_5',
            # for TOUCH use these:       animal_list='FK_10, FK_11, FK_12, FK_17, FK_18, FK_19, FK_20, MA_15, MA_16, MA_17, MA_23, MA_24, MA_25, MA_26, MA_27', #
            'MA_29epi, MA_30epi, MA_32epi, MA_33epi, MA_34epi,MA_35epi,MA_36epi,MA_37epi,MA_38epi,MA_39epi, MA_40epi' ,
            output_filename='C:\\Users\\acuna\\Documents\\Two_photon_imaging_data_analysis\\Figures_paper_SNI_anxiety\\_data\\confusion_matrices_all_cells.json',
            verbose=True)

