import sys
import pickle
import numpy as np
import pandas as pd

from Utilities import matlab_file


def run(classifier_filename, data_filename, output_filename, verbose=False):

    # Load classifier
    if verbose:
        print('Loading classifier')
    mcc = pickle.load(open(classifier_filename, 'rb'))

    # Load SP data
    if verbose:
        print('Loading data and computing posterior probabilities')
    data = matlab_file.load(data_filename)

    # Compute posterior probabilities
    posterior_probabilities = list()
    for clf in mcc.trained_classifiers:
        posterior_probabilities.append(clf.predict_proba(data))
    posterior_probabilities = np.dstack(posterior_probabilities)
    # Average probabilities across classifiers in the ensemble
    posterior_probabilities = np.mean(posterior_probabilities, axis=2)
    # Get label of predictions
    predicted_labels = np.argmax(posterior_probabilities, axis=1)

    # Make DataFrame to collect results
    predicted_class, n_predicted_class = np.unique(predicted_labels, return_counts=True)
    n_predicted_class = n_predicted_class / n_predicted_class.sum()
    all_n_predicted_class = np.zeros((mcc.n_classes, ), dtype=float)
    all_n_predicted_class[np.in1d(np.arange(mcc.n_classes), predicted_class)] = n_predicted_class

    df = pd.DataFrame(all_n_predicted_class, columns=['probability'])
    df['stim'] = mcc.class_names
    df['perf'] = mcc.sig_performance
    df['perc'] = mcc.correct_classifications_perc
    df['p'] = np.mean(posterior_probabilities, axis=0)
    df.sort_values(by='probability', ascending=False, inplace=True)

    # Write results in JSON format
    df.to_csv(output_filename, index=False)

    if verbose:
        print('Completed')


################################################################################
# Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        animal_ID = 'FK_17'
        run(#classifier_filename=r'V:\Ca_imaging_pain\6_data\ACC_CCI_anesth\response_decoding\%s\%s_cond1_decoding_classifier.p' % (animal_ID, animal_ID),
            classifier_filename=r'V:\Ca_imaging_pain\6_data\ACC_CCI_anesth\response_decoding\%s\selective\session-5__HPS\%s_cond1_decoding_classifier.p' % (animal_ID, animal_ID),
            data_filename= r'D:\_MATLAB_CaImaging\SP_decoding_data.mat',
            output_filename= r'D:\_MATLAB_CaImaging\SP_decoding_results.csv',
            verbose=True)



