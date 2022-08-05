import os
import numpy as np
import pickle
from matplotlib import pyplot as plt

from Utilities.IO_operations import log
from Utilities import matlab_file
from decoding.EnsembleTreeClassifier import EnsembleTreeClassifier as ETC

from decoding.treeexplainer.treeexplainer import TreeExplainer


input_filename = r'D:\_MATLAB_2PI\MA_3\MA_3_cond1_decoding_data.mat'

do_reload = False

filename = os.path.splitext(os.path.basename(input_filename))[0]
seed = int(np.clip(np.int64(np.sum([ord(char) for char in filename])), a_min=0, a_max=2 ** 32 - 1))
np.random.seed(seed)
log('Random number generator seed set to %i' % seed)

# Load parameter file
log('Analyzing \'%s\'' % input_filename)
PARAMETERS = matlab_file.load(input_filename)
# Unpack parameters
class_names = list(PARAMETERS['data'].keys())
n_classes = len(class_names)

if do_reload:
    # Make sure that data is organized in "tidy" format
    for cl in class_names:
        if PARAMETERS['data'][cl].ndim == 1:
            PARAMETERS['data'][cl] = np.atleast_2d(PARAMETERS['data'][cl]).transpose()

    # Prepare data and labels for multi-class classifier
    LABELS = []  # array of labels
    for ci, cl in enumerate(class_names):
        # Get number of samples and generate data to keep in the SAMPLE_IDX DataFrame
        n_samples = PARAMETERS['data'][cl].shape[0]
        # Make a list of labels for these data points
        these_labels = np.zeros((n_samples,), dtype=int) + ci
        LABELS.append(these_labels)
    # Concatenate the labels array, and make data array
    LABELS = np.hstack(LABELS)
    DATA = np.vstack(([PARAMETERS['data'][cl] for cl in class_names]))

    # Initialize classifier
    clf = ETC(ensemble_type='RandomForest',
              tree_split_criterion='entropy',
              verbose=True)
    # Train and test classifier
    clf.fit_predict(DATA, LABELS, class_names=class_names,
                    CV_scheme='hv-block', n_CV_splits=3,
                    fraction_training_for_tuning=.25, n_CV_splits_tuning=3,
                    calibration_method='sigmoid', n_CV_splits_calibration=3,
                    data_scaling='minmax', n_permutations=0)

    # Store classifier to disk
    pickle.dump((clf, TE, DATA), open(r'C:\Users\deluna\Downloads\clf.p', 'wb'))
    exit()

else:
    clf = pickle.load(open(r'C:\Users\deluna\Downloads\clf.p', 'rb'))

RF = clf.trained_classifiers[0].calibrated_classifiers_[0].base_estimator.named_steps['classifier']
estimator = RF.estimators_[0]

TE = TreeExplainer(RF, target_names=clf.class_names, n_jobs=-1, verbose=True)
TE.explain(DATA, compute_conditional_contribution=True, n_jobs=-1, verbose=True)


from decoding.third_party.treeinterpreter import treeinterpreter_old as ti
results = ti.compute_feature_contributions_ensemble(RF, clf.X[:2, :], compute_conditional_contribution=True, n_jobs=-1)

print()


# clf.compute_feature_contribution()
# clf.compute_feature_contribution()
print(clf.class_names)
clf.plot_feature_contribution_barplot(classes_to_show=['temp_48', 'HPS', 'SP'], show_errorbars=False, contribution_type='all', sort_features=True)

print()



from scipy.stats import ranksums
from itertools import combinations
import pandas as pd
pd.options.display.max_colwidth = 300
pd.options.display.max_columns = 50
pd.options.display.width = 1000

# correct_predictions = np.where((clf.y == i_class) & (clf.predicted_labels == i_class))[0]
# incorrect_predictions = np.where((clf.y == i_class) & (clf.predicted_labels != i_class))[0]
# data = np.zeros((clf.n_features, 3)) * np.nan

feature_pairs = np.vstack(list(combinations(np.arange(clf.n_features), 2)))
n_pairs = feature_pairs.shape[0]

percent_agreement = np.zeros((clf.n_features, clf.n_features, clf.n_classes)) * np.nan
percent_agreement_correct = percent_agreement.copy()
percent_agreement_incorrect = percent_agreement.copy()
for i_class in range(clf.n_classes):
    # Get samples in this class
    samples_this_class = np.where(clf.y == i_class)[0]
    n_samples_this_class = samples_this_class.shape[0]
    # Get samples correct and incorrect predictions
    correct_samples_this_class = np.array(np.intersect1d(samples_this_class, np.where(clf.predicted_labels == i_class)[0]))
    incorrect_samples_this_class = np.array(np.intersect1d(samples_this_class, np.where(clf.predicted_labels != i_class)[0]))
    n_correct = correct_samples_this_class.shape[0]
    n_incorrect = incorrect_samples_this_class.shape[0]

    # Get values of feature contributions
    all_contributions = (clf.relative_contribution[samples_this_class, :, i_class] > 0).astype(np.int8)
    all_contributions_correct = (clf.relative_contribution[correct_samples_this_class, :, i_class] > 0).astype(np.int8)
    all_contributions_incorrect = (clf.relative_contribution[incorrect_samples_this_class, :, i_class] > 0).astype(np.int8)

    for i_pair in range(n_pairs):
        f1 = feature_pairs[i_pair, 0]
        f2 = feature_pairs[i_pair, 1]

        x = np.sum((all_contributions[:, f1] - all_contributions[:, f2]) == 0) / n_samples_this_class * 100
        xc = np.sum((all_contributions_correct[:, f1] - all_contributions_correct[:, f2]) == 0) / n_correct * 100
        xi = np.sum((all_contributions_incorrect[:, f1] - all_contributions_incorrect[:, f2]) == 0) / n_incorrect * 100

        percent_agreement[f1, f2, i_class] = x
        percent_agreement[f2, f1, i_class] = x

        percent_agreement_correct[f1, f2, i_class] = xc
        percent_agreement_correct[f2, f1, i_class] = xc

        percent_agreement_incorrect[f1, f2, i_class] = xi
        percent_agreement_incorrect[f2, f1, i_class] = xi


from scipy.cluster import hierarchy
import seaborn as sns

sorting_orders = list()
for i_class in range(clf.n_classes):
    zc = percent_agreement_correct[:, :, i_class]
    zi = percent_agreement_incorrect[:, :, i_class]
    linkage = hierarchy.linkage(zc, method='weighted', metric='canberra')
    h = hierarchy.dendrogram(linkage, no_plot=True, color_threshold=-np.inf)
    sorting_orders.append(h['leaves'])

keep_order_from_class = None
plt.close('all'); fig, ax = plt.subplots(nrows=2, ncols=3); ax = ax.ravel()
for i_class in range(clf.n_classes):
    zc = percent_agreement_correct[:, :, i_class]
    zi = percent_agreement_incorrect[:, :, i_class]

    if keep_order_from_class is None:
        idx = i_class
    else:
        idx = keep_order_from_class
    order = sorting_orders[idx]
    zc = zc[:, order][order, :]
    zi = zi[:, order][order, :]

    maskc = np.zeros_like(zc); maskc[np.triu_indices_from(maskc)]=True
    maski = np.zeros_like(zi); maski[np.tril_indices_from(maskc)]=True
    this_ax = ax[i_class]
    sns.heatmap(zc, vmin=0, vmax=100, cmap='Reds', robust=True, square=True, mask=maskc, xticklabels=order, yticklabels=order, cbar=False, ax=this_ax)
    sns.heatmap(zi, vmin=0, vmax=100, cmap='Blues', robust=True, square=True, mask=maski, xticklabels=order, yticklabels=order, cbar=False, ax=this_ax)
    this_ax.set_title(clf.class_names[i_class])





