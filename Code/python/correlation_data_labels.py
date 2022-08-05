import matplotlib.pyplot as plt
import numpy as np
from Utilities import matlab_file
import pandas as pd
import seaborn as sns
sns.set(); np.random.seed(0)
from scipy.stats import norm

cc = list()

for (animal_ID, condition_no) in [('MA_13', 3), ('FK_2', 4), ('FK_1', 4)]:
    input_filename = r'D:\_MATLAB_2PI\%s\%s_cond%i_decoding_data.mat' % (animal_ID, animal_ID, condition_no)
    PARAMETERS = matlab_file.load(input_filename)
    class_names = list(PARAMETERS['data'].keys())
    n_classes = len(class_names)
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
        these_labels = np.zeros((n_samples, ), dtype=int) + ci
        LABELS.append(these_labels)
    # Concatenate the labels array, and make data array
    old_LABELS = np.hstack(LABELS)
    LABELS = old_LABELS.copy()
    # LABELS[np.in1d(old_LABELS, np.where(np.in1d(target_names, 'SP'))[0])] = 0
    # LABELS[np.in1d(old_LABELS, np.where(~np.in1d(target_names, 'SP'))[0])] = 1
    DATA = np.vstack(([PARAMETERS['data'][cl] for cl in class_names]))

    df = pd.DataFrame(np.hstack((LABELS.reshape(-1, 1), DATA)), columns=['label'] + ['cell_%i' % (i + 1) for i in range(DATA.shape[1])])
    df['label'] = df['label'].astype('category').cat.codes
    a = df.corr()['label'][:].values[1:]
    cc.append(a)

cc = np.hstack(cc)

print()

plt.clf(); sns.set_style('ticks'); ax = sns.distplot(cc, norm_hist=True, kde_kws={"color": "k", "lw": 3})






plt.clf(); sns.set_style('ticks')
ii = 2
ax = list()
for cl in range(3):
    plt.subplot(1, 3, cl+1)
    ax.append(sns.distplot(null_distribution[ii, cl, : ], norm_hist=True, kde_kws={"color": "k", "lw": 3}))
    plt.axvline(self.posterior_probabilities[ii, cl], color='r')
    plt.title(self.posterior_probabilities[ii, cl])
ax = np.hstack(ax)
xlim = np.vstack([ii.get_xlim() for ii in ax])
xlim = [np.min(xlim), np.max(xlim)]
[ii.set_xlim(xlim[0], xlim[1]) for ii in ax]
# ylim = np.vstack([ii.get_ylim() for ii in ax])
# ylim = [np.min(ylim), np.max(ylim)]
# [ii.set_ylim(ylim[0], ylim[1]) for ii in ax]











