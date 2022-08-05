# Authors: Jean-Remi King <jeanremi.king@gmail.com>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
#          Denis Engemann <denis.engemann@gmail.com>
#
# License: BSD (3-clause)

import matplotlib.pyplot as plt

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression

import mne
from mne.datasets import sample
from mne.decoding import GeneralizingEstimator

# Preprocess data
data_path = sample.data_path()
# Load and filter data, set up epochs
raw_fname = data_path + '/MEG/sample/sample_audvis_filt-0-40_raw.fif'
events_fname = data_path + '/MEG/sample/sample_audvis_filt-0-40_raw-eve.fif'
raw = mne.io.read_raw_fif(raw_fname, preload=True)
picks = mne.pick_types(raw.info, meg=True, exclude='bads')  # Pick MEG channels
raw.filter(1., 30., fir_design='firwin')  # Band pass filtering signals
events = mne.read_events(events_fname)
event_id = {'Auditory/Left': 1, 'Auditory/Right': 2,
            'Visual/Left': 3, 'Visual/Right': 4}
tmin = -0.050
tmax = 0.400
decim = 2  # decimate to make the example faster to run
epochs = mne.Epochs(raw, events, event_id=event_id, tmin=tmin, tmax=tmax,
                    proj=True, picks=picks, baseline=None, preload=True,
                    reject=dict(mag=5e-12), decim=decim)

clf = make_pipeline(StandardScaler(), LogisticRegression(solver='lbfgs'))
time_gen = GeneralizingEstimator(clf, scoring='roc_auc', n_jobs=1, verbose=True)

# Fit classifiers on the epochs where the stimulus was presented to the left.
# Note that the experimental condition y indicates auditory or visual
time_gen.fit(X=epochs['Left'].get_data(), y=epochs['Left'].events[:, 2] > 2)

X = epochs['Right'].get_data().copy()
y = epochs['Right'].events[:, 2] > 2
# import numpy as np
# X = np.tile(X, (1,1,2))
scores = time_gen.score(X=X, y=y)

fig, ax = plt.subplots(1)
im = ax.matshow(scores, vmin=0, vmax=1., cmap='RdBu_r', origin='lower', extent=epochs.times[[0, -1, 0, -1]])
ax.axhline(0., color='k'); ax.axvline(0., color='k')
ax.xaxis.set_ticks_position('bottom')
ax.set_xlabel('Testing Time (s)')
ax.set_ylabel('Training Time (s)')
ax.set_title('Generalization across time and condition')
plt.colorbar(im, ax=ax)
plt.show()
