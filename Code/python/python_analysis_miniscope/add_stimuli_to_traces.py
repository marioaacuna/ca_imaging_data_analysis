import os
from scipy.io.matlab import loadmat, savemat
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
# plt.ion()

animal_ID = 'I38_4'
pain_session = 'pain_2' # pain_ or pain_2
date = '190627'
signal_type = 'traces'  # 'spikes' or 'traces'
frame_rate = 5
pre_stimulus_window_s = 5
post_stimulus_window_s = 5


# Load traces
traces_root_path = r'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\%s\jointExtraction\sorted' % animal_ID
traces_filename = os.path.join(traces_root_path, '%s_%s.mat' % (signal_type, pain_session))
traces = loadmat(traces_filename).get('%s_%s' % (signal_type, pain_session))
# Make time axis
n_frames = traces.shape[1]
duration = n_frames / frame_rate
time = np.linspace(0, duration, traces.shape[1])

# Load timestamps
timestamps_root_path = r'T:\Mario\Behavior\ACC_pain\behavior cams\%s\%s\pain' % (animal_ID, date)
timestamps_filename = os.path.join(timestamps_root_path, 'stimuli_timestamps.csv')
timestamps = pd.read_csv(timestamps_filename)
n_trials = timestamps.shape[0]
# Add frame number to each stimulus
timestamps['frame'] = 0
for i_trial in range(n_trials):
    this_timestamp = timestamps.loc[i_trial, 'Time']
    timestamps.loc[i_trial, 'frame'] = np.abs(time - this_timestamp).argmin()

# Convert windows to number of frames
pre_stimulus_window = int(np.round(pre_stimulus_window_s * frame_rate))
post_stimulus_window = int(np.round(post_stimulus_window_s * frame_rate))

# Make new array with shape [n_cells x time x trials]
n_cells = traces.shape[0]
n_timepoints = pre_stimulus_window + post_stimulus_window
unique_stimuli, n_trials_per_stimulus = np.unique(timestamps['Stimulus_type'], return_counts=True)
n_stimuli = unique_stimuli.shape[0]
data = dict({unique_stimuli[i]: np.zeros((n_cells, n_timepoints, 0), dtype=traces.dtype) for i in range(n_stimuli)})
for i_trial in range(n_trials):
    stimulus_onset = timestamps.loc[i_trial, 'frame']
    if stimulus_onset < pre_stimulus_window:
        raise Exception('Check me!')
    if stimulus_onset > n_frames - post_stimulus_window:
        raise Exception('Check me!, trial :' + str(i_trial) )

    stimulus_type = timestamps.loc[i_trial, 'Stimulus_type']
    this_trial = traces[:, stimulus_onset-pre_stimulus_window:stimulus_onset+post_stimulus_window]
    data[stimulus_type] = np.dstack((data[stimulus_type], this_trial))

# Store file to disk
output_filename = os.path.join(traces_root_path, '%s_%s_stamped.mat' % (signal_type, pain_session))
savemat(output_filename, {'%s_pain_TS' % signal_type: data}, appendmat=False)

print('done')
