import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Utilities import matlab_file
from decoding.third_party.matplotlib_venn import venn2, venn3

input_filename = r'C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP\_data\count_selective_cells.mat'
stimuli_groups = [
                  ['nox. heat', ['temp_48']],
                  ['shock', ['HPS', 'FPS']],
                  ['inn. heat', ['temp_43', 'temp_38']],
                  ['aversive', ['puff', 'sound']],
                  ['heat 43', ['temp_43']],
                  ['heat 38', ['temp_38']],
                  ]
stimuli_to_take = ['nox. heat', 'shock', 'aversive']

# Set colors
color_nociceptive = [.42, .78, .76]
color_shock = [.55, .55, .55]
color_aversive = [.96, .56, .44]
color_heat_43 = [.45, .26, .59]
color_heat_38 = [.98, .93, .27]

# Load variables in Matlab file
variables = matlab_file.load(input_filename)
CELL_COUNT = variables['CELL_COUNT']
TIMEPOINTS = variables['TIMEPOINTS']
STIMULI = variables['STIMULI'].tolist()
STIMULI_table = variables['STIMULI_table']
del variables

# For each animal, take sessions before surgery
animals = list(CELL_COUNT.keys())
n_animals = len(animals)
n_STIMULI = len(STIMULI)
n_combinations = STIMULI_table.shape[0]
cell_count_CCI_pre = np.zeros((n_combinations, n_animals), dtype=np.float32) * np.nan
cell_count_CCI_post = np.zeros((n_combinations, n_animals), dtype=np.float32) * np.nan
cell_count_sham_pre = np.zeros((n_combinations, n_animals), dtype=np.float32) * np.nan
cell_count_sham_post = np.zeros((n_combinations, n_animals), dtype=np.float32) * np.nan
n_cells_per_stimuli_combination = np.zeros((n_combinations, ), dtype=int)
n_total_cells = sum([i[3] for i in TIMEPOINTS])

for iid, animal_ID in enumerate(animals):
    # Get indices of these sessions
    this_animal_idx = np.in1d(TIMEPOINTS[:, 0], animal_ID)
    timepoints = TIMEPOINTS[this_animal_idx, 2][0]
    experimental_group = TIMEPOINTS[this_animal_idx, 1][0]

    sessions_before_surgery = np.where(timepoints <= 0)[0]
    sessions_after_surgery = np.where(timepoints > 0)[0]
    if sessions_before_surgery.shape[0] == 0 or sessions_after_surgery.shape[0] == 0:
        continue
    n_cells = TIMEPOINTS[this_animal_idx, 3][0]

    # Get the corresponding cell counts
    this_cell_count = CELL_COUNT[animal_ID]

    if sessions_before_surgery.shape[0] == 1:
        idx_before = 0
    else:
        idx_before = -1
    if sessions_after_surgery.shape[0] == 1:
        idx_after = 0
    else:
        idx_after = 1

    this_cell_count_before_surgery = this_cell_count[:, sessions_before_surgery[idx_before]]
    this_cell_count_after_surgery = this_cell_count[:, sessions_after_surgery[idx_after]]
    if this_cell_count_before_surgery.ndim == 1:
        this_cell_count_before_surgery = this_cell_count_before_surgery.reshape(-1, 1)
    if this_cell_count_after_surgery.ndim == 1:
        this_cell_count_after_surgery = this_cell_count_after_surgery.reshape(-1, 1)
    # Average across sessions
    this_cell_count_before_surgery = np.nanmean(this_cell_count_before_surgery, axis=1)
    this_cell_count_after_surgery = np.nanmean(this_cell_count_after_surgery, axis=1)

    # Add number of cells that were stimulated with these stimuli
    n_cells_per_stimuli_combination[np.logical_not(np.isnan(this_cell_count_before_surgery))] += n_cells

    # Store values
    if experimental_group == 'CCI':
        cell_count_CCI_pre[:, iid] = this_cell_count_before_surgery
        cell_count_CCI_post[:, iid] = this_cell_count_after_surgery
    elif experimental_group == 'sham':
        cell_count_sham_pre[:, iid] = this_cell_count_before_surgery
        cell_count_sham_post[:, iid] = this_cell_count_after_surgery
    else:
        continue

plt.figure(facecolor='w')
for i_iter in range(4):
    if i_iter == 2:
        cell_count = cell_count_CCI_pre.copy()
        title_str = 'CCI before'
    elif i_iter == 3:
        cell_count = cell_count_CCI_post.copy()
        title_str = 'CCI after'
    elif i_iter == 0:
        cell_count = cell_count_sham_pre.copy()
        title_str = 'sham before'
    elif i_iter == 1:
        cell_count = cell_count_sham_post.copy()
        title_str = 'sham after'

    # Convert data for plotting function
    data = pd.DataFrame(columns=['value'] + STIMULI)
    for iid, animal_ID in enumerate(animals):
        new_data = pd.DataFrame(STIMULI_table, columns=STIMULI)
        new_data['value'] = cell_count[:, iid]
        data = pd.concat((data, new_data), ignore_index=True, sort=False)
    # Remove 0s
    data.loc[np.where(data['value'] == 0)[0], 'value'] = np.nan
    # Convert selectivity for stimulus to boolean index
    for stimulus_name in STIMULI:
        data[stimulus_name] = data[stimulus_name].map(np.bool)

    # Convert count of cells to percentages
    data.set_index(STIMULI, inplace=True)
    data.dropna(inplace=True)
    data = data.value.groupby(level = list(range(n_STIMULI))).sum()
    data = data.reset_index()
    for ii in range(data.shape[0]):
        idx = np.all(STIMULI_table == data.loc[ii, STIMULI].values, axis=1).nonzero()[0][0]
        n = n_cells_per_stimuli_combination[idx]
        data.loc[ii, 'value'] /= n

    # Get number of non-selective cells
    data_to_show = pd.DataFrame(data['value'].values, columns=['value'])
    n_stimuli = len(stimuli_groups)
    group_names = list()
    for i_group in range(n_stimuli):
        group_name = stimuli_groups[i_group][0]
        group_names.append(group_name)
        stimuli_to_group = stimuli_groups[i_group][1]
        if len(stimuli_to_group) == 1:
            data_to_show[group_name] = data[stimuli_to_group]
        else:
            data_to_show[group_name] = np.any(data[stimuli_to_group].values == True, axis=1)
    data = data_to_show.set_index(group_names)

    data_to_plot = data.value.groupby(level = list(range(n_stimuli))).sum()
    data_to_plot = data_to_plot / data_to_plot.sum() * 100
    data_to_plot = data_to_plot.reset_index()

    data_to_show = data_to_plot[stimuli_to_take + ['value']].values.astype(np.float32)
    labels = [''.join(['%i' % i for i in x]) for x in data_to_show[:, :len(stimuli_to_take)].astype(int)]
    subsets = {lab: 0 for lab in labels}
    for ii, lab in enumerate(labels):
        value = data_to_show[ii, -1]
        subsets[lab] += value
    n = [data_to_show[data_to_show[:, i] == 1, -1].sum() for i in range(len(stimuli_to_take))]
    set_labels = [i + ' (%i%%)' % np.round(j, decimals=0) for i, j in zip(stimuli_to_take, n)]
    plt.subplot(1, 4, i_iter + 1)
    V = venn3(subsets, set_labels=set_labels, set_colors=[color_nociceptive, color_shock, color_aversive], alpha=1, subset_label_formatter=lambda x: '%i%%' % np.round(x, decimals=0))
    plt.title(title_str)

plt.show()
