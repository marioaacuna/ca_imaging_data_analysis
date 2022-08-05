import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from Utilities import matlab_file
from decoding.third_party.upsetplot import UpSet

input_filename = r'C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP\_data\count_selective_cells.mat'
output_filename = r'V:\Ca_imaging_pain\6_plots\Figures_paper_FMP\fig1_barplot_overlap_v2.pdf'

# Each item in the list is: name of group to show, group of stimuli, bar plot color
# stimuli_groups = [
#     ['FPS', ['FPS'], 'red'],
#     ['HPS', ['HPS'], 'orange'],
#     ]
# stimuli_groups = [
#     ['pinprick', ['temp_48', 'HPS', 'HPS_short', 'HPS_50_percent', 'pinprick'], 'red'],
#     ['salient', ['FPS', 'puff', 'odor_plus'], 'orange'],
#     ['innocuous', ['temp_38', 'temp_43', 'touch', 'sound'], 'lightgray'],
#     ]


stimuli_groups = [


    ['Shock hindpaw', ['HPS'], 'deeppink'],
    ['Shock forepaw', ['FPS'], 'violet'],
    ['48°', ['temp_48'], 'red'],
    ['43°', ['temp_43'], 'mediumturquoise'],
    ['38°', ['temp_38'], 'mediumturquoise'],
    ['Pinprick', ['pinprick'], 'violet'],
    ['Airpuff', ['puff'], 'dodgerblue'],
    # ['Odor', ['odor_plus'], 'dodgerblue'],
    ['Sound', ['sound'], 'lightseagreen'],
    ['Touch', ['touch'], 'deepskyblue'],
                ]

# stimuli_groups = [
#     ['temp_48', ['temp_48'], 'red'],
#     ['temps', ['temp_38', 'temp_43'], 'mediumturquoise'],
#     ['HPS', ['HPS', 'HPS_short', 'HPS_50_percent'], 'deeppink'],
#     ['FPS', ['FPS'], 'violet'],
#     ['pin', ['pinprick'], 'violet'],
#     ['puff', ['puff'], 'dodgerblue'],
#     ['odor', ['odor_plus'], 'dodgerblue'],
#     ['Touch', ['touch'], 'deepskyblue'],
#     ['Sound', ['sound'], 'lightseagreen'],
#                 ]

timepoints_of_interest = [-2]
separate_groups = False

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
n_combinations = STIMULI_table.shape[0]
cell_count = {k: np.zeros((n_combinations, n_animals), dtype=np.float32) * np.nan for k in ['sham', 'CCI']}
n_cells_per_stimuli_combination = {k: np.zeros((n_combinations, ), dtype=int) for k in ['sham', 'CCI']}
animals_in_group = {k: list() for k in ['sham', 'CCI']}
n_STIMULI = len(STIMULI)
n_stimuli = len(stimuli_groups)


for iid, animal_ID in enumerate(animals):
    # Get indices of these sessions
    this_animal_idx = np.in1d(TIMEPOINTS[:, 0], animal_ID)
    experimental_group = TIMEPOINTS[this_animal_idx, 1][0]
    if not experimental_group in list(cell_count.keys()):
        continue
    timepoints = TIMEPOINTS[this_animal_idx, 2][0]

    # Get sessions of interest
    session_before_surgery_idx = np.where(timepoints < 0)[0]
    session_after_surgery_idx = np.where(timepoints > 0)[0]
    if session_before_surgery_idx.shape[0] >= 2:
        session_before_surgery_idx = session_before_surgery_idx[-2:]
    if session_after_surgery_idx.shape[0] >= 2:
        session_after_surgery_idx = session_after_surgery_idx[:2]

    sessions_of_interest = list()
    if -2 in timepoints_of_interest:
        sessions_of_interest.append(session_before_surgery_idx[0])
    if -1 in timepoints_of_interest:
        sessions_of_interest.append(session_before_surgery_idx[-1])
    if +1 in timepoints_of_interest:
        sessions_of_interest.append(session_after_surgery_idx[0])
    if +2 in timepoints_of_interest:
        sessions_of_interest.append(session_after_surgery_idx[-1])
    sessions_of_interest = np.unique(sessions_of_interest)

    print('%s - %s' % (animal_ID, experimental_group))
    # Get number of cells
    n_cells = TIMEPOINTS[this_animal_idx, 3][0]
    # Get the corresponding cell counts
    this_cell_count = CELL_COUNT[animal_ID]
    if this_cell_count.ndim == 1:
        this_cell_count = this_cell_count.reshape(-1, 1)
    this_cell_count = this_cell_count[:, sessions_of_interest]
    # Average across sessions
    if this_cell_count.ndim > 1:
        this_cell_count = np.nanmean(this_cell_count, axis=1)

    # Add number of cells that were stimulated with these stimuli
    n_cells_per_stimuli_combination[experimental_group][np.logical_not(np.isnan(this_cell_count))] += n_cells
    # Store values
    cell_count[experimental_group][:, iid] = this_cell_count
    animals_in_group[experimental_group].append(animal_ID)

if not separate_groups:
    # Combine all data in a field called 'all groups' and delete experimental groups
    old_groups = list(cell_count.keys())
    n_cells_per_stimuli_combination['all groups'] = np.sum(np.vstack(list(n_cells_per_stimuli_combination.values())), axis=0)
    [n_cells_per_stimuli_combination.pop(k) for k in old_groups]
    cell_count['all groups'] = np.nansum(np.dstack(list(cell_count.values())), axis=2)
    [cell_count.pop(k) for k in old_groups]
    groups_to_analyze = ['all groups']

else:
    groups_to_analyze = list(cell_count.keys())
    for this_group in groups_to_analyze:
        cell_count[this_group] = cell_count[this_group][:, np.where(np.logical_not(np.all(np.isnan(cell_count[this_group]), axis=0)))[0]]


for i_group, this_group in enumerate(groups_to_analyze):
    this_group = groups_to_analyze[i_group] #'all groups' # 'sham'
    # Convert data for plotting function
    data = pd.DataFrame(columns=['value'] + STIMULI)
    for ii in range(cell_count[this_group].shape[1]):
        new_data = pd.DataFrame(STIMULI_table, columns=STIMULI)
        new_data['value'] = cell_count[this_group][:, ii]
        data = pd.concat((data, new_data), ignore_index=True, sort=False)
    # Remove 0s
    data.loc[np.where(data['value'] == 0)[0], 'value'] = np.nan
    # Convert selectivity for stimulus to boolean index
    for stimulus_name in STIMULI:
        data[stimulus_name] = data[stimulus_name].map(np.bool)

    # Convert count of cells to percentages
    data.set_index(STIMULI, inplace=True)
    data.dropna(inplace=True)
    data = data['value'].groupby(level = list(range(n_STIMULI))).sum()
    data = data.reset_index()
    for ii in range(data.shape[0]):
        idx = np.all(STIMULI_table == data.loc[ii, STIMULI].values, axis=1).nonzero()[0][0]
        n = n_cells_per_stimuli_combination[this_group][idx]
        data.loc[ii, 'value'] /= n


    # Get number of non-selective cells
    data_to_show = pd.DataFrame(data['value'].values, columns=['value'])
    group_names = list()
    for i_stim in range(n_stimuli):
        group_name = stimuli_groups[i_stim][0]
        group_names.append(group_name)
        stimuli_to_group = stimuli_groups[i_stim][1]
        try:
            if len(stimuli_to_group) == 1:
                data_to_show[group_name] = data[stimuli_to_group]
            else:
                data_to_show[group_name] = np.any(data[stimuli_to_group].values == True, axis=1)

        except:
            print()
    data = data_to_show.set_index(group_names)

    data_to_plot = data.value.groupby(level = list(range(n_stimuli))).sum()
    data_to_plot = data_to_plot / data_to_plot.sum() * 100

    data_to_plot2 = data_to_plot.reset_index()
    print(data_to_plot2)
    columns = list(data_to_plot2.columns)
    columns.remove('value')

    idx = np.where(np.sum(data_to_plot2[columns].values.astype(np.int8), axis=1) != 1)[0]
    data_to_plot2.loc[idx, 'value'] = 0
    data_to_plot2 = data_to_plot2.set_index(columns).squeeze()

    # idx = np.where(np.sum(data_to_plot2[columns].values.astype(np.int8), axis=1) <= 1)[0]
    # data_to_plot2 = data_to_plot2.loc[idx].set_index(columns).squeeze()

    # Generate plot
    axes = UpSet(data_to_plot, facecolor=[i[2] for i in stimuli_groups], linecolor=(.5, .5, .5)).plot()

    # Add labels
    axes['intersections'].set_ylabel('Percent of total cells')
    axes['totals'].set_xlabel('Percent of selective cells')
    title_str = 'session %s' % (', '.join([str(i) for i in timepoints_of_interest]))
    if groups_to_analyze[0] == 'all groups':
        title_str = 'Average ' + title_str
    else:
        title_str = this_group + ' ' + title_str
    axes['fig'].suptitle(title_str)
    # Set tick positions
    y_step = 5
    x_step = 10
    ylims = list(axes['intersections'].get_ylim())
    ylims[1] = np.ceil(ylims[1] / y_step) * y_step
    axes['intersections'].set_yticks(np.arange(0, 101, y_step))
    axes['intersections'].set_ylim(ylims)
    xlims = list(axes['totals'].get_xlim())
    xlims[0] = np.ceil(xlims[0] / x_step) * x_step
    axes['totals'].set_xticks(np.arange(0, 101, x_step))
    axes['totals'].set_xlim(xlims)

    # Fix figure size
    axes['fig'].set_size_inches([40, 9])
    # Save figure to disk
    # plt.savefig(fname=output_filename, bbox_inches='tight', format='pdf')
    # plt.close(axes['fig'])

plt.show()
