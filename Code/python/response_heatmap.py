# System packages
import sys

# Numerical packages
import numpy as np

# Graphical packages
from matplotlib import pyplot as plt
import seaborn as sns

# Local repository
from Utilities import matlab_file
from Utilities.visualization import adjust_spines


def draw_response_heatmap(info_filename):
    """
    """
    # Set default parameters
    normalize_colormap = 'session_1'  # 'none', 'all', 'cell', 'session_1', 'log'
    threshold = 0
    split_colors = True

    # Initialize variables
    im = None
    im_pos = None
    im_neg = None
    im_neuter = None

    # Load data and unpack inputs
    PARAMETERS = matlab_file.load(info_filename, variable_names='PARAMETERS')
    data = PARAMETERS['data']
    selectivity = PARAMETERS['selectivity']
    is_selective = PARAMETERS['is_selective']
    labels = PARAMETERS['labels']
    title = PARAMETERS['title']
    output_filename = PARAMETERS['output_filename']
    stimulus_timestamps = PARAMETERS['stimulus_timestamps']
    time_axis = PARAMETERS['time_axis']
    font_size = PARAMETERS['font_size']

    # Fix type of some variables
    if not isinstance(stimulus_timestamps, (list, np.ndarray)):
        stimulus_timestamps = [stimulus_timestamps]

    # Get size and range of data to plot
    n_cells, trial_duration, n_sessions = data.shape
    if normalize_colormap == 'log':
        log_data = np.log2(data)
        log_data[np.logical_not(np.isfinite(log_data))] = np.nan
        data_range = [np.nanmin(log_data), np.nanmax(log_data)]
    elif normalize_colormap == 'session_1':
        data_range = [data[:, :, 0].min(), data[:, :, 0].max()]
    else:
        data_range = [data.min(), data.max()]

    # Open blank figure
    sns.set_style("white")
    fig = plt.figure(facecolor="w")
    # Set size
    fig.set_size_inches([5.0 * n_sessions, 7.0])
    ax = list()

    for i_sess in range(n_sessions):
        this_data = data[:, :, i_sess].copy()

        if normalize_colormap == 'cell':
            this_data /= np.max(this_data, axis=0)
        elif normalize_colormap == 'log':
            this_data = np.log2(this_data)
            this_data[np.logical_not(np.isfinite(this_data))] = np.nan

        if threshold is not None:
            if normalize_colormap != 'log':
                this_data[this_data <= threshold] = np.nan
            else:
                this_data[this_data <= np.log2(threshold)] = np.nan

        # Prepare data and plot heatmap
        ax.append(fig.add_subplot(1, n_sessions, i_sess + 1))
        if split_colors:
            data_selective_pos = np.ma.masked_where(np.tile(np.logical_not(np.logical_and(is_selective[:, i_sess] == 1, selectivity[:, i_sess] > 0.5)).reshape(-1, 1), (1, trial_duration)), this_data, copy=True)
            data_selective_neg = np.ma.masked_where(np.tile(np.logical_not(np.logical_and(is_selective[:, i_sess] == 1, selectivity[:, i_sess] < 0.5)).reshape(-1, 1), (1, trial_duration)), this_data, copy=True)
            data_non_selective = np.ma.masked_where(np.tile(np.logical_not(is_selective[:, i_sess] == 0).reshape(-1, 1), (1, trial_duration)), this_data, copy=True)

            im_pos = ax[-1].pcolormesh(time_axis-stimulus_timestamps[0], np.arange(1, n_cells + 1), data_selective_pos, edgecolors='none', cmap='Reds')
            im_neg = ax[-1].pcolormesh(time_axis-stimulus_timestamps[0], np.arange(1, n_cells + 1), data_selective_neg, edgecolors='none', cmap='Blues')
            im_neuter = ax[-1].pcolormesh(time_axis-stimulus_timestamps[0], np.arange(1, n_cells + 1), data_non_selective, edgecolors='none', cmap='gray_r')
            ax[-1].invert_yaxis()

        else:
            im = ax[-1].imshow(this_data, aspect='auto', extent=[time_axis[0] - stimulus_timestamps[0], time_axis[-1] - stimulus_timestamps[0], n_cells, 1], interpolation='none', origin='upper', cmap='gray_r')

        # Normalize colors
        if normalize_colormap in ['none', 'cell', 'log']:
            pass
        else:
            if split_colors:
                im_pos.set_clim(vmin=data_range[0], vmax=data_range[1])
                im_neg.set_clim(vmin=data_range[0], vmax=data_range[1])
                im_neuter.set_clim(vmin=data_range[0], vmax=data_range[1])
            else:
                im.set_clim(vmin=data_range[0], vmax=data_range[1])

        # Draw timestamps
        for t in stimulus_timestamps:
            ax[-1].axvline(t - stimulus_timestamps[0], color='k', lw=1)

        # Adjust labels and axes appearance
        if i_sess == 0:
            adjust_spines(ax[-1], ['bottom', 'left'], smart_bounds=True)
            ax[-1].set_ylabel(r'Cell $\it{no.}$', fontsize=font_size-2)
            y_ticks = np.unique(np.hstack(([1], np.arange(0, n_cells, 10), [n_cells])))
            ax[-1].set_yticks(y_ticks)
        else:
            adjust_spines(ax[-1], 'bottom', smart_bounds=True)
        ax[-1].set_xlabel('Time (s)', fontsize=font_size-2)
        ax[-1].set_title(labels[i_sess], fontsize=font_size)

    # Add title
    fig.suptitle(title, fontsize=font_size)
    # Adjust subplot spacing
    fig.tight_layout(rect=[0, 0, 1, .9])

    # Save figure to disk
    if output_filename != '':
        plt.savefig(fname=output_filename, bbox_inches='tight', format=output_filename.split('.')[-1])
        # Close figure
        plt.close(fig)
    else:
        # Return the figure handle so it can be printed, for example.
        return fig


################################################################################
# Direct call
################################################################################
if __name__ == '__main__':
    # Get user inputs
    if len(sys.argv) > 1:
        info_filename = sys.argv[1]
    else:
        info_filename = r'D:\_MATLAB_2PI\heatmap_plot_data.mat'

    # Make graph
    draw_response_heatmap(info_filename)

