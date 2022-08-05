# Imports
import sys
from collections import OrderedDict
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

from Utilities import matlab_file
from Utilities.statistics import confidence_interval
from Utilities.visualization import adjust_spines


def make_graph(info_filename):
    # Load plotting instructions
    PARAMETERS = matlab_file.load(info_filename, variable_names='PARAMETERS')
    # Unpack values
    DATA = PARAMETERS['data']
    GROUPS = PARAMETERS['groups']
    group_colors = OrderedDict()
    for i in PARAMETERS['group_colors']:
        group_colors[i[0]] = i[1]
    group_names = list(group_colors.keys())
    n_groups = len(group_names)
    stimuli_names = PARAMETERS['stimuli']
    labels = PARAMETERS['labels']
    n_timepoints = len(labels)
    output_filename = PARAMETERS['output_filename']
    title = PARAMETERS['title']
    scatterplot_size = PARAMETERS.pop('scatterplot_size', 60)
    alpha = PARAMETERS.pop('alpha', .7)
    add_n = PARAMETERS.pop('add_info', False)
    add_average = PARAMETERS.pop('add_average', False)
    # Set other parameters
    identity_line_color = (0, 0, 0)
    max_n_plots_per_column = 3
    n_rows_subplots = int(np.ceil(n_timepoints / float(max_n_plots_per_column)))
    n_columns_subplots = min(max_n_plots_per_column, n_timepoints)
    subplot_pos = np.unravel_index(np.arange(n_timepoints), (n_rows_subplots, n_columns_subplots))

    # Get data
    data = np.vstack(DATA).ravel()
    # Find data limits
    data_limits = np.array([np.nanmin(data), np.nanmax(data)])
    # Round data
    data_step = 0.1
    data_limits = [np.floor(data_limits[0]/data_step)*data_step, np.ceil(data_limits[1]/data_step)*data_step]
    data_ticks = np.arange(data_limits[0], data_limits[1]+data_step, data_step)

    # Open blank figure
    sns.set_style("white")
    fig = plt.figure(facecolor="w")
    # Set size
    fig.set_size_inches([12.0, 5.0 * n_rows_subplots])
    ax = list()

    for ii in range(n_timepoints):
        ax.append(plt.subplot2grid((n_rows_subplots, n_columns_subplots), (subplot_pos[0][ii], subplot_pos[1][ii])))

        # Get data and make color sequence
        data = DATA[ii]
        groups = GROUPS[ii]

        # Plot data
        data_for_average = list()
        for group_idx, group in enumerate(group_names):
            rows = np.where(groups == group_idx + 1)[0]
            if rows.size == 0:
                continue
            this_color = list(group_colors[group])
            ax[-1].scatter(data[rows, 0], data[rows, 1], s=scatterplot_size, c=this_color, marker='o', linewidths=.1, edgecolors='w', label=group, zorder=n_groups-group_idx, alpha=alpha)
            if add_average:
                data_for_average.append(data[rows, :])

        # Add average across all groups
        if add_average:
            data_for_average = np.vstack(data_for_average)
            mean = np.mean(data_for_average, axis=0)
            ci = confidence_interval(data_for_average, confidence=0.99, parametric=True)
            ax[-1].scatter(mean[0], mean[1], s=100, c='k', marker='o', linewidths=1, edgecolors='w', label='mean', zorder=n_groups+2)
            ax[-1].plot([mean[0]]*2, [ci[0][1], ci[1][1]], color='k', lw=1.5, zorder=n_groups+1)
            ax[-1].plot([ci[0][0], ci[1][0]], [mean[1]]*2, color='k', lw=1.5, zorder=n_groups + 1)

        # Add identity lines
        limits = np.vstack([ax[-1].get_xlim(), ax[-1].get_ylim()])
        ax[-1].plot([0, 1], [.5, .5], color=identity_line_color, lw=1, zorder=0)
        ax[-1].plot([.5, .5], [0, 1], color=identity_line_color, lw=1, zorder=0)
        ax[-1].plot([0, 1], [0, 1], color=identity_line_color, lw=1, zorder=0)
        # Make axis square
        plt.axis('square')
        # Restore limits
        limits = [np.min(limits[:, 0]), np.max(limits[:, 1])]
        ax[-1].set_xlim(limits)
        ax[-1].set_ylim(limits)

        # Add title
        axes_title = labels[ii]
        if add_n:
            axes_title += ' ($\itn$=%i)' % data.shape[0]
        ax[-1].set_title(axes_title, fontsize=14)

    [iax.set_ylim(data_limits[0], data_limits[1]) for iax in ax]
    [iax.set_xlim(data_limits[0], data_limits[1]) for iax in ax]
    # Draw ticks and labels
    [iax.set_xticks(data_ticks) for iax in ax]
    [iax.set_yticks(data_ticks) for iax in ax]
    [iax.set_xlabel(stimuli_names[0], fontsize=12) for iax in ax]
    [iax.set_ylabel(stimuli_names[1], fontsize=12) for iax in ax]
    fig.suptitle(title, fontsize=16)
    # Adjust spines
    [adjust_spines(iax, ["left", "bottom"]) for iax in ax]
    # Adjust subplot spacing
    fig.tight_layout()

    # Make legend
    pos0 = ax[0].get_position()
    pos1 = ax[-1].get_position()
    ax_bg = fig.add_axes((pos0.x0, pos0.y0, pos1.x1-pos0.x0, pos0.y1-pos0.y0), zorder=1)
    ax_bg.set_axis_off()
    handles = list()
    for group_idx, group in enumerate(group_names):
        handles.append(ax_bg.scatter(0, 0, s=scatterplot_size, c=group_colors[group], marker='o', alpha=alpha, linewidths=1, edgecolors='w', zorder=n_groups-group_idx))
    fig.legend(handles, group_names, bbox_to_anchor=[.5, .01], loc=8, borderaxespad=0, ncol=n_groups, fontsize=12, handlelength=.5)
    [h.set_visible(False) for h in handles]

    # Save figure to disk
    plt.savefig(fname=output_filename, bbox_inches="tight", format=output_filename.split(".")[-1])
    # Close figure
    plt.close(fig)


if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        info_filename = sys.argv[1]
    else:
        info_filename = r'D:\_MATLAB_2PI\scatterplot_data.mat'

    # Make graph
    make_graph(info_filename)
