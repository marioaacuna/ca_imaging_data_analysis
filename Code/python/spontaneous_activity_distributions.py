## Imports
import sys
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np

from Utilities import matlab_file
from Utilities.statistics import confidence_interval
from Utilities.visualization import adjust_spines
from Utilities.IO_operations import log


def make_graph(info_filename, data_filename, output_filename):
    # Initialize some variables
    x_limits = None
    x_label = None
    default_colors = None
    plot_handles = None

    # Load plotting instructions
    GRAPH_PARAMETERS = matlab_file.load(info_filename, variable_names='GRAPH_PARAMETERS')
    # Unpack values
    experimental_group_colors = GRAPH_PARAMETERS['group_colors']
    experimental_group_colors = {i: j for i, j in zip(experimental_group_colors[:, 0], experimental_group_colors[:, 1])}
    font_size = GRAPH_PARAMETERS['font_size']
    discard_unstable_baseline = bool(GRAPH_PARAMETERS['discard_unstable_baseline'])
    discard_unchanged_after = bool(GRAPH_PARAMETERS['discard_unchanged_after'])
    overlap_timepoints = GRAPH_PARAMETERS['overlap_timepoints'] == 1

    # Get list of alternative colors
    if overlap_timepoints:
        default_colors = [[11, 198, 255],
                          [255, 170, 0],
                          [255, 0, 4],
                          [0, 170, 0],
                          [0, 0, 255],
                          [255, 85, 255],
                          [190, 92, 35],
                          [255, 6, 114]]
        default_colors = (np.vstack(default_colors) / 255).tolist()

    # Load data
    DATA = matlab_file.load(data_filename, variable_names='DATA')
    attribute_names = sorted(list(DATA['attributes'].keys()))
    n_attributes = len(attribute_names)
    bins = {ii: DATA['attributes'][ii] for ii in attribute_names}
    if 'interval' in list(bins.keys()):
        bins['interval'] = bins['interval'] / 8.49
    experimental_groups = DATA['allowed_groups']
    if isinstance(experimental_groups, str):
        experimental_groups = [experimental_groups]
    n_experimental_groups = len(experimental_groups)
    DATA.pop('attributes', None)
    DATA.pop('allowed_groups', None)
    rows2remove = np.where(np.logical_not(np.in1d(DATA['experimental_group'], experimental_groups)))[0]
    if rows2remove.size > 0:
        DATA = {i: np.delete(DATA[i], rows2remove, axis=0) for i in list(DATA.keys())}

    # Make sure that some variables have the correct shape
    if DATA['change_after'].ndim == 1:
        DATA['change_after'] = np.atleast_2d(DATA['change_after']).transpose()

    if discard_unstable_baseline:
        n = DATA['stable_baseline'].shape[0]
        rows2remove = np.where(DATA['stable_baseline'] == 0)[0]
        if rows2remove.size > 0:
            DATA = {i: np.delete(DATA[i], rows2remove, axis=0) for i in list(DATA.keys())}
            log('%d/%d cells removed because of unstable baseline' % (rows2remove.shape[0], n))

    if discard_unchanged_after:
        n = DATA['change_after'].shape[0]
        rows2remove = np.where(np.any(DATA['change_after'] == 0, axis=1))[0]
        if rows2remove.size > 0:
            DATA = {i: np.delete(DATA[i], rows2remove, axis=0) for i in list(DATA.keys())}
            log('%d/%d cells removed because of unchanged attributes after baseline' % (rows2remove.shape[0], n))

    # Check whether there is a session before
    all_column_names = list(DATA.keys())
    # Get number of timepoints
    # n_timepoints_after = DATA['change_after'].shape[1]
    str_timepoints = list()
    title_timepoints = list()
    if attribute_names[0] + '_dist_before' in all_column_names:
        data = DATA[attribute_names[0] + '_dist_before']
        if not np.all(np.isnan(data.ravel())):
            str_timepoints += ['before']
            title_timepoints += ['before']
    max_n_timepoints_after = DATA['change_after'].shape[1]
    ordinal = lambda n: '%d%s' % (n, {1: 'st', 2: 'nd', 3: 'rd'}.get(n if n < 20 else n % 10, 'th'))
    for it in range(max_n_timepoints_after):
        if '%s_dist_after_%i' % (attribute_names[0], it + 1) in all_column_names:
            data = DATA[attribute_names[0] + '_dist_after_1']
            if not np.all(np.isnan(data.ravel())):
                str_timepoints += ['after_%i' % (it + 1)]
                title_timepoints += ['%s session after' % ordinal(it + 1)]
    n_timepoints = len(str_timepoints)

    # Open blank figure
    sns.set_style('white')
    fig = plt.figure(facecolor='w')
    # Set size
    fig_height = 13
    if overlap_timepoints:
        fig_width = 6
    else:
        fig_width = n_timepoints * 4
    fig.set_size_inches(fig_width, fig_height)
    legend_shown = False  # Initialize as 'no legend present'

    for irow in range(n_attributes):
        # Set specific parameters for each type of data
        current_attribute = attribute_names[irow]
        x_data = bins[current_attribute]
        if current_attribute == 'amplitude':
            x_limits = [0.01, 1]
            x_label = 'Amplitude ($\mathregular{log_{10}(\Delta F/F_0)}$)'
        elif current_attribute == 'interval':
            x_limits = [0.1, 100]
            x_label = 'Time between consecutive $\mathregular{Ca^{2+}}$ transients ($\mathregular{log_{10}}$(s))'
        if overlap_timepoints:
            ax_idx_offset = irow
        else:
            ax_idx_offset = irow * n_timepoints
        # Make names of keys specific to these data
        key_names = [current_attribute + '_dist_' + timepoint for timepoint in str_timepoints] + ['grouping_index']
        data = {key_name: DATA[key_name] for key_name in key_names}

        # Initialize axes list
        if overlap_timepoints:
            plot_handles = list()
        ax = list()
        for icol in range(n_timepoints):
            # Make axes
            if overlap_timepoints:
                if icol == 0:
                    ax.append(fig.add_subplot(2, 1, 1 + ax_idx_offset))
            else:
                ax.append(fig.add_subplot(2, n_timepoints, icol + 1 + ax_idx_offset))
            # Find column with data to show
            col2take = current_attribute + '_dist_' + str_timepoints[icol]

            # Initialize variable for graphic handles
            if not overlap_timepoints:
                plot_handles = list()
            for ii in range(n_experimental_groups):
                current_group = experimental_groups[ii]
                # Get distributions
                rows2take = np.where(data['grouping_index'] == ii+1)[0]
                y = data[col2take][rows2take, :]
                if y.size == 0 or np.all(np.isnan(y.ravel())):
                    continue
                # Compute mean and 95CI
                mean_y = np.nanmean(y, axis=0)
                err_y = confidence_interval(y, confidence=0.99)
                # Get the color
                if overlap_timepoints:
                    color = default_colors[icol]
                else:
                    color = experimental_group_colors[current_group]

                # Plot data
                ax[-1].fill_between(x=x_data, y1=err_y[0], y2=err_y[1], facecolor=color, alpha=0.2)
                h = ax[-1].plot(x_data, mean_y, color=color, zorder=10)[0]
                plot_handles.append(h)  # Store graphic handle

            if not legend_shown and not overlap_timepoints:
                ax[-1].legend(plot_handles, experimental_groups, loc=4, fontsize=font_size, frameon=True, framealpha=1.0, handlelength=1.0, handletextpad=0.5)
                legend_shown = True  # Switch variable to True so this code condition won't be met again

            # Add title to topmost row
            if irow == 0 and not overlap_timepoints:
                ax[-1].set_title(title_timepoints[icol], fontsize=font_size, fontweight='bold')

        # Add legend if timepoints are overlapped
        if not legend_shown and overlap_timepoints:
            ax[-1].legend(plot_handles, title_timepoints, loc=4, fontsize=font_size, frameon=True, framealpha=1.0, handlelength=1.0, handletextpad=0.5)
            legend_shown = True  # Switch variable to True so this code condition won't be met again
            ax[-1].set_title(experimental_groups[0], fontsize=font_size, fontweight='bold')

        # Increase range of values shown by transforming the parameter scale
        [iax.set_xscale('log') for iax in ax]
        # Adjust axes limits
        [iax.set_xlim(x_limits[0], x_limits[1]) for iax in ax]
        [iax.set_ylim(0, 1) for iax in ax]
        # Set x-label
        middle_panel_index = int(np.floor(n_timepoints / 2))
        ax[middle_panel_index].set_xlabel(x_label, fontsize=font_size)
        # Remove boxes all around
        ax[0].set_ylabel('Cumulative probability', fontsize=font_size)
        # Adjust spines
        [adjust_spines(iax, ['bottom','left'], offset=5) for iax in ax]
        [iax.tick_params(direction='out', length=3, width=1) for iax in ax]
        # Hide spines of all axes except the first one
        [iax.spines['left'].set_color('None') for iax in ax[1:]]
        [iax.tick_params(axis='y', colors='None') for iax in ax[1:]]
        # Add grid
        [iax.grid(True, which='major', axis='both') for iax in ax]

    # Fix axes appearance
    plt.tight_layout(rect=[0, 0, 1, 1])
    plt.subplots_adjust(wspace=0.25, hspace=0.25)

    # Save figure to disk
    plt.savefig(fname=output_filename, bbox_inches='tight', format=output_filename.split('.')[-1])
    # Close figure
    plt.close(fig)


if __name__ == '__main__':
    # Get user inputs
    if len(sys.argv) > 1:
        info_filename = sys.argv[1]
        data_filename = sys.argv[2]
        output_filename = sys.argv[3]
    else:
        info_filename = r'D:\_MATLAB_2PI\plotting_instructions_spontaneous_activity.mat'
        data_filename = r'D:\_MATLAB_2PI\spontaneous_activity.mat'
        output_filename = r'D:\_MATLAB_2PI\figure.pdf'
    # Make graph
    make_graph(info_filename, data_filename, output_filename)
