## Imports
import sys
from matplotlib import pyplot as plt
import matplotlib.collections as mpl_collections
import matplotlib.colors as mpl_colors
from matplotlib.patches import Rectangle
import pandas as pd
import seaborn as sns
import numpy as np

from Utilities import matlab_file
from Utilities.statistics import confidence_interval
from Utilities.visualization import adjust_spines


def make_graph(info_filename, data_filename, output_filename):
    # Load plotting instructions
    GRAPH_PARAMETERS = matlab_file.load(info_filename, variable_names='GRAPH_PARAMETERS')
    # Unpack values
    surgery_types_colors = GRAPH_PARAMETERS['group_colors']
    surgery_types_colors = {k:v for k,v in surgery_types_colors}
    dist_width = GRAPH_PARAMETERS['dist_width']
    font_size = GRAPH_PARAMETERS['font_size']
    scatterplot_size = GRAPH_PARAMETERS['scatterplot_size']

    # Load data
    DATAPOINTS = matlab_file.load(data_filename, variable_names='DATAPOINTS')

    # Open blank figure
    sns.set_style('white')
    fig = plt.figure(facecolor='w')
    # Set size
    fig.set_size_inches([12.0, 5.0])

    ############################################################################
    # Get data
    data = np.array(DATAPOINTS['n_cells_per_condition']['data'])
    # Find data limits
    data_limits = data[:,1]
    data_limits = np.array([np.nanmin(data_limits), np.nanmax(data_limits)])
    # Round data
    data_step = 25.
    data_limits = [0, np.ceil(data_limits[1]/data_step)*data_step]
    data_ticks = np.arange(data_limits[0], data_limits[1]+data_step, data_step)
    # Make row names
    surgery_names = DATAPOINTS['n_cells_per_condition']['group_names']

    # Make DataFrame with data
    DF = pd.DataFrame(columns=['surgery_type', 'n'])
    DF['surgery_type'] = surgery_names[data[:, 0] - 1]
    DF['n'] = data[:, 1]

    # Open axes
    ax_left = fig.add_subplot(1, 2, 1)
    # Set limits of y-axes
    ax_left.set_ylim(data_limits[0], data_limits[1])

    # Plot data
    ax_left = sns.swarmplot(x='surgery_type', y='n', data=DF, size=scatterplot_size, ax=ax_left, clip_on=False, dodge=True, palette=surgery_types_colors, linewidth=1)
    # Get handle of points
    h = [child for child in ax_left.get_children() if isinstance(child, mpl_collections.PathCollection)]

    # Loop through groups
    dist_centers = list()
    for icol in range(len(surgery_names)):
        # Plot errorbars
        d = h[icol].get_offsets()[:,0]
        if d.shape[0] == 0:
            continue
        center = np.mean(d)
        dist_centers.append(center)

        if d.shape[0] != 1:
            y = DF.loc[np.in1d(DF['surgery_type'], surgery_names[icol]), 'n'].values
            sem = confidence_interval(y)

            ax_left.add_patch(Rectangle((center-dist_width/2, sem[0]), dist_width, sem[1]-sem[0], edgecolor='None', facecolor=[.85, .85, .85]))
            ax_left.plot([center-dist_width/2, center+dist_width/2], [np.mean(y), np.mean(y)], linewidth=2, color='k', zorder=10)

        # Adjust colors
        color = mpl_colors.colorConverter.to_rgba(surgery_types_colors[surgery_names[icol]], alpha=.5)
        h[icol].set_facecolors(color)
        color = mpl_colors.colorConverter.to_rgba(surgery_types_colors[surgery_names[icol]], alpha=.9)
        h[icol].set_edgecolors(color)
        h[icol].set_linewidth(1)

    # Add y-axis to first plot and adjust spines
    ax_left.set_ylabel('Number of ROIs detected per mouse', fontsize=font_size)
    adjust_spines(ax_left, ['left', 'bottom'])
    ax_left.tick_params(direction='out', length=3, width=1)
    # Set ticks of leftmost axis
    ax_left.set_yticks(data_ticks)
    ax_left.set_ylim(data_limits[0], data_limits[1])

    # Remove x-ticks and label, and increase font of xticklabels
    ax_left.set_xticks(range(len(surgery_names)))
    ax_left.tick_params(axis='x', length=0, width=0, labelsize=font_size)
    ax_left.set_xlabel('')
    ax_left.spines['bottom'].set_color('None')

    # Move all axes to top
    ax_left.set_zorder(2)
    ax_left.set_alpha(0.0)
    ax_left.set_facecolor('None')

    ############################################################################
    # Get data
    data = np.array(DATAPOINTS['proportion_active_cells_per_condition']['data'])
    data = np.vstack((data[0,:], np.sum(data[1:,:],axis=0)))
    # Make row names
    group_names = DATAPOINTS['proportion_active_cells_per_condition']['group_names']
    row_names = DATAPOINTS['proportion_active_cells_per_condition']['row_names']
    row_names = row_names[:-1]
    row_names = np.array(row_names.tolist() * len(group_names))

    # Open axes
    ax_right = fig.add_subplot(1, 2, 2)

    ax_right.cla()
    donut_width = 0.3
    # Make outer donut
    colors = [surgery_types_colors[i].tolist() + [.6] for i in group_names]
    data_sum = np.sum(data, axis=0)
    _, hLabels, hPct = ax_right.pie(data_sum, radius=1, labels=group_names, colors=colors, autopct='%i', pctdistance=.85, startangle=90, counterclock=False, wedgeprops={'edgecolor':'None', 'width':donut_width}, textprops={'fontsize':font_size})
    # Change font size of labels
    [ih.set_fontsize(font_size) for ih in hLabels]
    # Change the percent text into raw count
    for ig in range(len(group_names)):
        hPct[ig].set_text('%i' % data_sum[ig])

    # Make inner donut
    data_all = data.T.ravel()
    hWedges, _, hPct = ax_right.pie(data_all, radius=1-donut_width, autopct='%i%%', pctdistance=.77, startangle=90, counterclock=False, wedgeprops={'edgecolor':[1,1,1], 'width':donut_width}, textprops={'fontsize': font_size})
    # Fix color of wedges
    colors = [surgery_types_colors[i].tolist() + [.4] for i in group_names]
    ig = 0
    for isubg in range(len(row_names)):
        # Remove 0% labels
        if hPct[isubg].get_text() == '0%':
            if data_all[isubg] > 0:
                hPct[isubg].set_text('%i' % data_all[isubg])
            else:
                hPct[isubg].set_text('')
        else:
            # Change the percent text into raw count
            hPct[isubg].set_text('%i' % data_all[isubg])

        if row_names[isubg] == 'active':
            hWedges[isubg].set_facecolor(colors[ig])
            ig += 1
        elif row_names[isubg] == 'inactive':
            hWedges[isubg].set_facecolor([0,0,0, .05])

    # Move all axes to top
    ax_right.set_zorder(2)
    ax_right.set_alpha(0.0)
    ax_right.set_facecolor('None')

    # Make title
    n_cells_kept = np.sum(data[0,:])
    n_cells_discarded = np.sum(data[1:,:])
    ax_right.set_title('Cells kept: %i | Cells discarded: %i' % (int(n_cells_kept), int(n_cells_discarded)), fontsize=font_size)
    # Make proportions equal
    ax_right.axis('equal')

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
        info_filename = r'D:\_MATLAB_2PI\plotting_instructions.mat'
        data_filename = r'V:\2p_imaging_pain\6_data\ACC_pain\summary.mat'
        output_filename = r'D:\_MATLAB_2PI\figure.pdf'
    # Make graph
    make_graph(info_filename, data_filename, output_filename)
