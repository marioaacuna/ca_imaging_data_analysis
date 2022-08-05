# System packages
import sys, os
from collections import OrderedDict

# Numerical packages
import numpy as np
import pandas as pd

# Graphical packages
from matplotlib import pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
# Printing figures
from matplotlib.backends.backend_pdf import PdfPages

# Local repository
#from Utilities import matlab_file


def draw_alluvial_plot(data, timepoint_names, delete_hidden=True, group_names=None,
                       split_flows=True, curvature=0.3, time_0=0,
                       draw_blocks=True, write_block_name='', colors_only_in_blocks=False, blocks_gap=0.05, block_width=0.1, add_info=None,
                       blocks_border=4, other_blocks_border=0,
                       label_offset=0.05, title=None, title_font=20, label_font=14, groups_font=16, output_filename='',
                       add_summary_pie=True):
    """
    # This function draws an alluvial plot.

    ### DATA PARAMETERS ###
    :param data: [pandas DataFrame] It contains flows of data on rows. Each column
        represents a time-point (or more generally, a data dimension). Additional
        columns are:
            'freq': count or frequencies of each flow. This column is mandatory.
            'hide': boolean flag to set whether to hide the corresponding
                flow. The parameter ``delete_hidden`` controls what to do with
                these flows (see below).
            'alpha': Transparency level of the flow, in range [0, 1].
    :param timepoint_names: [list] List of names of time-points. These should be
        columns in the data DataFrame.
    :param delete_hidden: [bool] Whether rows where the 'hide' flag is set to
        True should be simply ignored in the plot or deleted completely from the
        data.
    :param group_names: [list] List whose items are 2-elements lists containing
        the name of the group (as found in 'data') and its color (a tuple in either
        in RGB or RGBA form). The order of the names decides which flows will be
        sorted first and drawn on top of the others. If the second element is None
        or empty, gray color will be used.

    ### FLOW PARAMETERS ###
    :param split_flows: [bool] Whether to add a thin border around each flow to
        separate it from the others.
    :param curvature: [float] This parameter sets the curvature of the flows
        ribbons. The higher the value the more curved they will be.
            0: straight lines.
    :param time_0: [int] Index in ``timepoint_names`` on which
        colors will be aligned. Therefore, at this timepoint, groups will appear
        in one color only, and one will see colored flows spread from this point
        forward and/or backward in time to follow radiation of specific groups in
        time.

    ### GROUP PARAMETERS ###
    :param draw_blocks: [list, bool] Whether to draw a solid rectangles around
        the point where all the flows of the same group depart or arrive.
            bool: Whether to draw blocks for all groups or no group at all.
            sequence: Draw block only for the groups listed. These are listed in
                2-element lists containing the time-point label and the group name,
                as in x- and y-coordinate of the block to draw.
    :param write_block_name: [str] Whether to write the blocks name.
            'all': Write name on all blocks (those drawn and those not).
            'drawn': Write name only on drawn blocks, as set by ``draw_blocks``.
            'time_0': Only on the timepoint used to set colors.
            'none' or '': Do not write any names.
    :param colors_only_in_blocks: [bool] Whether color is used for the entire flow
            or only within blocks.
    :param blocks_gap: [float] Gap between blocks.
    :param block_width: [float] Width of blocks.
    :param add_info: [str] Whether to write the number of items in each flow. If
        None, it won't write anything, otherwise it adds the number of items with
        'n' or their percentage with '%';

    ### FIGURE PARAMETERS ###
    :param blocks_border: [int] Thickness of drawn blocks.
    :param other_blocks_border: [int] Thickness of blocks not drawn.
    :param label_offset: [float] Offset of label from graph.
    :param title: [str] Title to add to top of the plot.
    :param title_font: [int] Font-size of title.
    :param label_font: [int] Font-size of bottom labels.
    :param groups_font: [int] Font-size of group name in blocks.
    :param output_filename: [str] Filename where to save figure.
    """
    # Set default hard-coded properties
    default_transparency = 0.5
    neutral_color = np.array([.7, .7, .7], dtype=float)
    arrow_size = 0.03
    if split_flows:
        flow_edge_color = 'w'
        flow_edge_thickness = 1
    else:
        flow_edge_color = 'None'
        flow_edge_thickness = 0

    # Make sure that data has all the needed columns
    data_columns = list(data.columns)
    if 'freq' not in data_columns:
        raise Exception('data must contain a \'freq\' column')
    if 'alpha' not in data_columns:
        data['alpha'] = default_transparency
    if 'hide' not in data_columns:
        data['hide'] = False

    # Get size of data
    n_flows = data.shape[0]

    # Automatically hide frequencies of 0
    data.loc[data['freq'] == 0, 'hide'] = True

    # Check that there is a color for each group
    if group_names is None:
        group_names = [[i, None] for i in np.unique(data[timepoint_names].values.ravel())]
    for idx, g in enumerate(group_names):
        if g[1] is None or (isinstance(g[1], list) and len(g[1])==0):
            group_names[idx][1] = neutral_color
    # Add colors to each flow
    colors = np.zeros((n_flows, 3), dtype=float)
    for idx in range(n_flows):
        g = data.loc[idx, timepoint_names[time_0]]
        colors[idx, :] = [ii[1] for ii in group_names if ii[0] == g][0][0]
    # Store colors in data
    data['color'] = colors.tolist()

    # Convert data to categories (like R's factors)
    sorted_categories = [g[0] for g in group_names]
    for t in timepoint_names:
        data[t] = pd.Categorical(data[t], ordered=True, categories=sorted_categories)
    # Check that there are no unmatched groups
    if data[timepoint_names].isna().any().any():
        raise Exception('All groups must be included in variable \'group_names\'')
    # Transform colors to categorical values
    sorted_colors = [g[1] for g in group_names]
    data['color_order'] = None
    for c_idx, c in enumerate(sorted_colors):
        color_rows = np.where(np.all(np.equal(np.vstack(data['color']), np.tile(c, (n_flows, 1))), axis=1))[0]
        data.loc[color_rows, 'color_order'] = c_idx
    data['color_order'] = pd.Categorical(data['color_order'], ordered=True, categories=range(len(sorted_colors)))

    # Convert draw_blocks to list or boolean
    if isinstance(draw_blocks, bool):
        if draw_blocks is True:
            # Get names of all groups
            draw_blocks = [[x, y] for x in timepoint_names for y in np.unique(data[timepoint_names].values)]
        else:
            draw_blocks = list()
    elif not isinstance(draw_blocks, list):
        raise Exception('\'draw_blocks\' must be either boolean or a list')

    # Get number of observations kept or discarded
    summary_n_discarded = sum(data.loc[data['hide'], 'freq'])
    summary_n_kept = sum(data['freq']) - summary_n_discarded

    # Delete hidden data
    if delete_hidden:
        data.drop(np.where(data['hide'])[0], inplace=True, axis=0)
        data.reset_index(inplace=True, drop=True)
        # Update number of flows
        n_flows = data.shape[0]

    # Compute endpoints of flows (as polygons)
    freq = data['freq'] / data['freq'].sum()  # Normalize frequencies to have sum 1
    flows_pos = list()
    for t in range(n_time_points):
        # Move current dimension at the front of the list
        times_sorted = list(range(n_time_points))
        times_sorted.remove(t)
        times_sorted.insert(0, t)
        # Re-arrange columns
        column_names = [data_columns[i] for i in times_sorted]
        # Add color as the second grouping hierarchy (that is, within blocks)
        column_names.insert(1, 'color_order')
        # Order of rows of data starting from i-th dimension
        rows_order = data.sort_values(by=column_names, axis=0, ascending=True, inplace=False).index.values
        # Invert rows to plot from top to bottom
        rows_order = rows_order[::-1]

        # Breakpoints on a dimension
        x = np.hstack(([0], np.cumsum(freq[rows_order]))) * (1 - blocks_gap)
        # Stripe coordinates on a dimension
        x = np.vstack((x[:-1], x[1:])).transpose()
        # By how much stripes need to be shifted upwards
        _, as_numeric_d = np.unique(data.loc[rows_order, timepoint_names[t]].values, return_inverse=True)
        gap = np.cumsum(np.hstack(([0], (np.diff(as_numeric_d) != 0).astype(int))))
        # Get max cap and cap it to 1
        mx = gap.max()
        if mx == 0:
            mx = 1
        # Shifts
        gap = gap.astype(float) / mx * blocks_gap
        # Add gap-related shifts to stripe coordinates on dimension i
        x = (x + np.atleast_2d(gap).transpose())[rows_order.argsort(), :]
        # Store final values
        flows_pos.append(x)

    fig = plt.figure(figsize=(27.22, 15.47))
    # Make empty axes
    if add_summary_pie:
        height = 0.9
    else:
        height = 1.0
    ax = fig.add_axes([0, 0, 1, height])
    ax.set_axis_off()

    # Sort data's categories in ascending order at time 0
    sorted_rows = data.sort_values(by=timepoint_names[time_0], axis=0, ascending=True, inplace=False).index.values
    sorted_rows = sorted_rows[::-1]

    # Draw stripes
    for idx, row in enumerate(sorted_rows):
        # Skip flow if it has to be hidden
        if data.loc[row, 'hide']:
            continue

        # Get points
        x_all_top = np.empty((0, ))
        y_all_top = np.empty((0, ))
        y_all_bottom = np.empty((0, ))
        for it in range(n_time_points - 1):
            # Get points of the whole segment
            x = np.array([it, it, it + curvature, it + 1 - curvature, it + 1, it + 1, it + 1 - curvature, it + curvature, it]) + 1 + np.hstack(([block_width] * 3, [-block_width] * 4, [block_width] * 2))
            y = np.hstack((flows_pos[it][row, [0, 1, 1]], flows_pos[it + 1][row, [0, 0, 1, 1]][::-1], flows_pos[it][row, [0, 0]]))
            # Split top and bottom borders
            x_all_top = np.hstack((x_all_top, x[1:5], [it+2]))
            y_all_top = np.hstack((y_all_top, y[1:5], y[4]))
            y0 = np.hstack((y[-4], y[-4:]))[::-1]
            y_all_bottom = np.hstack((y_all_bottom, y0))

            # If only the center block will contain color, dodge it to avoid alpha-compositing of the colors
            if colors_only_in_blocks:
                # Replace last point
                x_all_top = np.hstack((x_all_top[:-1], x_all_top[-2]))
                y_all_top = np.hstack((y_all_top[:-1], y0[-1]))
                y_all_bottom = np.hstack((y_all_bottom[:-1], y0[-1]))
                # Add dodge-point until the last segment is drawn
                if it != n_time_points - 2:
                    x_all_top = np.hstack((x_all_top, x_all_top[-1] + block_width * 2))
                    y_all_top = np.hstack((y_all_top, y0[-1]))
                    y_all_bottom = np.hstack((y_all_bottom, y0[-1]))

        # Concatenate segment
        if colors_only_in_blocks:
            # Replicate first point at both ends of the ribbon
            x = np.hstack((x_all_top[0], x_all_top, x_all_top[::-1], x_all_top[0]))
        else:
            # Add the block width to the last point
            x_all_top[-1] += block_width
            # Extend ribbon
            x = np.hstack(([1 - block_width], x_all_top, x_all_top[::-1], [1-block_width]))
        # Combine all segments in one path
        y = np.hstack((y_all_top[0], y_all_top, y_all_bottom[::-1], y_all_bottom[0]))
        # Combine vertices
        vertices = np.column_stack((x, y)).tolist()

        # Make list of codes for matplotlib to draw spline curves and a straight segments in one go
        codes = list()
        # Top border
        for it in range(n_time_points - 1):
            codes.append([Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO])
            if colors_only_in_blocks:
                codes.append([Path.LINETO])
        if not colors_only_in_blocks:
            codes.append(Path.LINETO)
        # Bottom border (drawn backward)
        for _ in range(n_time_points - 1):
            codes.append([Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO])
            if colors_only_in_blocks:
                codes.append(Path.LINETO)
        if colors_only_in_blocks:
            codes.pop()
        # The first code must be a move-to action
        codes.insert(0, Path.MOVETO)
        # Concatenate all codes in one list
        codes = np.hstack(codes).tolist()

        # Make path object
        path = Path(vertices, codes, closed=True)
        # Get properties of patch to draw
        if colors_only_in_blocks:
            this_color = np.hstack((neutral_color, data.loc[row, 'alpha']))
        else:
            this_color = np.hstack((data.loc[row, 'color'], data.loc[row, 'alpha']))
        # Draw patch
        ax.add_patch(patches.PathPatch(path, facecolor=this_color, edgecolor=flow_edge_color, linewidth=flow_edge_thickness, capstyle='projecting', zorder=idx))

    # Draw category blocks and their label
    for t in range(n_time_points):
        # Get number of unique groups at this time point
        groups = np.unique(data[timepoint_names[t]])
        # Get blocks position
        new_d = pd.DataFrame()
        new_d['id'] = data[timepoint_names[t]]
        new_d['ind'] = np.arange(n_flows)
        grouped = new_d.groupby('id')
        blocks_size = list()
        for g in groups:
            ind = grouped.get_group(g)['ind'].values
            y = flows_pos[t][ind, :].ravel()
            blocks_size.append([y.min(), y.max()])

        # Draw blocks for groups
        for g in range(len(groups)):
            if write_block_name == 'all' or t == time_0:
                do_write_block_name = True
            else:
                do_write_block_name = False
            # Get coordinates of this block
            x = [t + 1 - block_width, t + 1 + block_width]
            y = blocks_size[g]
            # Check whether this block has to be drawn
            if any([i[0]==timepoint_names[t] and i[1]==groups[g] for i in draw_blocks]):
                border_thickness = blocks_border
                # Add name to block
                if write_block_name == 'drawn':
                    do_write_block_name = True
            else:
                border_thickness = other_blocks_border
            if border_thickness > 0:
                block_border_color = [group[1] for group in group_names if group[0]==groups[g]][0][0]
                # Draw bottom line
                ax.plot(x, [y[0]] * 2, color=block_border_color, linewidth=border_thickness, solid_capstyle='projecting', zorder=500)
                # Draw top line
                ax.plot(x, [y[1]] * 2, color=block_border_color, linewidth=border_thickness, solid_capstyle='projecting', zorder=500)
                # Draw left line
                ax.plot([x[0]] * 2, y, color=block_border_color, linewidth=border_thickness, solid_capstyle='projecting', zorder=500)
                # Draw right line
                ax.plot([x[1]] * 2, y, color=block_border_color, linewidth=border_thickness, solid_capstyle='projecting', zorder=500)

            # Add name of block
            if do_write_block_name:
                text_str = groups[g]
                text_str = text_str.replace('+', '$^+$')
                text_str = text_str.replace('-', '$^-$')
                text_str = text_str.replace(' ', '\n')
                ax.text(t + 1, np.mean(blocks_size[g]), text_str, horizontalalignment='center', verticalalignment='center', fontsize=groups_font, zorder=1000)

        # Color block inset
        if colors_only_in_blocks:
            x = [t + 1 - block_width, t + 1 + block_width]
            # Loop through segments
            for row in range(n_flows):
                # Get vertical range of this flow
                y = flows_pos[t][row, :]
                # Draw patch
                this_color = np.hstack((data.loc[row, 'color'], data.loc[row, 'alpha']))
                ax.add_patch(patches.Rectangle(xy=(x[0], y[0]), width=2*block_width, height=y[1]-y[0], facecolor=this_color, edgecolor=flow_edge_color, linewidth=flow_edge_thickness, capstyle='projecting', zorder=400))

    # Add number of elements in each group
    if add_info is not None:
        # Get x position
        x = 1 - block_width - block_width * label_offset * 2 + time_0
        for row in range(n_flows):
            # Only if low is shown
            if data.loc[row, 'hide']:
                continue
            # Get count or percentage
            count = data.loc[row, 'freq']
            if add_info == '%':
                value_to_match = data.loc[row, timepoint_names[time_0]]
                total_count = data.loc[data[timepoint_names[time_0]] == value_to_match, 'freq'].sum()
                count = float(count) / total_count * 100
                count = '%i%%' % count
            else:
                count = '%i' % count
            # Get y position
            y = flows_pos[time_0][row, :].mean()
            ax.text(x, y, count, horizontalalignment='right', verticalalignment='center', fontsize=label_font-2)

    # Draw labels at each time-point
    pos = flows_pos[time_0].ravel()
    min_y = float(pos.min())
    range_y = float(pos.max() - pos.min())
    for it in range(n_time_points):
        if timepoint_names[it] != '':
            ax.text(it + 1, min_y - range_y * label_offset, timepoint_names[it], horizontalalignment='center', verticalalignment='top', fontsize=label_font)
    # Add title
    if title is not None:
        title += ' ($\itn$=%i)' % (summary_n_kept + summary_n_discarded)
        fig.text(0.01, .99, title, horizontalalignment='left', verticalalignment='top', fontsize=title_font)
    # Add arrow pointing at time used to determine colors
    ax.arrow(x=time_0+1, y=1.01+arrow_size, dx=0, dy=-arrow_size, fc='k', ec='k', head_width=arrow_size, head_length=arrow_size, length_includes_head=True, zorder=10000)

    # Scale limits to contain everything
    ax.relim()
    ax.autoscale()
    # Adjust axes to contain labels and title
    y_lims = list(ax.get_ylim())
    if any([l != '' for l in timepoint_names]):
        y_lims[0] -= (y_lims[1]-y_lims[0]) * label_offset
    ax.set_ylim(y_lims[0], y_lims[1])

    # Draw pie-chart
    if add_summary_pie:
        ax_pie = fig.add_axes([.9, .9, .1, .1])
        ax_pie.set_axis_off()
        colors = [(0, 0, 0), neutral_color]
        h_patches, _ = ax_pie.pie([summary_n_kept, summary_n_discarded], radius=1, colors=colors, startangle=90, counterclock=False, wedgeprops={'edgecolor':'None', 'width':1})
        ax_pie.set_aspect('equal', 'box')
        ax_pie.relim()
        ax_pie.autoscale()
        ax_pie.legend(h_patches, ['modulated (%i)' % summary_n_kept, 'unmodulated (%i)' % summary_n_discarded], loc='center right', frameon=False, bbox_to_anchor=(-.05, 0.5), fontsize=label_font)

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
        info_filename = r'D:\_MATLAB_2PI\alluvial_plot_data.mat'

    # Load data and unpack inputs
    PARAMETERS = matlab_file.load(info_filename, variable_names='PARAMETERS')
    delete_hidden = PARAMETERS.pop('delete_hidden', True)
    add_info = PARAMETERS.pop('add_info', None)
    split_flows = PARAMETERS.pop('split_flows', True)
    do_increase_visibility_of_stable_flows = PARAMETERS.pop('do_increase_visibility_of_stable_flows', True)
    add_summary_pie = PARAMETERS.pop('add_summary_pie', True)

    # Read parameters that are common to all plots
    timepoint_names = PARAMETERS['labels']
    title = PARAMETERS['title']
    # Load colors of each group and convert them to a dictionary
    group_colors = OrderedDict()
    for g in PARAMETERS['group_colors'].tolist():
        group_colors[g[0]] = g[1]

    # Make sure that output folder exists
    output_dir = os.path.split(PARAMETERS['output_filename'])[0]
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    # Create pdf file
    PDF_file = PdfPages(PARAMETERS['output_filename'])

    # Get number of time-steps
    n_time_points = PARAMETERS['data'].shape[1]
    # Draw an alluvial plot for each time step
    for time_step in range(n_time_points):
        # Assemble data DataFrame
        data = pd.DataFrame(PARAMETERS['data'], columns=timepoint_names)
        data['freq'] = PARAMETERS['freq']
        # Hide non-selective cells throughout all sessions
        data['hide'] = False
        if delete_hidden:
            data.loc[np.where(np.all(data[timepoint_names].values == 'none', axis=1))[0], 'hide'] = True

        # Get all groups
        all_groups = np.unique(data[timepoint_names].values.ravel())
        n_groups = all_groups.shape[0]
        # Set colors according to the selectivity at the current time-step
        group_colors_keys = list(group_colors.keys())
        # Add single selectivity groups
        group_names = pd.DataFrame(columns=['name','color','order'])
        for group in all_groups:
            group_idx = np.array([i in group for i in group_colors_keys], dtype=int)
            if group_idx.sum() > 1:  # more than one match
                order = group_colors_keys.index('mixed')
                this_color = group_colors['mixed']
            else:
                order = np.where(group_idx == 1)[0][0]
                this_color = group_colors[group_colors_keys[order]]
            # Store name and color
            last_row = group_names.shape[0]
            group_names.at[last_row, 'name'] = group
            group_names.at[last_row, 'color'] = [this_color]
            group_names.at[last_row, 'order'] = order
        # Sort rows
        group_names.sort_values(by=['order','name'], inplace=True)
        # Reformat group_names variable in 2-element lists
        group_names = [[ig[1]['name'], ig[1]['color']] for ig in group_names.iterrows()]

        # Add alpha level to each flow
        data['alpha'] = .5
        if do_increase_visibility_of_stable_flows:
            # Increase visibility of cells that never change selectivity across all sessions
            stable_idx = np.where(np.all(np.equal(data[timepoint_names].values[:, 1:], np.tile(data[timepoint_names].values[:, 0], (n_time_points-1, 1)).transpose()), axis=1))[0]
            data.loc[stable_idx, 'alpha'] = 1
        # Reduce visibility of cells not selective at the current time step
        data.loc[data[timepoint_names[time_step]] == 'none', 'alpha'] = .15

        # Draw only blocks at time 0
        group_names_at_t0 = np.unique(data[timepoint_names[time_step]].values)
        blocks_to_draw = [[timepoint_names[time_step], y] for y in group_names_at_t0]

        # Make graph
        fig = draw_alluvial_plot(data, timepoint_names=timepoint_names, delete_hidden=delete_hidden, group_names=group_names,
                                 curvature=0.3, time_0=time_step, split_flows=split_flows,
                                 draw_blocks=blocks_to_draw, write_block_name='all',
                                 colors_only_in_blocks=True, blocks_gap=0.1, block_width=0.1, add_info=add_info,
                                 blocks_border=4, other_blocks_border=2,
                                 label_offset=0.05, title=title, title_font=20, label_font=14, groups_font=16,
                                 output_filename='',
                                 add_summary_pie=add_summary_pie)

        # Print file
        PDF_file.savefig(fig, bbox_inches='tight')

    # Flush images to disk
    PDF_file.close()
