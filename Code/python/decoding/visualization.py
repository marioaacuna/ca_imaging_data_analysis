# System packages
import os

# Numerical packages
import numpy as np
from scipy.stats import sem

# Graphical packages
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt

# Local repository
from Utilities.visualization import adjust_spines, center_diverging_colormap, colorbar


def classification_results(PARAMETERS, labels, SAMPLE_IDX, print_figure=True, force_neurons=None):
    n_classes = np.unique(labels).shape[0]
    class_names = list(PARAMETERS['data'].keys())
    n_trials_per_class = [PARAMETERS['data_size'][cl][1] for cl in class_names]

    if force_neurons is not None:
        if not isinstance(force_neurons, (list, np.ndarray)):
            force_neurons = [force_neurons]

    fig, ax = plt.subplots(nrows=3, ncols=n_classes, figsize=(4 * n_classes, 12))
    # class_idx = 1; cl=target_names[class_idx]
    for class_idx, cl in enumerate(class_names):
        # Get indices of "true positive" assignments relative to each class
        TP = np.where((SAMPLE_IDX['class'] == class_idx) & (SAMPLE_IDX['predicted_label'].values.astype(int) == class_idx) & (SAMPLE_IDX['significant_classified'] == 1))[0]
        sample_idx_in_class = SAMPLE_IDX.loc[TP, 'sample_idx_in_class'].values.astype(int)

        # Make heatmap of posterior probabilities of belonging to each class
        pp_heatmap = np.zeros((PARAMETERS['data'][cl].shape[0],), dtype=float) * np.nan
        pp_heatmap[sample_idx_in_class] = SAMPLE_IDX.loc[TP, 'posterior_probability']
        pp_heatmap = pp_heatmap.reshape(n_trials_per_class[class_idx], -1)
        # Estimate reliability of assignment as proportion of trials in which a correct assignment was made
        y = pp_heatmap.copy()
        y[np.isnan(y)] = 0  # Convert to binary
        y[y > 0] = 1
        reliability = np.nanmean(y, axis=0)
        # Get neural activity
        if force_neurons is not None:
            if len(force_neurons) == 1:
                activity = PARAMETERS['data'][cl][:, force_neurons]
            else:
                activity = np.mean(PARAMETERS['data'][cl][:, force_neurons], axis=1)
            activity = activity.reshape(n_trials_per_class[class_idx], -1)
        else:
            activity = np.mean(PARAMETERS['data'][cl], axis=1).reshape(n_trials_per_class[class_idx], -1)

        # Make plots
        # Posterior probability
        ax[0, class_idx].imshow(pp_heatmap, cmap='gray_r', aspect='auto', interpolation=None)
        y_ticklabels = np.unique([1] + np.arange(0, n_trials_per_class[class_idx], 5)[1:].tolist() + [n_trials_per_class[class_idx]])
        y_ticks = y_ticklabels - 1
        y_ticklabels = ['%i' % i for i in y_ticklabels]
        ax[0, class_idx].set_yticks(y_ticks)
        ax[0, class_idx].set_yticklabels(y_ticklabels)
        ax[0, class_idx].set_title(cl)
        # Reliability
        ax[1, class_idx].plot(reliability, '-ok', markersize=3)
        adjust_spines(ax[1, class_idx], ['bottom', 'left'], smart_bounds=False)
        # Neural activity
        ax[2, class_idx].fill_between(x=np.arange(activity.shape[1]), y1=np.mean(activity, axis=0) - sem(activity, axis=0), y2=np.mean(activity, axis=0) + sem(activity, axis=0), facecolor=(0.04, 0.78, 1), alpha=0.5)
        ax[2, class_idx].plot(np.mean(activity, axis=0), '-ok', markersize=3)
        adjust_spines(ax[2, class_idx], ['bottom', 'left'], smart_bounds=False)
        # Set axes labels
        if class_idx == 0:
            ax[0, class_idx].set_ylabel('Trial #')
            ax[1, class_idx].set_ylabel('Reliability')
            ax[2, class_idx].set_ylabel('Activity')
        ax[2, class_idx].set_xlabel('Time bin')
        # Fix x-limits
        [i.set_xlim(-.5, activity.shape[1]-.5) for i in ax[:, class_idx]]

        # Add timestamps
        if len(PARAMETERS['timestamps'][cl]) > 0:
            timestamps = PARAMETERS['timestamps'][cl] - 1  # These are 1-based as in Matlab
            timestamps = timestamps - .5  # So the timestamp marker will fall between two bins
            for i_ax in range(3):
                for t in timestamps:
                    ax[i_ax, class_idx].axvline(t, color='r')

    # Adjust subplot limits
    for row in range(1, 3):
        y_limits = np.vstack([i.get_ylim() for i in ax[row, :]])
        y_limits = [y_limits.min(), y_limits.max()]
        [i.set_ylim(y_limits) for i in ax[row, :]]
        [i.set_yticklabels([]) for i in ax[row, 1:]]

    # Reduce white space
    plt.tight_layout()
    # Print figure
    if print_figure:
        output_filename = os.path.join(PARAMETERS['output_figures']['folder'], PARAMETERS['output_figures']['base_filename'] + 'classification_correctly_assigned.pdf')
        if not os.path.exists(PARAMETERS['output_figures']['folder']):
            os.mkdir(PARAMETERS['output_figures']['folder'])
        plt.savefig(fname=output_filename, bbox_inches='tight', format='pdf')
        plt.close(fig)

    else:
        return fig


def heatmap_response_probability_per_trial(X, performance_single_cell, significant_activations, bins_to_cumulate, times, timestamps, title=''):
    # Set whether to show non-significant bins in a different colormap
    non_significant_monochrome = False

    # Get data shape
    n_trials, n_bins = X.shape
    # Make time axis and get the bin of stimulus onset
    dt = times[1] - times[0]
    stimulus_onset_bin = np.where(bins_to_cumulate > 0)[0][0]
    times_center = times + dt / 2

    # Make mask of significant bins
    mask = np.zeros_like(X, dtype=np.bool)
    for i_bin in range(performance_single_cell.shape[0]):
        if not performance_single_cell.loc[i_bin, 'is_selective']:
            continue
        significant_trials = performance_single_cell.loc[i_bin, 'significant_trials']
        mask[significant_trials, i_bin + stimulus_onset_bin] = True
    # Split mask in significant activations and non-significant ones
    mask_significant = np.zeros_like(mask)
    if significant_activations is not None:
        for i in range(significant_activations.shape[0]):
            mask_significant[:, significant_activations.loc[i, 'start']:significant_activations.loc[i, 'finish'] + 1] = mask[:, significant_activations.loc[i, 'start']:significant_activations.loc[i, 'finish'] + 1]
    mask_non_significant = (mask.astype(int) - mask_significant.astype(int)).astype(bool)

    # Get mean activity in baseline window
    baseline_window = np.where(bins_to_cumulate == 0)[0]
    baseline_X = X[:, baseline_window].mean()
    X -= baseline_X

    # Get ranges of image
    y_values = np.arange(n_trials)
    extent = [times[0], times[-1], y_values[0], y_values[-1] + 1]
    vmin = X.min()
    vmax = X.max()
    time_lims = times[[0, -1]]
    ylim = y_values[0], y_values[-1] + 1
    # Get colormap and aspect ratio
    cmap = center_diverging_colormap(vmin, vmax, cmap_name='seismic')
    im_args = dict(interpolation='none', aspect='auto', extent=extent, vmin=vmin, vmax=vmax)

    # Draw heatmap
    fig, ax = plt.subplots(facecolor='w', figsize=(27.22, 15.47))
    if non_significant_monochrome:
        ax.imshow(X, cmap='Greys', **im_args)
        im = ax.imshow(np.ma.masked_where(~mask, X), cmap=cmap, **im_args)
    else:
        im = ax.imshow(X, cmap=cmap, **im_args)

    # Draw contour
    if np.any(mask):
        # Mask of non-significant bins
        if mask_non_significant.any():
            big_mask = np.kron(mask_non_significant, np.ones((10, 10)))
            ax.contour(big_mask, colors=['k'], extent=extent, linestyles=['--'], linewidths=[0.5], corner_mask=False, antialiased=True, levels=[.5], extend='both', origin='upper')
        # Mask of non-significant bins
        if mask_significant.any():
            big_mask = np.kron(mask_significant, np.ones((10, 10)))
            ax.contour(big_mask, colors=['k'], extent=extent, linestyles=['-'], linewidths=[2.5], corner_mask=False, antialiased=True, levels=[.5], extend='both', origin='upper')
            # Draw limits of each significant activation, and mark peak of response probability
            for i in range(significant_activations.shape[0]):
                ax.plot([times[significant_activations.loc[i, 'start']], times[significant_activations.loc[i, 'finish'] + 1]], [0, 0], color='b', linewidth=5, solid_capstyle='butt')
                ax.plot(times[significant_activations.loc[i, 'peak_probability']] + dt / 2, [0], color='b', marker='^', markersize=15, linestyle='none')

    # Set limits
    ax.set_xlim(time_lims[0], time_lims[-1])
    ax.set_ylim(ylim)
    ax.set_xlabel('Time from stimulus onset (s)')
    ax.set_ylabel(r'Trial ${\itno.}$')
    # Add info to title
    if significant_activations is not None:
        idx = significant_activations['response_probability'].idxmax()
        max_t = times[significant_activations.loc[idx, 'peak_probability']] + dt / 2
        title_str = '%s (max ${\itp}$=%.1f%% at ${\itt}$=%.1f${\its}$ (AUC=%.2f))' % (title, significant_activations.loc[idx, 'response_probability'], max_t, significant_activations.loc[idx, 'AUC'])
    else:
        title_str = '%s (non-selective)' % title
    ax.set_title(title_str)
    # Set ticks
    x_ticks = np.unique(np.hstack((ax.get_xticks(), times[0], times[-1])))
    x_ticks = x_ticks[np.where(np.logical_and(x_ticks >= times[0], x_ticks <= times[-1]))[0]]
    y_ticks = [tick for tick in ax.get_yticks() if tick < n_trials]
    y_ticks = np.array(y_ticks).astype(int)
    ax.set(yticks=y_ticks + .5, yticklabels=n_trials - y_ticks, xticks=x_ticks)
    # Mark timestamps
    ts = timestamps - timestamps[0]
    for t in ts:
        ax.axvline(times[np.abs(t - times_center).argmin()], color='k', linestyle='--', linewidth=1)
    # Mark baseline period
    ax.plot([times[baseline_window[0]], times[np.abs(ts[0] - times_center).argmin()]], [0, 0], color=(.7, .7, .7), linewidth=4, solid_capstyle='butt')
    # Add colorbar
    cbar = colorbar(im)
    cbar.ax.set_ylabel(r'$\Delta$F/F$_0$ relative to mean of baseline')

    # Reduce white space
    plt.tight_layout(rect=[0, 0, 1, 1])

    return fig


def heatmap_response_probability_per_cell(performance_all_cells, response_probability_traces, times, timestamps, title=''):
    # Get data shape
    n_cells, n_bins = response_probability_traces.shape
    # Make time axis and get the bin of stimulus onset
    dt = times[1] - times[0]
    stimulus_onset_bin = np.abs(times).argmin()
    times_evoked = times[stimulus_onset_bin:]
    times_center = times_evoked + dt / 2

    # Prepare data
    X = response_probability_traces.copy()
    p_positive_response_mask = np.zeros_like(X, dtype=bool)
    p_negative_response_mask = np.zeros_like(X, dtype=bool)
    for i in range(performance_all_cells.shape[0]):
        i_cell = performance_all_cells.loc[i, 'cell'] - 1
        AUC = performance_all_cells.loc[i, 'AUC']
        beginning = int(performance_all_cells.loc[i, 'start'] - stimulus_onset_bin)
        end = int(performance_all_cells.loc[i, 'finish'] - stimulus_onset_bin)
        if AUC >= 0:
            p_positive_response_mask[i_cell, beginning:end] = True
        else:
            p_negative_response_mask[i_cell, beginning:end] = True
            X[i_cell, beginning:end] *= -1

    # Get ranges of image
    y_values = np.arange(n_cells)
    extent = [times_evoked[0], times_evoked[-1], y_values[0], y_values[-1] + 1]
    vmin = X.min()
    vmax = X.max()
    time_lims = times_evoked[[0, -1]]
    ylim = y_values[0], y_values[-1] + 1
    # Get colormap and aspect ratio
    cmap = center_diverging_colormap(vmin, vmax, cmap_name='seismic')
    im_args = dict(interpolation='none', aspect='auto', extent=extent, vmin=vmin, vmax=vmax)

    # Draw heatmap
    fig, ax = plt.subplots(facecolor='w', figsize=(27.22, 15.47))
    # Draw all traces in gray
    ax.imshow(X, cmap='Greys', alpha=.5, **im_args)
    # Overlay positive and negative contribution
    im = ax.imshow(np.ma.masked_where(~p_positive_response_mask, X), cmap=cmap, **im_args)
    ax.imshow(np.ma.masked_where(~p_negative_response_mask, X), cmap=cmap, **im_args)

    # Mark timestamps
    ts = timestamps - timestamps[0]
    for t in ts:
        ax.axvline(times_evoked[np.abs(t - times_center).argmin()], color='k', linestyle='--', linewidth=1)
    # Set limits
    ax.set_xlim(time_lims[0], time_lims[-1])
    ax.set_ylim(ylim)
    ax.set_xlabel('Time from stimulus onset (s)')
    ax.set_ylabel(r'Cell ${\itno.}$')
    ax.set_title(title)
    # Set ticks
    x_ticks = np.unique(np.hstack((ax.get_xticks(), times_evoked[0], times_evoked[-1])))
    x_ticks = x_ticks[np.where(np.logical_and(x_ticks >= times_evoked[0], x_ticks <= times_evoked[-1]))[0]]
    y_ticks = np.arange(n_cells)
    y_ticklabels = np.empty(y_ticks.shape[0], dtype=object)
    y_ticklabels[:] = ''
    y_ticklabels[np.arange(0, n_cells, 10)] = ['%i' % i for i in np.arange(0, n_cells, 10)]
    y_ticklabels[0] = '1'
    y_ticklabels[-1] = '%i' % n_cells
    ax.set(yticks=y_ticks + .5, yticklabels=y_ticklabels[::-1], xticks=x_ticks)

    # Add colorbar
    cbar = colorbar(im)
    cbar.ax.set_ylabel(r'Response probability (%)')

    # Maximize figure
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    # Reduce white space
    plt.tight_layout(rect=[0, 0, 1, 1])

    return fig
