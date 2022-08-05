# System packages
import sys
import os
import shutil
import tempfile

# Numerical packages
import numpy as np
import pandas as pd

# Visualization packages
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt
# Printing figures
from matplotlib.backends.backend_pdf import PdfPages

# Local repository
from decoding import visualization
from Utilities.IO_operations import log
from Utilities import matlab_file

# Set whether to show non-significant bins in a different colormap
non_significant_monochrome = True


def run(parameters_filename: str, result_filename: str, traces_filename: str, verbose: bool = True):
    # Load parameter file
    if verbose:
        log('Reading data')
    RESULTS = matlab_file.load(result_filename, variable_names='RESULTS')
    all_PARAMETERS = matlab_file.load(parameters_filename, variable_names='PARAMETERS')
    # Because output file might be on a remote location, it is better to append the file locally and then copy it
    temporary_filename = os.path.join(tempfile.gettempdir(), 'temp.pdf')
    # Create empty PDF
    PDF_file = PdfPages(temporary_filename, keep_empty=False)

    # Convert MATLAB's tables to pandas DataFrames, and keep track of stimuli and sessions
    sessions = list(RESULTS.keys())
    n_sessions = len(sessions)
    STIMULI = list()
    for i_sess in range(n_sessions):
        stimuli =  list(RESULTS[sessions[i_sess]].keys()) # stimuli = ['temp_48', 'HPS', 'FPS']
        STIMULI += stimuli
        n_stimuli = len(stimuli)
        for i_stim in range(n_stimuli):

            if len(RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance'])== 0:
                print('no cells in', str(sessions[i_sess]),'-', stimuli[i_stim])
                break

            if np.shape(RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance']) == np.shape(RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance_columns']):
                this_performance = [RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance']]
            else:
                this_performance = RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance']

            performance = pd.DataFrame(this_performance, columns=RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance_columns'])
            # Cast to correct data type
            performance['cell'] = performance['cell'].map(int)
            performance['start'] = performance['start'].map(int) - 1
            performance['finish'] = performance['finish'].map(int) - 1
            performance['response_probability'] = performance['response_probability'].map(float)
            performance['peak_probability'] = performance['peak_probability'].map(int) - 1
            performance['activity_difference'] = performance['activity_difference'].map(float)
            performance['AUC'] = performance['AUC'].map(float)
            performance['p'] = performance['p'].map(float)
            performance['effect_size'] = performance['effect_size'].map(float)
            # Store new variable
            RESULTS[sessions[i_sess]][stimuli[i_stim]]['performance'] = performance.copy()

    # Get unique list of stimuli
    STIMULI = np.unique(STIMULI)

    # Load denoised traces
    DENOISED_TRACES = matlab_file.load(traces_filename, variable_names='DENOISED_TRACES')

    # Loop over sessions
    for session_name in sessions:
        # Get the date of the session
        session_date = session_name.split('_')[1]
        # Get names of stimuli for this session
        stimuli_this_session = list(all_PARAMETERS[session_name].keys())
        # Get number of cells
        n_cells = RESULTS[session_name][stimuli_this_session[0]]['response_probability_traces'].shape[0]

        for stimulus_name in STIMULI:
            # Check whether this stimulus has been administered during this session
            if stimulus_name not in stimuli_this_session:
                continue
            # Extracts parameters for this session and stimulus
            PARAMETERS = all_PARAMETERS[session_name][stimulus_name]

            # Get timestamps
            n_bins = PARAMETERS['data_size'][0]
            timestamps_s = PARAMETERS['timestamps_s']
            if not isinstance(timestamps_s, (np.ndarray, list)):
                timestamps_s = [timestamps_s]
            timestamps_s = np.array(timestamps_s, dtype=int)
            # Make time axes for both plots
            n_frames = DENOISED_TRACES[session_name][stimulus_name][0, 0].shape[0]
            times_traces = np.linspace(PARAMETERS['trial_time'][0], PARAMETERS['trial_time'][1], n_frames) - timestamps_s[0]
            times_heatmap = np.linspace(PARAMETERS['trial_time'][0], PARAMETERS['trial_time'][1], n_bins + 1) - timestamps_s[0]

            # Log progress
            if verbose:
                log('Processing %s - %s' % (session_name, stimulus_name))

            # Get the observation markers
            bins_to_cumulate = PARAMETERS['bins_to_cumulate']
            bins_to_cumulate[np.isnan(bins_to_cumulate)] = -1
            bins_to_cumulate = bins_to_cumulate.astype(int)
            # Get the data and its original [trial x time] shape
            data = PARAMETERS['data']
            data_shape = PARAMETERS['data_size']

            # Compute selectivity of individual cells
            performance = RESULTS[session_name][stimulus_name]['performance']
            significant_trials = RESULTS[session_name][stimulus_name]['significant_trials']
            if len(performance) == 0:
                print('no significant cells for ' + stimulus_name + session_name)
                break
            for i_cell in range(n_cells):
                # Reshape data
                X = data[:, i_cell].reshape(-1, 1).reshape(data_shape[::-1]).copy()

                # Get data shape
                n_trials, n_bins = X.shape
                # Make time axis and get the bin of stimulus onset
                dt = times_heatmap[1] - times_heatmap[0]
                stimulus_onset_bin = np.where(bins_to_cumulate > 0)[0][0]
                times_center = times_heatmap + dt / 2

                # Make mask of significant bins
                mask = np.zeros_like(X, dtype=np.bool)

                if not np.shape(mask[:, np.where(bins_to_cumulate>0)[0]]) == np.shape(significant_trials[i_cell, :, :].transpose()):
                    log('something went wrong with the binning, check bins_to_cumulate')
                    print()

                else:
                    mask[:, np.where(bins_to_cumulate > 0)[0]] = significant_trials[i_cell, :,:].transpose()

                # if np.size(mask[:, stimulus_onset_bin: - stimulus_onset_bin],1) < np.size(significant_trials[i_cell, :, :].transpose(),1):
                #     mask[:, stimulus_onset_bin: - stimulus_onset_bin + 1 ] = significant_trials[i_cell, :, :].transpose()
                # else:
                #     mask[:,
                #     stimulus_onset_bin: - stimulus_onset_bin] = significant_trials[i_cell, :,:].transpose()

                # Get mean activity in baseline window
                baseline_window = np.where(bins_to_cumulate == 0)[0]
                baseline_X = X[:, baseline_window].mean()
                X -= baseline_X

                # Get ranges of image
                y_values = np.arange(n_trials)
                extent = [times_heatmap[0], times_heatmap[-1], y_values[0], y_values[-1] + 1]
                vmin = X.min()
                vmax = X.max()
                time_lims = times_heatmap[[0, -1]]
                ylim = y_values[0], y_values[-1] + 1
                # Get colormap and aspect ratio
                cmap = visualization.center_diverging_colormap(vmin, vmax, cmap_name='seismic')
                im_args = dict(interpolation='none', aspect='auto', extent=extent, vmin=vmin, vmax=vmax)

                # Open figure and rename axes
                fig, ax = plt.subplots(nrows=1, ncols=2, facecolor='w', figsize=(15, 10))
                ax_traces = ax[0]
                ax_heatmap = ax[1]

                # Draw traces
                traces = np.vstack(DENOISED_TRACES[session_name][stimulus_name][i_cell, :])
                max_value = traces.max()
                # Assign mock value to create offset between trials in case of no activity
                if max_value == 0:
                    max_value = 1
                y_limits = [np.PINF, np.NINF]
                for i_trial in range(n_trials):
                    offset = (n_trials - i_trial + 1) * (max_value / 3)
                    this_trace = traces[i_trial, :] + offset
                    # Adjust y-limits to accommodate this trace
                    if this_trace.min() < y_limits[0]:
                        y_limits[0] = this_trace.min()
                    if this_trace.max() > y_limits[1]:
                        y_limits[1] = this_trace.max()
                    # Set a line-style based on where there are spikes or not
                    if np.all(traces[i_trial, :] == 0):
                        ls = '--'
                    else:
                        ls = '-'
                    # Draw trace
                    ax_traces.plot(times_traces, this_trace, color='k', ls=ls)
                # Draw stimulus onset and all the other timestamps
                [ax_traces.axvline(ts - timestamps_s[0], color=(.7, .7, .7), lw=2) for ts in timestamps_s]
                # Fix axes appearance
                ax_traces.set_xlim(time_lims[0], time_lims[-1])
                ax_traces.set_ylim(*y_limits)
                ax_traces.set_xlabel('Time from stimulus onset (s)')
                visualization.adjust_spines(ax_traces, spines=['bottom'])

                # Draw heatmap
                if non_significant_monochrome:
                    ax_heatmap.imshow(X, cmap='Greys', **im_args)
                    im = ax_heatmap.imshow(np.ma.masked_where(~mask, X), cmap=cmap, **im_args)
                else:
                    im = ax_heatmap.imshow(X, cmap=cmap, **im_args)

                # Get performance of this cell
                performance_this_cell = performance.loc[np.where(performance['cell'] == (i_cell + 1))[0]]
                performance_this_cell.reset_index(drop=True, inplace=True)
                has_significant_activations = performance_this_cell.shape[0] > 0 # !! Check for inhibited cells!!

                # Draw contour of significant activations
                if has_significant_activations:
                    ax_heatmap.contour(np.kron(mask, np.ones((10, 10))), colors=['k'],
                               extent=extent, linestyles=['-'], linewidths=[2.0],
                               corner_mask=False, antialiased=True, levels=[0.5],
                               extend='both', origin='upper')
                    # Draw limits of each significant activation, and mark peak of response probability
                    for i_act in range(performance_this_cell.shape[0]):
                        # Line to mark period
                        ax_heatmap.plot([times_heatmap[performance_this_cell.loc[i_act, 'start']],
                                 times_heatmap[performance_this_cell.loc[i_act, 'finish']]],
                                [0, 0], color='b', linewidth=5, solid_capstyle='butt')
                        # Arrow head to mark peak
                        ax_heatmap.plot(times_heatmap[performance_this_cell.loc[i_act, 'peak_probability']] + dt / 2,
                                [0], color='b', marker='^', markersize=15, linestyle='none')

                # Set limits
                ax_heatmap.set_xlim(time_lims[0], time_lims[-1])
                ax_heatmap.set_ylim(ylim)
                ax_heatmap.set_xlabel('Time from stimulus onset (s)')
                ax_heatmap.set_ylabel(r'Trial ${\itno.}$')
                # Set ticks
                x_ticks = np.unique(np.hstack((ax_heatmap.get_xticks(), times_heatmap[0], times_heatmap[-1])))
                x_ticks = x_ticks[np.where(np.logical_and(x_ticks >= times_heatmap[0], x_ticks <= times_heatmap[-1]))[0]]
                y_ticks = [tick for tick in ax_heatmap.get_yticks() if tick < n_trials]
                y_ticks = np.array(y_ticks).astype(int)
                ax_heatmap.set(yticks=y_ticks + .5, yticklabels=n_trials - y_ticks, xticks=x_ticks)
                ax_traces.set(xticks=x_ticks)
                # Mark timestamps
                ts = timestamps_s - timestamps_s[0]
                for t in ts:
                    ax_heatmap.axvline(times_heatmap[np.abs(t-times_center).argmin() +1], color='k', linestyle='--', linewidth=1)
#                    ax_heatmap.axvline(times_heatmap[np.abs(t - times_center).argmin()], color='k', linestyle='--', linewidth=1) # This somehow was giving the wrong position
                # Mark baseline period
                ax_heatmap.plot([times_heatmap[baseline_window[0]],
                         times_heatmap[np.abs(ts[0] - times_center).argmin()]], [0, 0],
                         color=(.7, .7, .7), linewidth=4, solid_capstyle='butt')
                # Add colorbar
                cbar = visualization.colorbar(im)
                cbar.ax.set_ylabel(r'$\Delta$F/F$_0$ relative to mean of baseline')

                # Add info to title
                title = 'Session %s - Stimulus %s - Cell %i' % (session_date, stimulus_name, i_cell + 1)
                if has_significant_activations:
                    idx = performance_this_cell['response_probability'].idxmax()
                    max_t = times_heatmap[performance_this_cell.loc[idx, 'peak_probability']] + dt / 2
                    title_str = '%s (max ${\itp}$=%.1f%% at ${\itt}$=%.1f${\its}$, AUC=%.2f)' % (
                            title, performance_this_cell.loc[idx, 'response_probability'],
                            max_t, performance_this_cell.loc[idx, 'AUC'])
                else:
                    title_str = '%s (non-selective)' % title
                fig.suptitle(title_str)

                # Reduce white space
                plt.tight_layout(rect=[0, 0, 1, .95])

                # Save and close figure
                PDF_file.savefig(fig, bbox_inches='tight')
                # plt.close(fig)


    # Close PDF
    PDF_file.close()
    # Move PDF to final location
    shutil.copy2(temporary_filename, all_PARAMETERS['figure_filename'])

    if verbose:
        log('Finished')


################################################################################
### Direct call
################################################################################
if __name__ == "__main__":
    # Get user inputs
    if len(sys.argv) > 1:
        run(**dict(arg.split('=') for arg in sys.argv[1:]))

    else:
        # See more rows and columns of variables printed in console
        np.set_printoptions(suppress=True)
        pd.set_option('display.max_columns', 500)
        pd.set_option('display.width', 1000)

        run(parameters_filename=r'D:\_MATLAB_CaImaging\response_detection_parameters.mat',
            result_filename=r'D:\_MATLAB_CaImaging\response_detection_results.mat',
            traces_filename=r'D:\_MATLAB_CaImaging\response_detection_traces.mat',
            verbose=True)

