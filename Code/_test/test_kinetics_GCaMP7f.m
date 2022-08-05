%%
clc
clear
global GC
GC.experiment_name = 'CLA_pain';
animal_ID = 'NRN_5d'
cell_ID = 12;
use_smooth = 0;

METADATA = SQL_database.read_table_where('trials', {'+', 'stimulus', 'experimental_condition','calcium_indicator' }, animal_ID, 'animal_ID') ;
frame_rate = max(METADATA.frame_rate);
calcium_indicator = cell2mat(unique(METADATA.calcium_indicator));
decay_time = GC.decay_time_constant_calcium_reporter.(calcium_indicator);
decay_time_new = 4000;
rise_time_new = 1.2;
spikes_filename = get_filename_of('spikes', animal_ID);
spikes_new = load_variable(spikes_filename, 'spikes');
spikes_filename = 'V:\Ca_imaging_pain\5_deconvolved_traces\NRN_5d_inferred_spikes_old.mat';
spikes_old = load_variable(spikes_filename, 'spikes');

traces_filename = get_filename_of ('dFF', animal_ID);
traces = load_variable(traces_filename, 'dFF');
stimulus_to_take = 'HPSnl';
trials_idx = ismember(METADATA.stimulus, stimulus_to_take);

spikes_old = spikes_old(:, 2);
spikes_new = spikes_new(:, 2);
traces = traces(:, 2);

spikes_old_single_cell =spikes_old(cell_ID, :);
spikes_new_single_cell = spikes_new(cell_ID, :);
traces_single_cell = cell2mat(traces(cell_ID, :));
% Smooth trace
smoothing_width = floor(1000 / 1000 * frame_rate / 2) * 2 + 1;
[~, g] = sgolay(1, smoothing_width);
diff_filter = 1;
if use_smooth
    traces_single_cell_smooth = conv(traces_single_cell, factorial(diff_filter-1)/(-(1/frame_rate))^(diff_filter-1) * g(:,diff_filter), 'same');
else
    traces_single_cell_smooth = traces_single_cell;
end
denoised_old = convolve_spikes(spikes_old_single_cell, frame_rate, decay_time);
denoised_new = convolve_spikes(spikes_new_single_cell, frame_rate, decay_time);
% denoise the traces with wider decay time
denoised_old_wider = convolve_spikes(spikes_old_single_cell, frame_rate, decay_time_new, rise_time_new);
denoised_new_wider = convolve_spikes(spikes_new_single_cell, frame_rate, decay_time_new, rise_time_new);

disp('parameteres loaded')
%%
clc
close all
xlimits = [1, length(traces_single_cell)];%length(traces_single_cell)]; %[250,1000]
fig = figure('Color', 'w');
clf
title_this_fig =  ['spikes_new, decay: ', num2str(decay_time)];
ax1 = subplot(2,1,1);
plot(ax1,traces_single_cell_smooth, 'Color', 'r')
hold on
plot(ax1,cell2mat(spikes_new_single_cell), 'Color', 'k')
plot(ax1,plot(cell2mat(denoised_new), 'Color', 'bl', 'LineWidth', 2));
ylim([0, max(traces_single_cell(xlimits(1):xlimits(2)))]) % min(traces_single_cell(xlimits(1):xlimits(2)))
xlim(xlimits)
title(ax1,title_this_fig, 'Interpreter', 'none')
% fig old_spikes
ax2 = subplot(2,1,2);
title_this_fig =  ['spikes_old, decay: ', num2str(decay_time)];
plot(ax2,traces_single_cell_smooth, 'Color', 'r')
hold on
plot(ax2,cell2mat(spikes_old_single_cell), 'Color', 'k')
plot(ax2,plot(cell2mat(denoised_old), 'Color', 'bl', 'LineWidth', 2));
ylim([0, max(traces_single_cell(xlimits(1):xlimits(2)))]) 
xlim(xlimits)
title(ax2,title_this_fig, 'Interpreter', 'none')

% New figure with wider decay time
fig_2 = figure('Color', 'w');
title_this_fig =  ['spikes_new, decay: ', num2str(decay_time_new)];
ax1 = subplot(2,1,1);
plot(ax1,traces_single_cell_smooth, 'Color', 'r')
hold on
plot(ax1,cell2mat(spikes_new_single_cell), 'Color', 'k')
plot(ax1,plot(cell2mat(denoised_new_wider), 'Color', 'bl', 'LineWidth', 2));
ylim([0, max(traces_single_cell(xlimits(1):xlimits(2)))])
title(ax1,title_this_fig, 'Interpreter', 'none')
xlim(xlimits)
% fig old_spikes
ax2 = subplot(2,1,2);
title_this_fig =  ['spikes_old, decay: ', num2str(decay_time_new)];
plot(ax2,traces_single_cell_smooth, 'Color', 'r')
hold on
plot(ax2,cell2mat(spikes_old_single_cell), 'Color', 'k')
plot(ax2,plot(cell2mat(denoised_old_wider), 'Color', 'bl', 'LineWidth', 2));
ylim([0, max(traces_single_cell(xlimits(1):xlimits(2)))])
xlim(xlimits)
title(ax2,title_this_fig, 'Interpreter', 'none')

if strcmp(stimulus_to_take, 'HPSnl')
    fig3 = figure('Color', 'w');
    title_this_fig =  ['spikes_old, decay: ', num2str(decay_time_new)];
    plot(traces_single_cell_smooth, 'Color', 'r')
    hold on
    plot(cell2mat(spikes_old_single_cell), 'Color', 'k')
    plot(cell2mat(denoised_old_wider), 'Color', 'bl', 'LineWidth', 2);
    ylim([0, max(traces_single_cell(xlimits(1):xlimits(2)))])
    xlim(xlimits)
    title(title_this_fig, 'Interpreter', 'none')
    n_frames = length(traces_single_cell);
    STIMULUS = Metadata_workflow.load_stimuli(stimulus_to_take, 'frame_rate', frame_rate, 'n_frames', n_frames, 'ignore_if_missing',true);
    timestamp = round(STIMULUS.timestamps .* frame_rate); 
    plot(timestamp .* [1, 1], [-1, 1+1], 'k', 'LineStyle','--', 'YLimInclude','off')

    
end


