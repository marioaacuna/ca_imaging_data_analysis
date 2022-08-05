% Run this script after the the 'add_stimuli_to_traces' python code is
% done. the input is ['traces_',name_session,'_stamped.mat'] / or
% spikes_stamped

clc;clear all
cd ('M:\Mario\MATLAB_scripts\Others\MovieAnalysis_B_Grewe\MovieAnalysis\segmentation\')

animal_ID = 'I38_4';
name_session = 'pain_2';
signal_type = 'traces';  % 'spikes' or 'traces'
stimulus_onset = 25.5;
baseline = 1:25;
evoked = 26:50;
selectivity_threshold = 3;  % in SD units

root_path = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\jointExtraction\sorted\'];
switch signal_type
    case 'traces'
        filename = [root_path, 'traces_',name_session,'_stamped.mat'];
    case 'spikes'
        filename = [root_path, 'spikes_',name_session,'_stamped.mat'];
end
data_pain_TS = load_variable(filename);

stimuli = fieldnames(data_pain_TS);
n_stimuli = length(stimuli);
n_cells = size(data_pain_TS.(stimuli{1}), 1);

% Allocate variable to mark which cells are selective for each stimulus
selective_cells = false(n_cells, n_stimuli);

for i_stim = 1:n_stimuli
    data = data_pain_TS.(stimuli{i_stim});
    
    % Transform to z-score
    if strcmp(signal_type, 'traces')
        avg_baseline = mean(data(:, baseline, :), 2);
        std_baseline = std(data(:, baseline, :), 0, 2);
        data = (data - avg_baseline) ./ std_baseline;
    end
    
    % Check whether the trial-averaged activity crosses the threshold in
    % the evoked epoch
    selective_cells(:, i_stim) = any(mean(data(:, evoked, :), 3) > selectivity_threshold, 2);
    
    % Put back
    data_pain_TS.(stimuli{i_stim}) = data;
end


fig = figure('color', 'w');
%% Plot results
clf
ax = [];
for i_stim = 1:n_stimuli
    data = data_pain_TS.(stimuli{i_stim});
    
    ax(i_stim) = subplot(2, 2, i_stim);
    hold on
    non_selective_idx = find(selective_cells(:, i_stim) == false);
    plot(mean(data(non_selective_idx, :, :), 3).', 'color',[.7, .7, .7])
    selective_idx = find(selective_cells(:, i_stim) == true);
    plot(mean(data(selective_idx, :, :), 3).', 'linewidth',1.5)
    %suptitle ([animal_ID, '-' , name_session])
    plot([stimulus_onset, stimulus_onset], [-100, 100], 'color','k', 'YLimInclude','off')
    
    set(gca, 'Box','off', 'TickDir','out')
    title(stimuli{i_stim})
end
