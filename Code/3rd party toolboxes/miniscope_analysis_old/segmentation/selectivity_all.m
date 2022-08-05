% Run this script after the the 'add_stimuli_to_traces' python code is
% done. the input is ['traces_',name_session,'_stamped.mat'] / or
% spikes_stamped

clc;clear all
cd ('M:\Mario\MATLAB_scripts\Others\MovieAnalysis_B_Grewe\MovieAnalysis\segmentation\')
%% prepare data
animal_list = {'I32_1', 'I32_3', 'I37_2', 'I38_3', 'I38_4', 'I39_1'};
is_CCI = [false, true, false, false, true, true];
% animals_folder = 'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\';
n_animals = length(animal_list);
n_responding = single(false(4, n_animals));
session = '2'; % '1' ; '2';
data_form = 'raw_amp';% select : 'raw_amp', 'z_score'
prop_per_animal = single(false( 4,n_animals));
for n_animal = 1: length (animal_list)
    
    animal_ID =  str2mat(animal_list(n_animal)); %#ok<DSTRMT>
    name_session = ['pain_', session];
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
    switch data_form
        case 'raw_amp'
            ampl_selective_per_animal_raw = single(false(n_cells, n_stimuli));
        case 'z_score'
            ampl_selective_per_animal_zs = single(false(n_cells, n_stimuli));
    end
   
    for i_stim = 1:n_stimuli
        data = data_pain_TS.(stimuli{i_stim});
        
        % Transform to z-score
        if strcmp(signal_type, 'traces')
            avg_baseline = mean(data(:, baseline, :), 2);
            std_baseline = std(data(:, baseline, :), 0, 2);
            data_zs = (data - avg_baseline) ./ std_baseline;
        end
        
        % Check whether the trial-averaged activity crosses the threshold in
        % the evoked epoch
        selective_cells(:, i_stim) = any(mean(data_zs(:, evoked, :), 3) > selectivity_threshold, 2);
        selective_idx  = (find((selective_cells(:, i_stim) == true)));
        for selective_cell = 1:length(selective_idx)
            switch data_form
                case 'raw_amp'
                    peak_trials = single(false(length (data_pain_TS.(stimuli{i_stim})(1,1,:)), 1));
                    
                    for i_trial = 1: length (data_pain_TS.(stimuli{i_stim})(1,1,:))
                        peak_trial = max(data(selective_idx(selective_cell), evoked, i_trial)) ;
                        peak_trials (i_trial) = peak_trial;
                        %ampl_selective_per_animal(selective_idx(selective_cell), i_stim) = (max(mean(data(selective_idx(selective_cell), evoked, :), 3)')');
                    end
                    ampl_selective_per_animal_raw(selective_idx(selective_cell), i_stim) = mean(peak_trials);
                case 'z_score'
                    ampl_selective_per_animal_zs(selective_idx(selective_cell), i_stim) = (max(mean(data_zs(selective_idx(selective_cell), evoked, :), 3)')');
            end
            
        end
        switch data_form
            case 'raw_amp' 
                ampl_selective_per_animal_raw(ampl_selective_per_animal_raw==0)=NaN;
            case 'z_score'
                ampl_selective_per_animal_zs(ampl_selective_per_animal_zs==0)=NaN;
        end
        % Put back
        data_pain_TS.(stimuli{i_stim}) = data;
    end
    % proportion activated cells per animal
    prop_per_animal(:,n_animal) = sum(any(selective_cells,3))/n_cells;
    
    % Save file to disk per animal
    output_folder = 'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\';
    if ~exist(output_folder, 'dir'), mkdir(output_folder), end
    
    % Save data
    filename = [output_folder, ['results_analyis_', animal_ID,'_ses_',session,'.mat']];
    switch data_form
        case 'raw_amp'
            save(filename, 'selective_cells','prop_per_animal' ,'stimuli', 'ampl_selective_per_animal_raw')
        case 'z_score'
            save(filename,'selective_cells', 'prop_per_animal' , 'stimuli', 'ampl_selective_per_animal_zs')
    end
    disp(['analysis of ',animal_ID,' stored in ''', filename, ''''])
end

%% load data

all_data_cci = [];
all_data_sham = [];
for n_animal = 3:n_animals % excluding the first 2 animals due to expression of Gcamp under different promoter (hsyn instead of CamkII)
    animal_ID =  str2mat(animal_list(n_animal));
    if is_CCI(n_animal)
        data= load (['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\results_analyis_', animal_ID,'_ses_',session,'.mat']);
        switch data_form
            case 'raw_amp'
                all_data_cci = vertcat(all_data_cci,data.ampl_selective_per_animal_raw);
            case 'z_score'
                all_data_cci = vertcat(all_data_cci,data.ampl_selective_per_animal_zs);
        end
    else
        data= load (['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\results_analyis_', animal_ID,'_ses_',session,'.mat']);
        switch data_form
            case 'raw_amp'
                all_data_sham = vertcat(all_data_sham,data.ampl_selective_per_animal_raw);
            case 'z_score'
                all_data_sham = vertcat(all_data_sham,data.ampl_selective_per_animal_zs);
        end
    end

end

%% Stats
avg_evoked_cci = nanmean (all_data_cci,1);
avg_evoked_sham = nanmean (all_data_sham,1);
prop_resp_cci = sum(any(all_data_cci,3))./ length(all_data_cci);
prop_resp_sham = sum(any(all_data_sham,3))/ length(all_data_sham);

[~ ,p_cold]=ttest2(all_data_cci(:,1),all_data_sham(:,1), 'tail', 'both');
[~ ,p_heat]=ttest2(all_data_cci(:,2),all_data_sham(:,2), 'tail', 'both');
[~ ,p_pin]=ttest2(all_data_cci(:,3),all_data_sham(:,3), 'tail', 'both');
[~ ,p_touch]=ttest2(all_data_cci(:,4),all_data_sham(:,4), 'tail', 'both');
disp(['p_cold: ', mat2str(p_cold)])
disp(['p_heat: ', mat2str(p_heat)])
disp(['p_pinprick: ', mat2str(p_pin)])
disp(['p_touch: ', mat2str(p_touch)])
disp(['prop_responding_cci: ', mat2str(prop_resp_cci)])
disp(['prop_responding_sham: ', mat2str(prop_resp_sham)])

% Save file to disk
output_folder = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\'];
if ~exist(output_folder, 'dir'), mkdir(output_folder), end

% Save data
filename = [output_folder, ['all_results','_ses_',session,'.mat']];
save(filename, 'all_data_cci', 'all_data_sham', 'p_cold', 'p_heat', 'p_pin', 'p_touch', 'prop_resp_cci', 'prop_resp_sham')
disp(['Results stored in ''', filename, ''''])

 %% Plot results
% fig = figure('color', 'w');
% 
% clf
% ax = [];
% for i_stim = 1:n_stimuli
%     data = data_pain_TS.(stimuli{i_stim});
%     
%     ax(i_stim) = subplot(2, 2, i_stim);
%     hold on
%     non_selective_idx = find(selective_cells(:, i_stim) == false);
%     plot(mean(data(non_selective_idx, :, :), 3).', 'color',[.7, .7, .7])
%     selective_idx = find(selective_cells(:, i_stim) == true);
%     plot(mean(data(selective_idx, :, :), 3).', 'linewidth',1.5)
%     %suptitle ([animal_ID, '-' , name_session])
%     plot([stimulus_onset, stimulus_onset], [-100, 100], 'color','k', 'YLimInclude','off')
%     
%     set(gca, 'Box','off', 'TickDir','out')
%     title(stimuli{i_stim})
% end
