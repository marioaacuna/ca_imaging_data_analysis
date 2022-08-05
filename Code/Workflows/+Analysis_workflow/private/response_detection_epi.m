function response_detection_epi (animal_ID, traces, missing_only )
% Selectivity for Stimuli. Part of MS analysis pipeline

global GC LOGGER
GC.experiment_name ='ACC_SNI_anxiety';

% get Parameters

METADATA = SQL_database.read_epi_trials (animal_ID);
% filename = get_filename_of('dFF', animal_ID);
% traces = load_variable (filename, 'dFF');
% experimental_group = SQL_database.read_table_where('experiments', 'experimental_group', animal_ID,'animal_ID', 'return_as_table',false);
%METADATA = METADATA(ismemberCellRows(METADATA.experiment(:),{'pain'}), :);

% delete sessions in which no timestamps have been uploaded yet
%METADATA(ismemberCellRows(METADATA.type, {'other'}),:) = [];
% delete stimuli where the time interval is lower that 15 secs 
%diff_timestamps = abs(diff(METADATA.timestamps(:,1)));
%diff_idx = find(diff_timestamps < 0.015);
%METADATA(diff_idx,:) = [];


n_cells = size(traces, 1);
frame_rate =  GC.epifluorescence_downsample_to_frame_rate;

trials_idx = find(ismemberCellRows(METADATA.type(:), {'evoked'}) & ismemberCellRows(METADATA.experiment(:),{'pain'}));
dates = unique(METADATA(trials_idx, :).date);
sessions = (unique(METADATA.day_from_surgery(trials_idx), 'stable'));
n_sessions = length (dates);
stimuli = unique(METADATA.stimulus(trials_idx));
n_stim = length (stimuli);

type_threshold = 'p_val';


% get output variables
result_filename = get_filename_of('response_detection_epi', animal_ID, type_threshold, missing_only);
root_path = [GC.data_root_path ,'\',GC.aggregated_data_path,'\', GC.experiment_name, '\', 'response_detection\'];
if ~exist (root_path, 'dir')
    mkdir(root_path)
end

if missing_only && exist(result_filename, 'file')
    RESULTS = load_variable(result_filename, 'RESULTS');
else  % Overwrite previous analysis or simply initialize file
    RESULTS = struct();
end
% RESULTS.session.stim --> {selectivity per cell}, {mean zScore per cell}, {std Zscore per cell}


for i_sess = 1: n_sessions
    date = dates{i_sess}; %  '190821';

    trials_idx = find(ismemberCellRows(METADATA.date, {date}) & ismemberCellRows(METADATA.type(:), {'evoked'})) ;
    traces_this_sess = traces(:, trials_idx);

    name_session = sessions{i_sess};
    idx_sess = find(ismemberCellRows(METADATA.day_from_surgery, {name_session}));
    METADATA_this_sess =  METADATA(ismemberCellRows(METADATA.day_from_surgery(:),{name_session}), :);
    GC.session_column_prefix = name_session;
    field_name = make_variable_name(['session_', GC.session_column_prefix]);
    for i_stim = 1 : n_stim
        
        name_stimulus = stimuli {i_stim};
        stim_idx = find(ismemberCellRows(METADATA_this_sess.stimulus, {name_stimulus}));
        n_trials = length(stim_idx);
        
        is_selective = NaN(n_cells,1);
        is_excited = NaN(n_cells,1);
        is_inhibited = NaN(n_cells,1);
        mean_traces = cell (n_cells, 1);
        std_traces = cell (n_cells, 1);
        results = struct();

        for i_cell = 1 : n_cells
           
            data_trials = cell(n_trials, 1);
            for i_trial = 1: n_trials
                data_trials{i_trial} = traces_this_sess{i_cell, stim_idx(i_trial)};
                data_trials{i_trial} = data_trials{i_trial} - median(data_trials{i_trial});
            end
            
            data_trials = cell2mat(horzcat(data_trials));
            trial_duration = GC.miniscope_trial_duration.(name_stimulus);
            baseline = 1: (frame_rate * trial_duration(1));
            evoked = (1 + (frame_rate * trial_duration (1)))  : ((trial_duration(2) + trial_duration(1))  * frame_rate);
           
            avg_baseline = mean(data_trials(:, baseline), 2);
            std_baseline = std(data_trials(:, baseline, :), 0, 2);
            z_data = (data_trials - avg_baseline) ./ std_baseline;
            switch type_threshold
                case 'std' 
                    mean_trace = mean(z_data,1);
                    selectivity_threshold = (std(mean_trace(:, baseline) * 3)) ;%+ median(mean_trace,2));
                    is_selective(i_cell,:) = any(mean_trace(:, evoked, :) > (selectivity_threshold), 2);
                    mean_traces{i_cell,:} =  mean(z_data,1);
                    std_traces{i_cell,:} = std(z_data,1);
                case 'p_val'
                    data_bl = sum(z_data(:, baseline),2);
                    data_ev  = sum(z_data(:, evoked),2); % or 5s after stim
                    p = wilcoxon_signrank(data_bl, data_ev);
                    is_selective(i_cell,:) = p < 0.05;
                    is_excited(i_cell,:)   = (p < 0.05 & mean(data_ev) > mean(data_bl));
                    is_inhibited(i_cell,:) = p < 0.05 & mean(data_ev) < mean(data_bl);
                    mean_traces{i_cell,:}  = mean(z_data,1);
                    std_traces{i_cell,:}   = std(z_data,1);

            end
        end
        results.selectivity.selective = is_selective;
        results.selectivity.excited = is_excited;
        results.selectivity.inhibited = is_inhibited;
        results.mean_traces =  mean_traces;
        results.std_traces = std_traces;
        RESULTS.(field_name).(name_stimulus) = results;
     end
    RESULTS.(field_name).(name_stimulus) = results;
end
 
    

% Wri_sesse file to disk
LOGGER.info('Storing data to disk')
save(result_filename, 'RESULTS', '-v7.3')
LOGGER.info('(ok)', 'append',true)



%%
% for i_stim = 1 : length(stimuli)
%     figure
%     hold on
%     title (strcat (stimuli(i_stim)))
%     traces = RESULTS.session__minus58.(cell2mat(stimuli(i_stim))).mean_traces(:,1);
%     for i_cell = 1 : n_cells
%         this_avg  = traces{i_cell,:};
%         plot(this_avg)
%     end
% end


%%

% 
% filename_spikes = get_filename_of('spikes', animal_ID);
% spikes = load_variable (filename_spikes, 'spikes');
% 
% traces_cell1 = traces {1, 70};
% spikes_cell1 = spikes {1, 70};
% 
% figure, plot(traces_cell1)
% ay = gca;  ay.YAxisLocation  ='left';
% hold on
% figure, plot(spikes_cell1, '-o')
% ax = gca;  ax.YAxisLocation  ='right';
% hold off
