function response_detection(experiment_type, animal_ID, spikes, sessions_to_analyze, METADATA, STIMULI, missing_only)
% This funciotn detects significant trials based on Logistic Regression analysis and ROC, implemented in python.


% Implementation 0: old, 1:epi (threshold p < 0.05)
if endsWith(animal_ID, 'epi')
    implementation = 1;
else
    implementation = 0;
end


if implementation == 1
    python_script = 'detect_activations_old_p_0_05';
else
    python_script = 'detect_activations_old';
end


% Get general_configs
global GC LOGGER
% experiment_type = experiment_type;
% Get options from general_configs
downsampling_window_width = GC.response_detection_window_width;
n_permutations_significance = GC.detection_n_permutations_significance;

% Get the number of ROIs
n_cells = size(spikes, 1);
% Get frame rate
frame_rate = max(METADATA.frame_rate);
time_bin_activity_n_frames = round(downsampling_window_width * frame_rate);
baseline_window_n_bins = round(GC.detection_baseline_window / downsampling_window_width);
indicator = SQL_database.read_table_where('experiments', {'Calcium_indicator'}, animal_ID, 'animal_ID', 'return_as_table', false);
if strcmp(indicator, 'GCaMP7f')
    GC.evoked_activity_max_latency = 20;
end
evoked_window_n_bins = round(GC.evoked_activity_max_latency * frame_rate / time_bin_activity_n_frames);

% Get output variable
result_filename = get_filename_of('response_detection_p_0_01', animal_ID, experiment_type);
if missing_only && exist(result_filename, 'file')
    RESULTS = load_variable(result_filename, 'RESULTS');
else  % Overwrite previous analysis or simply initialize file
    RESULTS = struct();
    if exist(result_filename, 'dir')
        file_names = strsplit(result_filename, animal_ID);
        mkdir(file_names{1})
    end
end
% make dir if it doesnt exist
file_names = strsplit(result_filename, animal_ID);
dir_name = file_names{1};
if ~exist(dir_name, 'dir'); mkdir(dir_name), end

% Make temporary folder
temp_folder = [GC.temp_dir, '\',animal_ID];
if ~exist(temp_folder, 'dir')
    mkdir(temp_folder)
end

%% Decoding analysis
n_conditions = size(sessions_to_analyze, 1);
data_type = SQL_database.read_table_where('experiments', 'data_type', animal_ID,'animal_ID', 'return_as_table',false);
switch data_type
    case '2p'
        
        for i_cond = 1:n_conditions
            columns_to_read = GC.analysis.(experiment_type).group_trials_by_session;
            str_message = 'Processing ';
            for ifield = 1:length(columns_to_read)
                str_message = [str_message, columns_to_read{ifield}, ' ''', sessions_to_analyze{i_cond, ifield}, ''', '];
            end
            LOGGER.info(str_message(1:end-2))
            
            fieldname = '';
            for ifield = 1:length(columns_to_read)
                if strcmp(columns_to_read{ifield}, 'date')
                    fieldname = [fieldname, 'session', sessions_to_analyze{i_cond, ifield}, '_'];
                else
                    fieldname = [fieldname, sessions_to_analyze{i_cond, ifield}, '_'];
                end
            end
            fieldname = fieldname(1:end-1);
            stimulus_name = make_variable_name(sessions_to_analyze{i_cond, ismember(columns_to_read, 'stimulus')});
            if strcmp(stimulus_name, 'temp_42') || strcmp(stimulus_name, 'temp_43')
                keyboard;
            end
            result_fieldname = sessions_to_analyze{i_cond, ismember(columns_to_read, GC.analysis.(experiment_type).session_column)};
            if endsWith(result_fieldname, ';'), result_fieldname = result_fieldname(1:end-1);end
            if endsWith(fieldname, ';'), fieldname = fieldname(1:end-1);end
            % Initialize structure containing data and parameters for analysis
            PARAMETERS = struct();
            PARAMETERS.frame_rate = frame_rate;
            PARAMETERS.n_frames_per_bin = time_bin_activity_n_frames;
            PARAMETERS.data = [];
            PARAMETERS.data_size = [];
            PARAMETERS.stimulus_profile = [];
            PARAMETERS.trial_time = [];
            PARAMETERS.timestamps = [];
            PARAMETERS.timestamps_s = [];
            PARAMETERS.stimulus_onset = [];
            PARAMETERS.bins_to_cumulate = [];
            PARAMETERS.evoked_activity_max_latency = 100;
            PARAMETERS.evoked_activity_max_latency_bins = round(PARAMETERS.evoked_activity_max_latency / downsampling_window_width);
            PARAMETERS.CV_k_folds = GC.detection_k_fold_crossvalidation;
            PARAMETERS.significance_test_n_permutations = n_permutations_significance;
            PARAMETERS.output_filename = [temp_folder, filesep(), fieldname, '_detection_results.txt'];
            
            % Get data
            trial_idx = METADATA.keep;
            for ifield = 1:length(columns_to_read)
                trial_idx = trial_idx & ismember(table2cell(METADATA(:, columns_to_read{ifield})), sessions_to_analyze{i_cond, ifield});
            end
            data = spikes(:, trial_idx);
            % Get timestamps
            if startsWith(stimulus_name, 'temp_')
                timestamps = round(STIMULI.temp_48.timestamps * frame_rate);
            else
                timestamps = round(STIMULI.(stimulus_name).timestamps * frame_rate);
            end
            if startsWith(stimulus_name, 'temp_')
                stimulus_onset = timestamps(2);
            else
                stimulus_onset = timestamps(1);
            end
            PARAMETERS.stimulus_profile = STIMULI.(stimulus_name).stimulus_profile;
            % Bin data to speed-up computations
            [data_binned, ~, bin_edges] = make_tensor(data, time_bin_activity_n_frames, stimulus_onset);
            PARAMETERS.trial_time = STIMULI.(stimulus_name).time_axis([bin_edges(1, 1), bin_edges(end, 2)]);
            n_bins = size(bin_edges, 1);
            data_size = size(data_binned);
            PARAMETERS.data_size = data_size([2, 3]);
            PARAMETERS.data = reshape(data_binned, n_cells, []).';
            stimulus_onset_bin = find(bin_edges(:,1) == stimulus_onset);
            timestamps_bin = NaN(length(timestamps), 1);
            for it = 1:length(timestamps)
                timestamps_bin(it) = find(timestamps(it) >= bin_edges(:,1) & timestamps(it) <= bin_edges(:,2));
            end
            PARAMETERS.timestamps = timestamps_bin;
            PARAMETERS.timestamps_s = STIMULI.(stimulus_name).timestamps;
            
            % Mark bins to cumulate with a unique id
            bins_to_cumulate = NaN(n_bins, 1);
            bins_to_cumulate(stimulus_onset_bin:end) = 1:length(bins_to_cumulate(stimulus_onset_bin:end));
            bins_to_cumulate(stimulus_onset_bin - baseline_window_n_bins : stimulus_onset_bin - 1) = 0;
            if startsWith(stimulus_name, 'temp_')
                bins_to_cumulate(bins_to_cumulate > 2) = NaN;
            else
                last_bin = stimulus_onset_bin + round(STIMULI.(stimulus_name).duration * frame_rate / time_bin_activity_n_frames) + evoked_window_n_bins;
                bins_to_cumulate(last_bin + 1:end) = NaN;
            end
            PARAMETERS.bins_to_cumulate = bins_to_cumulate;
            
            % Store data to disk
            temp_data_filename = [temp_folder, filesep(), fieldname, '_detection_data.mat'];
            save(temp_data_filename, 'PARAMETERS', '-v6')
            % Run detection analysis in python
            run_in_python(os.path.join('decoding', python_script), ['input_filename=', temp_data_filename], 'verbose=True')
            
            % Reformat results
            results = jsondecode(fileread(PARAMETERS.output_filename));
            p = results.performance;
            % If no cell responded, allocate an empty array to create an empty
            % table
            if isempty(p{2})
                p{2} = NaN(0, length(p{1}));
            end
            results.performance = array2table(p{2}, 'VariableNames',p{1}(:).');
            if implementation == 0
                % Add 1 to timestamps because they have been 0-indexed in python
                results.performance{:, 'start'} = results.performance{:, 'start'} + 1;
                results.performance{:, 'finish'} = results.performance{:, 'finish'} + 1;
                results.performance{:, 'peak_probability'} = results.performance{:, 'peak_probability'} + 1;
                
                % Convert "significant_trials" to a 3D matrix (cell x bin x trial)
                n_trials = PARAMETERS.data_size(2);
                n_bins = size(results.response_probability_traces, 2);
                
                significant_trials = zeros(n_cells, n_bins, n_trials);
                for i_cell = 1:n_cells
                    % Make sure output is a cell array
                    if ~iscell(results.significant_trials{i_cell})
                        if all(results.significant_trials{i_cell} == 0)
                            results.significant_trials{i_cell} = cell(n_bins, 1);
                        else
                            results.significant_trials{i_cell} = num2cell(results.significant_trials{i_cell});
                        end
                    end
                    % Mark trials that are significant
                    for i_bin = 1:n_bins
                        if isempty(results.significant_trials{i_cell}{i_bin})
                            continue
                        else
                            significant_trials(i_cell, i_bin, results.significant_trials{i_cell}{i_bin} + 1) = 1;
                        end
                    end
                end
                % Replace the cell array of cell arrays with a matrix
                results.significant_trials = significant_trials;
            end
            
            RESULTS.([GC.analysis.(GC.experiment_name).session_column_prefix, result_fieldname]).(stimulus_name) = results;
            % Update file on disk
            save(result_filename, 'RESULTS', '-v7.3')
            
            % Clean up
            delete(temp_data_filename)
            delete(PARAMETERS.output_filename)
        end
        
        if ~exist(temp_folder, 'dir')
            delete(temp_folder)
        end
        
    case 'epi'
        downsampling_window_width = GC.response_detection_window_width;
        n_permutations_significance = GC.detection_n_permutations_significance;
        baseline_window_n_bins = round((GC.detection_baseline_window/2) / downsampling_window_width);

        for i_cond = 1:n_conditions
            columns_to_read = GC.analysis.(experiment_type).group_trials_by_session;
            str_message = 'Processing ';
            for ifield = 1:length(columns_to_read)
                str_message = [str_message, columns_to_read{ifield}, ' ''', sessions_to_analyze{i_cond, ifield}, ''', '];
            end
            LOGGER.info(str_message(1:end-2))
            
            fieldname = '';
            for ifield = 1:length(columns_to_read)
                if strcmp(columns_to_read{ifield}, 'date')
                    fieldname = [fieldname, 'session', sessions_to_analyze{i_cond, ifield}, '_'];
                else
                    fieldname = [fieldname, sessions_to_analyze{i_cond, ifield}, '_'];
                end
            end
            fieldname = fieldname(1:end-1);
            stimulus_name = make_variable_name(sessions_to_analyze{i_cond, ismember(columns_to_read, 'stimulus')});
            if strcmp(stimulus_name, 'temp_42') || strcmp(stimulus_name, 'temp_43')
                keyboard;
            end
            result_fieldname = sessions_to_analyze{i_cond, ismember(columns_to_read, GC.analysis.(GC.experiment_name).session_column)};
            
            % Initialize structure containing data and parameters for analysis
            PARAMETERS = struct();
            PARAMETERS.frame_rate = frame_rate;
            PARAMETERS.n_frames_per_bin = time_bin_activity_n_frames;
            PARAMETERS.data = [];
            PARAMETERS.data_size = [];
            PARAMETERS.stimulus_profile = [];
            PARAMETERS.trial_time = [];
            PARAMETERS.timestamps = [];
            PARAMETERS.timestamps_s = [];
            PARAMETERS.stimulus_onset = [];
            PARAMETERS.bins_to_cumulate = [];
            PARAMETERS.evoked_activity_max_latency = 100;
            PARAMETERS.evoked_activity_max_latency_bins = round(PARAMETERS.evoked_activity_max_latency / downsampling_window_width);
            PARAMETERS.CV_k_folds = GC.detection_k_fold_crossvalidation;
            PARAMETERS.significance_test_n_permutations = n_permutations_significance;
            PARAMETERS.output_filename = [temp_folder, filesep(), fieldname, '_detection_results.txt'];
            
            % Get data
%             trial_idx = METADATA.keep;
%             for ifield = 1:length(columns_to_read)
%                 trial_idx = trial_idx & ismember(table2cell(METADATA(:, columns_to_read{ifield})), sessions_to_analyze{i_cond, ifield});
%             end
            METADATA_epi = SQL_database.read_epi_trials(animal_ID);
            trial_idx = ismember(METADATA_epi.date, sessions_to_analyze{i_cond}) & ismember(METADATA_epi.stimulus, sessions_to_analyze{i_cond, 2});
            data = spikes(:, trial_idx); 
            % Get timestamps
            if startsWith(stimulus_name, 'temp_')
                timestamps = round(STIMULI.temp_48.timestamps * frame_rate);
            else
                if strcmp(experiment_type, 'Aversion_puff')
                    timestamps = round(STIMULI.(stimulus_name).timestamps * frame_rate);
                else
                    timestamps = round(STIMULI.cold.timestamps * frame_rate); % this does not work when puff or touch are there for the pain experiments; so, use 'cold'
                end
            end
            if startsWith(stimulus_name, 'temp_')
                stimulus_onset = timestamps(2);
            else
                stimulus_onset = timestamps(1);
            end
            PARAMETERS.stimulus_profile = STIMULI.(stimulus_name).stimulus_profile;
            % Bin data to speed-up computations
            [data_binned, ~, bin_edges] = make_tensor(data, time_bin_activity_n_frames, stimulus_onset);
            PARAMETERS.trial_time = STIMULI.(stimulus_name).time_axis([bin_edges(1, 1), bin_edges(end, 2)]);
            n_bins = size(bin_edges, 1);
            data_size = size(data_binned);
            PARAMETERS.data_size = data_size([2, 3]);
            PARAMETERS.data = reshape(data_binned, n_cells, []).';
            stimulus_onset_bin = find(bin_edges(:,1) == stimulus_onset);
            timestamps_bin = NaN(length(timestamps), 1);
            for it = 1:length(timestamps)
                timestamps_bin(it) = find(timestamps(it) >= bin_edges(:,1) & timestamps(it) <= bin_edges(:,2));
            end
            PARAMETERS.timestamps = timestamps_bin;
            if strcmp(experiment_type, 'Aversion_puff')
                PARAMETERS.timestamps_s = STIMULI.(stimulus_name).timestamps; % Set it at 5 seconds for all the stim
            else
                PARAMETERS.timestamps_s = 5; % Set it at 5 seconds for all the stim
            end
            
            % Mark bins to cumulate with a unique id
            bins_to_cumulate = NaN(n_bins, 1);
            bins_to_cumulate(stimulus_onset_bin:end) = 1:length(bins_to_cumulate(stimulus_onset_bin:end));
            bins_to_cumulate(stimulus_onset_bin - stimulus_onset_bin + 1 : stimulus_onset_bin - 1) = 0;
            last_bin = stimulus_onset_bin + round(0 * frame_rate / time_bin_activity_n_frames) + evoked_window_n_bins;
            bins_to_cumulate(last_bin + 1:end) = NaN;
            PARAMETERS.bins_to_cumulate = bins_to_cumulate;
            
            % Store data to disk
            temp_data_filename = [temp_folder, filesep(), fieldname, '_detection_data.mat'];
            save(temp_data_filename, 'PARAMETERS', '-v6')
            % Run detection analysis in python
            run_in_python(os.path.join('decoding', python_script), ['input_filename=', temp_data_filename], 'verbose=True')
            
            % Reformat results
            results = jsondecode(fileread(PARAMETERS.output_filename));
            p = results.performance;
            % If no cell responded, allocate an empty array to create an empty
            % table
            if isempty(p{2}) 
                p{2} = NaN(0, length(p{1}));
            end
            results.performance = array2table(p{2}, 'VariableNames',p{1}(:).');
             if implementation == 0 || implementation == 1
                % Add 1 to timestamps because they have been 0-indexed in python
                results.performance{:, 'start'} = results.performance{:, 'start'} + 1;
                results.performance{:, 'finish'} = results.performance{:, 'finish'} + 1;
                results.performance{:, 'peak_probability'} = results.performance{:, 'peak_probability'} + 1;
                
                % Convert "significant_trials" to a 3D matrix (cell x bin x trial)
                n_trials = PARAMETERS.data_size(2);
                n_bins = size(results.response_probability_traces, 2);
                
                significant_trials = zeros(n_cells, n_bins, n_trials);
                for i_cell = 1:n_cells
                    % Make sure output is a cell array
                    if ~iscell(results.significant_trials{i_cell})
                        if all(results.significant_trials{i_cell} == 0)
                            results.significant_trials{i_cell} = cell(n_bins, 1);
                        else
                            results.significant_trials{i_cell} = num2cell(results.significant_trials{i_cell});
                        end
                    end
                    % Mark trials that are significant
                    for i_bin = 1:n_bins
                        if isempty(results.significant_trials{i_cell}{i_bin})
                            continue
                        else
                            significant_trials(i_cell, i_bin, results.significant_trials{i_cell}{i_bin} + 1) = 1;
                        end
                    end
                end
                % Replace the cell array of cell arrays with a matrix
                results.significant_trials = significant_trials;
            end
            
            RESULTS.([GC.analysis.(GC.experiment_name).session_column_prefix, result_fieldname]).(stimulus_name) = results;
            % Update file on disk
            
            save(result_filename, 'RESULTS', '-v7.3')
            
            % Clean up
            delete(temp_data_filename)
            delete(PARAMETERS.output_filename)
        end
        
        if ~exist(temp_folder, 'dir')
            delete(temp_folder)
        end
        
end


%% MLint exceptions
%#ok<*PFBNS,*UNRCH,*PFOUS,*AGROW,*STRNU>
