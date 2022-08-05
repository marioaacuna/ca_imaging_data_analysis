%% notes: if experiment_name = CLA_pain, delete the last trial, line6 -> check

clear, clc
global GC
experiment_name =  'CLA_pain';
if strcmp (experiment_name, 'CLA_pain')
    delete_last_trial= 1;
else
    delete_last_trial= 0;
end

GC.experiment_name = experiment_name;
cd(os.path.join(GC.repository_root_path, 'Code', 'Workflows', '+Analysis_workflow', 'private'))

downsampling_window_width = GC.response_detection_window_width;
baseline_window_n_bins = round(GC.detection_baseline_window / downsampling_window_width);
if strcmp(experiment_name, 'ACC_SNI_anxiety')
    addpath('C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_SNI_anxiety')
end
    
overwrite = 1;
all_animals = {'NRN_8b'};
data_type = SQL_database.read_table_where('experiments', 'data_type', all_animals,'animal_ID', 'return_as_table',false);
n_mice = length(all_animals);
%% 
for i_mouse = 1:n_mice
    animal_ID = all_animals{i_mouse};
    disp(animal_ID)
%     cells_to_take = [1,4, 28];
    GC.experiment_name = experiment_name;

    output_filename =os.path.join(GC.data_root_path, GC.plots_path, experiment_name, 'response_detection', [animal_ID, '.pdf']);
    if exist(output_filename, 'file') && ~overwrite
        disp('Already processed. Skipped')
        keyboard
    end
    
    Calcium_indicator = SQL_database.read_table_where('experiments', 'Calcium_indicator', animal_ID,'animal_ID', 'return_as_table',false);
    tau_decay = GC.decay_time_constant_calcium_reporter.(Calcium_indicator);
    tau_rise = GC.rise_time_constant_calcium_reporter.(Calcium_indicator);
    
    % Convert results
    results_filename = get_filename_of('response_detection_p_0_01', animal_ID, GC.experiment_name);
    RESULTS_original = load_variable(results_filename, 'RESULTS');
    all_stimuli = [];
    sessions = fieldnames(RESULTS_original);
        
    for i_sess = 1:length(sessions)
        session_name = sessions{i_sess};
        stimuli = fieldnames(RESULTS_original.(session_name));
        if strcmp(experiment_name, 'ACC_SNI_anxiety')
            stimuli = natsort(stimuli);
        end
        all_stimuli = [all_stimuli; stimuli(:)];
        for i_stim = 1:length(stimuli)
            stimulus_name = stimuli{i_stim};
            performance = RESULTS_original.(session_name).(stimulus_name).performance;
            performance_columns = performance.Properties.VariableNames;
            performance = table2cell(performance);
            response_probability_traces = RESULTS_original.(session_name).(stimulus_name).response_probability_traces;
            significant_trials = RESULTS_original.(session_name).(stimulus_name).significant_trials;
            if delete_last_trial
                significant_trials = significant_trials(:,:,1:end-1);
            end
            RESULTS.(session_name).(stimulus_name) = struct();
            RESULTS.(session_name).(stimulus_name).performance_columns = performance_columns;
            RESULTS.(session_name).(stimulus_name).performance = performance;
            RESULTS.(session_name).(stimulus_name).response_probability_traces = response_probability_traces;
            RESULTS.(session_name).(stimulus_name).significant_trials = significant_trials;
        end
    end
    
    results_filename = os.path.join(GC.temp_dir, 'response_detection_results.mat');
    save(results_filename, 'RESULTS', '-v6')
    disp('Saved RESULTS')
    
    all_stimuli = unique(all_stimuli);
    n_stimuli = length(all_stimuli);
    
    
    %% Make parameters and compute denoised traces
    % Load spikes
    spikes_filename = get_filename_of('spikes', animal_ID);
    spikes = load_variable(spikes_filename, 'spikes');
    
    try
        Ca_events_filename = get_filename_of('Ca_events', animal_ID);
        keep_cells = load_variable(Ca_events_filename, 'keep_cell');
        keep_cells = logical(keep_cells);
        spikes= (spikes(keep_cells,:));
    catch
        spikes = (spikes);
    end
    
    n_cells = size(spikes, 1);

    switch data_type
        case '2p'
            % Load trials table
            METADATA = SQL_database.read_table_where('trials', {'+', 'date','stimulus','experimental_condition'}, {animal_ID, experiment_name}, {'animal_ID', 'experiment_name'});
            % Add column to mark unwanted trials
            METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
            % Remove trials of non-allowed stimuli
            all_stim_metadata = cell(n_stimuli,1);
            for i_stim = 1:n_stimuli
                stimulus_name = all_stimuli(i_stim);
                all_stim_metadata(i_stim) =  make_variable_name(stimulus_name, true, true);
            end
            METADATA.keep(~ismember(METADATA{:, 'stimulus'}, all_stim_metadata)) = false;
            % if sum(METADATA.keep) < 1, continue, end
            % Mark sessions where a compound was administered
            columns_to_read = {'date', 'stimulus'};%GC.session_name_prefix.(GC.experiment_name);
            METADATA = METADATA(:, [{'compound', 'frame_rate', 'n_frames', 'keep'}, columns_to_read]);
            % Get other info
            frame_rate = max(METADATA.frame_rate);
            time_bin_activity_n_frames = round(downsampling_window_width * frame_rate);
            evoked_window_n_bins = round(GC.evoked_activity_max_latency * frame_rate / time_bin_activity_n_frames);
            
        case 'epi'
            % Check that all stimuli of all sessions have been analyzed
            METADATA = SQL_database.read_table_where('sessions', {'date', 'stimulus', 'experimental_condition'}, {animal_ID, experiment_name}, {'animal_ID', 'experiment_name'});
            columns_to_read = GC.analysis.(experiment_name).group_trials_by_session;%GC.session_name_prefix.(GC.experiment_name);
            stimulus_column = find(ismember(columns_to_read, 'stimulus'));
            compound = GC.analysis.ACC_SNI_anxiety.analysis_allowed_compounds;
            compound = repmat((compound), height(METADATA),1);
            METADATA.compound = compound;
            METADATA = METADATA(:, [{'compound', 'day_from_surgery'}, columns_to_read]);
            % Remove sessions where a compound was administered
            METADATA(:, 'keep') = {true};
            METADATA(~ismember(METADATA.compound, compound), 'keep') = {false};
            METADATA(~METADATA.keep, :) = [];
            METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
            % Remove trials of non-allowed stimuli
            stimuli_this_dataset = unique(METADATA.stimulus);
            % Remove SPs
            stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, {'SP', 'SP_2', 'SPn', 'SP_2n', 'SPnl', 'SP_2nl'}), :);
            % Remove other stimuli that cannot be analyzed by this algorithm
            stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, GC.detection_no_timestamp_stimuli));
            METADATA.keep(~ismember(METADATA{:, 'stimulus'}, stimuli_this_dataset)) = false;
            % Remove SP and SP_2
            METADATA_epi = SQL_database.read_epi_trials(animal_ID);
%             analyzed_sessions = unique(METADATA_epi.date(ismemberCellRows(METADATA_epi.type, {'evoked'})), 'stable') ;
            frame_rate = GC.epifluorescence_downsample_to_frame_rate;
            frame_rate = repmat(frame_rate, height(METADATA),1);
            METADATA.frame_rate = frame_rate;
            frame_rate = max(METADATA.frame_rate);
            n_frames = size(spikes{1},2);
            n_frames = repmat(n_frames, height(METADATA), 1);
            METADATA.n_frames = n_frames;
            % Mark sessions where a compound was administered
            METADATA = METADATA(:, [{'compound', 'frame_rate', 'n_frames', 'keep','day_from_surgery' }, columns_to_read]);
            % Sessions to analyze
            sessions_to_analyze = table2cell(unique(METADATA(:, columns_to_read), 'rows'));
            analyzed_sessions = unique(METADATA_epi.date(ismemberCellRows(METADATA_epi.type, {'evoked'})), 'stable') ;
            % Remove other stimuli that cannot be analyzed by this algorithm
            sessions_to_analyze = sessions_to_analyze(ismember(sessions_to_analyze,analyzed_sessions),:);
            sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), {'SP', 'SP_2', 'SPnl', 'SP_2nl', 'SPn', 'SP_2n'}), :);
            sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), GC.detection_no_timestamp_stimuli), :);
            time_bin_activity_n_frames = round(downsampling_window_width * frame_rate);
%             GC.evoked_activity_max_latency = 10;
            evoked_window_n_bins = round(GC.evoked_activity_max_latency * frame_rate / time_bin_activity_n_frames);
    end
    

    % Load stimuli
    STIMULI = struct();
    for i_stim = 1:n_stimuli
        stimulus_name = all_stim_metadata{i_stim};
        stimulus_name_metadata = make_variable_name(stimulus_name, true, true);
        n_frames = unique(METADATA{ismember(table2cell(METADATA(:, 'stimulus')), stimulus_name_metadata), 'n_frames'});
        n_frames = n_frames(1);
        fieldname = make_variable_name(all_stim_metadata{i_stim});
        STIMULI.(fieldname) = Metadata_workflow.load_stimuli(all_stim_metadata{i_stim}, 'frame_rate',frame_rate, 'n_frames',n_frames, 'ignore_if_missing',true);
    end
    
    PARAMETERS = struct();
    PARAMETERS.figure_filename = output_filename;
    DENOISED_TRACES = struct();
    
    for i_sess = 1:length(sessions)
        session_name = sessions{i_sess};
        session_date = strsplit(session_name, '_');
        session_date = session_date{2};
        
        PARAMETERS.(session_name) = struct();
        for i_stim =1:n_stimuli
            stimulus_name = all_stimuli{i_stim};
            stimulus_name_metadata = make_variable_name(stimulus_name, true, true);
            
            % Get indices of trials
            switch data_type
                case '2p'
                    trial_idx = METADATA.keep & ...
                        ismember(METADATA.date, session_date) & ...
                        ismember(METADATA.stimulus, stimulus_name_metadata);
                case 'epi'
                    METADATA_epi = SQL_database.read_epi_trials(animal_ID);
                    sessions_to_analyze_ = unique(sessions_to_analyze(:,1));
                    trial_idx = ismember(METADATA_epi.date, sessions_to_analyze_{i_sess}) & ismember(METADATA_epi.stimulus, all_stimuli{i_stim});
                    
            end
            if ~any(trial_idx)
                disp([stimulus_name, ' not present in session ', session_date])
                continue
            end
            if strcmp(Calcium_indicator, 'GCaMP7f')
                GC.evoked_activity_max_latency = 20;
            end
            % Get data
            PARAMETERS.(session_name).(stimulus_name) = struct();
            PARAMETERS.(session_name).(stimulus_name).frame_rate = frame_rate;
            PARAMETERS.(session_name).(stimulus_name).n_frames_per_bin = time_bin_activity_n_frames;
            PARAMETERS.(session_name).(stimulus_name).evoked_activity_max_latency = GC.evoked_activity_max_latency;
            PARAMETERS.(session_name).(stimulus_name).evoked_activity_max_latency_bins = round(GC.evoked_activity_max_latency / downsampling_window_width);
            data = spikes(:, trial_idx);
            if delete_last_trial
                data = data(:,1:end-1);
            end
            % Compute denoised traces
            switch Calcium_indicator
                case 'GCaMP6f'
                    DENOISED_TRACES.(session_name).(stimulus_name) = convolve_spikes(data, frame_rate, tau_decay);
                case 'GCaMP7f'
                    DENOISED_TRACES.(session_name).(stimulus_name) = convolve_spikes(data, frame_rate, tau_decay, tau_rise);
            end
            % Get timestamps
            if strcmp(experiment_name, 'ACC_SNI_anxiety')
                trial_duration = GC.miniscope_trial_duration.(all_stimuli{i_stim});
                timestamps = trial_duration(1) * frame_rate+1;
            else
                timestamps = round(STIMULI.(stimulus_name).timestamps * frame_rate);
            end
            stimulus_onset = timestamps(1);
            PARAMETERS.(session_name).(stimulus_name).stimulus_profile = STIMULI.(stimulus_name).stimulus_profile;
            % Bin data to speed-up computations
            [data, ~, bin_edges] = make_tensor(data, time_bin_activity_n_frames, stimulus_onset);
            PARAMETERS.(session_name).(stimulus_name).trial_time = STIMULI.(stimulus_name).time_axis([bin_edges(1, 1), bin_edges(end, 2)]);
            n_bins = size(bin_edges, 1);
            data_size = size(data);
            PARAMETERS.(session_name).(stimulus_name).data_size = data_size([2, 3]);
            PARAMETERS.(session_name).(stimulus_name).data = reshape(data, n_cells, []).';
            stimulus_onset_bin = find(bin_edges(:,1) == stimulus_onset);
            timestamps_bin = NaN(length(timestamps), 1);
            for it = 1:length(timestamps)
                timestamps_bin(it) = find(timestamps(it) >= bin_edges(:,1) & timestamps(it) <= bin_edges(:,2));
            end
            PARAMETERS.(session_name).(stimulus_name).stimulus_onset = stimulus_onset_bin;
            PARAMETERS.(session_name).(stimulus_name).timestamps = timestamps_bin;
            if strcmp(experiment_name, 'ACC_SNI_anxiety')
                PARAMETERS.(session_name).(stimulus_name).timestamps_s = STIMULI.cold.timestamps;
            else
                PARAMETERS.(session_name).(stimulus_name).timestamps_s = STIMULI.(stimulus_name).timestamps;
            end
                
            
            % Mark bins to cumulate with a unique id
            switch data_type
                case '2p'
                    bins_to_cumulate = NaN(1, n_bins);
                    bins_to_cumulate(stimulus_onset_bin:end) = 1:length(bins_to_cumulate(stimulus_onset_bin:end));
                    bins_to_cumulate(stimulus_onset_bin-stimulus_onset_bin+1:stimulus_onset_bin-1) = 0;  % the baseline window
                    if startsWith(stimulus_name, 'temp_')
                        bins_to_cumulate(bins_to_cumulate > 2) = NaN;
                    else
                        last_bin = stimulus_onset_bin + round(STIMULI.(stimulus_name).duration * frame_rate / time_bin_activity_n_frames) + evoked_window_n_bins;
                        bins_to_cumulate(last_bin + 1:end) = NaN;
                    end
                    
                case 'epi'
                    bins_to_cumulate = NaN(n_bins, 1);
                    bins_to_cumulate(stimulus_onset_bin:end) = 1:length(bins_to_cumulate(stimulus_onset_bin:end));
                    bins_to_cumulate(stimulus_onset_bin - stimulus_onset_bin + 1 : stimulus_onset_bin - 1) = 0;
                    last_bin = stimulus_onset_bin + round(0 * frame_rate / time_bin_activity_n_frames) + evoked_window_n_bins;
                    bins_to_cumulate(last_bin + 1:end) = NaN;
            end
            bins_to_cumulate = bins_to_cumulate(:);
            PARAMETERS.(session_name).(stimulus_name).bins_to_cumulate = bins_to_cumulate;
        end
    end
    
    parameters_filename = os.path.join(GC.temp_dir, 'response_detection_parameters.mat');
    save(parameters_filename, 'PARAMETERS', '-v6')
    disp('Saved PARAMETERS')
    
    traces_filename = os.path.join(GC.temp_dir, 'response_detection_traces.mat');
    save(traces_filename, 'DENOISED_TRACES', '-v6')
    disp('Saved DENOISED_TRACES')
    
    
    %% Make figures
    run_in_python('detect_activations_plot.py', ['parameters_filename=', parameters_filename], ['result_filename=', results_filename], ['traces_filename=', traces_filename], 'verbose=True')
    
end


% Clean up and finish
delete(parameters_filename)
% delete(result_filename)
delete(traces_filename)
disp('Finished')


%#ok<*AGROW>
