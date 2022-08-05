%% PREAMBLE
% This function launches routines to analyze calcium traces or deconvolved
% spikes recorded during the experiments labeled as "ACC_CCI_anesth" and 'ACC_SNI_anxiety'.

%% MAIN
function done = main(INFO)
done = true;  % This is the output of the workflow, which is read by the launcher GUI and used to continue to the next workflow

% Get logger and general_configs
global LOGGER GC

% Unpack INFO to make all variables of this script explicit
all_animal_names = natsort(INFO.experiments);
% Unpack switches to decide what analyses to run
do_trace_extraction                                            = INFO.actions.analysis.ca_extraction.main;
detect_response_modulation                                     = INFO.actions.analysis.response_modulation.main;
response_detection_missing_only                                = INFO.actions.analysis.response_modulation.missing;
compute_response_decoding                                      = INFO.actions.analysis.response_decoding.main;
response_decoding_missing_only                                 = INFO.actions.analysis.response_decoding.missing;
compute_response_decoding_on_all_cells                         = INFO.actions.analysis.response_decoding.all;
compute_response_decoding_on_selective_cells                   = INFO.actions.analysis.response_decoding.selective;
compute_response_decoding_on_stable_cells                      = INFO.actions.analysis.response_decoding.stable;
compute_response_decoding_on_ensemble_cells                    = INFO.actions.analysis.response_decoding.ensembles;
compute_response_decoding_on_specific_cells                    = INFO.actions.analysis.response_decoding.specific;
compute_response_decoding_on_H_noci_cells                      = INFO.actions.analysis.response_decoding.H_noci;
compute_response_decoding_on_H_noci_cells_random               = INFO.actions.analysis.response_decoding.H_noci_random;
compute_response_decoding_on_H_sal_cells                       = INFO.actions.analysis.response_decoding.H_sal;
compute_response_decoding_on_H_sal_cells_random                = INFO.actions.analysis.response_decoding.H_sal_random;
compute_response_decoding_on_H_sal_cells_categories            = INFO.actions.analysis.response_decoding.H_sal_categories;
compute_response_decoding_on_H_sal_cells_random_categories     = INFO.actions.analysis.response_decoding.H_sal_random_categories;
compute_response_decoding_on_all_cells_categories              = INFO.actions.analysis.response_decoding.all_cell_categories;

compute_spontaneous_activity                                   = INFO.actions.analysis.spontaneous_activity.main;
spontaneous_activity_missing_only                              = INFO.actions.analysis.spontaneous_activity.missing;
compute_ensembles                                              = INFO.actions.analysis.ensemble_detection.main;
ensemble_detection_missing_only                                = INFO.actions.analysis.ensemble_detection.missing;
analyze_EPM_track                                              = INFO.actions.analysis.EPM_track.main;
EPM_missing_only                                               = INFO.actions.analysis.EPM_track.missing;
% Check what to do here
n_actions_to_perform = sum([do_trace_extraction*2, ...
                            detect_response_modulation, ...
                            compute_response_decoding_on_all_cells, ...
                            compute_response_decoding_on_selective_cells, ...
                            compute_response_decoding_on_stable_cells, ...
                            compute_response_decoding_on_ensemble_cells, ...
                            compute_response_decoding_on_specific_cells, ...
                            compute_response_decoding_on_H_noci_cells, ...
                            compute_response_decoding_on_H_noci_cells_random, ...                            
                            compute_response_decoding_on_H_sal_cells, ...   
                            compute_response_decoding_on_H_sal_cells_random, ...
                            compute_response_decoding_on_H_sal_cells_categories,...
                            compute_response_decoding_on_H_sal_cells_random_categories,...
                            compute_response_decoding_on_all_cells_categories,...
                            compute_spontaneous_activity, ...
                            compute_ensembles, ...
                            analyze_EPM_track, ...
                            ]);
if n_actions_to_perform < 1, return, end
% Initialize counter
current_action = 1;

% List allowed stimuli and compounds for the analysis
allowed_stimuli = GC.analysis.(INFO.selected_experiment).analysis_allowed_stimuli;
allowed_compounds = GC.analysis.(INFO.selected_experiment).analysis_allowed_compounds;


%% CALCIUM TRACES ANALYSIS
if do_trace_extraction
    % Enable toolboxes used by this workflow
    toolboxes_to_use = {'Toolbox_Romano', 'OASIS_matlab'};
    toggle_toolbox(toolboxes_to_use, 'on')

    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Computing ', char(916), 'F/F0'])
    for iid = 1:length(all_animal_names)
        % Get codename of this animal
        animal_ID = all_animal_names{iid};
        LOGGER.info(['Processing ', animal_ID])
        % Check if CNMFE was done
        neuron_filename =os.path.join(GC.data_root_path,GC.segmentation_path, [animal_ID, '_Neuron_Source2D.mat']);
        neuron_exist = exist(neuron_filename, 'file');
        % Extract calcium traces
        ca_trace_extraction(animal_ID, neuron_exist)
    end
    current_action = current_action + 1;

    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Deconvolving spiking events'])
    for iid = 1:length(all_animal_names)
        % Get codename of this animal
        animal_ID = all_animal_names{iid};
        LOGGER.info(['Processing ', animal_ID])
        % Get filename of output
        inferred_spikes_filename = get_filename_of('spikes', animal_ID);
        % Get filenames of input
        dFF_file = get_filename_of('dFF', animal_ID);
        % Load data
        dFF = load_variable(dFF_file, 'dFF');

        % Get trial frames indices
        data_type = SQL_database.read_table_where('experiments', 'data_type', animal_ID,'animal_ID', 'return_as_table',false);
        switch data_type
            case '2p'
                if any(strcmp(animal_ID, {'NRN_15a', 'NRN_15b' ,'NRN_16a', 'NRN_16b'}))
                    frames_idx = get_trial_frames(animal_ID, 'reset_index',true, 'return_all_sessions',1);
                else
                    frames_idx = get_trial_frames(animal_ID, 'reset_index',true, 'return_all_sessions', false);
                end
                frame_rate = SQL_database.read_table_where('trials', 'frame_rate', animal_ID, 'animal_ID', 'return_as_table',false);
                frame_rate = frame_rate{1};
                
            case 'epi'
                trial_duration = cellfun(@(x) length(x), dFF(1, :))';
                frames_idx = cumsum(trial_duration);
                frames_idx = [[1; frames_idx(1:end-1) + 1], frames_idx];
                frame_rate = GC.epifluorescence_downsample_to_frame_rate;
        end
        Ca_indicator = SQL_database.read_table_where('experiments', {'Calcium_indicator'}, animal_ID, 'animal_ID', 'return_as_table', false);      
        n_ROIs = size(dFF,1);
        % Deconcolve traces and infer spikes
        spikes = cell(size(dFF));
        Ca_events = cell(n_ROIs, 1);
        
        % Infer spikes on the whole trace
        % Load Neuron source2D
        % Check if CNMFE was done
        neuron_filename =os.path.join(GC.data_root_path,GC.segmentation_path, [animal_ID, '_Neuron_Source2D.mat']);
        neuron_exist = exist(neuron_filename, 'file');
        if neuron_exist
            neuron = load_variable(os.path.join(GC.data_root_path,GC.segmentation_path, [animal_ID, '_Neuron_Source2D.mat']));
            SPIKES = neuron.S;
        else
            SPIKES = ones(n_ROIs,1);
        end
        
        for iroi = 1:n_ROIs
            % Log progress
            LOGGER.trace(['ROI ', num2str(iroi), ' (out of ', num2str(n_ROIs), ')'], 'overwrite',iroi>1)
            % Get data of this ROI
            dFF_ROI = cell2mat(dFF(iroi, :)).';
            % Run spiking events deconvolution and compute event dynamics
            [spikes_ROI, Ca_events_ROI] = deconvolve_spiking_events(animal_ID, dFF_ROI, frame_rate, frames_idx, SPIKES(iroi,:));
            spikes(iroi,:) = spikes_ROI;
            Ca_events{iroi} = Ca_events_ROI;
        end
        % Store data on disk
        LOGGER.trace('Writing inferred spikes to disk')
        save(inferred_spikes_filename, 'spikes','Ca_events', '-v7.3')
        
        % Clean-up temporary files
        files = dir([os.path.join(GC.temp_dir, animal_ID), '_*']);
        files = [[animal_ID, '.mat']; {files.name}'];
        files = strcat(GC.temp_dir, files);
        if GC.preprocessing_delete_temp_file
            LOGGER.info('Deleting local temporary files')
            for ii = 1:length(files)
                f = files{ii};
                if exist(f, 'file')  % Make sure file exists before trying to delete it
                    delete(f)
                end
            end
        end
    end
    % Increase counter
    current_action = current_action + 1;

    % Old implementation
%     LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Detecting cell assemblies from the spontaneous activity'])
%     for iid = 1:length(all_animal_names)
%         animal_ID = all_animal_names{iid};
%         LOGGER.info(['Processing ', animal_ID])
%         
%         % Detect assemblies
%         filename_assemblies = get_filename_of('assemblies', animal_ID);
%         if exist(filename_assemblies, 'file')
%             ASSEMBLIES = load_variable(filename_assemblies, 'ASSEMBLIES');
%         else
%             ASSEMBLIES = struct();
%         end
%         % Load spike data
%         filename = get_filename_of('spikes', animal_ID);
%         spikes = load_variable(filename, 'spikes');
%         % Read metadata
%         METADATA = SQL_database.read_table_where('trials', {'date','stimulus'}, animal_ID, 'animal_ID');
%         [sessions, ~, session_idx] = unique(METADATA.date);
%         n_sessions = length(sessions);
%         for i_sess = 1:n_sessions
%             trial_idx = session_idx == i_sess & ismember(METADATA.stimulus, 'SP');
%             data = spikes(:, trial_idx);
%             
%             cond_fieldname = ['cond_', num2str(i_sess)];
%             ASSEMBLIES.(cond_fieldname) = estimate_cell_assemblies_PCApromax(data);
%         end
%         save(filename_assemblies, 'ASSEMBLIES', '-v7.3')
%     end
    
    % Disable toolboxes used by this workflow
    toggle_toolbox(toolboxes_to_use, 'off')

    % Increase counter
    current_action = current_action + 1;
end


%% ANALYSES
if detect_response_modulation || compute_spontaneous_activity
    % Log beginning of workflow
    LOGGER.info('Single-cell analyses workflow', 'decorate',true)
    
    %% ANALYZE RESPONSE MODULATION
    % --- Response detection ---
    if detect_response_modulation
        % Log action
        LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Detecting activity modulation'])
        experiment_type =  INFO.selected_experiment;
        data_type = unique(SQL_database.read_table_where('experiments', {'data_type'}, INFO.experiments, 'animal_ID', 'return_as_table', false), 'stable');
        if iscell(data_type), data_type = cell2mat(data_type);end
%         data_type = cell2mat(unique(data_type));
        switch data_type
%             case {'ACC_CCI_anesth', 'CLA_pain', 'Anaesth_check', 'CLA_ACC_dreadd'}
            case '2p'
                stats_filename = get_filename_of('spontaneous_activity_stats', experiment_type);
                
                for iid = 1:length(all_animal_names)
                    % Read previous analysis for this animal
                    animal_ID = all_animal_names{iid};
                    % Get filename of inputs
                    inferred_spikes_filename = get_filename_of('spikes', animal_ID);
                    if ~exist(inferred_spikes_filename, 'file')
                        LOGGER.info(['''', animal_ID, ''' was not fully processed. Skipping'])
                        continue
                    end
                    logger_msg = '';
                    to_skip = false;
                    
                    % Check that all stimuli of all sessions have been analyzed
                    METADATA = SQL_database.read_table_where('sessions', {'date', 'stimulus', 'experimental_condition'}, {animal_ID, experiment_type}, {'animal_ID', 'experiment_name'});
                    columns_to_read = GC.analysis.(experiment_type).group_trials_by_session;%GC.session_name_prefix.(GC.experiment_name);
                    stimulus_column = find(ismember(columns_to_read, 'stimulus'));
                    METADATA = METADATA(:, [{'compound'}, columns_to_read]);
                    % Remove sessions where a compound was administered
                    METADATA(:, 'keep') = {true};
                    METADATA(~ismember(METADATA.compound, allowed_compounds), 'keep') = {false};
                    METADATA(~METADATA.keep, :) = [];
                    
                    filename = get_filename_of('response_detection_p_0_01', animal_ID, experiment_type);
                    if exist(filename, 'file')
                        RESULTS = load_variable(filename, 'RESULTS');
                        already_analyzed = fieldnames(RESULTS);
                        session_column = columns_to_read{ismember(columns_to_read, GC.analysis.(experiment_type).session_column)};
                        session_column_prefix = GC.analysis.(experiment_type).session_column_prefix;
                        sessions = unique(METADATA{:, session_column});
                        % Either take only sessions left to analyze or not
                        if response_detection_missing_only
                            sessions_to_analyze = cell(0, length(columns_to_read));
                            for i_session = 1:length(sessions)
                                % Get stimuli to analyze
                                all_stimuli = METADATA{ismember(METADATA{:, session_column}, sessions{i_session}), 'stimulus'};
                                % Remove SP
                                all_stimuli = all_stimuli(~ismember(all_stimuli, 'SP'));
                                % Remove other stimuli that cannot be analyzed by this algorithm
                                all_stimuli = all_stimuli(~ismember(all_stimuli, GC.detection_no_timestamp_stimuli));
                                if isempty(all_stimuli)
                                    continue
                                end
                                
                                % Get name of session
                                session_name = [session_column_prefix, sessions{i_session}];
                                % If session hasn't been analyzed, add all stimuli
                                if ~ismember(already_analyzed, session_name)
                                    stimuli_to_analyze = all_stimuli;
                                else
                                    analyzed_stimuli = fieldnames(RESULTS.(session_name));
                                    stimuli_to_analyze = all_stimuli(~ismember(all_stimuli, analyzed_stimuli));
                                end
                                % Mark the stimuli left to analyze
                                if ~isempty(stimuli_to_analyze)
                                    switch GC.experiment_name
                                        case {'ACC_CCI_anesth', 'CLA_pain', 'Anaesth_check'}
                                            sessions_to_analyze = [sessions_to_analyze; [repmat(sessions(i_session), length(stimuli_to_analyze), 1), stimuli_to_analyze(:)]];
                                            if strcmp(GC.experiment_name, 'Anaesth_check'), keyboard, end % Check for stimuli already analyzed
                                        case 'ACC_pain_LP-211'
                                            session_date = METADATA.date{1};
                                            sessions_to_analyze = [sessions_to_analyze; [repmat(session_date, length(stimuli_to_analyze), 1), stimuli_to_analyze(:), repmat(sessions(i_session), length(stimuli_to_analyze), 1)]];
%                                         case 'CLA_pain'
%                                             sessions_to_analyze = [sessions_to_analyze; [repmat(sessions(i_session), length(stimuli_to_analyze), 1), stimuli_to_analyze(:)]];
% %                                             session_date = METADATA.date{1};
% %                                             sessions_to_analyze = [sessions_to_analyze; [repmat(session_date, length(stimuli_to_analyze), 1), stimuli_to_analyze(:), repmat(sessions(i_session), length(stimuli_to_analyze), 1)]];
                                    end
                                end
                            end
                        else
                            sessions_to_analyze = table2cell(METADATA(:, columns_to_read));
                        end
                        
                        % Skip if analyzed all stimuli
                        if isempty(sessions_to_analyze)
                            LOGGER.info(['All done for ', animal_ID, '. Skipping'])
                            continue
                        else
                            logger_msg = ['Will continue previous analysis on ', animal_ID];
                        end
                        
                    else
                        logger_msg = ['Processing ', animal_ID];
                        sessions_to_analyze = table2cell(unique(METADATA(:, columns_to_read), 'rows'));
                    end
                    
                    % Remove SP and SP_2
                    sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), {'SP', 'SP_2', 'SPnl', 'SP_2nl', 'SPn', 'SP_2n'}), :);
                    % Remove other stimuli that cannot be analyzed by this algorithm
                    sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), GC.detection_no_timestamp_stimuli), :);
                    if isempty(sessions_to_analyze)
                        logger_msg = ['All done for ', animal_ID, '. Skipping'];
                        to_skip = true;
                    end
                    LOGGER.info(logger_msg)
                    if to_skip, continue, end
                    
                    % Load spikes
                    LOGGER.trace('Loading spikes')
                    spikes = load_variable(inferred_spikes_filename, 'spikes');
                    try
                        Ca_events_filename = get_filename_of('Ca_events', animal_ID);
                        keep_cells = load_variable(Ca_events_filename, 'keep_cell');
                        keep_cells = logical(keep_cells);
                        spikes= (spikes(keep_cells,:));
                    catch
                        spikes = (spikes);
                    end
                    
                    % Load trials table
                    METADATA = SQL_database.read_table_where('trials', {'+', 'date','stimulus', 'experimental_condition'}, {animal_ID, experiment_type}, {'animal_ID', 'experiment_name'});
                    stimuli_this_dataset = unique(METADATA.stimulus);
                    % Remove SPs
                    stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, {'SP', 'SP_2', 'SPn', 'SP_2n', 'SPnl', 'SP_2nl'}), :);
                    % Remove other stimuli that cannot be analyzed by this algorithm
                    stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, GC.detection_no_timestamp_stimuli));
                    % Add column to mark unwanted trials
                    METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
                    % Remove trials of non-allowed stimuli
                    METADATA.keep(~ismember(METADATA{:, 'stimulus'}, stimuli_this_dataset)) = false;
                    if sum(METADATA.keep) < 1, continue, end
                    % Mark sessions where a compound was administered
                    METADATA = METADATA(:, [{'compound', 'frame_rate', 'n_frames', 'keep','day_from_surgery' }, columns_to_read]);
                    % Get other info
                    frame_rate = max(METADATA.frame_rate);
%                     stimuli_this_dataset = unique(sessions_to_analyze(:, 2));
                    
                    % Detect response to common set of stimuli
                    n_stimuli = length(stimuli_this_dataset);
                    if n_stimuli > 0
                        LOGGER.info(['Analyzing responses to stimuli: [', strjoin(stimuli_this_dataset, ', '), ']'])
                        % Load stimuli
                        STIMULI = struct();
                        for i_stim = 1:n_stimuli
                            n_frames = METADATA{ismember(table2cell(METADATA(:, 'stimulus')), stimuli_this_dataset{i_stim}), 'n_frames'}; n_frames = n_frames(1);
                            fieldname = make_variable_name(stimuli_this_dataset{i_stim});
                            STIMULI.(fieldname) = Metadata_workflow.load_stimuli(stimuli_this_dataset{i_stim}, 'frame_rate',frame_rate, 'n_frames',n_frames, 'ignore_if_missing',true);
                        end
                        response_detection(experiment_type, animal_ID, spikes, sessions_to_analyze, METADATA, STIMULI, response_detection_missing_only);
%                         response_detection_z_score(animal_ID, spikes, sessions_to_analyze, METADATA, STIMULI, response_detection_missing_only);
                    end
                end
%             case {'ACC_SNI_anxiety', 'GBP_pain'}
            case 'epi'
                type_detection = 'permutation';
                %                 type_detection = 'z_score';
                switch type_detection
                    case 'permutation'
                        for iid = 1:length(all_animal_names)
                            % Read previous analysis for this animal
                            animal_ID = all_animal_names{iid};
                            % Get filename of inputs
                            inferred_spikes_filename = get_filename_of('spikes', animal_ID);
                            if ~exist(inferred_spikes_filename, 'file')
                                LOGGER.info(['''', animal_ID, ''' was not fully processed. Skipping'])
                                continue
                            end
                            logger_msg = '';
                            to_skip = false;
                            
                            % Check that all stimuli of all sessions have been analyzed
                            METADATA = SQL_database.read_table_where('sessions', {'date', 'stimulus', 'experimental_condition'}, {animal_ID, experiment_type}, {'animal_ID', 'experiment_name'});
                            columns_to_read = GC.analysis.(experiment_type).group_trials_by_session;%GC.session_name_prefix.(GC.experiment_name);
                            stimulus_column = find(ismember(columns_to_read, 'stimulus'));
                            compound = GC.analysis.ACC_SNI_anxiety.analysis_allowed_compounds;
                            compound = repmat((compound), height(METADATA),1);
                            METADATA.compound = compound;
                            METADATA = METADATA(:, [{'compound', 'day_from_surgery'}, columns_to_read]);
                            % Remove sessions where a compound was administered
                            METADATA(:, 'keep') = {true};
                            METADATA(~ismember(METADATA.compound, allowed_compounds), 'keep') = {false};
                            METADATA(~METADATA.keep, :) = [];
                            
                            filename = get_filename_of('response_detection_p_0_01', animal_ID, experiment_type);
                            if exist(filename, 'file')
                                RESULTS = load_variable(filename, 'RESULTS');
                                already_analyzed = fieldnames(RESULTS);
                                session_column = columns_to_read{ismember(columns_to_read, GC.analysis.(experiment_type).session_column)};
                                session_column_prefix = GC.analysis.(experiment_type).session_column_prefix;
                                sessions = unique(METADATA{:, session_column});
                                % Either take only sessions left to analyze or not
                                if response_detection_missing_only
                                    sessions_to_analyze = cell(0, length(columns_to_read));
                                    for i_session = 1:length(sessions)
                                        % Get stimuli to analyze
                                        all_stimuli = METADATA{ismember(METADATA{:, session_column}, sessions{i_session}), 'stimulus'};
                                        % Remove SP
                                        all_stimuli = all_stimuli(~ismember(all_stimuli, 'SP'));
                                        % Remove other stimuli that cannot be analyzed by this algorithm
                                        all_stimuli = all_stimuli(~ismember(all_stimuli, GC.detection_no_timestamp_stimuli));
                                        if isempty(all_stimuli)
                                            continue
                                        end
                                        
                                        % Get name of session
                                        session_name = [session_column_prefix, sessions{i_session}];
                                        % If session hasn't been analyzed, add all stimuli
                                        if ~ismember(already_analyzed, session_name)
                                            stimuli_to_analyze = all_stimuli;
                                        else
                                            analyzed_stimuli = fieldnames(RESULTS.(session_name));
                                            stimuli_to_analyze = all_stimuli(~ismember(all_stimuli, analyzed_stimuli));
                                        end
                                        % Mark the stimuli left to analyze
                                        if ~isempty(stimuli_to_analyze)
                                            switch GC.experiment_name
                                                case 'ACC_SNI_anxiety'
                                                    sessions_to_analyze = [sessions_to_analyze; [repmat(sessions(i_session), length(stimuli_to_analyze), 1), stimuli_to_analyze(:)]];
                                                case 'ACC_pain_LP-211'
                                                    session_date = METADATA.date{1};
                                                    sessions_to_analyze = [sessions_to_analyze; [repmat(session_date, length(stimuli_to_analyze), 1), stimuli_to_analyze(:), repmat(sessions(i_session), length(stimuli_to_analyze), 1)]];
                                                case 'CLA_pain'
                                                    sessions_to_analyze = [sessions_to_analyze; [repmat(sessions(i_session), length(stimuli_to_analyze), 1), stimuli_to_analyze(:)]];
                                                    %                                             session_date = METADATA.date{1};
                                                    %                                             sessions_to_analyze = [sessions_to_analyze; [repmat(session_date, length(stimuli_to_analyze), 1), stimuli_to_analyze(:), repmat(sessions(i_session), length(stimuli_to_analyze), 1)]];
                                            end
                                        end
                                    end
                                else
                                    sessions_to_analyze = table2cell(METADATA(:, columns_to_read));
                                end
                                
                                % Skip if analyzed all stimuli
                                if isempty(sessions_to_analyze)
                                    LOGGER.info(['All done for ', animal_ID, '. Skipping'])
                                    continue
                                else
                                    logger_msg = ['Will continue previous analysis on ', animal_ID];
                                end
                                
                            else
                                logger_msg = ['Processing ', animal_ID];
                                sessions_to_analyze = table2cell(unique(METADATA(:, columns_to_read), 'rows'));
                            end
                            
                            % Remove SP and SP_2
                            METADATA_epi = SQL_database.read_epi_trials(animal_ID);
                            analyzed_sessions = unique(METADATA_epi.date(ismemberCellRows(METADATA_epi.type, {'evoked'})), 'stable') ;
                            % Remove other stimuli that cannot be analyzed by this algorithm
                            sessions_to_analyze = sessions_to_analyze(ismember(sessions_to_analyze,analyzed_sessions),:);
                            sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), {'SP', 'SP_2', 'SPnl', 'SP_2nl', 'SPn', 'SP_2n'}), :);
                            sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), GC.detection_no_timestamp_stimuli), :);
                            if isempty(sessions_to_analyze)
                                logger_msg = ['All done for ', animal_ID, '. Skipping'];
                                to_skip = true;
                            end
                            LOGGER.info(logger_msg)
                            if to_skip, continue, end
                            
                            % Load spikes
                            LOGGER.trace('Loading spikes')
                            spikes = load_variable(inferred_spikes_filename, 'spikes');
                            
                            % Load trials table
                            p_filename = get_filename_of('miniscope_movie_parameters', animal_ID);
                            miniscope_parameters = load_variable(p_filename, 'p');
                            %                     METADATA = SQL_database.read_table_where('trials', {'+', 'date','stimulus', 'experimental_condition'}, {animal_ID, experiment_type}, {'animal_ID', 'experiment_name'});
                            stimuli_this_dataset = unique(METADATA.stimulus);
                            % Remove SPs
                            stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, {'SP', 'SP_2', 'SPn', 'SP_2n', 'SPnl', 'SP_2nl'}), :);
                            % Remove other stimuli that cannot be analyzed by this algorithm
                            stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, GC.detection_no_timestamp_stimuli));
                            % Add column to mark unwanted trials
                            METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
                            % Remove trials of non-allowed stimuli
                            METADATA.keep(~ismember(METADATA{:, 'stimulus'}, stimuli_this_dataset)) = false;
                            if sum(METADATA.keep) < 1, continue, end
                            % Add frame_rate, n_frames
                            frame_rate = GC.epifluorescence_downsample_to_frame_rate;
                            frame_rate = repmat(frame_rate, height(METADATA),1);
                            METADATA.frame_rate = frame_rate;
                            n_frames = size(spikes{1},2);
                            n_frames = repmat(n_frames, height(METADATA), 1);
                            METADATA.n_frames = n_frames;
                            % Mark sessions where a compound was administered
                            METADATA = METADATA(:, [{'compound', 'frame_rate', 'n_frames', 'keep','day_from_surgery' }, columns_to_read]);
                            % Get other info
                            frame_rate = max(METADATA.frame_rate);
                            %                     stimuli_this_dataset = unique(sessions_to_analyze(:, 2));
                            % Detect response to common set of stimuli
                            n_stimuli = length(stimuli_this_dataset);
                            if n_stimuli > 0
                                LOGGER.info(['Analyzing responses to stimuli: [', strjoin(stimuli_this_dataset, ', '), ']'])
                                % Load stimuli
                                STIMULI = struct();
                                for i_stim = 1:n_stimuli
                                    n_frames = METADATA{ismember(table2cell(METADATA(:, 'stimulus')), stimuli_this_dataset{i_stim}), 'n_frames'}; n_frames = n_frames(1);
                                    fieldname = make_variable_name(stimuli_this_dataset{i_stim});
                                    STIMULI.(fieldname) = Metadata_workflow.load_stimuli(stimuli_this_dataset{i_stim}, 'frame_rate',frame_rate, 'n_frames',n_frames, 'ignore_if_missing',true);
                                end
                                response_detection(experiment_type, animal_ID, spikes, sessions_to_analyze, METADATA, STIMULI, response_detection_missing_only);
                            end
                        end
                    case 'z_score'
                        for iid = 1:length(all_animal_names)
                            if response_detection_missing_only
                                % Read previous analysis for this animal
                                animal_ID = all_animal_names{iid};
                                % Get filename of inputs
                                traces_filename = get_filename_of('dFF', animal_ID);
                                if ~exist(traces_filename, 'file')
                                    LOGGER.info(['''', animal_ID, ''' was not fully processed. Skipping'])
                                    continue
                                end
                                logger_msg = '';
                                to_skip = false;
                                LOGGER.trace('Loading traces')
                                traces = load_variable(traces_filename, 'dFF');
                                response_detection_epi(animal_ID, traces, response_detection_missing_only)
                            end
                        end
                end
        end
        % Increase counter
        current_action = current_action + 1;
    end
    
    %% ANALYZE SPONTANEOUS ACTIVITY
    if compute_spontaneous_activity
        % Log action
        LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Analyze spontaneous activity'])
        % Make filename of output file
        experiment_type =  INFO.selected_experiment;
        switch experiment_type
            case {'ACC_CCI_anesth', 'CLA_pain'}        
                stats_filename = get_filename_of('spontaneous_activity_stats', experiment_type);
            case 'ACC_SNI_anxiety'
                stats_filename = get_filename_of('spontaneous_activity_stats_epi', experiment_type);
        end
        % Make folder if it doesn't exist
        [folder, ~] = fileparts(stats_filename);
        if ~exist(folder,'dir'), mkdir(folder), end
        % Load file
        if exist(stats_filename, 'file')
            switch  experiment_type
                case 'ACC_CCI_anesth'
                    STATS_spontaneous_activity = load_variable(stats_filename, 'STATS_spontaneous_activi_sessy');
                case 'ACC_SNI_anxiety'
                    STATS_spontaneous_activity = load_variable(stats_filename, 'EVENT_STATS_spontaneous_activi_sessy');
            end
        else
            STATS_spontaneous_activity = struct();
        end
        if spontaneous_activity_missing_only
            already_processed = fieldnames(STATS_spontaneous_activity);
            already_processed_names = already_processed(ismember(all_animal_names, already_processed));
            if ~isempty(already_processed_names)
                LOGGER.info(['Will skip: ', strjoin(already_processed_names, ', ')])
            end
            not_processed_names = all_animal_names(~ismember(all_animal_names, already_processed));
        else
            not_processed_names = all_animal_names;
        end
        if isempty(not_processed_names)
            LOGGER.info('Nothing left to analyze')
        else
            switch experiment_type
                case {'ACC_CCI_anesth', 'CLA_pain'}
                    spontaneous_activity_stats(not_processed_names, stats_filename)
                case 'ACC_SNI_anxiety'
                    if ~endsWith(all_animal_names{1},'epi')
                        keyboard
                    else
                        spontaneous_activity_stats_epi(not_processed_names, stats_filename)
                    end
            end
            
        end
        
        
        
        % Make filename of output file
        stats_filename = get_filename_of('spontaneous_activity_stats_epi', GC.experiment_name);
        % Make folder if it doesn't exist
        [folder, ~] = fileparts(stats_filename);
        if ~exist(folder,'dir'), mkdir(folder), end
        if exist(stats_filename, 'file')
            STATS_spontaneous_activity = load_variable(stats_filename, 'EVENT_STATS_spontaneous_activi_sessy');
        else
            STATS_spontaneous_activity = struct();
        end
        
        if spontaneous_activity_missing_only
            already_processed = fieldnames(STATS_spontaneous_activity);
            already_processed_names = already_processed(ismember(all_animal_names, already_processed));
            if ~isempty(already_processed_names)
                LOGGER.info(['Will skip: ', strjoin(already_processed_names, ', ')])
            end
            not_processed_names = all_animal_names(~ismember(all_animal_names, already_processed));
        else
            not_processed_names = all_animal_names;
        end
        
        if isempty(not_processed_names)
            LOGGER.info('Nothing left to analyze')
        else
            spontaneous_activity_stats_epi(not_processed_names, stats_filename)
        end
           
        
        
        
        
        
        % Increase counter
        current_action = current_action + 1;
    end
end


%% ENSEMBLES
if compute_ensembles
    % Enable toolboxes used by this workflow
    toolboxes_to_use = {'SVDEnsemble'};
    toggle_toolbox(toolboxes_to_use, 'on')
    
    % Take the stimuli to analyze
    stimuli_combinations = GC.assemblies_stimuli_to_combine;
    n_stimuli_combinations = length(stimuli_combinations);
    % Define names of columns of sessions table
    sessions_table_columns = {'date', 'day_from_surgery', 'stimuli', 'session_fieldname', 'stimuli_fieldname'};

    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Detecting cell ensembles (SVD implementation)'])
    for iid = 1:length(all_animal_names)
        animal_ID = all_animal_names{iid};
        LOGGER.info(['Processing ', animal_ID])
        
        % Detect assemblies
        filename_assemblies = get_filename_of('assemblies', animal_ID);
        if exist(filename_assemblies, 'file') && ensemble_detection_missing_only
            ASSEMBLIES = load_variable(filename_assemblies, 'ASSEMBLIES');
            if  ~all(ismember(cellfun(@cell2mat,stimuli_combinations, 'UniformOutput', false), cellfun(@cell2mat, ASSEMBLIES.info.stimuli, 'UniformOutput', false)))
                % Initialize output structure
                
                ASSEMBLIES.info.stimuli = stimuli_combinations(:);
                
            end
%             for i_stim_comb = 1:n_stimuli_combinations
%                 stimuli_to_analyze = stimuli_combinations{i_stim_comb};
%                 for i_stim = 1:length(ASSEMBLIES.info.stimuli)
%                     if ~all(ismember(ASSEMBLIES.info.stimuli{i_stim}, stimuli_to_analyze{1}))
%                         keyboard  % add combination
%                     end
%                 end
%             end
            
        else
            % Initialize output structure
            ASSEMBLIES = struct();
            ASSEMBLIES.info = struct();
            ASSEMBLIES.info.sessions = cell2table(cell(0, length(sessions_table_columns)), 'VariableNames',sessions_table_columns);
            ASSEMBLIES.info.stimuli = stimuli_combinations(:);
        end
        
        % Flag that spikes are not loaded yet
        is_data_loaded = false;
        
        % Read metadata
        METADATA = SQL_database.read_table_where('trials', {'+','date','experimental_condition','stimulus'}, animal_ID, 'animal_ID');
        session_info = unique(METADATA(:, {'date', 'day_from_surgery'}), 'rows');
        sessions = session_info.date;
        n_sessions = length(sessions);
        
        % Make table with new data
        sessions_to_analyze = [repmat(session_info, n_stimuli_combinations, 1), ...
            reshape(repmat(stimuli_combinations, n_sessions, 1), [], 1), ...
            strcat('cond_', value2str(reshape(repmat(1:n_sessions, n_stimuli_combinations, 1)', [], 1), '%i', true)), ...
            ];
        sessions_to_analyze.stimuli_fieldname(:) = {''};  % Don't know yet
        sessions_to_analyze.Properties.VariableNames = sessions_table_columns;
        n_sessions_to_analyze = height(sessions_to_analyze);
        
        % Remove combinations previously analyzed, if any
        if ensemble_detection_missing_only && ~isempty(ASSEMBLIES.info.sessions)
            sessions_already_analyzed = false(n_sessions_to_analyze, 1);
            for i_sess = 1:n_sessions_to_analyze
                % Start checking date and experimental condition
                match_rows = ismember(ASSEMBLIES.info.sessions{:, 'date'}, sessions_to_analyze{i_sess, 'date'}) & ...
                    ismember(ASSEMBLIES.info.sessions{:, 'day_from_surgery'}, sessions_to_analyze{i_sess, 'day_from_surgery'});
                
                if any(match_rows)
                    match_rows_idx = find(match_rows);
                    % Check whether stimuli are the same
                    for i_match = 1:length(match_rows_idx)
                        if all(ismember(ASSEMBLIES.info.sessions.stimuli{i_match}, sessions_to_analyze.stimuli{i_sess}))
                            sessions_already_analyzed(i_sess) = true;
                            break
                        end
                    end
                end
            end
            
            % Remove rows of sessions already analyzed
            sessions_to_analyze(sessions_already_analyzed, :) = [];
            if isempty(sessions_to_analyze)
                LOGGER.info('No more sessions left to analyze. Skipped')
                continue
            end
        end
%         n_sessions_to_analyze = height(sessions_to_analyze);
        n_sessions_to_analyze = length(unique(sessions_to_analyze.date, 'stable'));

        for i_sess = 1:n_sessions_to_analyze
            session_date = sessions_to_analyze.date{i_sess};
            session_idx = find(ismember(session_info.date, session_date));
            
            % Create fieldname for ASSEMBLIES structure
            session_fieldname = ['cond_', num2str(session_idx)];
            LOGGER.info(['Session ''', sessions{session_idx}, ''''])
            
            % Check whether session has been already analyzed
            if ~isfield(ASSEMBLIES, session_fieldname) || ~ensemble_detection_missing_only || length(fieldnames(ASSEMBLIES.(session_fieldname))) ~= n_stimuli_combinations
                if ~is_data_loaded
                    % Load spike data
                    LOGGER.info('Loading data')
                    spikes = load_variable(get_filename_of('spikes', animal_ID), 'spikes');
                    Ca_events_filename = get_filename_of('Ca_events', animal_ID);
                    if exist(Ca_events_filename, 'file')
                        keep_cell = logical(load_variable(Ca_events_filename, 'keep_cell'));
                        spikes = spikes(keep_cell, :);
                    end
                    is_data_loaded = true;
                end
                
                % Initialize structure field
                if ~isfield(ASSEMBLIES, session_fieldname)
                    ASSEMBLIES.(session_fieldname) = struct();
                end
            end
            
            % Loop through stimulus combinations
            for i_stim_comb = 1:n_stimuli_combinations
                % Select stimuli to analyze and make a name for the structure field
%                 if strcmp(animal_ID, 'MA_15'), keyboard, end
                stimuli_to_analyze = stimuli_combinations{i_stim_comb};
                
                for ii = 1:length(ASSEMBLIES.info.stimuli)
                    if all(ismember(stimuli_to_analyze, ASSEMBLIES.info.stimuli{ii}))
                        this_stim_combination = ii;
                        break
                    end
                end
                stimuli_fieldname = ['stim_', num2str(this_stim_combination)];
                
                % Check whether this combination has been already analyzed
                if ~isfield(ASSEMBLIES.(session_fieldname), stimuli_fieldname) || ~ensemble_detection_missing_only
                    LOGGER.info(['\tStimuli: ', strjoin(stimuli_to_analyze, ' - ')])
                else
                    LOGGER.info('\tAlready analyzed. Skipped')
                    continue
                end
                    
                % Perform analysis
                ASSEMBLIES.(session_fieldname).(stimuli_fieldname) = ...
                    population_vector_analysis(animal_ID,spikes, METADATA, sessions{session_idx}, stimuli_to_analyze, 'only_evoked',false);
                
                % Store information on session
                sessions_info = sessions_to_analyze(i_sess, :);
                sessions_info.stimuli_fieldname = {stimuli_fieldname};
                ASSEMBLIES.info.sessions = [ASSEMBLIES.info.sessions; sessions_info];
                % Save file
                save(filename_assemblies, 'ASSEMBLIES', '-v7.3')
            end
        end
    end
    
    % Disable toolboxes used by this workflow
    toggle_toolbox(toolboxes_to_use, 'off')

    % Increase counter
    current_action = current_action + 1;
end


%% DECODING
if compute_response_decoding
    % Initialize flags
    is_selectivity_stats_loaded = false;

    for analysis_type = 1:12  % all, selective, stable, ensembles, specific, etc.
        error_msg = '';
        switch analysis_type
            case 1
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_all_cells, continue, end
                analysis_type_str = 'all';
                analysis_type_msg = 'all cells';

            case 2
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_selective_cells, continue, end
                analysis_type_str = 'selective';
                analysis_type_msg = 'only selective cells';

            case 3
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_stable_cells, continue, end
                analysis_type_str = 'stable';
                analysis_type_msg = 'only stable cells';

            case 4
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_ensemble_cells, continue, end
                analysis_type_str = 'ensembles';
                analysis_type_msg = 'only ensemble cells';
                
            case 5
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_specific_cells, continue, end
                analysis_type_str = 'specific';
                analysis_type_msg = 'only specific cells';
                
            case 6 
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_H_noci_cells, continue, end
                analysis_type_str = 'H_noci';
                analysis_type_msg = 'only non-H_noci cells';
                
            case 7
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_H_noci_cells_random, continue, end
                analysis_type_str = 'H_noci_random';
                analysis_type_msg = 'only non-H_noci random cells';
                
            case 8
                if ~compute_response_decoding_on_H_sal_cells, continue, end
                analysis_type_str = 'H_sal';
                analysis_type_msg = 'only non-H_sal cells';
                
            case 9
                if ~compute_response_decoding_on_H_sal_cells_random, continue, end
                analysis_type_str = 'H_sal_random';
                analysis_type_msg = 'random non-H_sal cells';
            
            case 10
                if ~compute_response_decoding_on_H_sal_cells_categories, continue, end
                analysis_type_str = 'H_sal_categories';
                analysis_type_msg = 'only non-H_sal cells in categories';
                
            case 11
                if ~compute_response_decoding_on_H_sal_cells_random_categories, continue, end
                analysis_type_str = 'H_sal_random_categories';
                analysis_type_msg = 'random non-H_sal cells in categories';
                
            case 12
                % Check whether this step can be performed, and set its message
                if ~compute_response_decoding_on_all_cells_categories, continue, end
                analysis_type_str = 'all_categories';
                analysis_type_msg = 'all cells';

        end

        % Load info about selectivity
        if ismember(analysis_type, [2, 3, 5, 6, 7, 8, 9, 10,11]) && ~is_selectivity_stats_loaded
            try
                switch INFO.selected_experiment
                    case 'ACC_CCI_pain'
                        data_path_str = 'Figures_paper_FMP';
                    case 'ACC_SNI_anxiety'
                         data_path_str = 'Figures_paper_SNI_anxiety';
                    otherwise
                        keyboard
                end
                SELECTIVITY_stats = load_variable(os.path.join(GC.repository_root_path, data_path_str, '_data', 'SELECTIVITY_stats__p.mat'), 'SELECTIVITY_stats');
                is_selectivity_stats_loaded = true;
            catch
                error_msg = 'The file ''SELECTIVITY_stats__p.mat'' does not exist, probably because you haven''t identified selective cells yet. Skipping this analysis step';
            end
        end

        % Log action
        if isempty(error_msg)
            LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Decode response with ', analysis_type_msg])
        else
            LOGGER.error(error_msg)
            continue
        end

        % Get filename of output file
        stats_filename = get_filename_of(['response_decoding_stats_', analysis_type_str], INFO.selected_experiment);

        % Make folder if it doesn't exist
        [folder, ~] = fileparts(stats_filename);
        if ~exist(folder,'dir'), mkdir(folder), end
        % Load file
        if exist(stats_filename, 'file')
            STATS_response_decoding = load_variable(stats_filename, 'STATS_response_decoding');
        else
            STATS_response_decoding = struct();
        end

        for iid = 1:length(all_animal_names)
            % Get codename of this animal
            animal_ID = all_animal_names{iid};
            if ismember(analysis_type, [2, 3, 5, 6,7, 8, 9, 10,11]) && ~is_selectivity_stats_loaded
                if ~ismember(animal_ID, fieldnames(SELECTIVITY_stats))
                    LOGGER.warn([animal_ID, ' not fully analyzed'])
                    continue
                end
            end
             experiment_type =  INFO.selected_experiment;
             data_type = SQL_database.read_table_where('experiments', 'data_type', animal_ID,'animal_ID', 'return_as_table',false);
             switch data_type
                 case '2p'
                   % Read metadata from database
                    METADATA = SQL_database.read_table_where('trials', {'+','stimulus','experimental_condition','date'}, {animal_ID, GC.experiment_name}, {'animal_ID', 'experiment_name'});
                    stimuli_this_dataset = unique(METADATA.stimulus);
                    % Remove other stimuli that cannot be analyzed by this algorithm
                    stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, GC.detection_no_timestamp_stimuli));
                    stimuli_this_dataset(ismember(stimuli_this_dataset, {'SP', 'SP_2', 'HPS_short', 'HPS_50%', 'odor-'})) = [];
                    % Add column to mark unwanted trials
                    METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
                    % Remove trials of non-allowed stimuli
                    METADATA.keep(~ismember(METADATA{:, 'stimulus'}, stimuli_this_dataset)) = false;
                    % Remove trials where a compound had been administered
                    METADATA.keep(~ismember(METADATA.compound, allowed_compounds)) = false;
                    % Select sessions to analyze
                    sessions = unique(str2double(METADATA.day_from_surgery));
                    sessions_to_analyze_idx = 1:length(sessions);
                    sessions_to_analyze = sessions(sessions_to_analyze_idx);
                    n_sessions_to_analyze = length(sessions_to_analyze);
                    METADATA.keep(~ismember(str2double(METADATA.day_from_surgery), sessions_to_analyze)) = false;
                    if sum(METADATA.keep) < 1, continue, end
                    % Get frame rate
                    frame_rate = max(METADATA.frame_rate);
                    % Load stimuli
                    STIMULI = struct();
                    for istim = 1:length(stimuli_this_dataset)
                        this_stimulus = stimuli_this_dataset{istim};
                        
                        n_frames = METADATA{ismember(table2cell(METADATA(:, 'stimulus')), this_stimulus), 'n_frames'};
                        n_frames = n_frames(1);
                        fname = make_variable_name(stimuli_this_dataset{istim});
                        % Use the stimulus profile of temp_48 for all temperature stimuli
                        if startsWith(this_stimulus, 'temp_'), this_stimulus = 'temp_48'; end
                        STIMULI.(fname) = Metadata_workflow.load_stimuli(this_stimulus, 'frame_rate',frame_rate, 'n_frames',n_frames);
                    end
                    
                    % Initialize action to true
                    do_analyze = true;
                    msg = '';
                    
                    % Check whether this animal has ever been analyzed
                    if isfield(STATS_response_decoding, animal_ID) && any(strcmp(analysis_type_str, {'all', 'all_categories'})) %|| strcmp(analysis_type_str, 'all_categories')
                        stimuli_analyzed = fieldnames(STATS_response_decoding.(animal_ID));
                        % If all stimuli have been analyzed, ...
                        if all(ismember(fieldnames(STIMULI), stimuli_analyzed))
                            if response_decoding_missing_only  % No need to continue
                                do_analyze = false;
                                msg = ['All done for ', animal_ID, '. Skipped'];
                            else
                                msg = ['Will overwrite previous analysis of ', animal_ID];
                            end
                            
                        else
                            msg = ['Not all stimuli have been analyzed, so it will overwrite previous incomplete analysis of ', animal_ID];
                        end
                    else
                        msg = ['Processing ', animal_ID];
                    end
                    LOGGER.info(msg)
                    if ~do_analyze, continue, end
                    
                    % Get filename of inputs
                    inferred_spikes_filename = get_filename_of('spikes', animal_ID);
                    Ca_events_filename = get_filename_of('Ca_events', animal_ID);
                    if ~exist(inferred_spikes_filename, 'file') || ~exist(Ca_events_filename, 'file')
                        LOGGER.info('Dataset was not fully processed. Skipping')
                        continue
                    end
                    keep_cells = load_variable(Ca_events_filename, 'keep_cell');
                    keep_cells = logical(keep_cells);
                    data_is_loaded = false;
                    % Get number of cells in this FOV
                    n_cells = sum(keep_cells);
                    % Initialize cell arrays to contain set of cells to analyze
                    cell_subset = {};
                    cell_subset_info = {};
                    
                    
                    
                 case 'epi'
                     % Get filename of inputs
                     inferred_spikes_filename = get_filename_of('spikes', animal_ID);
                     Ca_events_filename = get_filename_of('Ca_events', animal_ID);
                     if ~exist(inferred_spikes_filename, 'file')
                         LOGGER.info('Dataset was not fully processed. Skipping')
                         continue
                     end
                     n_cells  = SQL_database.read_table_where ('experiments', 'n_ROIs', animal_ID, 'animal_ID', 'return_as_table', false);
                     keep_cells = true(n_cells, 1);
                     data_is_loaded = false;
                     METADATA = SQL_database.read_epi_trials(animal_ID);
%                      METADATA = SQL_database.read_table_where('sessions', {'date', 'stimulus', 'experimental_condition'}, {animal_ID, experiment_type}, {'animal_ID', 'experiment_name'});
                     columns_to_read = GC.analysis.(experiment_type).group_trials_by_session;%GC.session_name_prefix.(GC.experiment_name);
                     stimulus_column = find(ismember(columns_to_read, 'stimulus'));
                     compound = GC.analysis.ACC_SNI_anxiety.analysis_allowed_compounds;
                     compound = repmat((compound), height(METADATA),1);
                     METADATA.compound = compound;
                     METADATA = METADATA(:, [{'compound', 'day_from_surgery', 'experiment', 'type'}, columns_to_read]);
                     % Remove sessions where a compound was administered
                     METADATA(:, 'keep') = {true};
                     METADATA(~ismember(METADATA.compound, allowed_compounds), 'keep') = {false};
                     METADATA(~METADATA.keep, :) = [];
                     % Remove SP and SP_2
                     METADATA_epi = SQL_database.read_epi_trials(animal_ID);
                     analyzed_sessions = unique(METADATA_epi.date(ismemberCellRows(METADATA_epi.type, {'evoked'})), 'stable') ;
                     % Remove other stimuli that cannot be analyzed by this algorithm
                     session_date = METADATA.date{1};
                     %                     sessions_to_analyze = [sessions_to_analyze; [repmat(session_date, length(stimuli_to_analyze), 1), stimuli_to_analyze(:), repmat(sessions(i_session), length(stimuli_to_analyze), 1)]];
                     %                     sessions_to_analyze = sessions_to_analyze(ismember(sessions_to_analyze,analyzed_sessions),:);
                     %                     sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), {'SP', 'SP_2', 'SPnl', 'SP_2nl', 'SPn', 'SP_2n'}), :);
                     %                     sessions_to_analyze = sessions_to_analyze(~ismember(sessions_to_analyze(:, stimulus_column), GC.detection_no_timestamp_stimuli), :);
                     stimuli_this_dataset = unique(METADATA.stimulus);
                     % Remove SPs
                     stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, {'SP', 'SP_2', 'SPn', 'SP_2n', 'SPnl', 'SP_2nl'}), :);
                     % Remove other stimuli that cannot be analyzed by this algorithm
                     stimuli_this_dataset = stimuli_this_dataset(~ismember(stimuli_this_dataset, GC.detection_no_timestamp_stimuli));
                     % Add column to mark unwanted trials
                     METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
                     % Remove trials of non-allowed stimuli
                     METADATA.keep(~ismember(METADATA{:, 'stimulus'}, stimuli_this_dataset) | ~ismember(METADATA{:, 'type'}, 'evoked')) = false;
                     if sum(METADATA.keep) < 1, continue, end
                     % Add frame_rate, n_frames
                     frame_rate = GC.epifluorescence_downsample_to_frame_rate;
                     frame_rate = repmat(frame_rate, height(METADATA),1);
                     METADATA.frame_rate = frame_rate;
                     n_frames = GC.analysis.(experiment_type).n_frames;
                     n_frames = repmat(n_frames, height(METADATA), 1);
                     METADATA.n_frames = n_frames;
                     % Mark sessions where a compound was administered
                     METADATA = METADATA(:, [{'compound', 'frame_rate', 'n_frames', 'keep','day_from_surgery', 'experiment' }, columns_to_read]);
                     cell_subset = {};
                     cell_subset_info = {};
                     
                     % Load stimuli
                    STIMULI = struct();
                    for istim = 1:length(stimuli_this_dataset)
                        this_stimulus = stimuli_this_dataset{istim};
                        n_frames = METADATA{ismember(table2cell(METADATA(:, 'stimulus')), this_stimulus), 'n_frames'};
                        n_frames = n_frames(1);
               
                        fname = make_variable_name(stimuli_this_dataset{istim});
                        STIMULI.(fname) = Metadata_workflow.load_stimuli(this_stimulus, 'frame_rate',GC.epifluorescence_downsample_to_frame_rate, 'n_frames',n_frames);
                        if sum(strcmp (fname, {'touch', 'pinprick'})) == 1
                            STIMULI.(fname).timestamps = 5;
                            STIMULI.(fname).duration = 0;
                        end
                    end
                    
                    % Initialize action to true
                    do_analyze = true;
                    msg = '';
                    
                    % Check whether this animal has ever been analyzed
                    if isfield(STATS_response_decoding, animal_ID) && any(strcmp(analysis_type_str, {'all', 'all_categories'}))
                        stimuli_analyzed = fieldnames(STATS_response_decoding.(animal_ID));
                        % If all stimuli have been analyzed, ...
                        if all(ismember(fieldnames(STIMULI), stimuli_analyzed))
                            if response_decoding_missing_only  % No need to continue
                                do_analyze = false;
                                msg = ['All done for ', animal_ID, '. Skipped'];
                            else
                                msg = ['Will overwrite previous analysis of ', animal_ID];
                            end
                            
                        else
                            msg = ['Not all stimuli have been analyzed, so it will overwrite previous incomplete analysis of ', animal_ID];
                        end
                    else
                        msg = ['Processing ', animal_ID];
                    end
                    LOGGER.info(msg)
                    if ~do_analyze, continue, end
             end

            % Select subset of cells to analyze. Each subset is enclosed in
            % a cell of the array `cell_subset`.
            switch analysis_type_str
                case {'all' , 'all_categories'}
                    cell_subset = {1:n_cells};

                    if  endsWith(analysis_type_str, 'categories') % Set variables to pass to decoding script
                        categories = {'noci', 'ave', 'inn'};
                        noci =  {'pinprick','temp_48', 'HPS', 'FPS'};
                        ave  =  {'puff','sound'};
                        inn =  {'touch','temp_43', 'temp_42'};
                        CATEGORIES = cell2table ({noci, ave, inn}, 'VariableNames', categories);
                        cell_subset_info = {'all_categories'};
                    else
                        cell_subset_info = {'all'};
                    end
                    
                case 'selective'
                    %                     keyboard
                    % Loop through each session and select only selective
                    % cells; loop through set of stimuli to which cells
                    % should be selective to
                    for i_cond = 1:n_sessions_to_analyze
                        %                         GC.response_decoding_selective_to         = {
                        %                                                  {'touch'}, ...
                        %                                                   ...
                        %                                                  {'touch'}, ...
                        %                                                  };
                        
                        for i_stimset = 1:length(GC.response_decoding_selective_to)
                            stimulus_set = GC.response_decoding_selective_to{i_stimset};
                            if ~ismember(fieldnames(SELECTIVITY_stats.(animal_ID)), stimulus_set)
                                continue
                            end
                            selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                            for i_stim = 1:length(stimulus_set)
                                selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                            end
                            selective_cells = selective_cells(keep_cells, :);
                            n_selective = sum(selective_cells);
                            % Keep only cells that are selective to all stimuli
                            selective_cells = find(all(selective_cells, 2));
                            if isempty(selective_cells), continue, end
                            cell_subset{end + 1, 1} = selective_cells;
                            cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                        end
                    end
                    
                case 'stable'
                    %                     keyboard
                    % Get pairs of consecutive sessions
                    condition_pairs = sessions_to_analyze([1:n_sessions_to_analyze-1; 2:n_sessions_to_analyze]');
                    n_pairs = size(condition_pairs, 1);
                    for i_pair = 1:n_pairs
                        session_1 = condition_pairs(i_pair, 1);
                        try
                            session_2 = condition_pairs(i_pair, 2);
                        catch
                            continue
                        end
                        session_columns = ismember(sessions, [session_1, session_2]);
                        stim_to_decode = cellfun(@cell2mat, GC.response_decoding_selective_to, 'UniformOutput', false);
                        for i_stimset = 1:length(GC.response_decoding_selective_to)
                            stimulus_set = GC.response_decoding_selective_to{i_stimset};
                            if ~ismember(stimulus_set,stimuli_this_dataset)
                                continue
                            end
                            selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                            for i_stim = 1:length(stimulus_set)
                                % Get only cells that are selective in both sessions
                                selective_cells(:, i_stim) = all(SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, session_columns) > .5, 2);
                            end
                            selective_cells = selective_cells(keep_cells, :);
                            % Keep only cells that are selective to all stimuli
                            selective_cells = find(all(selective_cells, 2));
                            if isempty(selective_cells), continue, end
                            cell_subset{end + 1, 1} = selective_cells;
                            cell_subset_info{end + 1, 1} = {[session_1, session_2], stimulus_set};
                        end
                    end
                    
                case 'ensembles'
                    
                    % Load info on assemblies
                    %                     assemblies_filename = get_filename_of('assemblies', animal_ID);
                    assemblies_filename = ['V:\Ca_imaging_pain\5_deconvolved_traces\Assemblies_SP_HPS_temp_48_touch\',animal_ID, '_assemblies.mat'];
                    
                    if ~exist(assemblies_filename, 'file')
                        LOGGER.info([animal_ID ' had no ensembles.'])
                        continue
                    end
                    ASSEMBLIES = load_variable(assemblies_filename, 'ASSEMBLIES');
                    stimuli_names = GC.assemblies_stimuli_to_combine{1};
                    % Loop through conditions and get info
                    conditions = fieldnames(ASSEMBLIES); conditions(ismember(conditions, 'info')) = [];
                    n_conditions = length(conditions);
                    % Set number of iterations for bootstrap
                    n_permutations = 1000;
                    days_from_surgery_this_animal = unique(ASSEMBLIES.info.sessions.day_from_surgery, 'stable');
                    for i_cond = 1:n_conditions
                        fn = conditions{i_cond};
                        stimuli_sets = fieldnames(ASSEMBLIES.(fn));
                        day_from_surgery =  days_from_surgery_this_animal{i_cond};
                        % select only the SP
                        %                         keyboard % (try to select the ensembles per stim that significantly are active during the evoked period)
                        %                         ensemble_to_take = 'SP';
                        pref_ens = cell(length(stimuli_sets), 1);
                        for i_stim = 1:length(stimuli_sets)
                            % Get session when this assembly was detected
                            %                             idx = find(ismember(ASSEMBLIES.info.sessions.session_fieldname, fn) & ismember(ASSEMBLIES.info.sessions.stimuli_fieldname, stimuli_sets{i_stim}), 1);
                            %                             day_from_surgery = ASSEMBLIES.info.sessions.day_from_surgery{idx};
                            % Get stimuli
                            %                             stimulus_set = ASSEMBLIES.info.sessions.stimuli(idx, :);
                            this_stim = stimuli_names{i_stim};
                            if isempty(ASSEMBLIES.(fn).(stimuli_sets{i_stim})), continue, end
                            % Get the ensembles
                            core_svd = ASSEMBLIES.(fn).(stimuli_sets{i_stim}).core_svd;
                            state_pks = ASSEMBLIES.(fn).(stimuli_sets{i_stim}).state_pks;
                            param =  ASSEMBLIES.(fn).(stimuli_sets{i_stim}).param;
                            info_ens =   ASSEMBLIES.(fn).(stimuli_sets{i_stim}).info;
                            n_ens = length(core_svd);
                            % Set parameters
                            a = NaN(1,n_ens);
                            b = a;
                            P = a;
                            if ~strcmp(this_stim, 'SP')
                                if ~startsWith(this_stim, 'temp_')
                                    time_stamp = info_ens.timestamps_bin;
                                    time_stamp_bs =  time_stamp;
                                    
                                else
                                    trials_idx = ismember(METADATA.stimulus, this_stim);
                                    frame_rate = unique(METADATA{trials_idx, 'frame_rate'});
                                    if length(frame_rate) > 1, error('frameRate:MoreThanOneFrameRate', '!'), end
                                    n_frames = unique(METADATA{trials_idx, 'n_frames'});
                                    STIMULUS_properties = Metadata_workflow.load_stimuli(this_stim, 'frame_rate',frame_rate, 'n_frames',n_frames);
                                    %                     load_variable(['/Volumes/GroupNevian4/Ca_imaging_pain/0_stimuli/', this_stim{1},'.mat'], 'STIMULUS');
                                    ts = round(STIMULUS_properties.timestamps * frame_rate);
                                    time_stamp = ts(2);
                                    time_stamp_bs= ts(1);
                                end
                                % evoked window
                                n_frames_trial = info_ens.n_bins_per_trial;
                                evoked_frames = [(time_stamp +1) : (time_stamp + time_stamp_bs )];
                                baseline_frames = [1:time_stamp_bs];
                                %                 this_stim_bsl_idx = ismember(param.all_activations.stimulus, this_stim) & ismember(param.all_activations.bin,baseline_frames );
                                length_stim = n_frames_trial;
                                
                                %                 total_length_stim = sum(ismember(param.all_activations.stimulus, this_stim));
                                if ~startsWith(this_stim, 'temp_')
                                    sweeps_stim = height(info_ens.bin_edges)/length_stim;
                                else
                                    sweeps_stim = height(info_ens.bin_edges)/n_frames;
                                end
                                % Loop through ensembles
                                for i_ens = 1: n_ens
                                    activation_this_ens = state_pks == i_ens;
                                    % Run permutation test - bootstrap
                                    activity_this_stim = activation_this_ens;
                                    act_reshaped= reshape(activity_this_stim', [length_stim,sweeps_stim]);
                                    act_reshaped = act_reshaped';
                                    count_activations_evoked_stim = sum(sum(act_reshaped(:,evoked_frames)));
                                    count_activations_bsl_stim    = sum(sum(act_reshaped(:,baseline_frames)));
                                    b(1,i_ens) =  count_activations_evoked_stim / (count_activations_evoked_stim +count_activations_bsl_stim );
                                    
                                    m = length(activity_this_stim); n = 1; d = count_activations_evoked_stim + count_activations_bsl_stim;
                                    
                                    count_n_per = NaN(n_permutations,1);
                                    for i_per = 1 :n_permutations
                                        activity_shifted = zeros(m, n);
                                        for column = 1:n
                                            activity_shifted(:, column) = randperm(m) <= d;
                                        end
                                        activity_shifted_reshaped_ctrl = reshape(activity_shifted', [length_stim,sweeps_stim]);
                                        activity_shifted_reshaped = activity_shifted_reshaped_ctrl';
                                        count_bsl_BS = sum(sum(activity_shifted_reshaped(:,baseline_frames)));
                                        count_evk_BS = sum(sum(activity_shifted_reshaped(:,evoked_frames)));
                                        count_n_per(i_per, 1) = count_evk_BS / (count_evk_BS + count_bsl_BS);
                                        
                                    end
                                    null_distribution = (count_n_per);
                                    P(1,i_ens) =  (sum(null_distribution(:) >= b(1,i_ens)) +1) / (n_permutations+1);
                                end
                                pref_ens{i_stim} = find(P(1,:)<0.05);
                                
                            else
                                % Count activations in SP and take the one with the highest
                                act_counts = NaN(n_ens,1);
                                for i_ens = 1: n_ens
                                    %                     this_ensemble = ['ensemble_', num2str(i_ens)];
                                    act_counts(i_ens,1) = sum(state_pks == i_ens);
                                end
                                [~,pref_ens_idx] = max(act_counts);
                                %                                 a(1,:) =  NaN;
                                %                                 b(1,:) =  NaN;
                                
                                %                                 P(1,i_ens) =  0.05;
                                pref_ens{i_stim} = pref_ens_idx;
                            end
                            %                             if ~strcmp(cell2mat(stimulus_set), ensemble_to_take)
                            %                                 continue
                            %                             else
                            %                                 % Get indices of cells
                            %                                 try
                            %                                     ensembles = ASSEMBLIES.(fn).(stimuli_sets{i_stim_set}).core_svd;
                            %                                 catch
                            %                                     continue
                            %                                 end
                            %                                 n_ensembles = length(ensembles);
                            %                                 for i_ens = 1:n_ensembles
                            %                                     if length(ensembles{i_ens}) < 2, continue, end
                            %                                     cell_subset{end + 1, 1} = ensembles{i_ens};
                            %                                     cell_subset_info{end + 1, 1} = {day_from_surgery, stimulus_set};
                            %                                 end
                            %                             end
                            if isempty(pref_ens{i_stim})
                                cell_subset{end + 1, 1} = [];
                            else
                                cell_subset{end + 1, 1} = core_svd{pref_ens{i_stim}};
                            end
                            cell_subset_info{end + 1, 1} = {day_from_surgery, this_stim};
                        end
                    end
                    
                case 'specific'
                    keyboard
                    % Loop through each session and select only selective
                    % cells; loop through set of stimuli to which cells
                    % should be selective to
                    for i_cond = 1:n_sessions_to_analyze
                        for i_stimset = 1:length(GC.response_decoding_selective_to)
                            stimulus_set = GC.response_decoding_selective_to{i_stimset};
                            selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                            for i_stim = 1:length(stimulus_set)
                                selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                            end
                            selective_cells = selective_cells(keep_cells, :);
                            % Keep only cells that are selective to all stimuli
                            selective_cells = find(all(selective_cells, 2));
                            
                            % Discard cells that were selective to other stimuli
                            stimuli_this_condition = unique(METADATA.stimulus(ismember(str2double(METADATA.day_from_surgery), sessions_to_analyze(i_cond))));
                            stimuli_this_condition(ismember(stimuli_this_condition, stimulus_set)) = [];
                            stimuli_this_condition(ismember(stimuli_this_condition, {'SP', 'SP_2', 'HPS_short', 'HPS_50%', 'odor-'})) = [];
                            stimuli_this_condition(ismember(stimuli_this_condition, GC.detection_no_timestamp_stimuli)) = [];
                            cells_to_discard = false(height(SELECTIVITY_stats.(animal_ID)), 1);
                            for i_stim = 1:length(stimuli_this_condition)
                                this_stimulus = make_variable_name(stimuli_this_condition{i_stim});
                                if ismember(this_stimulus, SELECTIVITY_stats.(animal_ID).Properties.VariableNames)
                                    cells_to_discard = cells_to_discard | SELECTIVITY_stats.(animal_ID).(this_stimulus)(:, sessions_to_analyze_idx(i_cond)) > .5;
                                end
                            end
                            cells_to_discard = find(cells_to_discard);
                            selective_cells = setdiff(selective_cells, cells_to_discard);
                            if isempty(selective_cells), continue, end
                            cell_subset{end + 1, 1} = selective_cells;
                            cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                        end
                    end
                case {'H_noci', 'H_noci_random'}
                    %                     keyboard
                    % Loop through each session and select only selective
                    % cells; loop through set of stimuli to which cells
                    % should be selective to
                    noci_stim ={'HPS', 'FPS', 'temp_48', 'pinprick'};
                    H_class = 3; % At least 3 stim present;
                    % If lower than 3 noci stim are present, continue
                    if sum(ismember(noci_stim, fieldnames(SELECTIVITY_stats.(animal_ID))))< H_class
                        LOGGER.info('No subset of H_noci stimuli is available. Skipped')
                        continue
                    end
                    % set the selective cells for at least 3 stim
                    stimulus_set = noci_stim(ismember(noci_stim, fieldnames(SELECTIVITY_stats.(animal_ID))));
                    for i_cond = 1:n_sessions_to_analyze
                        day_from_surgery = sessions(i_cond);
                        if day_from_surgery>0
                            day_from_surgery = ['+', num2str(day_from_surgery)];
                        else
                            day_from_surgery = ['-', num2str(abs(day_from_surgery))];
                        end
                        % Check if at least H_class nr stim are present in this session, otherwise, skip
                        stim_this_cond = unique(METADATA.stimulus(ismember(METADATA.day_from_surgery, day_from_surgery)));
                        if sum(ismember(noci_stim, stim_this_cond))< H_class
                            LOGGER.info('Skiping session, not enough noci stimuli')
                            continue
                        end
                        selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                        for i_stim = 1: length(stimulus_set)
                            selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                        end
                        selective_cells = selective_cells(keep_cells, :);
                        H_cells = find(sum(selective_cells,2)>=H_class);
                        if isempty(H_cells)
                            LOGGER.info(['No H_cells present in session ',day_from_surgery ,  '. Skipped']),
                            continue
                        end
                        
                        if strcmp(analysis_type_str, 'H_noci_random')
                            n_H_cells = length(H_cells);
                            init_cell_vector = 1:sum(keep_cells);
                            init_cell_vector= init_cell_vector';
                            cell_vector = init_cell_vector;
                            % Delete H_cells from vector
                            cell_vector(cell_vector(H_cells,1)) = [];
                            % Pick a random subset from the cell_vector that are not H_cells
                            random_H_cells = randi(numel(cell_vector), n_H_cells, 1);
                            % Select the pseudo H_cells
%                             cell_subset{end + 1, 1} = [1:n_cells];(Try to do decoding with ALL neurosn
                            cell_subset{end + 1, 1} = init_cell_vector(~ismember(init_cell_vector, random_H_cells)); % Original
%                             cell_subset{end + 1, 1} = random_H_cells;% (Try  to do decoding ONLY with random neurons
                            cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                        else
                            cell_vector = 1:sum(keep_cells);
                            cell_vector = cell_vector';
                            cell_subset{end + 1, 1} = cell_vector(~ismember(cell_vector, H_cells)); % Original
%                             cell_subset{end + 1, 1} = H_cells; % (Try  to do decoding ONLY with H_noci - 17.12.20)
                            cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                        end
                    end
                    %                     for i_cond = 1:n_sessions_to_analyze
                    %                         for i_stimset = 1:length(GC.response_decoding_selective_to)
                    %                             stimulus_set = GC.response_decoding_selective_to{i_stimset};
                    %                             if ~ismember(fieldnames(SELECTIVITY_stats.(animal_ID)), stimulus_set)
                    %                                 continue
                    %                             end
                    %                             selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                    %                             for i_stim = 1:length(stimulus_set)
                    %                                 selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                    %                             end
                    %                             selective_cells = selective_cells(keep_cells, :);
                    %                             n_selective = sum(selective_cells);
                    %                             % Keep only cells that are selective to all stimuli
                    %                             % Remove only H_noci cells where n_remove = n_selective
                    %                             rand_cells = randperm(sum(keep_cells), n_cells - n_selective);
                    %                             selective_cells = rand_cells';
                    % %                             selective_cells = find(all(selective_cells, 2));
                    %                             if isempty(selective_cells), continue, end
                    %                             cell_subset{end + 1, 1} = selective_cells;
                    %                             cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                    %                         end
                    %                     end
                case 'H_sal'
                    %                     keyboard
                    % Loop through each session and select only selective
                    % cells; loop through set of stimuli to which cells
                    % should be selective to
                    % This takes H-R saliency neurons
                    sal_stim ={'HPS', 'FPS', 'temp_48', 'pinprick', 'puff', 'sound', 'temp_43', 'temp_42', 'touch'};
                    H_class = 3; % At least 3 stim present;
                    % If lower than 3 noci stim are present, continue
                    if sum(ismember(sal_stim, fieldnames(SELECTIVITY_stats.(animal_ID))))< H_class
                        LOGGER.info('No subset of H_sal stimuli is available. Skipped')
                        continue
                    end
                    % set the selective cells for at least 3 stim
                    stimulus_set = sal_stim(ismember(sal_stim, fieldnames(SELECTIVITY_stats.(animal_ID))));
                    for i_cond = 1:n_sessions_to_analyze
                        day_from_surgery = sessions(i_cond);
                        if day_from_surgery>0
                            day_from_surgery = ['+', num2str(day_from_surgery)];
                        else
                            day_from_surgery = ['-', num2str(abs(day_from_surgery))];
                        end
                        % Check if at least H_class nr stim are present in this session, otherwise, skip
                        stim_this_cond = unique(METADATA.stimulus(ismember(METADATA.day_from_surgery, day_from_surgery)));
                        if sum(ismember(sal_stim, stim_this_cond))< H_class
                            LOGGER.info(['Skiping session ', day_from_surgery, ', not all stim were present'])
                            continue
                        end
                        selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                        for i_stim = 1: length(stimulus_set)
                            selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                        end
                        selective_cells = selective_cells(keep_cells, :);
                        H_cells = find(sum(selective_cells,2)>=H_class);
                        if isempty(H_cells)
                            LOGGER.info(['No H_cells present in session ',day_from_surgery ,  '. Skipped']),
                            continue
                        end

                        cell_vector = 1:sum(keep_cells);
                        cell_vector = cell_vector';
                        cell_subset{end + 1, 1} = cell_vector(~ismember(cell_vector, H_cells));
                        cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                    end
                    
                case 'H_sal_random'
                    %                     keyboard
                    % Loop through each session and select only selective
                    % cells; loop through set of stimuli to which cells
                    % should be selective to
                    % This takes H-R saliency neurons
                    sal_stim ={'HPS', 'FPS', 'temp_48', 'pinprick', 'puff', 'sound', 'temp_43', 'temp_42', 'touch'};
                    H_class = 3; % At least 3 stim present;
                    % If lower than 3 sal stim are present, continue
                    if sum(ismember(sal_stim, fieldnames(SELECTIVITY_stats.(animal_ID))))< H_class
                        LOGGER.info('No subset of H_sal stimuli is available. Skipped')
                        continue
                    end
                    % set the selective cells for at least 3 stim
                    stimulus_set = sal_stim(ismember(sal_stim, fieldnames(SELECTIVITY_stats.(animal_ID))));
                    msg = '';
                    % check if already analyzed
                    if isfield(STATS_response_decoding, animal_ID)
                        if response_decoding_missing_only  % No need to continue
                            msg = ['All done for ', animal_ID, '. Skipped'];
                            do_analyze = false;
                        else
                            msg = ['Will overwrite previous analysis of ', animal_ID];
                            do_analyze = true;
                        end
                        LOGGER.info(msg)
                        if ~do_analyze, continue, end
                    end
                    % Loop through conditions
                    for i_cond = 1:n_sessions_to_analyze
                        day_from_surgery = sessions(i_cond);
                        if day_from_surgery>0
                            day_from_surgery = ['+', num2str(day_from_surgery)];
                        else
                            day_from_surgery = ['-', num2str(abs(day_from_surgery))];
                        end
                        % Check if at least H_class nr stim are present in this session, otherwise, skip
                        stim_this_cond = unique(METADATA.stimulus(ismember(METADATA.day_from_surgery, day_from_surgery)));
                        if sum(ismember(sal_stim, stim_this_cond))< H_class
                            LOGGER.info(['Skiping session ', day_from_surgery, ', not all stim were present'])
                            continue
                        end
                        selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                        for i_stim = 1: length(stimulus_set)
                            selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                        end
                        selective_cells = selective_cells(keep_cells, :);
                        H_cells = find(sum(selective_cells,2)>=H_class);
                        if isempty(H_cells)
                            LOGGER.info(['No H_cells present in session ',day_from_surgery ,  '. Session skipped']),
                            continue
                        end
                        % Select random cells that are not part of the H-cells
                        n_H_cells = length(H_cells);
                        init_cell_vector = 1:sum(keep_cells);
                        init_cell_vector= init_cell_vector';
                        cell_vector = init_cell_vector;
                        % Delete H_cells from vector
                        cell_vector(cell_vector(H_cells,1)) = [];
                        % Pick a random subset from the cell_vector that are not H_cells
                        random_H_cells = randi(numel(cell_vector), n_H_cells, 1);
                        % Select the pseudo H_cells
                        cell_subset{end + 1, 1} = init_cell_vector(~ismember(init_cell_vector, random_H_cells));
                        cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                    end
            case {'H_sal_categories' , 'H_sal_random_categories'}
                % So far, we are taking cells that are not H-R for saliency
%                     keyboard
                    categories  = {'noci', 'ave', 'inn'};
                    sal_stim ={'HPS', 'FPS', 'temp_48', 'pinprick', 'puff', 'sound', 'temp_43', 'temp_42', 'touch'};
                    H_class = 3;
                    stimulus_set = sal_stim(ismember(sal_stim, fieldnames(SELECTIVITY_stats.(animal_ID))));
                    % Check the categories for which each stimulus belongs to
                    cat_this_animal = cell(0);
                    for i_cat = 1 : length(categories)
                        this_cat = categories{i_cat};
                        switch this_cat
                            case 'noci'
                                noci =  {'pinprick','temp_48', 'HPS', 'FPS'};
                                does_belong = ismember(noci, stimulus_set);
                                if any (does_belong)
                                    cat_this_animal(end+1) = {this_cat};
                                end
                            case 'ave'
                                ave = {'puff', 'sound'};
                                does_belong = ismember(ave, stimulus_set);
                                if any (does_belong)
                                    cat_this_animal(end+1) = {this_cat};
                                end
                            case 'inn'
                                inn = {'touch','temp_43', 'temp_42'};
                                does_belong = ismember(inn, stimulus_set);
                                if any (does_belong)
                                    cat_this_animal(end+1) = {this_cat};
                                end
                        end                    
                    end
                    CATEGORIES = cell2table ({noci, ave, inn}, 'VariableNames', categories);
                    msg = '';
                    % check if already analyzed
                    if isfield(STATS_response_decoding, animal_ID)
                        if response_decoding_missing_only  % No need to continue
                            msg = ['All done for ', animal_ID, '. Skipped'];
                            do_analyze = false;
                        else
                            msg = ['Will overwrite previous analysis of ', animal_ID];
                            do_analyze = true;
                        end
                        LOGGER.info(msg)
                        if ~do_analyze, continue, end
                    end
                    % Loop through conditions
                    for i_cond = 1:n_sessions_to_analyze
                        day_from_surgery = sessions(i_cond);
                        if day_from_surgery>0
                            day_from_surgery = ['+', num2str(day_from_surgery)];
                        else
                            day_from_surgery = ['-', num2str(abs(day_from_surgery))];
                        end
                        % Check if at least H_class nr stim are present in this session, otherwise, skip
                        stim_this_cond = unique(METADATA.stimulus(ismember(METADATA.day_from_surgery, day_from_surgery)));
                        if sum(ismember(sal_stim, stim_this_cond))< H_class
                            LOGGER.info(['Skiping session ', day_from_surgery, ', not all stim were present'])
                            continue
                        end
                        selective_cells = false(height(SELECTIVITY_stats.(animal_ID)), length(stimulus_set));
                        for i_stim = 1: length(stimulus_set)
                            selective_cells(:, i_stim) = SELECTIVITY_stats.(animal_ID).(stimulus_set{i_stim})(:, sessions_to_analyze_idx(i_cond)) > .5;
                        end
                        selective_cells = selective_cells(keep_cells, :);
                        H_cells = find(sum(selective_cells,2)>=H_class);
                        if isempty(H_cells)
                            LOGGER.info(['No H_cells present in session ',day_from_surgery ,  '. Session skipped']),
                            continue
                        end
                        n_H_cells = length(H_cells);
                        
                        if strcmp (analysis_type_str, 'H_sal_random_categories')
                            % Select random cells that are not part of the H-cells
                            init_cell_vector = 1:sum(keep_cells);
                            init_cell_vector= init_cell_vector';
                            cell_vector = init_cell_vector;
                            % Delete H_cells from vector
                            cell_vector(cell_vector(H_cells,1)) = [];
                            % Pick a random subset from the cell_vector that are not H_cells
                            random_H_cells = randi(numel(cell_vector), n_H_cells, 1);
                            % Select the pseudo H_cells
                            cell_subset{end + 1, 1} = init_cell_vector(~ismember(init_cell_vector, random_H_cells));
                            cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), cat_this_animal};
                        else
                            cell_vector = 1:sum(keep_cells);
                            cell_vector = cell_vector';
                            cell_subset{end + 1, 1} = cell_vector(~ismember(cell_vector, H_cells));
                            cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), cat_this_animal};
                        end
                    end
            end

            % Run analysis on each cell subset
            n_subsets = length(cell_subset);
            if n_subsets == 0
                LOGGER.info('No subset of cells is available. Skipped')
                continue
            end
            
            for i_cellset = 1:n_subsets
                % check wether cell subset is empty
                if isempty(cell_subset{i_cellset}), continue, end
                % Check whether this dataset needs to be analyzed
                if response_decoding_missing_only && ~any(strcmp(analysis_type_str, {'all','all_categories'}))
                    try
                        switch analysis_type_str
                            case 'stable'
                                same_session = ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}, 'rows');
                                
                            case 'ensembles'
                                same_session = ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
                                    cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding.(animal_ID).info.cell_subset);
                                
                            otherwise
                                same_session = ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1});
                        end
                        row = find(same_session & ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));
                        if ~isempty(row)
                            if strcmp(analysis_type_str, 'stable')
                                msg = ['sessions ', strjoin(value2str(cell_subset_info{i_cellset}{1}), ' and ')];
                            else
                                msg = ['session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                            end
                            if ismember(analysis_type_str, {'stable', 'selective', 'specific'})
                                msg2 = [' cells selective to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            elseif ismember(analysis_type_str,'H_noci')
                                msg2 = [' cells H_noci taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            elseif ismember(analysis_type_str,'H_noci')
                                msg2 = [' cells H_noci random taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            elseif ismember(analysis_type_str,'H_sal')
                                msg2 = [' cells H_saliency taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            elseif ismember(analysis_type_str,'H_sal_random')
                                msg2 = [' random cells H_saliency taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            elseif ismember(analysis_type_str,'H_sal_categories')
                                msg2 = [' cells H_saliency in categories taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            elseif ismember(analysis_type_str,'H_sal_random_categories')
                                msg2 = [' random cells H_saliency in categories taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and ')];
                            else
                                msg2 = ' cells in an ensemble present';
                            end
                            LOGGER.info(['Already analyzed the ', num2str(length(cell_subset{i_cellset})), msg2, ' in ', msg])
                            continue
                        end
                    end
                end

                % Make folder, if it does not exist
                data_dir = os.path.join(GC.data_root_path, GC.aggregated_data_path, experiment_type, filesep,'response_decoding', filesep, animal_ID, filesep);
                if ~strcmp(analysis_type_str,'all')
                    if strcmp(analysis_type_str, 'all_categories')
                        data_dir = os.path.join(data_dir, analysis_type_str);
                    else
                        data_dir = os.path.join(data_dir, analysis_type_str, ['session', num2str(cell_subset_info{i_cellset}{1}, '%+i'), '__', cell2mat(cell_subset_info{i_cellset}{2})]);
%                     data_dir = os.path.join(data_dir, analysis_type_str, ['session', num2str(cell_subset_info{i_cellset}{1}, '%+i'), '__', cell_subset_info{i_cellset}{2}]);
                    end
                end
                if ~exist(data_dir, 'dir'), mkdir(data_dir), end

                % Log action
                switch analysis_type_str
                    case 'all'
                        msg = 'Analyzing all cells';

                    case 'selective'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells selective to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];

                    case 'stable'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells stably selective to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in sessions ', strjoin(value2str(cell_subset_info{i_cellset}{1}), ' and ')];

                    case 'ensembles'
%                         msg = ['Analyzing an ensemble of ', num2str(length(cell_subset{i_cellset})), ' cells present in session ', cell_subset_info{i_cellset}{1}];
                        msg = ['Analyzing an ensemble of ', num2str(length(cell_subset{i_cellset})), ' cells present in session ', cell_subset_info{i_cellset}{1}, ' for stimulus ', cell_subset_info{i_cellset}{2}];

                    
                    case 'specific'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells specific to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                    case 'H_noci'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells not H_noci taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
 
                    case 'H_noci_random'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells random for not H_noci taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                        
                    case 'H_sal'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells not-H_sal taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                        
                    case 'H_sal_random'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells random not-H_sal taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                        
                    case 'H_sal_categories'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells not-H_sal in categories taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                        
                    case 'H_sal_random_categories'
                        msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells random not-H_sal in categories taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
                         
                    case 'all_categories'
                        msg = 'Analyzing all cells in categories';        
                end
                LOGGER.info(msg)

                % Load data
                if ~data_is_loaded
                    LOGGER.trace('Loading spikes')
                    spikes = load_variable(inferred_spikes_filename, 'spikes');
                    spikes = spikes(keep_cells, :);
                    data_is_loaded = true;
                end

                % When analyzing ensembles, only consider the session of interest
                switch analysis_type_str
                    case 'stable'
                        trials_idx = ismember(str2double(METADATA.day_from_surgery), cell_subset_info{i_cellset}{1});
                        [result, cm] = response_decoding(animal_ID, spikes(cell_subset{i_cellset}, trials_idx), METADATA(trials_idx, :), STIMULI, data_dir, experiment_type);
                        
                    case 'ensembles'
                        [result, cm] = response_decoding(animal_ID, spikes(cell_subset{i_cellset}, :), METADATA, STIMULI, data_dir, experiment_type);
                        
                    case {'H_sal_categories', 'H_sal_random_categories', 'all_categories'}
                        [result, cm] = response_decoding_categories(animal_ID, spikes(cell_subset{i_cellset}, :), METADATA, STIMULI, CATEGORIES, data_dir, experiment_type, cell_subset_info{i_cellset}); 
                        
                    otherwise
                        [result, cm] = response_decoding(animal_ID, spikes(cell_subset{i_cellset}, :), METADATA, STIMULI, data_dir, experiment_type);
                end
                
                % Reload file if there are more datasets to analyze
                if exist(stats_filename, 'file')
                    STATS_response_decoding = load_variable(stats_filename, 'STATS_response_decoding');
                    STATS_response_decoding_cm = load_variable(stats_filename, 'STATS_response_decoding_cm');
                end

                % Create variables if they don't exist
                if ~isfield(STATS_response_decoding, animal_ID)
                    STATS_response_decoding.(animal_ID) = struct();
                    STATS_response_decoding_cm.(animal_ID) = struct();
                end
                if ~isfield(STATS_response_decoding.(animal_ID), 'info')
                    STATS_response_decoding.(animal_ID).info = cell2table(cell(0, 3), 'VariableNames',{'stimulus', 'session', 'cell_subset'});
                    STATS_response_decoding_cm.(animal_ID).info = STATS_response_decoding.(animal_ID).info;
                end

                % Store results
                switch analysis_type_str
                    case {'all', 'all_categories'}
                        STATS_response_decoding.(animal_ID) = result;
                        STATS_response_decoding_cm.(animal_ID) = cm;

                    otherwise
                        % Find row in info table
                        switch analysis_type_str
                            case 'stable'
                                try
                                    row = find(ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}, 'rows') & ...
                                         ismember(STATS_response_decoding.(animal_ID).info.stimulus, cell_subset_info{i_cellset}{2}));
%                                      ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));

                                catch
                                    row = 1;
                                end
                                
                            case 'ensembles'
                                try
%                                     row = find(ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                         ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')) & ...
%                                         cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding.(animal_ID).info.cell_subset));
                                    row = find(ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
                                        ismember(STATS_response_decoding.(animal_ID).info.stimulus, cell_subset_info{i_cellset}{2}) & ...
                                        cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding.(animal_ID).info.cell_subset));
%                                     ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')) & ...
                                catch
                                    row = 1;
                                end
                            case {'selective', 'H_noci', 'H_noci_random', 'H_sal','H_sal_random', 'H_sal_categories', 'H_sal_random_categories'}
                                try
                                    row = find(ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
                                          ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));
                                
%                                     row = find(ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                         ismember(STATS_response_decoding.(animal_ID).info.stimulus, cell_subset_info{i_cellset}{2}) & ...
%                                         cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding.(animal_ID).info.cell_subset));
%                                     ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')) & ...
                                catch
%                                     keyboard
                                    row = 1;
                                end
                                
                            otherwise
                                row = find(ismember(STATS_response_decoding.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
                                    ismember(STATS_response_decoding.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));
                                
                        end
                        if isempty(row)
                            row = height(STATS_response_decoding.(animal_ID).info) + 1;
                        end
                        if row ~= 1
                            STATS_response_decoding.(animal_ID).info.stimulus{row} = strjoin(cell_subset_info{i_cellset}{2}, '_');
%                             STATS_response_decoding.(animal_ID).info.stimulus{row}
%                             = cell_subset_info{i_cellset}{2}; % Use this
%                             for ensembles (fix code)

                            if strcmp(analysis_type_str, 'stable')
                                STATS_response_decoding.(animal_ID).info.session(row, :) = cell_subset_info{i_cellset}{1};
                            else
                                try
                                    STATS_response_decoding.(animal_ID).info.session{row} = cell_subset_info{i_cellset}{1};
                                catch ME
                                    if strcmp(ME.identifier, 'MATLAB:cellAssToNonCell')
                                        STATS_response_decoding.(animal_ID).info.session(row) = cell_subset_info{i_cellset}{1};
                                    else
                                        rethrow(ME)
                                    end
                                end
                            end
                            STATS_response_decoding.(animal_ID).info.cell_subset{row} = cell_subset{i_cellset};
                        else
                            STATS_response_decoding.(animal_ID).info = cell2table({cell_subset_info{i_cellset}{2}, cell_subset_info{i_cellset}{1}, cell_subset(i_cellset)}, 'VariableNames',STATS_response_decoding.(animal_ID).info.Properties.VariableNames);
%                             STATS_response_decoding.(animal_ID).info = cell2table({strjoin(cell_subset_info{i_cellset}{2}, '_'), cell_subset_info{i_cellset}{1}, cell_subset(i_cellset)}, 'VariableNames',STATS_response_decoding.(animal_ID).info.Properties.VariableNames);
                        end
                        % Copy info table
                        STATS_response_decoding_cm.(animal_ID).info = STATS_response_decoding.(animal_ID).info;
                        % Update output structure
                        fn = ['cell_set_', num2str(i_cellset)];
                        STATS_response_decoding.(animal_ID).(fn) = result;
                        STATS_response_decoding_cm.(animal_ID).(fn) = cm;
                end
                % Overwrite file on disk
                save(stats_filename, 'STATS_response_decoding','STATS_response_decoding_cm', '-v7.3')
            end
        end
    end

    % Increase counter
    current_action = current_action + 1;
end

%% EPM
if analyze_EPM_track
   
    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Computing ', char(916), 'EPM_tracking'])
    for iid = 1:length(all_animal_names)
        % Get codename of this animal
        animal_ID = all_animal_names{iid};
        LOGGER.info(['Processing EPM ', animal_ID])
        % Extract calcium traces
        EPM_track(animal_ID, EPM_missing_only)
    end
    current_action = current_action + 1;
    disp('done with EPM analysis')
end




%% MLint exceptions
%#ok<*AGROW,*NASGU,*CTCH,*EFIND>
