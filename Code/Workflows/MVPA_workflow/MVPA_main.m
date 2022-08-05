%% INITIAL CONDITIONS
% This code will gather data from 2p imaging and run the MVPA analysis,
% using RUN_Decoding_MVPA.m
% This script takes some parts form main.m file from the analysis workflow
% For any update, copy-paste from that script
% Version 1.0 - update : so far it can do all and all_categories
% 
clc
clear
close all
% add to path where the foler for the analysis is
% addpath('C:\Users\acuna\Documents') % Local folder. (To be changed in future releases)
CFG_MVPA = general_configs_MVPA;
INFO.selected_experiment = 'ACC_CCI_anesth';
global GC LOGGER CFG_MVPA
GC.experiment_name = 'ACC_CCI_anesth';
% addpath('/Users/marioacuna/ca_imaging_data_analysis/Figures_paper_FMP')
all_animal_names = animal_list();
do_random = 0; % random labels
classifier = 'multicalss_lda';

current_action = 1;
n_actions_to_perform = 1;
analysis_type = 1; % all cells
% analysis_type = 7; % excluding H-noci and random cells 
% Select the time of analysis to do on
compute_response_decoding_on_all_cells                         = true;
compute_response_decoding_on_selective_cells                   = false;
compute_response_decoding_on_stable_cells                      = false;
compute_response_decoding_on_ensemble_cells                    = false;
compute_response_decoding_on_specific_cells                    = false;
compute_response_decoding_on_H_noci_cells                      = false;
compute_response_decoding_on_H_noci_cells_random               = false;
compute_response_decoding_on_H_sal_cells                       = false;
compute_response_decoding_on_H_sal_cells_random                = false;
compute_response_decoding_on_H_sal_cells_categories            = false;
compute_response_decoding_on_H_sal_cells_random_categories     = false;
compute_response_decoding_on_all_cells_categories              = false;


response_decoding_missing_only = 0;
allowed_stimuli = {'HPS', 'FPS' ,'temp_48', 'pinprick', 'puff', 'sound', 'temp_43', 'touch'};
allowed_compounds = {''};
%% Start loading data
is_selectivity_stats_loaded = false;

for analysis_type = 1:12  % all, selective, stable, ensembles, specific
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
            if ~compute_response_decoding_on_H_noci_cells_random, continue, end
            analysis_type_str = 'H_noci_random';
            analysis_type_msg = 'random non-H_noci cells';
            
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
    if ismember(analysis_type, [2, 3, 5, 6, 7, 8, 9, 10, 11]) && ~is_selectivity_stats_loaded
        try
            SELECTIVITY_stats = load_variable(os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', 'SELECTIVITY_stats__p.mat'), 'SELECTIVITY_stats');
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
    %     stats_filename = get_filename_of(['response_decoding_stats_', analysis_type_str], INFO.selected_experiment);
    if ~do_random
        stats_filename = os.path.join(GC.data_root_path,  GC.aggregated_data_path, INFO.selected_experiment, ['response_decoding_stats_MVPA_',classifier,'_',analysis_type_str,'.mat']);
    else
        stats_filename = os.path.join(GC.data_root_path,  GC.aggregated_data_path, INFO.selected_experiment, ['response_decoding_stats_MVPA_',classifier,'_',analysis_type_str,'_RAN.mat']);
    end
    % Make folder if it doesn't exist
    [folder, ~] = fileparts(stats_filename);
    if ~exist(folder,'dir'), mkdir(folder), end
    % Load file
    if exist(stats_filename, 'file')
        STATS_response_decoding_MVPA_cm = load_variable(stats_filename, 'STATS_response_decoding_MVPA_cm');
        STATS_response_decoding_MVPA_days = load_variable(stats_filename, 'STATS_response_decoding_MVPA_days');
        STATS_response_decoding_MVPA_info =  load_variable(stats_filename, 'STATS_response_decoding_MVPA_info');
    else
        STATS_response_decoding_MVPA_cm = struct();
        STATS_response_decoding_MVPA_days = struct();
        STATS_response_decoding_MVPA_info = struct();
    end
    
    for iid = 16:length(all_animal_names)
        % Get codename of this animal
        animal_ID = all_animal_names{iid};
        if ismember(analysis_type, [2, 3, 5, 6,7, 8, 9, 10]) && ~is_selectivity_stats_loaded
            if ~ismember(animal_ID, fieldnames(SELECTIVITY_stats))
                LOGGER.warn([animal_ID, ' not fully analyzed'])
                continue
            end
        end
        experiment_type =  INFO.selected_experiment;
        data_type = SQL_database.read_table_where('experiments', 'data_type', animal_ID,'animal_ID', 'return_as_table',false);
        switch data_type % So far only works for 2p-data
            case '2p'
                % Read metadata from database
                METADATA = SQL_database.read_table_where('trials',{'+','stimulus','experimental_condition','date'}, {animal_ID, 'ACC_CCI_anesth'}, {'animal_ID', 'experiment_name'});
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
                if isfield(STATS_response_decoding_MVPA_cm, animal_ID) && any(strcmp(analysis_type_str, {'all', 'all_categories'}))
                    %                     stimuli_analyzed = fieldnames(STATS_response_decoding_MVPA.(animal_ID));
%                     stimuli_analyzed = STATS_response_decoding_MVPA_classes.(animal_ID);
                    
                    % If all stimuli have been analyzed, ...
                    
                    if response_decoding_missing_only  % No need to continue
                        do_analyze = false;
                        msg = ['All done for ', animal_ID, '. Skipped'];
                    else
                        msg = ['Will overwrite previous analysis of ', animal_ID];
                    end
                    
                    %                     if all(ismember(fieldnames(STIMULI), stimuli_analyzed))
                    %                         if response_decoding_missing_only  % No need to continue
                    %                             do_analyze = false;
                    %                             msg = ['All done for ', animal_ID, '. Skipped'];
                    %                         else
                    %                             msg = ['Will overwrite previous analysis of ', animal_ID];
                    %                         end
                    %
                    %                     else
                    %                         msg = ['Not all stimuli have been analyzed, so it will overwrite previous incomplete analysis of ', animal_ID];
                    %                     end
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
                    stimuli_analyzed = fieldnames(STATS_response_decoding_MVPA.(animal_ID));
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
        CATEGORIES = []; % Pass an empty variable in case if there's no categorization
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
                keyboard
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
                    
                    for i_stimset = 1:length(GC.response_decoding_selective_to)
                        stimulus_set = GC.response_decoding_selective_to{i_stimset};
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
            case 'random'
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
                        % Remove only random cells where n_remove = n_selective
                        rand_cells = randperm(sum(keep_cells), n_cells - n_selective);
                        selective_cells = rand_cells';
                        %                             selective_cells = find(all(selective_cells, 2));
                        if isempty(selective_cells), continue, end
                        cell_subset{end + 1, 1} = selective_cells;
                        cell_subset_info{end + 1, 1} = {sessions_to_analyze(i_cond), stimulus_set};
                    end
                end
        end
        
        % Run analysis on each cell subset
        n_subsets = length(cell_subset);
        if n_subsets == 0
            LOGGER.info('No subset of cells is available. Skipped')
            continue
        end
        
%         for i_cellset = 1:n_subsets
            % check wether cell subset is empty
%             if isempty(cell_subset{i_cellset}), continue, end
            % Check whether this dataset needs to be analyzed
           
            
            % Make folder, if it does not exist
            %                 data_dir = os.path.join(GC.data_root_path, GC.aggregated_data_path, experiment_type, filesep,'response_decoding', filesep, animal_ID, filesep);
%             data_dir = os.path.join(GC.data_root_path, GC.aggregated_data_path, experiment_type,'response_decoding_MVPA', animal_ID);
            
%             if ~strcmp(analysis_type_str, 'all')
%                 if strcmp(analysis_type_str,'all_categories')
%                      data_dir = os.path.join(data_dir, analysis_type_str);
%                 else
%                     data_dir = os.path.join(data_dir, analysis_type_str, ['session', num2str(cell_subset_info{i_cellset}{1}, '%+i'), '__', cell2mat(cell_subset_info{i_cellset}{2})]);
%                 %                     data_dir = os.path.join(data_dir, analysis_type_str, ['session', num2str(cell_subset_info{i_cellset}{1}, '%+i'), '__', cell_subset_info{i_cellset}{2}]);
%                 end
%             end
%             if ~exist(data_dir, 'dir'), mkdir(data_dir), end
            
            % Log action
%             switch analysis_type_str
%                 case 'all'
%                     msg = 'Analyzing all cells';
%                     
%                 case 'selective'
%                     msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells selective to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
%                     
%                 case 'stable'
%                     msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells stably selective to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in sessions ', strjoin(value2str(cell_subset_info{i_cellset}{1}), ' and ')];
%                     
%                 case 'ensembles'
%                     %                         msg = ['Analyzing an ensemble of ', num2str(length(cell_subset{i_cellset})), ' cells present in session ', cell_subset_info{i_cellset}{1}];
%                     msg = ['Analyzing an ensemble of ', num2str(length(cell_subset{i_cellset})), ' cells present in session ', cell_subset_info{i_cellset}{1}, ' for stimulus ', cell_subset_info{i_cellset}{2}];
%                     
%                     
%                 case 'specific'
%                     msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells specific to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
%                 case 'random'
%                     msg = ['Analyzing the ', num2str(length(cell_subset{i_cellset})), ' cells randomly taken to ', strjoin(cell_subset_info{i_cellset}{2}, ' and '), ' in session ', num2str(cell_subset_info{i_cellset}{1}, '%+i')];
%                     
%                  case 'all_categories'
%                         msg = 'Analyzing all cells in categories';   
%             end
%             LOGGER.info(msg)
            
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
                    [result, cm] = response_decoding(animal_ID, spikes(cell_subset{i_cellset}, trials_idx), METADATA(trials_idx, :), STIMULI, data_dir);
                    
                case 'ensembles'
                    [result, cm] = response_decoding(animal_ID, spikes(cell_subset{i_cellset}, :), METADATA, STIMULI, data_dir, experiment_type);
                    
                case {'H_sal_categories', 'H_sal_random_categories', 'all_categories'}
                    [result, cm, classes, info] = RUN_Decoding_MVPA(animal_ID, spikes(cell_subset{i_cellset}, :), METADATA, STIMULI, CATEGORIES, data_dir, experiment_type);
                    
                otherwise
                    %                     [result, cm] = response_decoding(animal_ID, spikes(cell_subset{i_cellset}, :), METADATA, STIMULI, data_dir, experiment_type);
                    [~, cm, days, info] = RUN_Decoding_MVPA(animal_ID, spikes, METADATA, STIMULI, allowed_stimuli, CATEGORIES, experiment_type, do_random, analysis_type_str ,cell_subset_info, cell_subset, classifier);
            end
            
            % Reload file if there are more datasets to analyze
            if exist(stats_filename, 'file')
%                 STATS_response_decoding_MVPA = load_variable(stats_filename, 'STATS_response_decoding_MVPA');
                STATS_response_decoding_MVPA_cm = load_variable(stats_filename, 'STATS_response_decoding_MVPA_cm');
                STATS_response_decoding_MVPA_info =load_variable(stats_filename, 'STATS_response_decoding_MVPA_info');
                STATS_response_decoding_MVPA_days =load_variable(stats_filename, 'STATS_response_decoding_MVPA_days');
            end
            
            % Create variables if they don't exist
            if ~isfield(STATS_response_decoding_MVPA_cm, animal_ID)
%                 STATS_response_decoding_MVPA.(animal_ID) = struct();
                STATS_response_decoding_MVPA_cm.(animal_ID) = struct();
                STATS_response_decoding_MVPA_info.(animal_ID) = struct();
                STATS_response_decoding_MVPA_days.(animal_ID) = struct();
            end
            %             if ~isfield(STATS_response_decoding_MVPA.(animal_ID), 'info') && ~isfield(STATS_response_decoding_MVPA.(animal_ID), 'classes')
            %                 STATS_response_decoding_MVPA.(animal_ID).info = cell2table(cell(0, 3), 'VariableNames',{'stimulus', 'session', 'cell_subset'});
            %                 STATS_response_decoding_MVPA_cm.(animal_ID).info = STATS_response_decoding_MVPA.(animal_ID).info;
            %             end
            
            % Store results
%             switch analysis_type_str
%                 case { 'all', 'all_categories'}
%                     STATS_response_decoding_MVPA.(animal_ID) = result;
                    STATS_response_decoding_MVPA_cm.(animal_ID) = cm;
                    STATS_response_decoding_MVPA_info.(animal_ID) = info;
                    STATS_response_decoding_MVPA_days.(animal_ID) = days;
                    
%                 otherwise
%                     keyboard
%                     % Find row in info table
%                     switch analysis_type_str
%                         case 'stable'
%                             try
%                                 row = find(ismember(STATS_response_decoding_MVPA.(animal_ID).info.session, cell_subset_info{i_cellset}{1}, 'rows') & ...
%                                     ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, cell_subset_info{i_cellset}{2}));
%                                 %                                      ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));
%                                 
%                             catch
%                                 row = 1;
%                             end
%                             
%                         case 'ensembles'
%                             try
%                                 %                                     row = find(ismember(STATS_response_decoding_MVPA.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                 %                                         ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')) & ...
%                                 %                                         cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding_MVPA.(animal_ID).info.cell_subset));
%                                 row = find(ismember(STATS_response_decoding_MVPA.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                     ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, cell_subset_info{i_cellset}{2}) & ...
%                                     cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding_MVPA.(animal_ID).info.cell_subset));
%                                 %                                     ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')) & ...
%                             catch
%                                 row = 1;
%                             end
%                         case {'selective', 'random'}
%                             try
%                                 row = find(ismember(STATS_response_decoding_MVPA.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                     ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));
%                                 
%                                 %                                     row = find(ismember(STATS_response_decoding_MVPA.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                 %                                         ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, cell_subset_info{i_cellset}{2}) & ...
%                                 %                                         cellfun(@(x) isequal(cell_subset{i_cellset}, x), STATS_response_decoding_MVPA.(animal_ID).info.cell_subset));
%                                 %                                     ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')) & ...
%                             catch
%                                 keyboard
%                                 row = 1;
%                             end
%                             
%                         otherwise
%                             row = find(ismember(STATS_response_decoding_MVPA.(animal_ID).info.session, cell_subset_info{i_cellset}{1}) & ...
%                                 ismember(STATS_response_decoding_MVPA.(animal_ID).info.stimulus, strjoin(cell_subset_info{i_cellset}{2}, '_')));
%                             
%                     end
% %                     if isempty(row)
% %                         row = height(STATS_response_decoding_MVPA.(animal_ID).info) + 1;
% %                     end
% %                     if row ~= 1
% %                         STATS_response_decoding_MVPA.(animal_ID).info.stimulus{row} = strjoin(cell_subset_info{i_cellset}{2}, '_');
% %                         %                             STATS_response_decoding_MVPA.(animal_ID).info.stimulus{row}
% %                         %                             = cell_subset_info{i_cellset}{2}; % Use this
% %                         %                             for ensembles (fix code)
% %                         
% %                         if strcmp(analysis_type_str, 'stable')
% %                             STATS_response_decoding_MVPA.(animal_ID).info.session(row, :) = cell_subset_info{i_cellset}{1};
% %                         else
% %                             try
% %                                 STATS_response_decoding_MVPA.(animal_ID).info.session{row} = cell_subset_info{i_cellset}{1};
% %                             catch ME
% %                                 if strcmp(ME.identifier, 'MATLAB:cellAssToNonCell')
% %                                     STATS_response_decoding_MVPA.(animal_ID).info.session(row) = cell_subset_info{i_cellset}{1};
% %                                 else
% %                                     rethrow(ME)
% %                                 end
% %                             end
% %                         end
% %                         STATS_response_decoding_MVPA.(animal_ID).info.cell_subset{row} = cell_subset{i_cellset};
% %                     else
% %                         STATS_response_decoding_MVPA.(animal_ID).info = cell2table({cell_subset_info{i_cellset}{2}, cell_subset_info{i_cellset}{1}, cell_subset(i_cellset)}, 'VariableNames',STATS_response_decoding_MVPA.(animal_ID).info.Properties.VariableNames);
% %                         %                             STATS_response_decoding_MVPA.(animal_ID).info = cell2table({strjoin(cell_subset_info{i_cellset}{2}, '_'), cell_subset_info{i_cellset}{1}, cell_subset(i_cellset)}, 'VariableNames',STATS_response_decoding_MVPA.(animal_ID).info.Properties.VariableNames);
% %                     end
%                     % Copy info table
% %                     STATS_response_decoding_MVPA_cm.(animal_ID).info = STATS_response_decoding_MVPA.(animal_ID).info;
% %                     % Update output structure
% %                     fn = ['cell_set_', num2str(i_cellset)];
% %                     STATS_response_decoding_MVPA.(animal_ID).(fn) = result;
% %                     STATS_response_decoding_MVPA_cm.(animal_ID).(fn) = cm;
%             end
            % Overwrite file on disk
            save(stats_filename, 'STATS_response_decoding_MVPA_cm','STATS_response_decoding_MVPA_info' ,'STATS_response_decoding_MVPA_days' , '-v7.3')
%         end
    end
end
current_action = current_action + 1;
disp('Finished')
