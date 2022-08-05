%% PREAMBLE
% This function launches the Metadata workflow, which updates the metadata for
% each trial of two-photon calcium imaging data.

%% MAIN
function done = main(INFO, action)   
done = true;  % This is the output of the workflow, which is read by the launcher GUI and used to continue to the next workflow

% Get logger and general_configs
global LOGGER GC
% Log beginning of workflow
LOGGER.info('Metadata workflow', 'decorate',true)

% Check whether user wants to add or remove data from database
switch action
    case 'add'
        % 1. Write experiments
        % Log action
        LOGGER.info('Updating metadata')
        toggle_toolbox('ScanImageTiffReader', 'on')

        % Get animal ID
        animal_ID = INFO.selected_animal_ID;
        
        % Read from database once more to get the experiment id that was automatically assigned.
        METADATA_experiments = SQL_database.read_table_where('experiments', [], animal_ID, 'animal_ID');

        % 2. Write sessions and trials
        % Get sessions
        sessions = INFO.experiments;

        METADATA_sessions = [];
        METADATA_trials = [];
        % Loop through animals
        for isess = 1:size(sessions,1)
            % Get date
            experiment_date = sessions.date{isess};

            % Fill in metadata for this session
            metadata_session = struct();
            exp_row = ismember(table2cell(METADATA_experiments(:, 'animal_ID')), animal_ID);
            experiment_id                  = METADATA_experiments{exp_row, 'experiment_id'};
            metadata_session.experiment_id = experiment_id;
            metadata_session.date          = experiment_date;
            metadata_session.stimulus      = sessions.stimulus{isess};
            metadata_session.experimental_condition = sessions.experimental_condition{isess};

            % Read metadata from csv file
            switch INFO.selected_data_type
                case '2p'
                    metadata_file = os.path.join(GC.data_root_path, GC.tiff_metadata_path, [animal_ID, '_', experiment_date, '.csv']);
                    session_metadata = Metadata_workflow.read_preprocessing_metadata(metadata_file);
                    % Keep only metadata for the current session
                    session_rows = true(height(session_metadata), 1);
                    for icol = 1:length(GC.filename_match_metadata)
                        session_rows = session_rows & ismember(session_metadata.(GC.filename_match_metadata{icol}), sessions.(GC.filename_match_metadata{icol}){isess});
                    end
                    session_metadata = session_metadata(session_rows, :);
                    % Get number of trials and store it
                    n_trials = height(session_metadata);
                    metadata_session.n_trials = n_trials;

                    % Loop through trials
                    for itrial = 1:n_trials
                        % Log action
                        LOGGER.trace([animal_ID, ' ', experiment_date, ' ', sessions.experimental_condition{isess}, ' ', sessions.stimulus{isess}, ': trial ', num2str(itrial), '/', num2str(n_trials)])

                        % Get the name of the tif file
                        tiff_name_to_store = [animal_ID, filesep(), experiment_date, filesep(), session_metadata.filename{itrial}, '.tif'];
                        tiff_name = os.path.join(GC.data_root_path, GC.tiff_raw_path, tiff_name_to_store);

                        % Read tif metadata
                        [frame_rate, frame_width, frame_height, n_frames, bidirectional_scanning, ~, tiff_date, n_channels] = read_tiff_info(tiff_name, false);

                        % Fill in metadata for trials
                        metadata_trials = struct();
                        metadata_trials.experiment_id  = experiment_id;
                        metadata_trials.session_id     = isess;
                        metadata_trials.tiff_path      = tiff_name_to_store;
                        metadata_trials.date           = experiment_date;
                        metadata_trials.recording_time = datestr(tiff_date, 'HHMMSS');
                        metadata_trials.frame_width    = frame_width;
                        metadata_trials.frame_height   = frame_height;
                        metadata_trials.frame_rate     = frame_rate;
                        metadata_trials.n_frames       = n_frames;
                        metadata_trials.needs_fix_for_bidirectional_scanning = bidirectional_scanning;
                        metadata_trials.has_structural_channel = n_channels > 1;

                        % Concatenate metadata of trials
                        METADATA_trials = [METADATA_trials; metadata_trials];
                    end
                    
                case 'epi'
                    metadata_session.n_trials = -1;
            end
            
            % Concatenate metadata of sessions
            METADATA_sessions = [METADATA_sessions; metadata_session];
        end

        % Convert to table
        if strcmp(INFO.selected_data_type, '2p')
            METADATA_trials = struct2table(METADATA_trials);
        else
            METADATA_trials = cell2table(cell(0, 1));
        end
        METADATA_sessions = struct2table(METADATA_sessions);

        % Find out whether these sessions exist already in the database
        in_DB = SQL_database.contains('sessions', METADATA_sessions);
        % Add experiments not in database
        if any(~in_DB)
            % Log progress
            LOGGER.trace(['Updating ', GC.database_table_prefix, '_sessions table'])
            % Update SQL
            sessions_not_in_DB = METADATA_sessions(~in_DB, :);
            SQL_database.update('sessions', [], sessions_not_in_DB)
        end

        if strcmp(INFO.selected_data_type, '2p')
            % Read metadata once more to get the correct session id
            METADATA_sessions_id = SQL_database.read_table_where('sessions', 'session_id', METADATA_sessions);
            METADATA_sessions_id = METADATA_sessions_id.session_id;
            % Fix session id in trials table
            METADATA_trials(:, 'session_id') = num2cell(METADATA_sessions_id(METADATA_trials{:, 'session_id'}));

            % Sort trials by date
            METADATA_trials = sortrows(METADATA_trials, {'date','recording_time'});
            METADATA_trials.date = [];  % This column was only used for sorting

            % Find out whether these trials exist already in the database
            in_DB = SQL_database.contains('trials', METADATA_trials(:, 'tiff_path'));
            % Add experiments not in database
            if any(~in_DB)
                % Log progress
                LOGGER.trace(['Updating ', GC.database_table_prefix, '_trials table'])
                % Update SQL
                trials_not_in_DB = METADATA_trials(~in_DB, :);
                SQL_database.update('trials', [], trials_not_in_DB)
            end
        end
        
        LOGGER.info('done')
        toggle_toolbox('ScanImageTiffReader', 'off')
        
    case 'remove'
        % Log action
        LOGGER.info('Deleting metadata')
        
        % Make list of experiments to delete
        sessions_to_delete = INFO.experiments;
        
        % Find all IDs to delete
        IDs_to_delete = SQL_database.read_table_where('trials', {'trial_id','session_id','experiment_id'}, sessions_to_delete, GC.experiments_columns);
        if ~isempty(IDs_to_delete)  % There are trials associated to these sessions
            LOGGER.trace(['Deleting ', num2str(height(IDs_to_delete)), ' trials'])
            SQL_database.delete_table_where('trials', num2cell(unique(IDs_to_delete.trial_id)), 'trial_id')
        else % If no trials for these sessions, look for sessions records    
            LOGGER.trace('No trials left for these experiments. Looking for sessions records')
            IDs_to_delete = SQL_database.read_table_where('sessions', {'session_id','experiment_id'}, sessions_to_delete, GC.experiments_columns);
        end
        if isempty(IDs_to_delete)  % There are no sessions associated to these experiments
            LOGGER.trace('No sessions left for these experiments')
        else
            LOGGER.trace(['Deleting sessions of ', num2str(height(IDs_to_delete)), ' trials'])
            SQL_database.delete_table_where('sessions', num2cell(unique(IDs_to_delete.session_id)), 'session_id')
        end
        
        % Find empty experiments records
        experiments_ID = SQL_database.read_table_where('experiments', 'experiment_id', unique(sessions_to_delete(:,2)), 'animal_ID');
        experiments_ID = num2cell(unique(experiments_ID.experiment_id));
        still_sessions = SQL_database.read_table_where('sessions', 'session_id', experiments_ID, 'experiment_id');
        still_trials   = SQL_database.read_table_where('trials',   'trial_id',   experiments_ID, 'experiment_id');
        % If there are no trials or sessions associated to this experiment, delete it
        if isempty(still_sessions) && isempty(still_trials)
            LOGGER.trace(['No sessions or trials left for these experiments. Deleting ', num2str(length(experiments_ID)), ' experiments'])
            SQL_database.delete_table_where('experiments', experiments_ID, 'experiment_id')
        end
        
        LOGGER.info('done')
end
        

%% MLint exceptions
%#ok<*AGROW>
