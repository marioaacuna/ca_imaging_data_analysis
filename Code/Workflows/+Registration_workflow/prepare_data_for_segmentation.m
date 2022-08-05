function prepare_data_for_segmentation(animal_ID, data_type, experiment_name)

% clc
% animal_ID = 'MA_28epi'
% data_type = 'epi'

if ~exist('data_type', 'var'), data_type = '2p'; end
% Get global variables
global GC LOGGER

filename_params = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_params.mat']);

switch data_type
    case '2p'
         % Read metadata from SQL database
        switch experiment_name
            case 'CLA_pain'
                METADATA = SQL_database.read_table_where('trials', {'+','animal_ID','date','stimulus','experimental_condition','has_structural_channel'}, animal_ID, 'animal_ID', 'return_all_sessions',true, 'split_experimental_condition', false);
            case {'ACC_CCI_anesth', 'CLA_ACC_dreadd' }
                METADATA = SQL_database.read_table_where('trials', {'+','animal_ID','date','stimulus','experimental_condition','has_structural_channel'}, animal_ID, 'animal_ID', 'return_all_sessions',true);
        end
        
        has_structural_channel = logical(unique(METADATA{:, 'has_structural_channel'}));
        % Get performed experiments
        experiments_info = uniqueCellRows(METADATA{:, {'date', 'experimental_condition'}}, 'rows');
        
    case 'epi'
        METADATA = SQL_database.read_table_where('sessions', {}, animal_ID, 'animal_ID', 'return_all_sessions',true);
        has_structural_channel = false;

        % Get performed experiments
        experiments_info = uniqueCellRows(METADATA{:, {'date', 'experiment'}}, 'rows');

        % Check whether to rerun conversion
        filename_local = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_joint.h5']);
        if ~exist(filename_local, 'file')
            % Check whether it was deleted by this same script after conversion
            if ~exist(filename_params, 'file')
                rerun_conversion = true;
            else
                PARAMETERS = load_variable(filename_params, 'PARAMETERS');
                condition_names = PARAMETERS.condition_names;
                condition_names = cellfun(@(x) regexp(x, ' - frames', 'split'), condition_names, 'UniformOutput',false);
                condition_names = cellfun(@(x) x{1}, condition_names, 'UniformOutput',false);
                condition_names = cellfun(@(x) regexp(x, ' ', 'split'), condition_names, 'UniformOutput',false);
                condition_dates = cellfun(@(x) x{1}, condition_names, 'UniformOutput',false);
                condition_names = cellfun(@(x) x{2}, condition_names, 'UniformOutput',false);
                condition_names = [condition_dates(:), condition_names(:)];
                conditions = uniqueCellRows(condition_names);
                full_experiments_info = uniqueCellRows(METADATA{:, {'date', 'experimental_condition'}}, 'rows');
                is_condition_present = ismemberCellRows(full_experiments_info, conditions);
                rerun_conversion = ~all(is_condition_present);
            end
            
            if rerun_conversion
                % Align cell maps
                runAlignment(animal_ID, experiments_info)
                % Concatenate and normalize movies
                extractSignalsJoint(animal_ID, experiments_info)
            end
        end
        GC.experiment_name =experiment_name;
end

% Convert data to be used with python GUI for ROI segmentation
local_directory = os.path.join(GC.temp_dir, animal_ID);
if ~exist(local_directory, 'dir'), mkdir(local_directory), end
filename_projections = os.path.join(local_directory, [animal_ID, '_projections.mat']);

if ~exist(filename_projections, 'file')
    switch data_type
        case '2p'
            LOGGER.info('Computing image stack projections to pass to ROI segmentation GUI')
            make_trial_projections(animal_ID, filename_projections,experiment_name,  'data_type',data_type, 'reload_data',false);
            if has_structural_channel
                make_trial_projections(animal_ID, filename_projections, 'data_type',data_type, 'reload_data',false, 'only_structural_channel',true);
            end
            
        case 'epi'
            LOGGER.info('Computing image stack projections to pass to ROI segmentation GUI')
            make_trial_projections(animal_ID, filename_projections,experiment_name, 'data_type',data_type, 'reload_data',false);
               
            % Add cell maps from preprocessing
            PROJECTIONS = load_variable(filename_projections, 'PROJECTIONS');
            METADATA = SQL_database.read_table_where('sessions', {'+','animal_ID'}, animal_ID, 'animal_ID', 'return_all_sessions',true);
            METADATA{:, 'has_structural_channel'} = false;
            group_sessions_by = {'date', 'experiment'};
            % Load preprocessing parameters
            p_filename = get_filename_of('miniscope_movie_parameters', animal_ID);
            miniscope_parameters = load_variable(p_filename, 'p');
            % Read sessions. In case there's no timestamps, take only
            % sessions that were preprocessed, otherwise it's done already
            % (in the read_epi_trials script)
            try
                sessions = uniqueCellRows(METADATA{logical(METADATA.done_preprocessing), group_sessions_by}, 'rows');
            catch
                sessions = uniqueCellRows(METADATA{:, group_sessions_by}, 'rows');
            end
            n_sessions = size(sessions, 1);

            for i_sess = 1:n_sessions
                session_fieldname = ['session_', strjoin(sessions(i_sess, :), '_')];
                filename_cellMap = get_filename_of('miniscope_cell_map', animal_ID, sessions{i_sess, 1}, sessions{i_sess, 2});
                cellMap = load_variable(filename_cellMap, 'cellMap');
                n_repeats = size(PROJECTIONS.(session_fieldname).mean, 3);
                PROJECTIONS.(session_fieldname).cell_map = repmat(double(cellMap), [1, 1, n_repeats]);
            end
            save(filename_projections, 'PROJECTIONS', '-v7.3')
    end
end

if ~exist(filename_params, 'file')
    % Get order of sessions
    PROJECTIONS = load_variable(filename_projections, 'PROJECTIONS');
    session_order = fieldnames(PROJECTIONS);
    
    % Set parameters for each data type
    switch data_type
        case '2p'
            variable_name = 'Y';
            filename_local = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '.mat']);
            filename_remote = get_filename_of('motion_corrected', animal_ID);
            condition_names = cellfun(@(x,y) [x,y], METADATA.experimental_condition, METADATA.stimulus, 'UniformOutput',false);
            condition_names = strrep(condition_names, ';', ' ');
            frames_idx = get_trial_frames(animal_ID, 'reset_index',false, 'return_all_sessions',true);
            frame_size = [max(METADATA.frame_height), max(METADATA.frame_width)];
            frame_size = max(frame_size);
            frame_size = [frame_size, frame_size];
            is_stimulus_evoked = true;
            [~, idx, ~] = unique(METADATA.date, 'last');
            sessions_last_frame = frames_idx(idx, 2);
            
        case 'epi'
            variable_name = '1';
            filename_local = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_joint.h5']);
            filename_remote = '';
            % Load miniscope experiment parameters
            miniscope_parameters = load_variable(get_filename_of('miniscope_movie_parameters', animal_ID), 'p');
            % Set names and indices of each trial
            preprocessed_frame_edges = miniscope_parameters.user.jointExtraction.preprocessed_frame_edges;
%             sessions = uniqueCellRows(METADATA{:, {'date', 'experimental_condition'}});
            try
                sessions = uniqueCellRows(METADATA{logical(METADATA.done_preprocessing), {'date', 'experiment'}});
            catch
                sessions = uniqueCellRows(METADATA{:, {'date', 'experiment'}});
            end
             n_sessions = size(sessions, 1);
            frames_idx = [];
            condition_names = {};
            
            for i_sess = 1:n_sessions%miniscope_parameters.nSessions
                n_frames = preprocessed_frame_edges(i_sess, 2) - preprocessed_frame_edges(i_sess, 1) + 1;
                n_intervals = size(PROJECTIONS.(session_order{i_sess}).mean, 3);
                % Set beginning and end of each interval
                interval_start = linspace(1, n_frames, n_intervals + 1);
                interval_start = round(interval_start(1:end-1));
                interval_end = [interval_start(2:end) - 1, n_frames];
                this_session_frames_idx = [interval_start(:), interval_end(:)];
                frames_idx = [frames_idx; this_session_frames_idx + preprocessed_frame_edges(i_sess, 1) - 1];
                % Set name of conditions as beginning and end of each interval in seconds
                this_frame_rate = 5;% miniscope_parameters.user.frameRate(i_sess);
                this_session_condition_names = strcat(strjoin(sessions(i_sess, :), ' '), strcat({' - '}, cellfun(@(x,y) [x, '-', y], value2str(this_session_frames_idx(:, 1) / this_frame_rate, '%.0f'), value2str(this_session_frames_idx(:, 2) / this_frame_rate, '%.0f'), 'UniformOutput',false)), 's');
                condition_names = [condition_names; this_session_condition_names];
            end
            
            info = h5info(filename_local);
            frame_size = info.Datasets.Dataspace.Size([1, 2]);
            is_stimulus_evoked = false;
            
            sessions_last_frame = preprocessed_frame_edges(:, 2);
    end

    % Make structure containing parameters to pass to GUI
    PARAMETERS = struct();
    
    PARAMETERS.variable_name = variable_name;
    % Files and folders
    PARAMETERS.temp_dir = GC.temp_dir;
    PARAMETERS.animal_ID = animal_ID;
    PARAMETERS.filename_local = filename_local;
    PARAMETERS.filename_remote = filename_remote;
    PARAMETERS.filename_projections = filename_projections;
    PARAMETERS.filename_data_frame = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_data_frame.npy']);
    % PARAMETERS.filename_data_time = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_data_time.npy']);
    PARAMETERS.filename_data_time = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_data_time.hdf5']);
    PARAMETERS.filename_output = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_ROI_info.mat']);
    % Conditions and frames corresponding to trials
    PARAMETERS.condition_names = condition_names;
    PARAMETERS.frames_idx = frames_idx;
    PARAMETERS.session_order = session_order;
    % Store info regarding the imaging stack
    switch data_type
        case '2p', PARAMETERS.frame_rate = unique(METADATA.frame_rate);
        case 'epi', PARAMETERS.frame_rate = GC.epifluorescence_downsample_to_frame_rate;
    end
    PARAMETERS.n_frames = max(frames_idx(:));
    PARAMETERS.frame_height = frame_size(1);
    PARAMETERS.frame_width = frame_size(2);
    PARAMETERS.n_pixels = PARAMETERS.frame_height * PARAMETERS.frame_width;        
    % Add info on stimuli
    PARAMETERS.stimulus_evoked = is_stimulus_evoked;
    if is_stimulus_evoked
        % Load stimulus profile of each stimulus
        stimuli = unique(METADATA.stimulus, 'stable');
        stimulus_profile = cell(length(stimuli), 2);
        condition_has_stimulus = zeros(size(PARAMETERS.condition_names,1), 1);

        for istim = 1:length(stimuli)
            stimulus_name = stimuli{istim};
            % Load file
            filename = [GC.data_root_path, GC.stimuli_path, stimulus_name, '.mat'];
            if ~exist(filename, 'file')
                n_frames = max(table2array(SQL_database.read_table_where('trials', 'n_frames', stimulus_name, 'stimulus')));
                this_profile = zeros(n_frames, 1);
                LOGGER.error(['Stimulus ''', stimulus_name, ''' has not been added to the database yet. A blank profile will be passed to the GUI.'])
            else
                STIMULUS = load_variable(filename, 'STIMULUS');
                frame_rate = unique(METADATA{ismember(table2cell(METADATA(:, 'stimulus')), stimulus_name), 'frame_rate'});
                n_frames = unique(METADATA{ismember(table2cell(METADATA(:, 'stimulus')), stimulus_name), 'n_frames'});
                row = ismember(STIMULUS.stimulus_profile{:, 'frame_rate'}, frame_rate) & ismember(STIMULUS.stimulus_profile{:, 'n_frames'}, n_frames);
                time_axis = STIMULUS.stimulus_profile{row, 'time_axis'};
                stimulus_profile = STIMULUS.stimulus_profile{row, 'stimulus_profile'};
                STIMULUS.time_axis = time_axis{1};
                STIMULUS.stimulus_profile = stimulus_profile{1};
                this_profile = STIMULUS.stimulus_profile(:);
            end
            % Store info on this stimulus
            stimulus_profile{istim, 1} = stimulus_name;
            stimulus_profile{istim, 2} = this_profile;
            % Get where this stimulus was shown
            condition_has_stimulus(ismember(METADATA.stimulus, stimulus_name)) = istim;
        end
        
    else
        stimulus_profile = [];
        condition_has_stimulus = [];
    end
    % Store info in output variable
    PARAMETERS.stimuli = stimulus_profile;
    PARAMETERS.condition_has_stimulus = condition_has_stimulus;
    % Get index of frame that separates sessions
    PARAMETERS.sessions_last_frame = sessions_last_frame;
    % Store structure to disk
    save(filename_params, 'PARAMETERS', '-v6')
    
    filename_ROI_info = get_filename_of('ROI_info', animal_ID);
    if strcmp(data_type, 'epi') && ~exist(filename_ROI_info, 'file')
        LOGGER.info('Storing alignment of cell maps into GUI data')
        
        % Use cell maps alignment as suggested starting point for translating the FOV in the GUI
        filename = os.path.join(GC.registered_images, animal_ID, [animal_ID, '_alignedCellMaps.mat']);
        tforms = load_variable(filename, 'tforms');
        
        % Make TRANSFORMATION_IDX to mark which TRANSFORMATION_MATRICES apply to
        % which frame
        n_conditions = size(PARAMETERS.condition_names, 1);
        TRANSFORMATION_IDX = cell(n_conditions, 3);
        TRANSFORMATION_IDX(:, 1) = num2cell(0:n_conditions - 1);  % 0-indexed for python
        TRANSFORMATION_IDX(:, 3) = PARAMETERS.condition_names;
        condition_idx = [];
        for i_sess = 1:miniscope_parameters.nSessions
            n_frames = preprocessed_frame_edges(i_sess, 2) - preprocessed_frame_edges(i_sess, 1) + 1;
            n_intervals = size(PROJECTIONS.(session_order{i_sess}).mean, 3);
            condition_idx = [condition_idx; ones(n_intervals, 1) * i_sess - 1];
        end
        TRANSFORMATION_IDX(:, 2) = num2cell(condition_idx);

        % https://uk.mathworks.com/help/images/matrix-representation-of-geometric-transformations.html
        TRANSFORMATION_MATRICES = cell(miniscope_parameters.nSessions, 3);
        for i_sess = 1:miniscope_parameters.nSessions
            this_tform = tforms{i_sess}.T;

            % Compute rotation angles (and account for rounding errors)
            q1 = acosd(this_tform(1, 1));
            q2 = asind(this_tform(1, 2));
            % rotation_matrix = int32(round(mean([q1, q2]), 1));
            rotation_matrix = int32(0);  % Set it to 0 until I find a way to import it correctly into openCV

            % Compute translation offsets
            offset = [round(this_tform(3, 2), 0), round(this_tform(3, 1), 0)];
            original_canvas_size = miniscope_parameters.user.position{i_sess}([4, 3]);
            odd_number_of_pixels = original_canvas_size / 2 > floor(original_canvas_size / 2);
            % Compute canvas size
            canvas_size = original_canvas_size;
            if rotation_matrix ~= 0, canvas_size = canvas_size .* sqrt(2); end
            if any(offset ~= 0), canvas_size = canvas_size + (max(abs(offset)) .* 2); end
            canvas_size = ceil(canvas_size);
            canvas_size = floor(canvas_size ./ 2) .* 2;
            if odd_number_of_pixels(1), canvas_size(1) = canvas_size(1) + 1; end
            if odd_number_of_pixels(2), canvas_size(2) = canvas_size(2) + 1; end
            canvas_size = canvas_size + 4;
            % Compute margins
            center = floor(canvas_size / 2);
            new_center = center - round(offset);
            half_height = floor(original_canvas_size(1) / 2);
            half_width = floor(original_canvas_size(2) / 2);
            top = new_center(1) - half_height;
            left = new_center(2) - half_width;
            translation_matrix = {int32([top; top + original_canvas_size(1); left; left + original_canvas_size(2)]); int32(canvas_size(:))};

            % Leave scaling unaltered
            scaling_matrix = int32([1; 1]);

            % Store data
            TRANSFORMATION_MATRICES{i_sess, 1} = scaling_matrix;
            TRANSFORMATION_MATRICES{i_sess, 2} = rotation_matrix;
            TRANSFORMATION_MATRICES{i_sess, 3} = translation_matrix;
        end
        % Test if Tr Matx when empty runs well on GUI
        TRANSFORMATION_MATRICES = [];
        
        % Create other variables expected by GUI
        ROIs = [];
        ROI_info = [];
        MAPPING_TABLE = [];
        modified_MAPPING = [];
        ROI_fluorescence = [];
        
        % Store variables in local file
        ROI_filename_local = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_ROI_info.mat']);
        save(ROI_filename_local, 'ROIs','ROI_info','TRANSFORMATION_MATRICES','TRANSFORMATION_IDX','MAPPING_TABLE','modified_MAPPING','ROI_fluorescence', '-v6')
    end
    
else
    PARAMETERS = load_variable(filename_params, 'PARAMETERS');
end

% Convert data to be used in GUI
if exist(PARAMETERS.filename_local, 'file')
    run_in_python('GUI_ROI_segmentation', filename_params, true)
end

% If files exist, delete local file
% if exist(PARAMETERS.filename_data_frame, 'file') && exist(PARAMETERS.filename_data_time, 'file') && exist(PARAMETERS.filename_local, 'file')
%     keyboard
% %     delete(PARAMETERS.filename_local) % this is going to be used for
% %     CNMF-e 
% end

% Check whether there are ROI info already in the local folder. Otherwise, copy
% the file from the server
filename_local = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_ROI_info.mat']);
if ~exist(filename_local, 'file')
    filename_remote = get_filename_of('ROI_info', animal_ID);
    if exist(filename_remote, 'file')
        LOGGER.info(['Copying previous ROI segmentation from ''', filename_remote, ''''], 'contains_path',true)
        copyfile(filename_remote, filename_local)
        LOGGER.info('done', 'append',true)
    end
end


%% MLint exceptions
%#ok<*STRNU,*AGROW,*NASGU>
