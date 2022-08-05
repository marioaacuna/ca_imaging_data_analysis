function varargout = GUI_mark_stimuli(varargin)
    gui_Singleton = 1;
    gui_State = struct('gui_Name',mfilename, 'gui_Singleton', gui_Singleton, 'gui_OpeningFcn',@GUI_mark_stimuli_OpeningFcn, 'gui_OutputFcn',@GUI_mark_stimuli_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
    if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
    if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else, gui_mainfcn(gui_State, varargin{:}); end
function varargout = GUI_mark_stimuli_OutputFcn(~, ~, handles), varargout{1} = handles.output;


%% INIT
function GUI_mark_stimuli_OpeningFcn(hObject, ~, handles, ~)
    global GC

    % Assign GUI object to handles
    handles.output = hObject;
    % Link callback function when window has to close
    handles.figure.CloseRequestFcn = @quit_GUI;
    % Make sure window is not docked
    handles.figure.WindowStyle = 'normal';

    % Set maximum number of cameras
    n_cameras = 2;
    
    % Store axes and other items
    cams_axes = cell(1, n_cameras);
    sliders = cell(1, n_cameras);
    video_times = cell(1, n_cameras);
    for i_cam = 1:n_cameras
        cams_axes{i_cam} = handles.(['cam', num2str(i_cam), '_video']);
        sliders{i_cam} = handles.(['cam', num2str(i_cam), '_slider']);
        video_times{i_cam} = handles.(['cam', num2str(i_cam), '_time']);
    end

    % Make axes equal
    for i_cam = 1:n_cameras
        set(handles.figure, 'CurrentAxes',handles.(['cam', num2str(i_cam), '_video']))
        axis equal
    end
    
    % Make INFO structure
    INFO = struct();
    INFO.animal_ID = '';
    INFO.session_date = '';
    INFO.experiment = '';
    INFO.max_n_cameras = n_cameras;
    INFO.n_cameras = NaN;
    INFO.cams = [];
    INFO.cams_axes = cams_axes;
    INFO.sliders = sliders;
    INFO.video_times = video_times;
    INFO.image_handles = [];
    INFO.cameras_offset = [0, 0];
    INFO.videos_in_sync = false;
    INFO.last_sp_behav_timestamp = [];
    INFO.is_data_loaded = false;
    INFO.is_database_updated = true;
    
    % Make sure that GUI is empty from previous runs
    delete(handles.menu_choose_animal.Children);
    % Clear axes and table from previous info
    for i_cam = 1:length(INFO.cams_axes)
        cla(INFO.cams_axes{i_cam})
    end
    handles.table_stimuli.Data = cell(0, 5);

    % Store names of buttons
    fn = fieldnames(handles);
    INFO.stimuli_button_names = fn(cellfun(@(x) startsWith(x, 'btn_stim_'), fn));
    % Delete callbacks assigned to each one, and make sure to assign the correct one
    for i_btn = 1:length(INFO.stimuli_button_names)
        this_btn = handles.(INFO.stimuli_button_names{i_btn});
        stimulus = strrep(INFO.stimuli_button_names{i_btn}, 'btn_stim_', '');
        this_btn.Callback = {@btn_stim_add_new_Callback, handles, stimulus};
    end
    
    % Populate menu
    METADATA = SQL_database.read_table_where('sessions', {'experiment_name', 'animal_ID', 'date', 'experimental_condition'}, 1,'has_behavior');
    % Convert table to nested structure
    METADATA = METADATA(:, {'experiment_name', 'animal_ID', 'date', 'experiment'});
    % Rename 'date' column
    METADATA.Properties.VariableNames{3} = 'session_date';
    INFO.metadata_columns = METADATA.Properties.VariableNames;
    setappdata(handles.figure, 'INFO',INFO);
    % Prepend a prefix to each value so that it can be used as field name in the nested structure
    for i_col = 1:length(INFO.metadata_columns)
        this_col = INFO.metadata_columns{i_col};
        METADATA(:, this_col) = strcat('f_', METADATA{:, this_col});
    end    
    nested_METADATA = table2nested_struct(METADATA);
    % Make nested menus
    handles.menus = create_nested_menu(struct(), nested_METADATA, 'menu_choose_animal_experiments', handles.menu_choose_animal, hObject, true);
    
    % Fill in options for table_stimuli
    column_idx = ismember(handles.table_stimuli.ColumnName, 'Event');
    handles.table_stimuli.ColumnFormat{column_idx} = {' '};
    column_idx = ismember(handles.table_stimuli.ColumnName, 'Affective_response');
    handles.table_stimuli.ColumnFormat{column_idx} = [{' '}, GC.freely_moving_affective_responses(:)'];
    
    % Assign callback to other graphical elements
    set(handles.cam1_slider, 'Callback',{@cam_slider_Callback, hObject, 1})
    set(handles.cam2_slider, 'Callback',{@cam_slider_Callback, hObject, 2})
    set(handles.btn_play_cam1, 'Callback',{@btn_play_cam_Callback, hObject, 1})
    set(handles.btn_play_cam2, 'Callback',{@btn_play_cam_Callback, hObject, 2})
    
    % Turn off the GUI apart from the menu to choose the animal ID
    setappdata(handles.figure, 'INFO',INFO);
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    handles.figure.CloseRequestFcn = @quit_GUI;
    toggle_menus(handles, 'on')

    % Update handles structure
    guidata(hObject, handles);

    log('Ready')
    
    
%% CALLBACKS
function menu_choose_animal_Callback(~, ~, data_info, hObject)
    global GC
    % Retrieve handles
    handles = guidata(hObject);
    % Toggle off all menus while data is loaded
    toggle_menus(handles, 'off')
    % Get current info
    INFO = getappdata(handles.figure, 'INFO');
    
    % If a sessions is already loaded, alert the user
    if INFO.is_data_loaded && ~INFO.is_database_updated
        % Ask user whether to continue without saving
        selection = MessageBox({'Do you want to load another session'; 'without saving the latest changes'; 'to the current one?'}, 'Load new session', 'Yes', 'No', 'No');
        switch selection
            case 'Yes'
                % Simply continue
            case 'No'
                toggle_menus(handles, 'on')
                return
        end
    end
    
    % Unpack 'data_info'
    for i_col = 1:length(INFO.metadata_columns)
        eval([INFO.metadata_columns{i_col}, ' = data_info{i_col};']);
    end
    log(['Checking ', animal_ID, ' - ', session_date, ' - ', experiment])

    % Check content of video folder
    video_folder = os.path.join(GC.temp_dir, animal_ID);
%     GC.behavior_raw_video_root_path = 'T:\Mario\Behavior\';
%     video_folder = os.path.join(GC.behavior_raw_video_root_path, 'ACC_pain', GC.behavior_video_subfolder, animal_ID, session_date, experiment);
    content = dir(video_folder);
    files = {content.name}';
    files = files(cellfun(@(x) endsWith(x, '.mp4'), files));
    files = files(cellfun(@(x) startsWith(x, ['20',session_date]), files));
    files = strcat(video_folder, filesep, files);
    n_cameras = length(files);
    if n_cameras < 1
        log('Did not find any file. Make sure that raw videos have been concatenated')
        handles.menu_choose_animal.Enable = 'on';
        toggle_GUI(handles, 'on')
        handles.figure.CloseRequestFcn = @quit_GUI;
        toggle_menus(handles, 'on')
        return
    end
    INFO.animal_ID = animal_ID;
    INFO.session_date = session_date;
    INFO.experiment = experiment;
    INFO.n_cameras = n_cameras;
    if INFO.n_cameras > 1, suffix = 's'; else, suffix = ''; end
    log(['Found ', num2str(n_cameras), ' file', suffix, '. Loading ...'])

    % Clear axes and table from previous info
    for i_cam = 1:length(INFO.cams_axes)
        cla(INFO.cams_axes{i_cam})
    end
    handles.table_stimuli.Data = cell(0, 5);
    
    % Open files for reading
    INFO.cams = cell(n_cameras, 1);
    for i_cam = 1:n_cameras
        INFO.cams{i_cam} = VideoReader(files{i_cam});
    end

    % Load first frame
    current_frame = cell(n_cameras, 1);
    for i_cam = 1:n_cameras
        current_frame{i_cam} = readFrame(INFO.cams{i_cam});
    end

    % Show first frame
    INFO.image_handles = cell(n_cameras, 1);
    for i_cam = 1:n_cameras
        INFO.image_handles{i_cam} = imshow(current_frame{i_cam}, 'Parent',INFO.cams_axes{i_cam});
    end

    % Get info on movies
    movie_duration = zeros(n_cameras, 1);
    movie_frame_rate = zeros(n_cameras, 1);
    for i_cam = 1:n_cameras
        movie_duration(i_cam) = INFO.cams{i_cam}.Duration;
        movie_frame_rate(i_cam) = INFO.cams{i_cam}.FrameRate;
    end
    
    % Set ranges for sliders
    for i_cam = 1:n_cameras
        INFO.sliders{i_cam}.Min = 0;
        INFO.sliders{i_cam}.Max = movie_duration(i_cam);
        INFO.sliders{i_cam}.Value = 0;
        % The 2 SliderStep values are:
        % 1) when pressing arrows
        % 2) when clicking on the slider
        INFO.sliders{i_cam}.SliderStep = [1 / movie_frame_rate(i_cam), 5] ./ movie_duration(i_cam);
    end
    
    % Disable controls for second camera, if it doesn't exist
    if n_cameras < 2
        set(handles.cam2_video, 'Visible','off')
        set(handles.cam2_slider, 'Visible','off')
        set(handles.btn_play_cam2, 'Visible','off')
        set(handles.cam2_title, 'Visible','off')
        set(handles.cam2_time, 'Visible','off')
    end

    % Find out whether this dataset has been already analyzed
    METADATA = SQL_database.read_table_where('stimulations', {'+', 'experimental_condition', 'stimulus'}, {experiment_name, animal_ID, session_date},{'experiment_name', 'animal_ID', 'date'});
    if ~isempty(METADATA)
        METADATA = METADATA(ismember(METADATA.experiment, experiment), :);
    end
    if ~isempty(METADATA)
        log('Loading results of previous analysis ...')
        % If this dataset has been already analyzed, get the events already detected
        % and assemble them in a table
        table_data = Metadata_workflow.flatten_stimulations_table(METADATA);

        % Make 'table_stimuli' in GUI
        handles.table_stimuli.Data = prepare_data_for_GUI_table(table_data);

        % Set callback for when cell selection changes in table
        set(handles.table_stimuli, 'CellEditCallback', {@table_stimuli_Callback, handles})
        resize_table_columns(handles)
    end
    
    % Fill in options for table_stimuli
    METADATA = SQL_database.read_table_where('sessions', {'experimental_condition', 'stimulus'}, {experiment_name, animal_ID, session_date},{'experiment_name', 'animal_ID', 'date'});
    METADATA = METADATA(ismember(METADATA.experiment, experiment), :);
    stimuli = natsort(METADATA.stimulus);
    column_idx = ismember(handles.table_stimuli.ColumnName, 'Event');
    handles.table_stimuli.ColumnFormat{column_idx} = stimuli(:)';

    % Hide buttons that will not be used here
    switch experiment
        case {'SP', 'SP_2'}
            state_SP = 'on';
            state_stim = 'off';
        case 'pain'    
            state_SP = 'off';
            state_stim = 'on';
    end
    handles.btn_stim_cold.Visible = state_stim;
    handles.btn_stim_heat.Visible = state_stim;
    handles.btn_stim_pinprick.Visible = state_stim;
    handles.btn_stim_touch.Visible = state_stim;
    handles.btn_stim_puff.Visible = state_stim;
    handles.btn_sp_behav_affective.Visible = state_SP;
    handles.text_info_affective_behaviors.Visible = state_SP;
    
    % Store INFO
    INFO.is_data_loaded = true;
    INFO.is_database_updated = true;
    setappdata(handles.figure, 'INFO',INFO);

    % Enable GUI for interaction
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;
    toggle_menus(handles, 'on')
    log('Ready')


function btn_save_Callback(~, ~, handles)
    % Disable GUI
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')    

    global GC
    INFO = getappdata(handles.figure, 'INFO');
    
    % Convert data to a table
    data = cell2table(handles.table_stimuli.Data, 'VariableNames', handles.table_stimuli.ColumnName);
    % Convert Time_start and Time_end to seconds
    try
        data.Time_start = cellfun(@(x) convert_time_from_str_to_num(x), data.Time_start);
        data.Time_end = cellfun(@(x) convert_time_from_str_to_num(x), data.Time_end);
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:cellfun:NotACell')
            rethrow(ME)
        end
    end
    
    % Get session id
    METADATA = SQL_database.read_table_where('sessions', {'experiment_id', 'session_id', 'experimental_condition', 'stimulus'}, {INFO.animal_ID, INFO.session_date}, {'animal_ID', 'date'});
    METADATA = METADATA(ismember(METADATA.experiment, INFO.experiment), :);
    
    % Fix values for non-evoked pain sessions
    if ~strcmp(INFO.experiment, 'pain')
        data(:, 'Event') = {'SP'};
        data(:, 'Withdraw') = {false};
    end
    
    % Loop through events
    stimuli_detected = unique(data.Event);
    if ~isempty(stimuli_detected)
        log('Storing data to database ...')
        
        for i_stim = 1:length(stimuli_detected)
            this_stimulus = stimuli_detected{i_stim};
            rows = find(ismember(data.Event, this_stimulus));

            % --- Timestamps ----------
            timestamps = value2str(data{rows, {'Time_start', 'Time_end'}}, '%.3f');
            % Replace 'NaNs' with ''
            timestamps(ismember(timestamps, 'NaN')) = {''};
            % Add double quotes
            timestamps = cellfun(@(x) ['"', x, '"'], timestamps, 'UniformOutput',false);
            % Convert to a long string
            timestamps_ts = cell(size(timestamps, 1), 1);
            for i_ts = 1:size(timestamps, 1)
                timestamps_ts{i_ts} = ['[', strjoin(timestamps(i_ts, :), ', '), ']'];
            end
            timestamps_ts = ['[', strjoin(timestamps_ts, ', '), ']'];

            % --- Withdrawal response ----------
            response = strjoin(value2str(data.Withdraw(rows), '%i'), '');

            % --- Affective response ----------
            all_affective_responses = [{''}, GC.freely_moving_affective_responses(:)'];
            affective_responses = data.Affective_response(rows);
            affective_responses(ismember(affective_responses, ' ')) = {''};
            affective_response = strjoin(value2str(cellfun(@(x) find(ismember(all_affective_responses, x)), affective_responses), '%i'), '');

            % --- Valid stimuli ----------
            valid = strjoin(value2str(true(length(rows), 1), '%i'), '');

            % --- Stimulus type ----------
            if ismember(this_stimulus, GC.freely_moving_stimuli)
                stimulus_type = 'evoked';
            else
                stimulus_type = 'spontaneous';
            end

            % Get stimulus id
            session_id = METADATA.session_id(ismember(METADATA.stimulus, this_stimulus));
            experiment_id = METADATA.experiment_id(ismember(METADATA.stimulus, this_stimulus));
            stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id', {experiment_id, session_id}, {'experiment_id', 'session_id'}, 'return_as_table',false);

            if isempty(stimulus_id)
                % Add experiment and session id, and get newly assigned stimulus_id
                SQL_database.update('stimulations', {'experiment_id', 'session_id'}, {experiment_id, session_id})
                stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id', {experiment_id, session_id}, {'experiment_id', 'session_id'}, 'return_as_table',false);
            end
            SQL_database.update('stimulations', {'timestamps', 'response', 'type', 'valid', 'affective'}, {timestamps_ts, response, stimulus_type, valid, affective_response}, stimulus_id,'stimulus_id')
        end
    
    else  % Pass empty values in case there are no data in the stimuli_table
        log('Deleting data from database ...')
        
        stimuli_in_database = unique(METADATA.stimulus);
        for i_stim = 1:length(stimuli_in_database)
            this_stimulus = stimuli_in_database{i_stim};
            % Get stimulus id
            session_id = METADATA.session_id(ismember(METADATA.stimulus, this_stimulus));
            experiment_id = METADATA.experiment_id(ismember(METADATA.stimulus, this_stimulus));
            stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id', {experiment_id, session_id}, {'experiment_id', 'session_id'}, 'return_as_table',false);
            if isempty(stimulus_id)
                continue
            end
            SQL_database.delete_table_where('stimulations', stimulus_id,'stimulus_id')
        end
    end
    
    % Log outcome
    t = now;
    d = datetime(t,'ConvertFrom','datenum');
    log(['done', ' - ', datestr(d)])
    % Restore GUI
    toggle_GUI(handles, 'on')
    

function cam_slider_Callback(~, ~, hObject, which_camera)
    handles = guidata(hObject);
    if which_camera == 1
        value = handles.cam1_slider.Value;
    else
        value = handles.cam2_slider.Value;
    end
    update_video_frame(handles.figure, value, which_camera)


function btn_play_cam_Callback(~, ~, hObject, which_camera)
    handles = guidata(hObject);
    % If paused already, don't do anything
    if which_camera == 1
        if handles.btn_play_cam1.Value == 0, return, end
    else
        if handles.btn_play_cam2.Value == 0, return, end
    end
    % Play movie
    play_movie(handles)

    
function btn_sync_cameras_Callback(~, ~, handles)
    % Get INFO
    INFO = getappdata(handles.figure, 'INFO');
    INFO.videos_in_sync = handles.btn_sync_cameras.Value == 1;

    if handles.btn_sync_cameras.Value == 1
        % Store offset of both cameras
        for from_camera = 1:INFO.n_cameras
            INFO.cameras_offset(from_camera) = INFO.cams{from_camera}.CurrentTime;
        end
        % Reset time in both cameras to 0, and change color
        for from_camera = 1:INFO.n_cameras
            update_time_string(0, 0, INFO.video_times{from_camera})
            INFO.video_times{from_camera}.ForegroundColor = [0, .8, 0];
        end
        handles.btn_sync_cameras.BackgroundColor = [0, .8, 0];
        handles.btn_sync_cameras.String = 'Videos synchronized';
        
    else
        % Reset offset of both cameras
        for from_camera = 1:INFO.n_cameras
            INFO.cameras_offset(from_camera) = 0;
        end
        % Reset time in both cameras, and change color
        for from_camera = 1:INFO.n_cameras
            update_time_string(INFO.cams{from_camera}.CurrentTime, 0, INFO.video_times{from_camera})
            INFO.video_times{from_camera}.ForegroundColor = [0, 0, 0];
        end
        handles.btn_sync_cameras.BackgroundColor = [0.6510, 0.6510, 0.6510];
        handles.btn_sync_cameras.String = 'Videos not synchronized';        
    end
    
    setappdata(handles.figure, 'INFO', INFO);
  
    
function btn_stim_add_new_Callback(~, ~, handles, stim_type)
    % Get INFO
    INFO = getappdata(handles.figure, 'INFO');

    % Get current time from camera 1 + offset
    timestamp = INFO.cams{1}.CurrentTime - INFO.cameras_offset(1);
    % Convert to string in mm:ss.ms
    time_minutes = floor(timestamp / 60);
    time_seconds = floor(timestamp - time_minutes * 60);
    time_milliseconds = round((timestamp - time_minutes * 60 - time_seconds) * 1000);
    timestamp_str = [num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];
    % Update table in GUI
    handles.table_stimuli.Data(end+1, :) = {timestamp_str, '', stim_type, false, ''};
    % Sort timestamps
    handles.table_stimuli.Data = sortrows(handles.table_stimuli.Data, 1);
    resize_table_columns(handles)
    
    % Mark that the new state is not updated
    INFO.is_database_updated = false;
    setappdata(handles.figure, 'INFO',INFO)
    

function table_stimuli_Callback(~, event, handles)
    % Delete rows with empty timestamp
    rows_to_delete = ismember(handles.table_stimuli.Data(:, 1), '');
    if any(rows_to_delete)
        handles.table_stimuli.Data(rows_to_delete, :) = [];
        resize_table_columns(handles)
    end
        
    INFO = getappdata(handles.figure, 'INFO');
    INFO.is_database_updated = false;
    setappdata(handles.figure, 'INFO', INFO);

 
function btn_sp_behav_affective_Callback(~, ~, handles)
    % Get INFO
    INFO = getappdata(handles.figure, 'INFO');

    % Choose what to do based on whether the button is toggled or not
    switch handles.btn_sp_behav_affective.Value
        case 0
            % Change background color
            set(handles.btn_sp_behav_affective, 'BackgroundColor',[.94, .94, .94])
            % Turn on other buttons
            set([handles.btn_stim_cold, handles.btn_stim_heat, handles.btn_stim_pinprick, handles.btn_stim_touch, handles.btn_stim_puff], 'Enable', 'on')
            
            % Get current time from camera 1 + offset
            timestamp = INFO.cams{1}.CurrentTime - INFO.cameras_offset(1);
            % Convert to string in mm:ss.ms
            time_minutes = floor(timestamp / 60);
            time_seconds = floor(timestamp - time_minutes * 60);
            time_milliseconds = round((timestamp - time_minutes * 60 - time_seconds) * 1000);
            timestamp_str = [num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];

            % Update table
            handles.table_stimuli.Data(end+1, :) = {INFO.last_sp_behav_timestamp, timestamp_str, '', false, ' '};
            % Sort timestamps
            handles.table_stimuli.Data = sortrows(handles.table_stimuli.Data, 1);
            resize_table_columns(handles)
            
            % Reset variable in memory
            INFO.last_sp_behav_timestamp = [];
            
        case 1
            % Change background color
            set(handles.btn_sp_behav_affective, 'BackgroundColor',[0, 1, 0])
            % Turn off other buttons
            set([handles.btn_stim_cold, handles.btn_stim_heat, handles.btn_stim_pinprick, handles.btn_stim_touch, handles.btn_stim_puff], 'Enable', 'off')

            % Get current time from camera 1 + offset
            timestamp = INFO.cams{1}.CurrentTime - INFO.cameras_offset(1);
            % Convert to string in mm:ss.ms
            time_minutes = floor(timestamp / 60);
            time_seconds = floor(timestamp - time_minutes * 60);
            time_milliseconds = round((timestamp - time_minutes * 60 - time_seconds) * 1000);
            timestamp_str = [num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];
            INFO.last_sp_behav_timestamp = timestamp_str;
    end

    % Store INFO
    INFO.is_database_updated = false;
    setappdata(handles.figure, 'INFO', INFO);

    
%% UTILITY FUNCTIONS
function log(msg, varargin)
    global LOGGER
    if isempty(LOGGER)
        disp(msg)
    else
        LOGGER.info(msg, varargin{:})
    end
       

function output_menu = create_nested_menu(output_menu, nested_struct, fieldname, parent_menu, hObject, is_top_level)
    % Make sure that 'fieldname' is a cell array
    if ~iscell(fieldname)
        fieldname = {fieldname};
    end

    % Check whether this is not the deepest level
    if isstruct(nested_struct)
        % Make parent menu
        if ~is_top_level
            fieldname_info = regexprep(fieldname{end}, '^f_', '');
            parent_menu = uimenu(parent_menu, 'Text',fieldname_info);
        end

        % Loop through the fields of this structure
        fn = fieldnames(nested_struct);
        for i_fn = 1:length(fn)
            if is_top_level
                this_fieldname = fn(i_fn);
            else
                this_fieldname = [fieldname, fn(i_fn)];
            end
            output_menu = create_nested_menu(output_menu, nested_struct.(fn{i_fn}), this_fieldname, parent_menu, hObject, false);
        end
        % Stpre parent menu
        if is_top_level
            fieldname{end} = 'main';
        else
            fieldname{end + 1} = 'main';
        end
        output_menu = setfield(output_menu, fieldname{:}, parent_menu);

    % This is the deepest level. Make a menu and store it in the output menu
    else
        % Remove prefix that was used to convert each value to a valid field name
        fieldname_info = regexprep(fieldname, '^f_', '');
        submenu = uimenu(parent_menu, 'Text',fieldname_info{end}, ...
            'MenuSelectedFcn',{@menu_choose_animal_Callback, fieldname_info, hObject}, ...
            'Interruptible','off');
        output_menu = setfield(output_menu, fieldname{:}, submenu);
    end

   
function toggle_menus(handles, action)
    handles.menu_choose_animal.Enable = action;
    menus = flatten_struct(handles.menus);
    set(menus, {'Enable'}, {action})

    
function table_data = prepare_data_for_GUI_table(data)
    table_columns = {'Time_start', 'Time_end', 'Event', 'Withdrawal', 'Affective_response'};
    table_data = cell(size(data, 1), length(table_columns));
    % Convert timestamps
    time_s = cellfun(@(x) convert_time_from_num_to_str(x), num2cell(data{:, 'timestamps'}), 'UniformOutput',false);
    table_data(:, ismember(table_columns, 'Time_start')) = time_s(:, 1);
    table_data(:, ismember(table_columns, 'Time_end')) = time_s(:, 2);
    % Add event type
    table_data(:, ismember(table_columns, 'Event')) = data.stimulus;
    % Add response
    table_data(:, ismember(table_columns, 'Withdrawal')) = num2cell(data.response);
    % Add affective responses
    table_data(:, ismember(table_columns, 'Affective_response')) = data.affective;
    
    % Remove invalid trials
    table_data(~data.valid, :) = [];
    % Fill in empty cells
    table_data(cellfun(@(x) isempty(x), table_data)) = {''};
    
    % Sort rows by Time_start
    table_data = sortrows(table_data, find(ismember(table_columns, 'Time_start')));

    
function s = convert_time_from_num_to_str(n)
    if isnan(n)
        s = '';
    else
        time_minutes = floor(n / 60);
        time_seconds = floor(n - time_minutes * 60);
        time_milliseconds = round((n - time_minutes * 60 - time_seconds) * 1000);
        s = [num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];
    end
    
    
function n = convert_time_from_str_to_num(s)
    if isempty(s)
        n = NaN;
    else
        output = strsplit(s, ':');
        time_minutes = str2double(output{1}) * 60;
        time_seconds = str2double(output{2});
        n = time_minutes + time_seconds;
    end

    
function update_video_frame(figure_handle, slider_value, from_camera)
    try
        % Get INFO
        INFO = getappdata(figure_handle, 'INFO');
        % Update time in VideoReader
        INFO.cams{from_camera}.CurrentTime = slider_value;
        % Grab new frame and show it
        INFO.image_handles{from_camera}.CData = readFrame(INFO.cams{from_camera});
        % Update time counter
        update_time_string(slider_value, INFO.cameras_offset(from_camera), INFO.video_times{from_camera})

        % Also update other camera
        if INFO.videos_in_sync
            other_camera = setdiff(1:INFO.n_cameras, from_camera);
            if ~isempty(other_camera)
                INFO.cams{other_camera}.CurrentTime = INFO.cams{from_camera}.CurrentTime - INFO.cameras_offset(from_camera) + INFO.cameras_offset(other_camera);
                INFO.sliders{other_camera}.Value = INFO.cams{other_camera}.CurrentTime;
                INFO.image_handles{other_camera}.CData = readFrame(INFO.cams{other_camera});
                update_time_string(INFO.cams{other_camera}.CurrentTime, INFO.cameras_offset(other_camera), INFO.video_times{other_camera})
            end
        end
        drawnow
    
    catch ME
        if strcmp(ME.identifier, 'MATLAB:audiovideo:VideoReader:EndOfFile')
            return
        else
            rethrow(ME)
        end
    end
 
    
function play_movie(handles)
    % Do not allow sync button to work while videos are playing
    handles.btn_sync_cameras.Enable = 'off';

    % Get INFO
    INFO = getappdata(handles.figure, 'INFO');

    if INFO.videos_in_sync
        are_cams_playing = [true, true];
        handles.btn_play_cam1.Value = 1;
        handles.btn_play_cam2.Value = 1;
    else
        are_cams_playing = [handles.btn_play_cam1.Value == 1, handles.btn_play_cam2.Value == 1];
    end
    
    % Change text on button
    if are_cams_playing(1), handles.btn_play_cam1.String = '||'; end
    if are_cams_playing(2), handles.btn_play_cam2.String = '||'; end
    % Loop through videos
    while (are_cams_playing(1) && hasFrame(INFO.cams{1})) || ...
          (are_cams_playing(2) && hasFrame(INFO.cams{2}))
        
        idx_playing_cams = find(are_cams_playing);
        for from_camera = idx_playing_cams
            INFO.image_handles{from_camera}.CData = readFrame(INFO.cams{from_camera});
            % Set position of slider and time on clock
            current_time = INFO.cams{from_camera}.CurrentTime;
            INFO.sliders{from_camera}.Value = current_time;
            update_time_string(current_time, INFO.cameras_offset(from_camera), INFO.video_times{from_camera})
        end
        % Update figure
        drawnow
        
        % Check if camera are still playing
        are_cams_playing = [handles.btn_play_cam1.Value == 1, handles.btn_play_cam2.Value];
        % If videos have to stay in sync, but one is not playing anymore, stop both
        if INFO.videos_in_sync && ~all(are_cams_playing)
            break
        end
    end
    % Reset text on button
    handles.btn_play_cam1.String = '>';
    handles.btn_play_cam2.String = '>';

    % Store INFO
    setappdata(handles.figure, 'INFO', INFO);

    % Restore sync button so it works
    handles.btn_sync_cameras.Enable = 'on';


function resize_table_columns(handles)
    % Set size of table columns
    set(handles.table_stimuli, 'ColumnWidth', {100, 100, 90, 40, 200})
    
    
function update_time_string(slider_value, offset, video_times)
    % Subtract offset
    slider_value = slider_value - offset;
    slider_value_abs = abs(slider_value);
    slider_value_sign = sign(slider_value);
    if slider_value_sign < 0
        slider_value_sign = '-';
    else
        slider_value_sign = '+';
    end
    % Calculate time
    time_minutes = floor(slider_value_abs / 60);
    time_seconds = floor(slider_value_abs - time_minutes * 60);
    time_milliseconds = round((slider_value_abs - time_minutes * 60 - time_seconds) * 1000);
    % Update time counter
    video_times.String = [slider_value_sign, num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];


function quit_GUI(~, ~)
    warning('off')
    % Close request function to display a question dialog box
    selection = MessageBox('Stop analysis?', 'Close Request Function', 'Yes','No', 'No');
    switch selection
        case 'Yes'
%             close all force hidden
             close
        case 'No'
            return
    end

    
%#ok<*DEFNU,*TRYNC,*TNMLP,*INUSL>>
