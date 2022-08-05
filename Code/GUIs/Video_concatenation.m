function varargout = Video_concatenation(varargin)
    gui_Singleton = 1;
    gui_State = struct('gui_Name',mfilename, 'gui_Singleton',gui_Singleton, 'gui_OpeningFcn',@Video_concatenation_OpeningFcn, 'gui_OutputFcn',@Video_concatenation_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
    if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
    if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else gui_mainfcn(gui_State, varargin{:}); end

function varargout = Video_concatenation_OutputFcn(~, ~, handles), varargout{1} = handles.output;

    
function Video_concatenation_OpeningFcn(hObject, ~, handles, ~)
    % Set default values for inputs
    handles.txt_animal_ID.String = '';
    handles.check_overwrite.Value = 0;
    handles.table_experiments.Data = cell(0, size(handles.table_experiments.Data, 2));
    handles.btn_run.Enable = 'off';

    % Store handles
    handles.output = hObject;
    guidata(hObject, handles);


function txt_animal_ID_Callback(~, ~, handles)
    % Reset table
    handles.table_experiments.Data = cell(0, size(handles.table_experiments.Data, 2));
    handles.btn_run.Enable = 'off';

    % Check whether this animal is associate to behavioral videos
    animal_ID = handles.txt_animal_ID.String;
    has_behavior = SQL_database.read_table_where('experiments', 'has_behavior', animal_ID,'animal_ID', 'return_as_table',false);
    if isempty(has_behavior)
        log([animal_ID, ' does not exist in the database'])
        return
    end    
    has_behavior = logical(has_behavior);
    if ~has_behavior
        log([animal_ID, ' does not have behavioral videos'])
        return
    end
    
    % Update list of sessions
    METADATA = SQL_database.read_table_where('sessions', {'date', 'experimental_condition'}, animal_ID,'animal_ID');
    experiments = uniqueCellRows(METADATA{:, {'date', 'experimental_condition'}});
    data = cell(size(experiments, 1), size(handles.table_experiments.Data, 2));
    data(:, 1) = num2cell(true(size(experiments, 1), 1));
    data(:, [2, 3]) = experiments;
    handles.table_experiments.Data = data;
    handles.btn_run.Enable = 'on';
    
    
function btn_run_Callback(~, ~, handles)
    global GC
    
    % Turn off GUI
    handles_to_toggle = {'txt_animal_ID', 'table_experiments', 'check_overwrite', 'check_remove_raw', 'check_copy', 'btn_run'};
    for i_handle = 1:length(handles_to_toggle)
        handles.(handles_to_toggle{i_handle}).Enable = 'off';
    end
    
    % Get state of checkboxes
    overwrite = logical(handles.check_overwrite.Value);
    remove_raw = logical(handles.check_remove_raw.Value);
    do_copy_files = logical(handles.check_copy.Value);

    try
        % Get info on animal and sessions to analyze
        animal_ID = handles.txt_animal_ID.String;
        sessions = handles.table_experiments.Data;
        selected_sessions = cell2mat(sessions(:, 1));
        sessions = sessions(selected_sessions, :);
        n_sessions = size(sessions, 1);

        for i_sess = 1:n_sessions
            % Get info on this session
            this_session_date = sessions{i_sess, 2};
            this_session_experiment = strsplit(sessions{i_sess, 3}, ';');
            this_session_experiment = this_session_experiment{2};

            % Set path to folder where videos are located
            video_path = os.path.join(GC.behavior_raw_video_root_path, GC.behavior_video_subfolder, animal_ID, this_session_date, this_session_experiment);
            video_folder_destination = os.path.join(GC.temp_dir, animal_ID);
            if ~exist(video_folder_destination, 'dir'), mkdir(video_folder_destination), end
            cameras = dir(video_path);
            cameras = {cameras.name}';
            cameras = cameras(cellfun(@(x) startsWith(x, 'cam'), cameras));
            n_cameras = length(cameras);
            if n_cameras == 0
                log(sprintf('No videos were found for %s - %s - %s', animal_ID, this_session_date, this_session_experiment))
                continue
            else
                log(sprintf('%s - %s - %s: found videos from %i cameras', animal_ID, this_session_date, this_session_experiment, n_cameras))
            end
%             folder_items = dir(video_path);
%             if sum(endsWith({folder_items.name}, '.mp4')) == 2 && ~overwrite
%                 log(sprintf('Videos already concatenated for %s - %s - %s', animal_ID, this_session_date, this_session_experiment))
%                 continue
%             end
            % Loop through cameras
            for i_cam = 1:n_cameras
                % Get video
                camera_folder = os.path.join(video_path, cameras{i_cam});
                videos = dir(camera_folder);
                videos = {videos.name}';
                videos = videos(cellfun(@(x) endsWith(x, '.mp4'), videos));
                n_videos = length(videos);
                if n_videos == 0
                    log(sprintf('\tNo videos from camera %i', i_cam))
                    continue
                end

                % Make new filename
                f = strsplit(videos{1}, '_');
                filename = os.path.join(video_path, [f{1}, '_', cameras{i_cam}, '.mp4']);

                if ~exist(filename, 'file') || overwrite
                    % Prepend folder
                    videos = strcat(camera_folder, filesep, videos);

                    if n_videos == 1
                        log(sprintf('Copying ''%s''', filename), 'contains_path',true)
                        copyfile(videos{1}, filename)

                    else
                        videos_list_filename = [tempname(), '.txt'];
                        % Write videos list to disk
                        fileID = fopen(videos_list_filename, 'w');
                        for i_video = 1:n_videos
                            fprintf(fileID, 'file ''%s''\n', videos{i_video});
                        end
                        fclose(fileID);

                        % Run conversion in FFmpeg
                        cmd = sprintf('%s -safe 0 -f concat -i %s -c copy "%s"', GC.FFmpeg_exe, videos_list_filename, filename);
                        system(cmd)
                        delete(videos_list_filename)
                    end

                else
                    log(sprintf('\tAlready processed videos from camera %i', i_cam))
                end
                
                % Delete folder with raw videos
                if remove_raw
                    log('\tRemoving raw files ...')
                    rmdir(camera_folder, 's')
                    log('\t\tdone')
                end
                
                % Copy concatenated file to local folder
                if do_copy_files
                    log(sprintf('\tCopying video to ''%s''', video_folder_destination), 'contains_path',true)
                    copyfile(filename, video_folder_destination)
                    log('\t\tdone')
                end
            end
        end
        
    catch ME  % Grab exception
        log([ME.identifier, ': ', ME.stack(1).name, '\n', ME.message], 'contains_path',true)
    end
    
    % Re-enable GUI
    for i_handle = 1:length(handles_to_toggle)
        handles.(handles_to_toggle{i_handle}).Enable = 'on';
    end
    log('done')
    
    
function log(msg, varargin)
    % Print message
    global LOGGER
    if isempty(LOGGER)
        disp(msg)
    else
        LOGGER.info(msg, varargin{:})
    end

    
%% MLint exceptions
%#ok<*SEPEX,*DEFNU>
