function varargout = Metadata_GUI_sessions(varargin)
% Last Modified by GUIDE v2.5 27-Jun-2018 09:01:11
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, 'gui_Singleton',gui_Singleton, 'gui_OpeningFcn',@Metadata_GUI_sessions_OpeningFcn, 'gui_OutputFcn',@Metadata_GUI_sessions_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else, gui_mainfcn(gui_State, varargin{:}); end
% End initialization code - DO NOT EDIT

%% INIT
function varargout = Metadata_GUI_sessions_OutputFcn(hObject, ~, ~), varargout{1}=hObject;

function Metadata_GUI_sessions_OpeningFcn(hObject, ~, handles, varargin)
    % Unpack inputs
    INFO = varargin{1};

    % Update dropdown menu in GUI
    set(handles.menu_experiment_name, 'String', INFO.experiments_list(:, 1))
    set(handles.menu_experiment_name, 'Value', 1)
    % Store selected experiment in memory
    INFO.selected_experiment = INFO.experiments_list{1, 1};
    setappdata(handles.figure, 'INFO', INFO);

    % Make dropdown menu with list of animals
    update_list_of_animals(handles)
    
    % Update table data
    update_list_of_sessions(handles)
    
    setappdata(handles.figure, 'INFO', INFO);
    % Update GUI data
    handles.output = hObject;
    % Store handles in object
    guidata(hObject, handles);

    % Bring figure to front
    figure(handles.figure)

    % Update list of graphic handles that can be turned on and off
    toggle_GUI(handles, 'list')    


%% RUN
function btn_add_to_database_Callback(~, ~, handles)
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Read INFO and METADATA
    INFO = getappdata(handles.figure, 'INFO');
    
    % Read only selected sessions
    selected_experiments = handles.table_experiments.Data(:,2);
    selected_experiments = cellfun(@(x) ~isempty(x)&&x==1, selected_experiments);
    if any(selected_experiments)
        % Edit the info structure
        INFO.experiments = INFO.sessions_table(selected_experiments, :);

        try
            % Create progress bar and start
            done = Metadata_workflow.main(INFO, 'add');
            waitfor(done);  % Stop execution until workflow is done
            
            % Update list of sessions
            update_list_of_sessions(handles)
            
            % Re-enable GUI controls
            toggle_GUI(handles, 'on')
            handles.figure.CloseRequestFcn = @quit_GUI;
            figure(handles.figure)

        catch ME % any exception
            % Delete figure
            delete(handles.figure)
            % Rethrow exception
            disp([ME.identifier, ': ', ME.stack(1).name])
            disp(ME.message)
            rethrow(ME)
        end
    end    
    
%% BUTTONS
function btn_select_all_Callback(~, ~, handles)
    handles.table_experiments.Data(:,2) = {true};

function btn_select_none_Callback(~, ~, handles)
    handles.table_experiments.Data(:,2) = {false};

function menu_experiment_name_Callback(~, ~, handles)
    % Turn GUI off
    toggle_GUI(handles, 'off')
    % Update GUI
    update_list_of_animals(handles)
    update_list_of_sessions(handles)
    % Turn GUI on
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;

function menu_animal_ID_Callback(~, ~, handles)
    % Turn GUI off
    toggle_GUI(handles, 'off')
    % Update GUI
    update_list_of_sessions(handles)
    % Turn GUI on
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;

%% UTILITIES
function update_list_of_animals(handles)
    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    
    % Update name of selected experiment
    list_exps = get(handles.menu_experiment_name, 'String');
    INFO.selected_experiment = list_exps{get(handles.menu_experiment_name, 'Value')};
    % Get list of animals for this experiment
    list_of_animal_IDs = SQL_database.read_table_where('experiments', 'animal_ID', INFO.selected_experiment,'experiment_name', 'return_as_table',false);
    if ~iscell(list_of_animal_IDs)
        list_of_animal_IDs = {list_of_animal_IDs};
    end
    if length(list_of_animal_IDs) > 1
        list_of_animal_IDs = natsort(list_of_animal_IDs);
    end
    % Update dropdown menu in GUI
    set(handles.menu_animal_ID, 'Value', 1)  % Automatically select first one
    set(handles.menu_animal_ID, 'String', list_of_animal_IDs)
    % Store selected animal name in memory
    INFO.selected_animal_ID = list_of_animal_IDs{end};
    setappdata(handles.figure, 'INFO', INFO);
    
function update_list_of_sessions(handles)
    % Read INFO and general_configs
    global GC
    INFO = getappdata(handles.figure, 'INFO');
    
    % Update name of selected animal
    list_animals = get(handles.menu_animal_ID, 'String');
    INFO.selected_animal_ID = list_animals{get(handles.menu_animal_ID, 'Value')};
    
    % Reset table
    handles.table_experiments.Data = cell(0, 5);
    % Get list of sessions for this animal that are already in the database
    sessions_in_DB = SQL_database.read_table_where('sessions', {'date', 'experimental_condition', 'stimulus'}, INFO.selected_animal_ID, 'animal_ID', 'split_experimental_condition',false, 'return_all_sessions',true);
    sessions_in_DB = unique(sessions_in_DB, 'rows', 'stable');
    if height(sessions_in_DB) > 0
        % Assign values to column 'in_DB'
        sessions_in_DB(:, 'in_DB') = num2cell(true(height(sessions_in_DB),1));
        % Sort columns
        sessions_in_DB = sessions_in_DB(:, {'in_DB','date','experimental_condition','stimulus'});
    end
    
    % Add sessions from metadata
    this_experiment_row = find(ismember(INFO.experiments_list(:, 1), INFO.selected_experiment));
    this_mouse_row = ismember(INFO.experiments_list{this_experiment_row, 2}, INFO.selected_animal_ID);
    data_type = INFO.experiments_list{this_experiment_row, 3}{this_mouse_row};
    INFO.selected_data_type = data_type;
        
    switch INFO.selected_data_type
        case '2p'
            list_of_files = dir([os.path.join(GC.data_root_path, GC.tiff_metadata_path, INFO.selected_animal_ID), '*.csv']);
            list_of_files = {list_of_files.name}';
            % Keep only files that start with the animal ID
            list_of_files = list_of_files(cellfun(@(x) startsWith(x,[INFO.selected_animal_ID, '_']), list_of_files));
            % Sort sessions 'naturally', which should sort them chronologically
            list_of_files = natsort(list_of_files);

            % Make another table for sessions not in the database
            sessions_not_in_DB = cell(0, size(handles.table_experiments.Data, 2) - 1);
            % Read each file and get info out of them
            for icsv = 1:length(list_of_files)
                % Read file
                metadata = Metadata_workflow.read_preprocessing_metadata(os.path.join(GC.data_root_path, GC.tiff_metadata_path, list_of_files{icsv}));

                % Get experimental_condition
                experimental_conditions = table2cell(unique(metadata(:, {'experimental_condition', 'stimulus'})));

                % Get date
                experiment_date = strrep(strrep(list_of_files{icsv}, INFO.selected_animal_ID,''),'_','');
                experiment_date = experiment_date(1:end-4);  % Remove extension

                % Fill table
                for irow = 1:size(experimental_conditions, 1)
                    sessions_not_in_DB(end+1, :) = {false, experiment_date, experimental_conditions{irow,:}};
                end
            end
            sessions_not_in_DB = cell2table(sessions_not_in_DB, 'VariableNames',{'in_DB', 'date', 'experimental_condition', 'stimulus'});
        
        case 'epi'
            list_of_files = dir([os.path.join(GC.data_root_path, GC.movie_metadata_path, INFO.selected_animal_ID), '*.csv']);
            list_of_files = {list_of_files.name}';
            % Keep only files that start with the animal ID
            list_of_files = list_of_files(cellfun(@(x) startsWith(x, INFO.selected_animal_ID), list_of_files));
            % Sort sessions 'naturally', which should sort them chronologically
            list_of_files = natsort(list_of_files);
            
            % Make another table for sessions not in the database
            sessions_not_in_DB = cell(0, size(handles.table_experiments.Data, 2) - 1);
            % Read each file and get info out of them
            for icsv = 1:length(list_of_files)
                % Read file
                this_filename = os.path.join(GC.data_root_path, GC.movie_metadata_path, list_of_files{icsv});
                sessions_not_in_DB = Metadata_workflow.read_preprocessing_metadata(this_filename, true);
                % Add column for 'in_DB'
                sessions_not_in_DB = [num2cell(false(height(sessions_not_in_DB), 1)), table2cell(sessions_not_in_DB)];
            end
            sessions_not_in_DB = cell2table(sessions_not_in_DB, 'VariableNames',{'in_DB', 'date', 'stimulus', 'experimental_condition'});
    end
    
    sessions_not_in_DB = sessions_not_in_DB(:, {'in_DB', 'date', 'experimental_condition', 'stimulus'});
    % Remove sessions that are already in the database
    if height(sessions_in_DB) > 0
        sessions_not_in_DB(ismember(sessions_not_in_DB.date, sessions_in_DB.date) & ismember(sessions_not_in_DB.stimulus,sessions_in_DB.stimulus) & ismember(sessions_not_in_DB.experimental_condition, sessions_in_DB.experimental_condition), :) = [];
    end

    % Concatenate the two tables
    sessions_table = [sessions_in_DB; sessions_not_in_DB];
    % Sort by date
    sessions_table = sortrows(sessions_table, 'date');
    % Store data in memory
    INFO.sessions_table = sessions_table;
    setappdata(handles.figure, 'INFO', INFO)
    % Update table
    data = cell(height(sessions_table), 5);
    data(:, [1, 3:end]) = table2cell(sessions_table);
    data(:,2) = {false};  % Default is not selected
    data(~ismember(sessions_table.in_DB, true),2) = {true};  % Select sessions not in database
    handles.table_experiments.Data = data;

    
%% CLOSE FIGURE
function quit_GUI(handle_figure, ~)
    % Close request function to display a question dialog box
    selection = MessageBox('Are you sure?', 'Close Request Function', 'Yes','No', 'No');
    switch selection
        case 'Yes'
            delete(handle_figure)
        case 'No'
            return
    end
    

%% MLint exceptions
%#ok<*DEFNU,*TRYNC,*AGROW>
