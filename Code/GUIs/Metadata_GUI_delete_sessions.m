function varargout = Metadata_GUI_delete_sessions(varargin)
% Last Modified by GUIDE v2.5 16-Aug-2018 14:28:39
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, 'gui_Singleton',gui_Singleton, 'gui_OpeningFcn',@Metadata_GUI_delete_sessions_OpeningFcn, 'gui_OutputFcn',@Metadata_GUI_delete_sessions_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else, gui_mainfcn(gui_State, varargin{:}); end
% End initialization code - DO NOT EDIT

%% INIT
function varargout = Metadata_GUI_delete_sessions_OutputFcn(hObject, ~, ~), varargout{1}=hObject;

function Metadata_GUI_delete_sessions_OpeningFcn(hObject, ~, handles)
    % Make global structure
    INFO = struct();

    % Make dropdown menu with list of experiments
    experiment_names = SQL_database.read_table_where('experiments', 'experiment_name');
    experiment_names = unique(table2cell(experiment_names));
    % Get list of animals for each experiment
    experiments_list = cell(length(experiment_names), 2);
    for iexp = 1:length(experiment_names)
        animal_IDs = table2cell(SQL_database.read_table_where('experiments', 'animal_ID', experiment_names{iexp}, 'experiment_name'));
        experiments_list{iexp, 1} = experiment_names{iexp};
        experiments_list{iexp, 2} = animal_IDs;
    end
    INFO.experiments_list = experiments_list;
    % Update dropdown menu in GUI
    set(handles.menu_experiment_name, 'String', experiment_names)
    set(handles.menu_experiment_name, 'Value', 1)
    % Store selected experiment in memory
    INFO.selected_experiment = experiment_names{1};
    setappdata(handles.figure, 'INFO', INFO);
    setappdata(handles.figure,'INFO',INFO);

    % Update GUI data
    handles.output = hObject;
    handles.figure.CloseRequestFcn = @quit_GUI;
    % Store handles in object
    guidata(hObject, handles);

    % Make dropdown menu with list of animals
    update_list_of_animals(handles)
    % Update table data
    update_list_of_sessions(handles)
    
    % Bring figure to front
    figure(handles.figure)
    
    % Update list of graphic handles that can be turned on and off
    toggle_GUI(handles, 'list')
    figure(handles.figure)
    
%% RUN
function btn_save_Callback(~, ~, handles)
    % Get LOGGER
    global LOGGER

    % Turn GUI off
    toggle_GUI(handles, 'off')
    
    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    % Get current table
    keep = cell2mat(handles.table_experiments.Data(:, 1));
    % Log action
    LOGGER.info([INFO.selected_animal_ID, ': ', num2str(sum(keep)), ' sessions to keep, ', num2str(sum(~keep)), ' to discard'])
    % Update database
    SQL_database.update('sessions', 'keep',keep, INFO.sessions_table(:, 'session_id'),'session_id')
    % Log outcome
    LOGGER.info('Updated database')
    
    % Re-enable GUI controls
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;
    figure(handles.figure)
    
%% BUTTONS
function menu_experiment_name_Callback(~, ~, handles)
    % Turn GUI off
    toggle_GUI(handles, 'off')
    % Update GUI
    update_list_of_animals(handles)
    update_list_of_sessions(handles)
    % Turn GUI on
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;
    figure(handles.figure)    

function menu_animal_ID_Callback(~, ~, handles)
    % Turn GUI off
    toggle_GUI(handles, 'off')
    % Update GUI
    update_list_of_sessions(handles)
    % Turn GUI on
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;
    figure(handles.figure)

%% UTILITIES
function update_list_of_animals(handles)
    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    
    % Update name of selected experiment
    list_exps = get(handles.menu_experiment_name, 'String');
    INFO.selected_experiment = list_exps{get(handles.menu_experiment_name, 'Value')};
    % Get list of animals for this experiment
    list_of_animal_IDs = natsort(INFO.experiments_list{ismember(INFO.experiments_list(:,1), INFO.selected_experiment), 2});
    % Update dropdown menu in GUI
    set(handles.menu_animal_ID, 'Value', 1)  % Automatically select first one
    set(handles.menu_animal_ID, 'String', list_of_animal_IDs)
    % Store selected animal name in memory
    INFO.selected_animal_ID = list_of_animal_IDs{end};
    setappdata(handles.figure, 'INFO', INFO);
    
function update_list_of_sessions(handles)
    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    
    % Update name of selected animal
    list_animals = get(handles.menu_animal_ID, 'String');
    INFO.selected_animal_ID = list_animals{get(handles.menu_animal_ID, 'Value')};
    
    % Reset table
    handles.table_experiments.Data = cell(0,4);
    % Get state of each session from database
    kept_table = SQL_database.read_table_where('sessions', {'session_id', 'keep', 'date', 'experimental_condition', 'stimulus'}, INFO.selected_animal_ID, 'animal_ID', 'split_experimental_condition',false, 'return_all_sessions',true);
    % Sort by date
    kept_table = sortrows(kept_table, 'date');
    % Add columns that are in uitable
    kept_table.keep = logical(kept_table.keep);
    % Fill in data
    sessions_table = kept_table(:, {'keep', 'date', 'experimental_condition', 'stimulus'});
    handles.table_experiments.Data = table2cell(sessions_table);
    INFO.sessions_table = kept_table;
    setappdata(handles.figure, 'INFO', INFO);

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
