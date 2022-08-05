function varargout = launcher_GUI(varargin)
% Last Modified by GUIDE v2.5 04-Dec-2020 14:18:18
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, 'gui_Singleton',gui_Singleton, 'gui_OpeningFcn',@launcher_GUI_OpeningFcn, 'gui_OutputFcn',@launcher_GUI_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else, gui_mainfcn(gui_State, varargin{:}); end
% End initialization code - DO NOT EDIT

%% INIT
function varargout = launcher_GUI_OutputFcn(hObject, ~, handles)
    varargout = {hObject, handles};

function launcher_GUI_OpeningFcn(hObject, ~, handles)
    % Disable warnings
    warning('off','MATLAB:table:RowsAddedExistingVars')

    % Get LOGGER and general_configs
    global GC LOGGER
    % Log launch of the GUI
    LOGGER.info(['Starting (', char(datetime), ')'])
    
    % Store the INFO structure in memory
    setappdata(handles.figure, 'INFO',struct())
    
    % Set figure renderer to faster OpenGL
    set(handles.figure, 'Renderer','OpenGL');
    % Set tag of window
    set(handles.figure, 'Tag','Ca_imaging_analysis_launcher_GUI');

    % Get names of experiments
    update_list_of_experiments(handles)
    % Update table
    update_list_of_animals(handles)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add new actions (checkboxes and buttons) below
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Copy handles that point to checkboxes in the GUI
    handles.actions = struct();
    % Preprocessing and segmentation section
    handles.actions.registration.main     = handles.check_registration_main;
    handles.actions.registration.skip     = handles.check_skip_image_registration;
    handles.actions.registration.missing  = handles.check_missing_registration;
    % Analysis section
    handles.actions.analysis.main                                       = handles.check_analysis_main;
    handles.actions.analysis.ca_extraction.main                         = handles.check_calcium_extraction_main;
    handles.actions.analysis.spontaneous_activity.main                  = handles.check_analysis_spontaneous_activity;
    handles.actions.analysis.spontaneous_activity.missing               = handles.check_missing_spontaneous_activity;
    handles.actions.analysis.response_modulation.main                   = handles.check_analysis_response_modulation;
    handles.actions.analysis.response_modulation.missing                = handles.check_missing_response_modulation;
    handles.actions.analysis.ensemble_detection.main                    = handles.check_analysis_ensemble_detection;
    handles.actions.analysis.ensemble_detection.missing                 = handles.check_missing_ensemble_detection;
    handles.actions.analysis.response_decoding.main                     = handles.check_analysis_response_decoding;
    handles.actions.analysis.response_decoding.missing                  = handles.check_missing_analysis_response_decoding;
    handles.actions.analysis.response_decoding.all                      = handles.check_analysis_response_decoding_all;
    handles.actions.analysis.response_decoding.selective                = handles.check_analysis_response_decoding_selective;
    handles.actions.analysis.response_decoding.stable                   = handles.check_analysis_response_decoding_stable;
    handles.actions.analysis.response_decoding.ensembles                = handles.check_analysis_response_decoding_ensembles;
    handles.actions.analysis.response_decoding.specific                 = handles.check_analysis_response_decoding_specific;
    handles.actions.analysis.response_decoding.H_noci                   = handles.check_analysis_response_decoding_H_noci;
    handles.actions.analysis.response_decoding.H_noci_random            = handles.check_analysis_response_decoding_H_noci_random;
    handles.actions.analysis.response_decoding.H_sal                    = handles.check_analysis_response_decoding_H_sal;
    handles.actions.analysis.response_decoding.H_sal_random             = handles.check_analysis_response_decoding_H_sal_random;
    handles.actions.analysis.response_decoding.H_sal_categories         = handles.check_analysis_response_decoding_H_sal_categories;
    handles.actions.analysis.response_decoding.H_sal_random_categories  = handles.check_analysis_response_decoding_H_sal_random_categories;
    handles.actions.analysis.response_decoding.all_cell_categories      = handles.check_analysis_response_decoding_all_categories;
    
    handles.actions.analysis.EPM_track.main                             = handles.check_analysis_EPM_track;
    handles.actions.analysis.EPM_track.missing                          = handles.check_missing_EPM_sessions;
    % Visualization section
    handles.actions.visualization.main                                  = handles.check_visualization_main;
    handles.actions.visualization.response_detection.main               = handles.check_plot_response_detection;
    handles.actions.visualization.response_detection.selected           = handles.check_plot_response_detection_only_selected;
    handles.actions.visualization.response_detection.summary            = handles.check_plot_response_detection_summary;
    handles.actions.visualization.response_detection.overwrite          = handles.check_plot_response_detection_overwrite;
    handles.actions.visualization.response_decoding.main                = handles.check_plot_response_decoding;
    handles.actions.visualization.response_decoding.selected            = handles.check_plot_response_decoding_only_selected;
    handles.actions.visualization.response_decoding.summary             = handles.check_plot_response_decoding_summary;
    handles.actions.visualization.response_decoding.overwrite           = handles.check_plot_response_decoding_overwrite;
    handles.actions.visualization.spontaneous_activity.main             = handles.check_plot_spontaneous_activity;
    handles.actions.visualization.spontaneous_activity.selected         = handles.check_plot_spontaneous_activity_only_selected;
    handles.actions.visualization.spontaneous_activity.summary          = handles.check_plot_spontaneous_activity_summary;
    handles.actions.visualization.spontaneous_activity.overwrite        = handles.check_plot_spontaneous_activity_overwrite;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    % Copy values of checkboxes in INFO structure
    read_state_of_checkboxes(handles)
    
    % Load last state of GUI
    load_GUI_settings(handles);

    % Adapt buttons state
    check_section_callback([], [], handles)
    
    % Fix window appearance
    handles = realign_GUI_controls(handles);

    % Update GUI data
    handles.output = hObject;
    % Link callback function when window has to close
    handles.figure.CloseRequestFcn = {@stop_ML_execution, handles};
    % Make sure window is not docked
    handles.figure.WindowStyle = 'normal';
    % Store handles in object
    guidata(hObject, handles);

    % Set callback for when cell selection changes
    set(handles.table_experiments, 'CellEditCallback',{@table_experiments_Callback, handles}, ...
                                   'CellSelectionCallback',{@table_experiments_Callback, handles})

    % Update list of graphic handles that can be turned on and off
    toggle_GUI(handles, 'list')
    
    % Show version of repository in GUI
    handles.txt_version.String = ['v. ' GC.version];
    
%     % Move figure
%     handles.figure.Position([1, 2]) = [50, -200];
    
    % Log end of init
    LOGGER.info('Launcher GUI is initialized')
    

%% Main callback
function btn_run_Callback(~, ~, handles)
    global LOGGER
    
    % Get list of selected experiment
    selected_experiments = handles.table_experiments.Data(:,1);
    selected_experiments = cellfun(@(x) ~isempty(x)&x==1, selected_experiments);
    % If no experiments selected, abort run
    if ~any(selected_experiments), return, end

    % Turn GUI off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Save last state of the GUI
    save_GUI_settings(handles);

    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    INFO.experiments = INFO.sessions_table{selected_experiments, 'animal_ID'};

    try
        % Run workflows
        if INFO.actions.registration.main
            done = Registration_workflow.main(INFO);
            waitfor(done)  % Stop execution until workflow is done                
        end
        if INFO.actions.analysis.main
            done = Analysis_workflow.main(INFO);
            waitfor(done)  % Stop execution until workflow is done
        end
        if INFO.actions.visualization.main
            done = Visualization_workflow.main(INFO);
            waitfor(done)  % Stop execution until workflow is done
        end
        
    catch MException  % Grab exception
        % Log exception
        LOGGER.critical([MException.identifier, ': ', MException.stack(1).name, '\n', MException.message], 'contains_path',true, 'print_on_screen',false)
        close('force')
        rethrow(MException)
    end

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)
 
    
%% GUI options and settings
function save_GUI_settings(handles)
    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    % Read general_configs
    global GC

    % Get user temporary folder location and make filename for GUI's settings
    filename = os.path.join(GC.temp_dir, 'Ca_imaging_data_analysis_launcherGUI_state.mat');
    
    % Get the state of each GUI checkbox    
    GUI_state = [];
    % Get all fields from handles
    gui_fields = fields(handles);
    % Loop through each field
    for ifield = 1:length(gui_fields)
        obj = handles.(gui_fields{ifield});
        % If the object is a scalar and a graphics handle
        if isscalar(obj) && ishghandle(obj) && isfield(get(obj),'Style')
            style = get(obj,'Style');
            if ismember(style, {'radiobutton', 'checkbox'})
                GUI_state.(gui_fields{ifield}) = get(obj, 'Value');
            end
        end
    end
    % Save state of table
    TABLE_state = struct();
    TABLE_state.selected_experiment = INFO.selected_experiment;
    TABLE_state.selected_animals = handles.table_experiments.Data(cell2mat(handles.table_experiments.Data(:,1)),2);
    
    % Get hash of last git commit to know to which version we were running on at
    % the time of save
    [~, git_hash] = system('git rev-parse HEAD');
    git_hash = git_hash(1:end-1);
    
    % Write variables to file
    save(filename, 'GUI_state','TABLE_state','git_hash', '-v7.3', '-nocompression');

    
function load_GUI_settings(handles)
    % Get LOGGER and general_configs
    global LOGGER GC

    % Get temporary filename
    filename = os.path.join(GC.temp_dir, 'Ca_imaging_data_analysis_launcherGUI_state.mat');
    % If file doesn't exist, do nothing
    if ~exist(filename, 'file'), return, end
    
    % Compare saved git hash with current one
    [~, current_git_hash] = system('git rev-parse HEAD');
    current_git_hash = current_git_hash(1:end-1);
    previous_git_hash = load_variable(filename, 'git_hash');
    if ~strcmp(current_git_hash, previous_git_hash)
        % Delete file and return
        delete(filename)
        return
    end
    
    % Log action
    LOGGER.info(['Reading last save from "', filename, '"'], 'contains_path',true)
    
    % Include file loading and setting variables in a try block to ignore errors
    try
        % Load state from file
        GUI_state = load_variable(filename, 'GUI_state');
        % Get all fields from state
        state_fields = fields(GUI_state);
        % Loop through each field
        for ifield = 1:length(state_fields)
            set(handles.(state_fields{ifield}), 'Value', GUI_state.(state_fields{ifield}));
        end
    end
    try
        TABLE_state = load_variable(filename, 'TABLE_state');
        % Read list of experiments
        experiment_names = SQL_database.read_table_where('experiments', 'experiment_name');
        experiment_names = unique(table2cell(experiment_names));
        exp_value = find(ismember(experiment_names, TABLE_state.selected_experiment));
        if ~isempty(exp_value)
            % Switch to previous experiment
            set(handles.menu_experiment_name, 'Value', exp_value);
            update_list_of_animals(handles)
            % Unselect all rows
            handles.table_experiments.Data(:,1) = {false};
            animals_in_table = handles.table_experiments.Data(:,2);
            % Loop through animals in table and check all previously saved ones
            for iid = 1:length(TABLE_state.selected_animals)
                if ismember(TABLE_state.selected_animals{iid}, animals_in_table)
                    handles.table_experiments.Data(ismember(animals_in_table, TABLE_state.selected_animals{iid}),1) = {true};
                end
            end
        end
        % Re-align memory with GUI data
        INFO = getappdata(handles.figure, 'INFO');
        INFO.sessions_table = cell2table(handles.table_experiments.Data, 'VariableNames',INFO.sessions_table.Properties.VariableNames);
        setappdata(handles.figure, 'INFO',INFO);
    end
    
    % Adapt button state
    check_section_callback([],[],handles)
  
    
%% BUTTONS
function btn_select_all_Callback(~, ~, handles)
    handles.table_experiments.Data(:,1) = {true};
    is_run_allowed(handles)

    
function btn_select_none_Callback(~, ~, handles)
    handles.table_experiments.Data(:,1) = {false};
    is_run_allowed(handles)

    
function menu_experiment_name_Callback(~, ~, handles)
    update_list_of_animals(handles)
    is_run_allowed(handles)

    
function menu_metadata_stimuli_description_Callback(~, ~, handles)
    % Turn GUI off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    % Save last state of the GUI
    save_GUI_settings(handles);
    try
        % Run workflow
        done = Stimuli_workflow.main();
        waitfor(done)  % Stop execution until workflow is done
        
    catch MException  % Grab exception
        % Log exception
        global LOGGER
        LOGGER.critical([MException.identifier, ': ', MException.stack(1).name, '\n', MException.message], 'contains_path',true, 'print_on_screen',false)
        close('force')
        rethrow(MException)
    end

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)

    
function table_experiments_Callback(~, ~, handles)
    is_run_allowed(handles)

    
function menu_metadata_update_database_Callback(~, ~, handles)
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Save last state of the GUI
    save_GUI_settings(handles);

    try
        % Read INFO
        INFO = getappdata(handles.figure, 'INFO');
        % Run workflow
        waitfor(Metadata_GUI_launcher(INFO));  % Stop execution until workflow is done

        % Read list of experiments
        experiment_names = SQL_database.read_table_where('experiments', 'experiment_name');
        experiment_names = unique(table2cell(experiment_names));

        % If there are more experiments than before, the user has added a new one.
        % Select the new experiment
        if length(experiment_names) > size(INFO.experiments_list, 1)
            new_experiment = setdiff(experiment_names, INFO.experiments_list(:,1));
            INFO.selected_experiment = new_experiment{1};
        end
        % Update list of experiments
        update_list_of_experiments(handles)
        % Update dropdown menu in GUI
        set(handles.menu_experiment_name, 'Value', find(ismember(experiment_names, INFO.selected_experiment)))
        % Update list of animals
        update_list_of_animals(handles)
                
    catch MException  % Grab exception
        % Log exception
        global LOGGER
        LOGGER.critical([MException.identifier, ': ', MException.stack(1).name, '\n', MException.message], 'contains_path',any(strfind(MException.message, filesep())), 'print_on_screen',false)
        close('force')
        rethrow(MException)
    end

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)

  
function menu_metadata_discard_sessions_Callback(~, ~, handles)
    % Get LOGGER
    global LOGGER
        
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
        
    try
        % Run workflow
        LOGGER.info('GUI to keep / discard sessions from analysis has started')
        done = Metadata_GUI_delete_sessions();
        waitfor(done)  % Stop execution until workflow is done
        LOGGER.info('GUI to keep / discard sessions from analysis has quit')
        
    catch MException  % Grab exception
        % Log exception
        LOGGER.critical([MException.identifier, ': ', MException.stack(1).name, '\n', MException.message], 'contains_path',any(strfind(MException.message, filesep())), 'print_on_screen',false)
        close('force')
        rethrow(MException)
    end

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)

    
function btn_run_segmentation_Callback(~, ~, handles)
    % Get LOGGER
    global LOGGER
    
    % Get list of selected experiment
    selected_experiments = handles.table_experiments.Data(:,1);
    selected_experiments = cellfun(@(x) ~isempty(x)&x==1, selected_experiments);
    % If no experiments selected, abort run
    if ~any(selected_experiments), return, end
    
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Save last state of the GUI
    save_GUI_settings(handles);

    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    experiments = INFO.sessions_table{selected_experiments, 'animal_ID'};
    % Update INFO
    INFO.experiments = experiments;
    
    try
        % Run workflow
        done = ROI_selection_workflow.main(INFO);
        waitfor(done)  % Stop execution until workflow is done
        
    catch MException  % Grab exception
        % Log exception
        LOGGER.critical([MException.identifier, ': ', MException.stack(1).name, '\n', MException.message], 'contains_path',any(strfind(MException.message, filesep())), 'print_on_screen',false)
        close('force')
        rethrow(MException)
    end

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)

    
function btn_ca_event_threshold_Callback(~, ~, handles)
    % Get LOGGER
    global LOGGER
    
    % Get list of selected experiment
    selected_experiments = handles.table_experiments.Data(:,1);
    selected_experiments = cellfun(@(x) ~isempty(x)&x==1, selected_experiments);
    % If no experiments selected, abort run
    if ~any(selected_experiments), return, end
    
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Save last state of the GUI
    save_GUI_settings(handles);

    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');
    experiments = INFO.sessions_table{selected_experiments, 'animal_ID'};
    % Update INFO
    INFO.experiments = experiments;
    
    try
        % Run workflow
        for i_exp = 1:length(INFO.experiments)
            animal_ID = INFO.experiments{i_exp};
            done = Ca_event_detection(animal_ID);
            waitfor(done)  % Stop execution until workflow is done
        end
        
    catch MException  % Grab exception
        % Log exception
        LOGGER.critical([MException.identifier, ': ', MException.stack(1).name, '\n', MException.message], 'contains_path',any(strfind(MException.message, filesep())), 'print_on_screen',false)
        close('force')
        rethrow(MException)
    end

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)

   
function check_section_callback(~, ~, handles)
    % Read state of checkboxes
    read_state_of_checkboxes(handles)
    % Adapt state of buttons
    adapt_button_state(handles)
    % Enable Run button
    is_run_allowed(handles)


function check_subsection_callback(~, ~, handles)
    % Read state of checkboxes
    read_state_of_checkboxes(handles)
    % Adapt state of buttons
    adapt_button_state(handles)
    
    
%% CALLBACK MENUs
function menu_concatenate_videos_Callback(~, ~, handles)
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Save last state of the GUI
    save_GUI_settings(handles);

    % Run GUI
    waitfor(Video_concatenation());  % Stop execution until workflow is done

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)


function menu_score_behavior_Callback(~, ~, handles)    
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Save last state of the GUI
    save_GUI_settings(handles);

    % Run GUI
    waitfor(GUI_mark_stimuli());  % Stop execution until workflow is done

    % Re-enable GUI controls
    toggle_GUI(handles, 'on', {@stop_ML_execution, handles})
    figure(handles.figure)
    

%% UTILITY FUNCTIONS
function read_state_of_checkboxes(handles)
    % Read INFO from memory
    INFO = getappdata(handles.figure, 'INFO');
    % Add actions field if not present
    if ~isstruct(INFO)
        INFO = struct('actions', struct());
    else
        if ~isfield(INFO, 'actions')
            INFO.actions = struct();
        end
    end
    
    % Loop through the checkboxes handles and store their value
    flds = fieldnames(handles.actions);
    for ii = 1:length(flds)
        INFO.actions.(flds{ii}) = operate_on_checkboxes(handles.actions.(flds{ii}), 'get', 'Value');
    end
    % Store data in memory
    setappdata(handles.figure, 'INFO',INFO);
  
    
function output = operate_on_checkboxes(sub_handles, get_set, field)
    if isstruct(sub_handles)
        flds = fieldnames(sub_handles);
        for ii = 1:length(flds)
            output.(flds{ii}) = operate_on_checkboxes(sub_handles.(flds{ii}), get_set, field);
        end
    else
        if strcmp(get_set, 'get')
            output = get(sub_handles, field);
        else
            set(sub_handles, field{1}, field{2});
            output = [];
        end
    end
   
    
function adapt_button_state(handles)
    % Read INFO from memory
    INFO = getappdata(handles.figure, 'INFO');
    % Get names of all handles
    all_handle_names = fieldnames(handles);
    % Get action sections
    sections = fieldnames(INFO.actions);
    for ii = 1:length(sections)
        current_section = sections{ii};
        % Toggle accessories of main checkbox
        fn = all_handle_names(~cellfun(@isempty, regexp(all_handle_names, [current_section, '_main_accessory'])));
        if INFO.actions.(current_section).main
            status = 'on';
        else
            status = 'off';
        end
        for ifield = 1:length(fn)
            set(handles.(fn{ifield}), 'Enable', status)
        end

        if INFO.actions.(current_section).main
            % Enable all checkboxes under this section
            operate_on_checkboxes(handles.actions.(current_section), 'set', {'Enable','on'});
        else
            % Disable all checkboxes under this section
            operate_on_checkboxes(handles.actions.(current_section), 'set', {'Enable','off'});
            % Re-enable main checkbox
            set(handles.actions.(current_section).main, 'Enable','on')
        end
        is_enabled = INFO.actions.(current_section).main;
        if is_enabled
            is_enabled_str = 'on';
        else
            is_enabled_str = 'off';
        end
        % Toggle subsections (apart from 'main')
        subsections = fieldnames(INFO.actions.(current_section));
        subsections(ismember(subsections, 'main')) = [];
        if isempty(subsections)
            continue
        end
        for jj = 1:length(subsections)
            % Toggle value of main checkbox in this subsection
            current_subsection = subsections{jj};
            if isstruct(INFO.actions.(current_section).(current_subsection))
                subsection_fields = fieldnames(INFO.actions.(current_section).(current_subsection));
                if ismember('main', subsection_fields)
                    value_main = INFO.actions.(current_section).(current_subsection).main;
                    if value_main
                        value_main_str = 'on';
                    else
                        value_main_str = 'off';
                    end
                    set(handles.actions.(current_section).(current_subsection).main, 'Value',value_main)
                    % Toggle other checkboxes
                    subsection_fields(ismember(subsection_fields, 'main')) = [];
                    if isempty(subsection_fields)
                        continue
                    end
                    for kk = 1:length(subsection_fields)
                        if is_enabled
                            value = INFO.actions.(current_section).(current_subsection).(subsection_fields{kk});
                            set(handles.actions.(current_section).(current_subsection).(subsection_fields{kk}), 'Value',value, 'Enable',value_main_str)
                        else
                            value_main_str = is_enabled_str;
                        end
                        % Toggle accessories
                        fn = all_handle_names(~cellfun(@isempty, regexp(all_handle_names, [current_section, '_', current_subsection, '_', subsection_fields{kk}, '_accessory'])));
                        for ifield = 1:length(fn)
                            set(handles.(fn{ifield}), 'Enable', value_main_str)
                        end
                    end
                end
            end
        end
        
        % Apply exceptions
        if strcmp(current_section, 'registration') && is_enabled
            statuses = {'on', 'on'};
            if INFO.actions.registration.skip
                statuses = {'on', 'off'};
            elseif INFO.actions.registration.missing
                statuses = {'off', 'on'};
            end
            set(handles.actions.registration.skip, 'Enable',statuses{1})
            set(handles.actions.registration.missing, 'Enable',statuses{2})
        end
        
        % At least one subsection needs to be enabled
        if strcmp(current_section, 'analysis') && is_enabled
            value_main = INFO.actions.analysis.response_decoding.main;
            if value_main
                fn = {'all', 'selective', 'stable', 'ensembles', 'specific', 'H_noci','H_noci_random', 'H_sal', 'H_sal_random'};
                states = false(1, 4);
                for ifield = 1:length(fn)
                    states(ifield) = INFO.actions.analysis.response_decoding.(fn{ifield});
                end
                if all(~states)
                    % Turn on "All cells"
%                     INFO.actions.analysis.response_decoding.all = 1;
%                     set(handles.actions.analysis.response_decoding.all, 'Value',1)
                end
            end            
        end        
    end
  
    
function is_run_allowed(handles)
    % Get INFO
    INFO = getappdata(handles.figure, 'INFO');
    % Find out whether we have all the elements to run the analysis:
    % At least one experiment selected
    selected_experiments = cell2mat(handles.table_experiments.Data(:,1));
    if isempty(selected_experiments)
        set(handles.btn_run_segmentation, 'Enable', 'off')
        return
    end
    n_sessions = INFO.sessions_table.n_sessions(selected_experiments);
    if isempty(n_sessions)
        n_sessions = 0;
    end
    any_experiment = sum(n_sessions) > 0;
    
    % At least one workflow selected
    fldnms = fieldnames(INFO.actions);
    workflows_status = false(length(fldnms),1);
    for ifld = 1:length(fldnms)
        workflows_status(ifld) = INFO.actions.(fldnms{ifld}).main;
    end
    % Check whether at least 1 step of a selected workflow is selected
    workflows_idx = find(workflows_status);
    any_workflow_step = false;
    for iwf = 1:length(workflows_idx)
        sub_fldnms = fieldnames(INFO.actions.(fldnms{workflows_idx(iwf)}));
        content = struct2cell(INFO.actions.(fldnms{workflows_idx(iwf)}));
        % Apply exceptions
        if ~strcmp(fldnms{workflows_idx(iwf)}, 'registration')
            % Ignore 'main' button alone; here we are checking for sub-sections
            content(ismember(sub_fldnms, 'main'), :) = [];
        end
        for ii = 1:length(content)
            if isstruct(content{ii})
                content{ii} = content{ii}.main;
            end
        end
        content = cell2mat(content);
        any_workflow_step = any(content);
        if any_workflow_step
            break
        end
    end
    % Set whether the Run button is enabled
    if any(workflows_status) && any_workflow_step && any_experiment
        new_state = 'on';
    else
        new_state = 'off';
    end
    % Toggle button accordingly
    set(handles.btn_run, 'Enable', new_state)
    
    % Set whether the Run Segmentation button is enabled
    if any_experiment
        new_state = 'on';
    else
        new_state = 'off';
    end
    set(handles.btn_run_segmentation, 'Enable', new_state)

    
function stop_ML_execution(~, ~, handles)
    % Close request function to display a question dialog box
    selection = MessageBox('Stop analysis?', 'Close Request Function', 'Yes','No', 'No');
    switch selection
        case 'Yes'
            % Save state of GUI
            save_GUI_settings(handles)
            % Close GUI
            close('force')
            % Get LOGGER
            global LOGGER
            % Log when GUI was closed
            LOGGER.info(['Finished (', char(datetime), ')'])
        case 'No'
            return
    end
   
    
function update_list_of_experiments(handles)
    % Read INFO
    INFO = getappdata(handles.figure, 'INFO');

    % Read list of experiments
    experiment_names = SQL_database.read_table_where('experiments', 'experiment_name');
    experiment_names = unique(table2cell(experiment_names));
    % Get list of animals for each experiment
    experiments_list = cell(length(experiment_names), 3);
    for iexp = 1:length(experiment_names)
        animal_IDs = table2cell(SQL_database.read_table_where('experiments', {'animal_ID', 'data_type'}, experiment_names{iexp}, 'experiment_name'));
        experiments_list{iexp, 1} = experiment_names{iexp};
        experiments_list{iexp, 2} = animal_IDs(:, 1);
        experiments_list{iexp, 3} = animal_IDs(:, 2);
    end
    INFO.experiments_list = experiments_list;
    
    % Update dropdown menu in GUI
    set(handles.menu_experiment_name, 'String', experiment_names)
    set(handles.menu_experiment_name, 'Value', 1)
    
    % Store selected experiment in memory
    if isempty(experiment_names)
        INFO.selected_experiment = {''};
    else
        INFO.selected_experiment = experiment_names{1};
    end
    setappdata(handles.figure, 'INFO', INFO);
    
    % Store experiment name in general configs
    global GC
    GC.experiment_name = INFO.selected_experiment;

    
function update_list_of_animals(handles)
    % Read INFO and general configs
    INFO = getappdata(handles.figure, 'INFO');
    global GC
    
    % Empty experiments table
    column_names = handles.table_experiments.ColumnName;
    n_columns = length(column_names);
    handles.table_experiments.Data = cell(0, n_columns);

    % Update name of selected experiment
    list_exps = get(handles.menu_experiment_name, 'String');
    if isempty(list_exps)
        INFO.selected_experiment = {''};
        INFO.sessions_table = cell(0, n_columns);
    else
        INFO.selected_experiment = list_exps{get(handles.menu_experiment_name, 'Value')};
        % Get list of animals for this experiment
        list_of_animal_IDs = INFO.experiments_list{ismember(INFO.experiments_list(:,1), INFO.selected_experiment), 2};
        % Fill table
        data = cell(0, n_columns);
        for iid = 1:length(list_of_animal_IDs)
            % Get info on experiment
            animal_ID = list_of_animal_IDs{iid};
            exp_metadata = SQL_database.read_table_where('experiments', [], animal_ID, 'animal_ID');

            % Get number of sessions
            n_sessions = SQL_database.count_sessions(animal_ID);
            the_selected = n_sessions > 0;  % Automatically select if there is at least one session
            
            % Fill in table
            data{iid, 1} = the_selected;
            data{iid, 2} = animal_ID;
            data{iid, 3} = exp_metadata.experimental_group{1};
            data{iid, 4} = n_sessions;
            data{iid, 5} = exp_metadata.brain_area{1};
            data{iid, 6} = exp_metadata.genotype{1};
            data{iid, 7} = exp_metadata.Calcium_indicator{1};
            data{iid, 8} = exp_metadata.promoter{1};
            data{iid, 9} = exp_metadata.data_type{1};
        end
        % Sort by animal ID
        data = natsortrows(data, 2);
        % Store data in memory as table
        INFO.sessions_table = cell2table(data, 'VariableNames',make_variable_name(column_names));
        % Update table data
        handles.table_experiments.Data = data;
    end
    
    % If there is at least 1 row and colors have been provided for experimental groups, enable row striping
    n_rows = size(handles.table_experiments.Data, 1);
    % Get color of each group, replacing white where no color was provided
    BackgroundColor = ones(n_rows, 3);

    if n_rows > 1 && ~isempty(GC.experimental_groups_colors.(INFO.selected_experiment))
        % Get row where colors are
        all_colors = GC.experimental_groups_colors.(INFO.selected_experiment);
        for irow = 1:n_rows
            idx = ismember(all_colors(:,1), INFO.sessions_table.group{irow});
            if any(idx)
                % Convert original color to HSV
                color = rgb2hsv(all_colors{idx,2});
                % Desaturate color (the second number)
                color(2) = color(2) * 0.5;
                % Store value in original RGB form
                BackgroundColor(irow, :) = hsv2rgb(color);
            end
        end
    end
    % Set colors for groups
    col_idx = find(ismember(handles.table_experiments.ColumnName, 'group'));
    for irow = 1:n_rows
        handles.table_experiments.Data{irow, col_idx} = add_background_color_to_table_cell(handles.table_experiments.Data{irow, col_idx}, BackgroundColor(irow,:));
    end
    % Set width of all columns
    handles.table_experiments.ColumnWidth = {50, 100, 75, 70, 70, 80, 100, 75, 60};
    setappdata(handles.figure, 'INFO', INFO);

    
function new_text = add_background_color_to_table_cell(text, color)
    if ~exist('color','var') || isempty(color)
        color = [1, 1, 1];
    end
    % Transform to hex
    color = rgb2hex(color);
    % Make HTML code
    new_text = ['<html><table border=0 width=400 bgcolor=', color, '><TR><TD>', text, '</TD></TR></table></html>'];
    
    
function handles = realign_GUI_controls(handles)
    % Set parameters
    unit = 'pixels';
    title_top_padding = 0;
    title_bottom_padding = 10;
    gui_panels_padding = 20;
    gui_panels_left = 20;
    gui_panels_right = 20;
    subpanel_padding_top = 50;
    panel_title_padding_left = 15;
    panel_experiments_vertical_offset = 11;
    subpanel_padding_left = 20;
    subpanel_padding_right = 20;
    subpanel_submenu_bottom = 20;
    subpanel_submenu_left_column = 320;
    subpanel_submenu_submenu_left = 50;
    subpanel_submenu_padding = 10;
    subpanel_accessory_padding = 0;
    subpanel_accessory_padding_comma = -10;
    
    % Set unit for all graphical elements
    fn = fieldnames(handles);
    for i_fn = 1:length(fn)
        if isprop(handles.(fn{i_fn}), 'Units')
            set(handles.(fn{i_fn}), 'Units',unit)
        end
    end
    
    % Make figure large enough
    handles.figure.Position(3) = 1500;
    handles.figure.Position(4) = 900;
    drawnow;
    
    % Adjust position of title
    handles.txt_main_title.Position([1, 2, 3]) = [0, handles.figure.Position(4) - handles.txt_main_title.Position(4) - title_top_padding, handles.figure.Position(3)];
    handles.txt_version.Position([1, 2]) = [handles.figure.Position(3) - handles.txt_version.Position(3), handles.figure.Position(4) - handles.txt_version.Position(4)];
    
    % --- Preprocessing ---
    % Set position of inner elements
    handles.check_registration_main.Position([1, 2]) = [subpanel_padding_left, handles.panel_preprocessing.Position(4) - subpanel_padding_top];
    handles = sort_on_same_row(handles, 'check_registration_main', {'registration_main_accessory1', 'check_missing_registration', 'registration_main_accessory2', 'check_skip_image_registration', 'registration_main_accessory3'}, subpanel_submenu_left_column, [5, subpanel_accessory_padding, subpanel_accessory_padding, subpanel_submenu_padding, -5]);
    handles.btn_run_segmentation.Position = [subpanel_padding_left, subpanel_submenu_padding, handles.panel_preprocessing.Position(3) - subpanel_padding_left - subpanel_padding_right, handles.check_registration_main.Position(2) - subpanel_submenu_padding * 2];
    
    % --- Analysis ---
    % Set position of inner elements
    handles.check_calcium_extraction_main.Position([1, 2])        = [subpanel_padding_left, handles.panel_analysis.Position(4) - subpanel_padding_top];
    handles.check_analysis_spontaneous_activity.Position([1, 2])  = [subpanel_padding_left, handles.check_calcium_extraction_main.Position(2) - handles.check_calcium_extraction_main.Position(4)];
    handles = sort_on_same_row(handles, 'check_analysis_spontaneous_activity', {'analysis_spontaneous_activity_missing_accessory1', 'check_missing_spontaneous_activity', 'analysis_spontaneous_activity_missing_accessory2'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, -10]);
    handles.check_analysis_response_modulation.Position([1, 2])   = [subpanel_padding_left, handles.check_analysis_spontaneous_activity.Position(2) - handles.check_analysis_spontaneous_activity.Position(4)];
    handles = sort_on_same_row(handles, 'check_analysis_response_modulation', {'analysis_response_modulation_missing_accessory1', 'check_missing_response_modulation', 'analysis_response_modulation_missing_accessory2'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, -10]);
    handles.check_analysis_ensemble_detection.Position([1, 2])    = [subpanel_padding_left, handles.check_analysis_response_modulation.Position(2) - handles.check_analysis_response_modulation.Position(4)];
    handles = sort_on_same_row(handles, 'check_analysis_ensemble_detection', {'analysis_ensemble_detection_missing_accessory1', 'check_missing_ensemble_detection', 'analysis_ensemble_detection_missing_accessory2'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, -10]);
    handles.check_analysis_response_decoding.Position([1, 2])     = [subpanel_padding_left, handles.check_analysis_ensemble_detection.Position(2) - handles.check_analysis_ensemble_detection.Position(4)];
    handles = sort_on_same_row(handles, 'check_analysis_response_decoding', {'analysis_response_decoding_missing_accessory1', 'check_missing_analysis_response_decoding', 'analysis_response_decoding_missing_accessory2'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, -10]);
    handles.check_analysis_response_decoding_all.Position([1, 2]) = [subpanel_padding_left + subpanel_submenu_submenu_left, handles.check_analysis_response_decoding.Position(2) - handles.check_analysis_response_decoding.Position(4)];
    handles.check_analysis_response_decoding_selective.Position([1, 2]) = [subpanel_padding_left + subpanel_submenu_submenu_left, handles.check_analysis_response_decoding_all.Position(2) - handles.check_analysis_response_decoding_all.Position(4)];
    handles.check_analysis_response_decoding_stable.Position([1, 2]) = [subpanel_padding_left + subpanel_submenu_submenu_left, handles.check_analysis_response_decoding_selective.Position(2) - handles.check_analysis_response_decoding_selective.Position(4)];
    handles.check_analysis_response_decoding_ensembles.Position([1, 2]) = [subpanel_padding_left + subpanel_submenu_submenu_left, handles.check_analysis_response_decoding_stable.Position(2) - handles.check_analysis_response_decoding_stable.Position(4)];
    handles.check_analysis_response_decoding_specific.Position([1, 2]) = [subpanel_padding_left + subpanel_submenu_submenu_left, handles.check_analysis_response_decoding_ensembles.Position(2) - handles.check_analysis_response_decoding_ensembles.Position(4)];
    handles.check_analysis_EPM_track.Position([1, 2])             = [subpanel_padding_left, handles.check_analysis_response_decoding_specific.Position(2) - handles.check_analysis_response_decoding_specific.Position(4)];
    % Adjust panel height
    handles = adjust_panel_height(handles, 'panel_analysis', {'check_calcium_extraction_main', 'check_analysis_spontaneous_activity', 'check_analysis_response_modulation', 'check_analysis_ensemble_detection', 'check_analysis_response_decoding', 'check_analysis_response_decoding_all', 'check_analysis_response_decoding_selective', 'check_analysis_response_decoding_stable', 'check_analysis_response_decoding_ensembles', 'check_analysis_response_decoding_specific', 'check_analysis_EPM_track'}, subpanel_padding_top);
    
    % --- Visualization ---
    % Set position of inner elements
    handles.check_plot_response_detection.Position([1, 2]) = [subpanel_padding_left, handles.panel_visualization.Position(4) - subpanel_padding_top];
    handles = sort_on_same_row(handles, 'check_plot_response_detection', {'visualization_response_detection_selected_accessory1', 'check_plot_response_detection_only_selected', 'visualization_response_detection_selected_accessory2', 'check_plot_response_detection_summary', 'visualization_response_detection_selected_accessory3', 'check_plot_response_detection_overwrite', 'visualization_response_detection_selected_accessory4'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, subpanel_accessory_padding_comma, subpanel_submenu_padding, subpanel_accessory_padding_comma, subpanel_submenu_padding, -10]);
    handles.check_plot_response_decoding.Position([1, 2]) = [subpanel_padding_left, handles.check_plot_response_detection.Position(2) - handles.check_plot_response_detection.Position(4)];
    handles = sort_on_same_row(handles, 'check_plot_response_decoding', {'visualization_response_decoding_selected_accessory1', 'check_plot_response_decoding_only_selected', 'visualization_response_decoding_selected_accessory2', 'check_plot_response_decoding_summary', 'visualization_response_decoding_selected_accessory3', 'check_plot_response_decoding_overwrite', 'visualization_response_decoding_selected_accessory4'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, subpanel_accessory_padding_comma, subpanel_submenu_padding, subpanel_accessory_padding_comma, subpanel_submenu_padding, -10]);
    handles.check_plot_spontaneous_activity.Position([1, 2]) = [subpanel_padding_left, handles.check_plot_response_decoding.Position(2) - handles.check_plot_response_decoding.Position(4)];
    handles = sort_on_same_row(handles, 'check_plot_spontaneous_activity', {'visualization_spontaneous_activity_selected_accessory1', 'check_plot_spontaneous_activity_only_selected', 'visualization_spontaneous_activity_selected_accessory2', 'check_plot_spontaneous_activity_summary', 'visualization_spontaneous_activity_selected_accessory3', 'check_plot_spontaneous_activity_overwrite', 'visualization_spontaneous_activity_selected_accessory4'}, subpanel_submenu_left_column, [0, subpanel_accessory_padding, subpanel_accessory_padding_comma, subpanel_submenu_padding, subpanel_accessory_padding_comma, subpanel_submenu_padding, -10]);
    % Adjust panel height
    handles = adjust_panel_height(handles, 'panel_visualization', {'check_plot_response_detection', 'check_plot_response_decoding', 'check_plot_spontaneous_activity'}, subpanel_padding_top);

    % --- Adjust position of panels ---
    top_gui = handles.txt_main_title.Position(2) - title_bottom_padding;
    left_panel_experiments = handles.panel_experiments.Position(1);
    % Adjust horizontal position
    handles.panel_preprocessing.Position([1, 3]) = [gui_panels_left, left_panel_experiments - gui_panels_left - gui_panels_padding];
    handles.panel_analysis.Position([1, 3])      = [gui_panels_left, left_panel_experiments - gui_panels_left - gui_panels_padding];
    handles.panel_visualization.Position([1, 3]) = [gui_panels_left, left_panel_experiments - gui_panels_left - gui_panels_padding];
    % Adjust vertical position
    handles.panel_preprocessing.Position(2) = top_gui - handles.panel_preprocessing.Position(4);
    handles.panel_analysis.Position(2)      = handles.panel_preprocessing.Position(2) - handles.panel_analysis.Position(4) - gui_panels_padding;
    handles.panel_visualization.Position(2) = handles.panel_analysis.Position(2) - handles.panel_visualization.Position(4) - gui_panels_padding;
    % Adjust position of panels title
    handles = fix_position_panel_title(handles, 'panel_preprocessing', 'panel_title_preprocessing', panel_title_padding_left);
    handles = fix_position_panel_title(handles, 'panel_analysis', 'check_analysis_main', panel_title_padding_left);
    handles = fix_position_panel_title(handles, 'panel_visualization', 'check_visualization_main', panel_title_padding_left);

    % --- Visualization ---
    % Adjust position of panel
    handles.panel_experiments.Position(3) = handles.figure.Position(3) - handles.panel_experiments.Position(1) - gui_panels_right;
    handles.panel_experiments.Position(2) = gui_panels_padding;
    handles.panel_experiments.Position(4) = handles.figure.Position(4) - gui_panels_padding - (handles.figure.Position(4) - top_gui) + panel_experiments_vertical_offset;
    % Adjust position of buttons
    handles.btn_select_none.Position(1) = handles.panel_experiments.Position(3) - handles.btn_select_none.Position(3) - subpanel_padding_right;
    handles.btn_select_none.Position(2) = handles.panel_experiments.Position(4) - subpanel_padding_top - gui_panels_padding;
    handles.btn_select_all.Position(1) = handles.btn_select_none.Position(1) - handles.btn_select_all.Position(3) - subpanel_submenu_padding;
    handles.btn_select_all.Position(2) = handles.btn_select_none.Position(2);
    handles.menu_experiment_name.Position(1) = subpanel_padding_left;
    handles.menu_experiment_name.Position(2) = handles.btn_select_all.Position(2) - 5;
    % Adjust position of table
    handles.table_experiments.Position(1) = subpanel_padding_left;
    handles.table_experiments.Position(2) = subpanel_submenu_bottom;
    handles.table_experiments.Position(3) = handles.panel_experiments.Position(3) - subpanel_padding_left - subpanel_padding_right;
    handles.table_experiments.Position(4) = handles.panel_experiments.Position(4) - subpanel_padding_top - gui_panels_padding - handles.btn_select_all.Position(4);

    % --- Run button ---
    uistack(handles.btn_run, 'top')
    handles.btn_run.Position(4) = 40;
    handles.btn_run.Position([1, 2]) = [5, handles.figure.Position(4) - handles.btn_run.Position(4)];
    handles.btn_run.Position(3) = 150;
    
    % Set unit for all graphical elements
    fn = fieldnames(handles);
    for i_fn = 1:length(fn)
        if isprop(handles.(fn{i_fn}), 'Units')
            set(handles.(fn{i_fn}), 'Units','normalized')
        end
    end

    
function handles = sort_on_same_row(handles, parent_item, items, left_first_item, padding)
    center_height = handles.(parent_item).Position(2) + handles.(parent_item).Position(4) / 2;
    % Copy left position of first item
    last_element_right = left_first_item;
    for ii = 1:length(items)
        item = items{ii};
        % Set units in pixels
        handles.(item).Units = 'pixels';
        % Reposition item to be vertically aligned to center
        handles.(item).Position([1, 2]) = [last_element_right + padding(ii), center_height - handles.(item).Position(4) / 2];
        % Set position for next item
        last_element_right = handles.(item).Position(1) + handles.(item).Position(3);
    end
   
    
function handles = fix_position_panel_title(handles, panel_handle, panel_title_handle, panel_title_padding_left)
    handles.(panel_title_handle).Position([1, 2]) = [handles.(panel_handle).Position(1) + panel_title_padding_left, handles.(panel_handle).Position(2) + handles.(panel_handle).Position(4) - handles.(panel_title_handle).Position(4) / 2];

    
function handles = adjust_panel_height(handles, panel_handle, panel_items, subpanel_padding_top)
    % Calculate height of each item
    height = subpanel_padding_top / 2;
    for item = 1:length(panel_items)
        height = height + handles.(panel_items{item}).Position(4);
    end
    % Calculate offset on y-axis
    y_diff = height - handles.(panel_handle).Position(4);
    % Adjust panel height
    handles.(panel_handle).Position(2) = handles.(panel_handle).Position(2) - y_diff;
    handles.(panel_handle).Position(4) = height;
    % Re-adjust position of children
    children = handles.(panel_handle).Children;
    for i_child = 1:length(children)
        children(i_child).Position(2) = children(i_child).Position(2) + y_diff;
    end

    
%% MLint exceptions
%#ok<*DEFNU,*AGROW,*TLEV,*NASGU,*STRNU>
