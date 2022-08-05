function varargout = Metadata_GUI_animals(varargin)
% Last Modified by GUIDE v2.5 27-Aug-2019 15:18:13
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, 'gui_Singleton',gui_Singleton, 'gui_OpeningFcn',@Metadata_GUI_animals_OpeningFcn, 'gui_OutputFcn',@Metadata_GUI_animals_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else, gui_mainfcn(gui_State, varargin{:}); end
% End initialization code - DO NOT EDIT

%% INIT
function varargout = Metadata_GUI_animals_OutputFcn(hObject, ~, ~), varargout{1}=hObject;

function Metadata_GUI_animals_OpeningFcn(hObject, ~, handles)
    % Set tag of window
    set(handles.figure, 'Tag','Ca_imaging_analysis_Metadata_animals_GUI');

    % Initialize metadata structure
    METADATA = table();
    METADATA(:, 'field_name') = {'experiment_name';'animal_ID';'experimental_group';'brain_area';'calcium_indicator';'imaged_neural_compartment';'genotype';'promoter';'data_type';'has_behavior'};
    METADATA(:, 'default_value') = {'';'';'';'ACC';'GCaMP6f';'soma';'';'hSyn';'2p';'0'};
    METADATA(:, 'user_value') = METADATA.default_value(:);
    METADATA(:, 'can_be_empty') = {false;false;true;false;false;false;true;false;true;true};
    
    % Find names of experiments already in memory
    EXPERIMENT_NAMES = table2cell(unique(SQL_database.read_table_where('experiments', 'experiment_name')));
    set(handles.menu_experiment_name, 'String',[{''}; EXPERIMENT_NAMES])
    
    % Set default values
    set(handles.menu_experiment_name, 'Value',1)  % 1 is empty string
    set(handles.txt_animal_name, 'String',METADATA.default_value{ismember(METADATA.field_name,'animal_ID')})
    set(handles.txt_experiment_group, 'String',METADATA.default_value{ismember(METADATA.field_name,'experimental_group')})
    set(handles.txt_brain_area, 'String',METADATA.default_value{ismember(METADATA.field_name,'brain_area')})
    set(handles.menu_calcium_indicator, 'Value',find(ismember(get(handles.menu_calcium_indicator,'String'), METADATA.default_value{ismember(METADATA.field_name,'calcium_indicator')})))
    set(handles.menu_cell_compartment, 'Value',find(ismember(get(handles.menu_cell_compartment,'String'), METADATA.default_value{ismember(METADATA.field_name,'imaged_neural_compartment')})))
    set(handles.menu_genotype, 'Value',find(ismember(get(handles.menu_genotype,'String'), METADATA.default_value{ismember(METADATA.field_name,'genotype')})))
    set(handles.menu_promoter, 'Value',find(ismember(get(handles.menu_promoter,'String'), METADATA.default_value{ismember(METADATA.field_name,'promoter')})))
    set(handles.menu_data_type, 'Value',find(ismember(get(handles.menu_data_type,'String'), METADATA.default_value{ismember(METADATA.field_name,'data_type')})))
    set(handles.menu_behavior, 'Value',find(ismember(get(handles.menu_behavior,'String'), METADATA.default_value{ismember(METADATA.field_name,'has_behavior')})))
    
    % Update GUI data
    setappdata(handles.figure,'METADATA',METADATA);
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
    
    % Get LOGGER
    global LOGGER
    
    % Read INFO and METADATA
    METADATA = getappdata(handles.figure, 'METADATA');

    % Check that no field has been left empty
    can_write = true;
    for irow = 1:height(METADATA)
        % If field cannot be empty, check whether it isn't
        if ~METADATA.can_be_empty(irow)
            value = METADATA.user_value{irow};
            if isempty(value)
                can_write = false;
                break  % 1 empty value is enough
            end
        end
    end
    
    caught_error = false;
    if can_write
        % Update database
        columns_to_update = METADATA.field_name';
        values_to_insert = METADATA.user_value';
        try
            % Check whether this animal already exists
            experiment_name = values_to_insert{ismember(columns_to_update, 'experiment_name')};
            animal_ID = values_to_insert{ismember(columns_to_update, 'animal_ID')};
            result = SQL_database.read_table_where('experiments', {}, {experiment_name, animal_ID},{'experiment_name', 'animal_ID'});
            if isempty(result)  % Insert new values
                SQL_database.update('experiments', columns_to_update, values_to_insert)
            else  % Update existing values
                SQL_database.update('experiments', columns_to_update, values_to_insert, {experiment_name, animal_ID},{'experiment_name', 'animal_ID'})
                
            end
                
        catch MException  % Grab exception
            caught_error = true;
            if strcmp(MException.identifier, 'execute_SQL_query:SQL_error')
                LOGGER.warn(MException.message)
            else
                MException.rethrow()
            end
        end
        % Log outcome
        animal_ID = METADATA.user_value{ismember(METADATA.field_name,'animal_ID')};
        LOGGER.info(['Added ', animal_ID, ' to database'])
    else
        LOGGER.warn('One or more fields have been left empty. Will not be able to continue until these fields are not empty')
    end
    
    % Turn GUI on
    toggle_GUI(handles, 'on')
    if ~caught_error
        delete(handles.figure)
    end
    
%% TEXT EDITORS
function menu_experiment_name_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_experiment_name, 'String');
    experiment_name = menu_values{get(handles.menu_experiment_name, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'experiment_name')} = experiment_name;
    setappdata(handles.figure, 'METADATA',METADATA);    

function btn_add_new_experiment_Callback(~, ~, handles)
    % Ask user for input
    answer = inputdlg('Insert new name', 'Add new experiment name', 1, {''}, struct('FontSize',12));
    answer = answer{1};
    % Remove all white spaces at the borders and replace the inner ones with underscores
    answer = strrep(strtrim(answer), ' ', '_');
    if ~isempty(answer)
        % Add new name
        menu_values = get(handles.menu_experiment_name, 'String');
        set(handles.menu_experiment_name, 'String',[{answer}; menu_values])
        % Select it
        set(handles.menu_experiment_name, 'Value',1)
        % Store it into memory
        METADATA = getappdata(handles.figure, 'METADATA');
        METADATA.user_value{ismember(METADATA.field_name,'experiment_name')} = answer;
        setappdata(handles.figure, 'METADATA',METADATA);
    end
    
function txt_animal_name_Callback(~, ~, handles)
    % Save state of this text editor
    METADATA = getappdata(handles.figure, 'METADATA');
    animal_ID = get(handles.txt_animal_name, 'String');
    % Remove all white spaces at the borders and replace the inner ones with underscores
    animal_ID = strrep(strtrim(animal_ID), ' ', '_');
    % Update data in memory and the GUI
    set(handles.txt_animal_name, 'String', animal_ID)
    METADATA.user_value{ismember(METADATA.field_name,'animal_ID')} = animal_ID;
    setappdata(handles.figure, 'METADATA',METADATA);

function txt_experiment_group_Callback(~, ~, handles)
    % Save state of this text editor
    METADATA = getappdata(handles.figure, 'METADATA');
    experiment_group = get(handles.txt_experiment_group, 'String');
    % Remove all white spaces at the borders and replace the inner ones with underscores
    experiment_group = strrep(strtrim(experiment_group), ' ', '_');
    % Update data in memory and the GUI
    set(handles.txt_experiment_group, 'String', experiment_group)
    METADATA.user_value{ismember(METADATA.field_name,'experimental_group')} = experiment_group;
    setappdata(handles.figure, 'METADATA',METADATA);

function txt_brain_area_Callback(~, ~, handles)
    % Save state of this text editor
    METADATA = getappdata(handles.figure, 'METADATA');
    brain_area = get(handles.txt_brain_area, 'String');
    % Remove all white spaces at the borders and replace the inner ones with underscores
    brain_area = strrep(strtrim(brain_area), ' ', '_');
    % Update data in memory and the GUI
    set(handles.txt_brain_area, 'String', brain_area)
    METADATA.user_value{ismember(METADATA.field_name,'brain_area')} = brain_area;
    setappdata(handles.figure, 'METADATA',METADATA);

%% MENUS
function menu_cell_compartment_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_cell_compartment, 'String');
    cell_compartment = menu_values{get(handles.menu_cell_compartment, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'imaged_neural_compartment')} = cell_compartment;
    setappdata(handles.figure, 'METADATA',METADATA);    
    
function menu_calcium_indicator_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_calcium_indicator, 'String');
    calcium_indicator = menu_values{get(handles.menu_cell_compartment, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'calcium_indicator')} = calcium_indicator;
    setappdata(handles.figure, 'METADATA',METADATA);    
    
function menu_genotype_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_genotype, 'String');
    genotype = menu_values{get(handles.menu_genotype, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'genotype')} = genotype;
    setappdata(handles.figure, 'METADATA',METADATA);    
    
function menu_promoter_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_genotype, 'String');
    promoter = menu_values{get(handles.menu_promoter, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'promoter')} = promoter;
    setappdata(handles.figure, 'METADATA',METADATA);    
    
function menu_data_type_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_data_type, 'String');
    data_type = menu_values{get(handles.menu_data_type, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'data_type')} = data_type;
    setappdata(handles.figure, 'METADATA',METADATA);    

function menu_behavior_Callback(~, ~, handles)
    METADATA = getappdata(handles.figure, 'METADATA');
    menu_values = get(handles.menu_behavior, 'String');
    has_behavior = menu_values{get(handles.menu_behavior, 'Value')};
    METADATA.user_value{ismember(METADATA.field_name,'has_behavior')} = has_behavior;
    setappdata(handles.figure, 'METADATA',METADATA);    

    
%% CLOSE FIGURE
function quit_GUI(handle_figure)
    % Close request function to display a question dialog box
    selection = MessageBox('Are you sure?', 'Close Request Function', 'Yes','No', 'No');
    switch selection
        case 'Yes'
            delete(handle_figure)
        case 'No'
            return
    end 

%% MLint exceptions
%#ok<*DEFNU,*TRYNC>
