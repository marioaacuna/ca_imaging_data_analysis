function varargout = Metadata_GUI_launcher(varargin)
% Last Modified by GUIDE v2.5 15-Jun-2018 16:17:21
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, 'gui_Singleton',gui_Singleton, 'gui_OpeningFcn',@Metadata_GUI_launcher_OpeningFcn, 'gui_OutputFcn',@Metadata_GUI_launcher_OutputFcn, 'gui_LayoutFcn',[], 'gui_Callback',[]);
if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end
if nargout, [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:}); else, gui_mainfcn(gui_State, varargin{:}); end
% End initialization code - DO NOT EDIT

%% INIT
function varargout = Metadata_GUI_launcher_OutputFcn(hObject, ~, ~), varargout{1}=hObject;

function Metadata_GUI_launcher_OpeningFcn(hObject, ~, handles, varargin)
    % Get logger and general_configs
    global LOGGER
    % Log beginning of workflow
    LOGGER.info('Metadata workflow', 'decorate',true)

    % Unpack inputs
    INFO = varargin{1};
        
    % Set tag of window
    set(handles.figure, 'Tag','Ca_imaging_analysis_Metadata_launcher_GUI');
    
    % Update GUI data
    setappdata(handles.figure,'INFO',INFO);
    handles.output = true;
    % Store handles in object
    guidata(hObject, handles);

    % Bring figure to front
    figure(handles.figure)


%% RUN
function btn_register_animals_Callback(~, ~, handles)
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')

    % Launch another GUI
    waitfor(Metadata_GUI_animals());  % Stop execution until workflow is done

    % Re-enable GUI controls
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;
    figure(handles.figure)
    
function btn_register_sessions_Callback(~, ~, handles)
    % Update list of objects that can be turned on and off, and then turn everything off
    toggle_GUI(handles, 'list')
    toggle_GUI(handles, 'off')
    
    % Launch another GUI
    INFO = getappdata(handles.figure,'INFO');
    waitfor(Metadata_GUI_sessions(INFO));  % Stop execution until workflow is done

    % Re-enable GUI controls
    toggle_GUI(handles, 'on')
    handles.figure.CloseRequestFcn = @quit_GUI;
    figure(handles.figure)

    
%% CLOSE FIGURE
function quit_GUI(handle_figure, ~)
    delete(handle_figure)

    
%% MLint exceptions
%#ok<*DEFNU,*TRYNC>
