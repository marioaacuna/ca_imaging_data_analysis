function varargout = Ca_event_detection(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, 'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Ca_event_detection_OpeningFcn, ...
                   'gui_OutputFcn',  @Ca_event_detection_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}), gui_State.gui_Callback = str2func(varargin{1}); end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
function varargout = Ca_event_detection_OutputFcn(~, ~, handles) 
varargout{1} = handles.output;

% End initialization code - DO NOT EDIT

%% INIT
function Ca_event_detection_OpeningFcn(hObject, ~, handles, varargin)
    if isempty(varargin)
        disp('Do not launch GUI directly')
        return
    end
    animal_ID = varargin{1};
    
    % Make figure immediately visible
    figure(handles.figure)
    
    % Add animal_Id to figure title
    set(handles.figure, 'Name', ['Ca event detection - ', animal_ID])
    
    % Fix axes appearance
    set(handles.axes_traces, 'Box','off', 'TickLength',[0, 0])
    set(handles.axes_events, 'Box','off', 'TickLength',[0, 0])
    handles.btn_normalize.BackgroundColor = [.7,.7,.7];
    handles.btn_normalize.Value = 0;
    
    % Make INFO structure
    INFO = struct();
    INFO.animal_ID = animal_ID;
    
    % Load data
    disp('Loading data')
    dff_filename = get_filename_of('dFF', animal_ID);
    INFO.data = load_variable(dff_filename, 'dFF');
    INFO.n_cells = size(INFO.data, 1);
    
    % Load previous work on this animal
    filename = get_filename_of('Ca_events', INFO.animal_ID);
    if exist(filename, 'file')
        disp('Loading previous state')
        load(filename);
        INFO.thresholds = thresholds;
        INFO.keep_cell = keep_cell;
        try
            INFO.normalized_trace = normalized_trace;
        catch
            INFO.normalized_trace = zeros(INFO.n_cells, 1);
        end
        INFO.Ca_events = Ca_events;
        INFO.current_cell = find(isnan(INFO.keep_cell), 1, 'first');
        if isempty(INFO.current_cell), INFO.current_cell = 1; end
    else
        INFO.current_cell = 1;
        INFO.thresholds = NaN(INFO.n_cells, 1);
        INFO.keep_cell = NaN(INFO.n_cells, 1);
        INFO.normalized_trace = zeros(INFO.n_cells, 1);
        INFO.Ca_events = cell(INFO.n_cells, 1);
    end
    
    % Set current cell to the first one
    if INFO.current_cell == 1
        set(handles.btn_previous, 'Enable','off')
    else
        set(handles.btn_previous, 'Enable','on')
    end
    if INFO.current_cell == INFO.n_cells
        set(handles.btn_next, 'Enable','off')
    else
        set(handles.btn_next, 'Enable','on')
    end
    setappdata(handles.figure, 'INFO',INFO);
    % Update handles structure
    guidata(hObject, handles);
    
    % Load data of first cell
    handles = update_cell_trace(handles);
    
    % Detect peaks
    handles = detect_peaks(handles);

    % Position GUI in top left corner of screen
    movegui(handles.figure, [100, -1])

    % Link callback function when window has to close
    handles.figure.CloseRequestFcn = @stop_ML_execution;
    % Make sure window is not docked
    handles.figure.WindowStyle = 'normal';

    % Choose default command line output for GUI_mark_stimuli
    handles.output = hObject;

    % Update handles structure
    guidata(hObject, handles);
    
    disp('GUI is ready')

    
%% Callbacks
function menu_save_Callback(~, ~, handles)    
    % Load info
    INFO = getappdata(handles.figure, 'INFO');
    
    % Convert Ca_events to table
    Ca_events = INFO.Ca_events;
    thresholds = INFO.thresholds;
    keep_cell = INFO.keep_cell;
    normalized_trace = INFO.normalized_trace;
    
    % Save data to disk
    filename = get_filename_of('Ca_events', INFO.animal_ID);
    save(filename, 'Ca_events', 'thresholds', 'keep_cell','normalized_trace',  '-v7.3')
    disp(['Results saved in ''', filename, ''''])
    

function btn_accept_Callback(~, ~, handles)
    INFO = getappdata(handles.figure, 'INFO');
    INFO.keep_cell(INFO.current_cell) = logical(get(handles.btn_accept, 'Value'));
    setappdata(handles.figure, 'INFO',INFO)

    
function change_cell(hObject, ~, handles, move_to)
    % Load INFO
    INFO = getappdata(handles.figure, 'INFO');
    % Check if we can move to the another cell
    if move_to == -1 && INFO.current_cell == 1, return, end
    if move_to == 1 && INFO.current_cell == INFO.n_cells, return, end

    % Update current cell
    INFO.current_cell = INFO.current_cell + move_to;
    setappdata(handles.figure, 'INFO',INFO)
    
    % Load data of first cell
    handles = update_cell_trace(handles);
    
    % Detect peaks
    handles = detect_peaks(handles);
    INFO = getappdata(handles.figure, 'INFO');
    handles = plot_peaks(INFO, handles);

    if INFO.current_cell == 1
        set(handles.btn_previous, 'Enable','off')
    else
        set(handles.btn_previous, 'Enable','on')
    end
    if INFO.current_cell == INFO.n_cells
        set(handles.btn_next, 'Enable','off')
    else
        set(handles.btn_next, 'Enable','on')
    end

    % Store info
    guidata(hObject, handles);

    
function slider_threshold_Callback(hObject, ~, handles)
    % Update value of threshold in memory
    INFO = getappdata(handles.figure, 'INFO');
    INFO.thresholds(INFO.current_cell) = get(handles.slider_threshold, 'Value');
    
    if handles.btn_normalize.Value == 1
        data_to_analyze = 'smoothed_trace_norm';
    else
        data_to_analyze = 'smoothed_trace';
    end
    set(handles.slider_threshold, 'Min',min(INFO.(data_to_analyze)), 'Max',max(INFO.(data_to_analyze)), 'Value',INFO.thresholds(INFO.current_cell));
    setappdata(handles.figure, 'INFO',INFO)

    % Detect peaks
    handles = detect_peaks(handles);
    INFO = getappdata(handles.figure, 'INFO');
    handles = plot_peaks(INFO, handles);
    % Show value
    set(handles.txt_threshold_value, 'String',num2str(get(handles.slider_threshold, 'Value'), '%.3f'))
    guidata(hObject, handles);
    
    
function btn_normalize_Callback(hObject, ~, handles)
    INFO = getappdata(handles.figure, 'INFO');

    if handles.btn_normalize.Value == 1
        data_to_analyze = 'smoothed_trace_norm';
        handles.btn_normalize.BackgroundColor = 'g';
        INFO.normalized_trace(INFO.current_cell) = 1;
    else
        data_to_analyze = 'smoothed_trace';
        handles.btn_normalize.BackgroundColor = [.7,.7,.7];
        INFO.normalized_trace(INFO.current_cell) = 0;
    end
    slider_range = [min(INFO.(data_to_analyze)), max(INFO.(data_to_analyze))];
    current_threshold_value = handles.slider_threshold.Value;
    if current_threshold_value < slider_range(1) || current_threshold_value > slider_range(2)
        set(handles.slider_threshold, 'Min',slider_range(1), 'Max',slider_range(2), 'Value',mad(INFO.(data_to_analyze)) * 3);
        INFO.thresholds(INFO.current_cell) = get(handles.slider_threshold, 'Value');
        set(handles.txt_threshold_value, 'String',num2str(get(handles.slider_threshold, 'Value'), '%.3f'))
        handles = detect_peaks(handles);
    end
    
    setappdata(handles.figure, 'INFO',INFO)
    handles = plot_cell_trace(INFO, handles);
    handles = plot_peaks(INFO, handles);
    % Store info
    guidata(hObject, handles);
 

function ch_show_raw_Callback(~, ~, handles)
    INFO = getappdata(handles.figure, 'INFO');
    y_limits = [];
    if handles.ch_show_raw.Value == 1
        visiblity = 'on';
        y_limits = [y_limits; [min(INFO.trace), max(INFO.trace)]];
    else
        visiblity = 'off';
    end
    set(handles.h_raw_trace, 'Visible',visiblity)

    if handles.ch_show_smoothed.Value == 1
        y_limits = [y_limits; [min(INFO.smoothed_trace), max(INFO.smoothed_trace)]];
    end
    if ~isempty(y_limits)
        set(handles.axes_traces, 'YLim',[min(y_limits(:)), max(y_limits(:))])
    end
    

function ch_show_smoothed_Callback(~, ~, handles)
    INFO = getappdata(handles.figure, 'INFO');
    y_limits = [];
    if handles.ch_show_smoothed.Value == 1
        visiblity = 'on';
        y_limits = [y_limits; [min(INFO.smoothed_trace), max(INFO.smoothed_trace)]];
    else
        visiblity = 'off';
    end
    set(handles.h_smoothed_trace, 'Visible',visiblity)

    if handles.ch_show_raw.Value == 1
        y_limits = [y_limits; [min(INFO.trace), max(INFO.trace)]];
    end
    if ~isempty(y_limits)
        set(handles.axes_traces, 'YLim',[min(y_limits(:)), max(y_limits(:))])
    end
    
    
%% Utilities
function handles = plot_cell_trace(INFO, handles)
    % Clear current axes
    set(handles.figure, 'CurrentAxes', handles.axes_traces)
    cla, hold on
    zoom on
    
    % Gather y-limits
    y_limits = [];
    
    % Plot dFF
    if handles.btn_normalize.Value == 1
        data_to_plot = INFO.trace_norm;
    else
        data_to_plot = INFO.trace;
    end
    if handles.ch_show_raw.Value == 1
        visiblity = 'on';
        y_limits = [y_limits; [min(data_to_plot), max(data_to_plot)]];
    else
        visiblity = 'off';
    end
    handles.h_raw_trace = plot(handles.axes_traces, 1:length(INFO.trace), data_to_plot, 'Color',[.7,.7,.7], 'Visible',visiblity, 'XLimInclude','off', 'YLimInclude','off');
    
    % Plot smoothed trace
    if handles.btn_normalize.Value == 1
        data_to_plot = INFO.smoothed_trace_norm;
    else
        data_to_plot = INFO.smoothed_trace;
    end
    if handles.ch_show_smoothed.Value == 1
        visiblity = 'on';
        y_limits = [y_limits; [min(data_to_plot), max(data_to_plot)]];
    else
        visiblity = 'off';
    end
    handles.h_smoothed_trace = plot(handles.axes_traces, 1:length(INFO.trace), data_to_plot, 'Color','k', 'Visible',visiblity, 'XLimInclude','off', 'YLimInclude','off');

    % Set x-lim
    set(handles.axes_traces, 'XLim',[1, length(INFO.trace)], 'YLim',[min(y_limits(:)), max(y_limits(:))])
    
    
function handles = update_cell_trace(handles)
    global GC
    INFO = getappdata(handles.figure, 'INFO');
    
    % Update value of btn_normalize
    handles.btn_normalize.Value = INFO.normalized_trace(INFO.current_cell);
    if handles.btn_normalize.Value == 1
        handles.btn_normalize.BackgroundColor = 'g';
    else
        handles.btn_normalize.BackgroundColor = [.7,.7,.7];
    end

    % Update info on screen
    set(handles.txt_info, 'String',['cell ', num2str(INFO.current_cell), '/', num2str(INFO.n_cells)])
    
    % Concatenate trials
    INFO.trace = INFO.data(INFO.current_cell, :);
    INFO.trace = horzcat(INFO.trace{:}).';
    
    % Normalize trace: subtract median and divide by the MAD (similar to z-score, but not affected by outliers, i.e., Ca events)
    data_type = SQL_database.read_table_where('experiments', 'data_type', INFO.animal_ID,'animal_ID', 'return_as_table',false);
    switch data_type
        case '2p'
            METADATA = SQL_database.read_table_where('trials', {'frames_idx', 'date'}, INFO.animal_ID, 'animal_ID');
            frames_idx = cellfun(@(x) strsplit(x, ','), METADATA.frames_idx, 'UniformOutput',false);
            frames_idx = vertcat(frames_idx{:});
            frames_idx = str2double(frames_idx);

        case 'epi'
            METADATA = SQL_database.read_epi_trials(INFO.animal_ID);
            trial_duration = cellfun(@(x) length(x), INFO.data(1, :))';
            frames_idx = cumsum(trial_duration);
            frames_idx = [[1; frames_idx(1:end-1) + 1], frames_idx];
            frame_rate = GC.epifluorescence_downsample_to_frame_rate;
    end
    [sessions, ~, session_idx] = unique(METADATA.date);
    session_of_frame_idx = NaN(length(INFO.trace), 1);
    for i_trial = 1:size(frames_idx)
        session_of_frame_idx(frames_idx(i_trial, 1):frames_idx(i_trial, 2)) = session_idx(i_trial);
    end
    % Remove frames that have been discarded
    session_of_frame_idx(isnan(session_of_frame_idx)) = [];
    % Normalize trace
    INFO.trace_norm = INFO.trace;
    for i_sess = 1:length(sessions)
        F = INFO.trace(session_of_frame_idx == i_sess);
        F = F - median(F);
        F = F ./ mad(F);
        INFO.trace_norm(session_of_frame_idx == i_sess) = F;
    end

    % Filter trace
    INFO = filter_trace(INFO.trace, INFO, handles);
    INFO = filter_trace(INFO.trace_norm, INFO, handles, 'smoothed_trace_norm', false);
    
    % Plot trace
    handles = plot_cell_trace(INFO, handles);
    
    if handles.btn_normalize.Value == 1
        data_to_analyze = 'smoothed_trace_norm';
        INFO.normalized_trace(INFO.current_cell) = 1;
    else
        data_to_analyze = 'smoothed_trace';
        INFO.normalized_trace(INFO.current_cell) = 0;
    end
    % Set range of thresholds in slider
    if isnan(INFO.thresholds(INFO.current_cell))
        INFO.thresholds(INFO.current_cell) = (mad(INFO.(data_to_analyze)) * 4) + median(INFO.(data_to_analyze));
    end
    min_slider = min(INFO.(data_to_analyze));
    max_slider = max(INFO.(data_to_analyze));
    slider_value = INFO.thresholds(INFO.current_cell);
    if slider_value < min_slider
        slider_value = min_slider;
    end
    if slider_value > max_slider
        slider_value = max_slider;
    end
    set(handles.slider_threshold, 'Min', min_slider, 'Max', max_slider, 'Value',slider_value);
    % Show value
    set(handles.txt_threshold_value, 'String',num2str(get(handles.slider_threshold, 'Value'), '%.3f'))
    
    % Update checkbox to accept
    if isnan(INFO.keep_cell(INFO.current_cell))
        INFO.keep_cell(INFO.current_cell) = 1;
    end
    set(handles.btn_accept, 'Value', INFO.keep_cell(INFO.current_cell))
    
    % Store info
    setappdata(handles.figure, 'INFO',INFO)

    
function INFO = filter_trace(dFF_cell, INFO, handles, output_name, compute_derivative)
    if ~exist('output_name','var'), output_name = 'smoothed_trace'; end
    if ~exist('compute_derivative','var'), compute_derivative = true; end

    % Filter trace
    frame_rate = 8.49;
    smoothing_width = floor(1000 / 1000 * frame_rate / 3) * 2 + 1;
    % Smooth fluorescence trace and compute its first-order derivative
    [~, g] = sgolay(1, smoothing_width);
    % Smooth the data
    diff_filter = 1;
    INFO.(output_name) = conv(dFF_cell, factorial(diff_filter-1)/(-(1/frame_rate))^(diff_filter-1) * g(:,diff_filter), 'same');
            
    if compute_derivative
        % Apply 1st order derivative to smoothed data
        diff_filter = 2;
        INFO.smoothed_derivative = conv(INFO.(output_name), factorial(diff_filter-1)/(-(1/frame_rate))^(diff_filter-1) * g(:,diff_filter), 'same');

        % Find 0-crossings in 1st derivative (i.e., sign of product of consecutive
        % samples is negative)
        zx = find(sign(INFO.smoothed_derivative(1:end-1).*INFO.smoothed_derivative(2:end)) < 0);
        % Remove spurious points
        zx(zx<1) = [];
        zx(zx>=length(INFO.smoothed_derivative)) = [];
        % Get the sign of points around 0-crossings
        yx = [INFO.smoothed_derivative(zx) INFO.smoothed_derivative(zx+1)];
        % Keep transitions from rising to falling
        INFO.positive_crossings = zx(yx(:,1)>=0 & yx(:,2)<0);
        % Keep transitions from falling to rising
        INFO.negative_crossings = zx(yx(:,1)<0 & yx(:,2)>=0);
    end
    
    
function handles = plot_peaks(INFO, handles)
    % Draw peak location
    try delete(handles.peaks), end
    if handles.btn_normalize.Value == 1
        data_to_plot = INFO.smoothed_trace_norm;
    else
        data_to_plot = INFO.smoothed_trace;
    end
    handles.peaks = plot(handles.axes_traces, INFO.peaks, data_to_plot(INFO.peaks), 'or', 'markerfacecolor','r', 'markersize',8, 'XLimInclude','off', 'YLimInclude','off');

    
function handles = detect_peaks(handles)
    INFO = getappdata(handles.figure, 'INFO');
    
    % Detect and store peaks
    if handles.btn_normalize.Value == 1
        data_to_analyze = 'smoothed_trace_norm';
    else
        data_to_analyze = 'smoothed_trace';
    end
    [~, INFO.peaks] = findpeaks(INFO.(data_to_analyze), 'MinPeakHeight',get(handles.slider_threshold, 'Value') , 'MinPeakWidth',5);
    % Draw peak location
    handles = plot_peaks(INFO, handles);
    
    % Extract peaks
    INFO = extract_events(INFO);
    if ~isempty(INFO.Ca_events{INFO.current_cell})
        % Draw peaks
        set(handles.figure, 'CurrentAxes', handles.axes_events)
        cla, hold on
        for i_peak = 1:size(INFO.Ca_events{INFO.current_cell}, 1)
            event = INFO.(data_to_analyze)(INFO.Ca_events{INFO.current_cell}(i_peak, 1):INFO.Ca_events{INFO.current_cell}(i_peak, 3));
            plot(event, 'k')
        end
        % Plot histogram of amplitudes
        bins = unique(linspace(min(INFO.Ca_events{INFO.current_cell}(:, 4)), max(INFO.Ca_events{INFO.current_cell}(:, 4)), 100));
        if length(bins) == 1
            bins = [bins - .5, bins, bins + .5];
        end
        h = hist(INFO.Ca_events{INFO.current_cell}(:, 4), bins);
        set(handles.figure, 'CurrentAxes', handles.axes_hist)
        cla
        bar(bins, h, 'BarWidth',1, 'FaceColor','k', 'EdgeColor','none')
        set(handles.axes_hist, 'Box','off', 'TickLength',[0, 0])
    end
    
    % Store info
    setappdata(handles.figure, 'INFO',INFO)
    

function INFO = extract_events(INFO)
    n_peaks = length(INFO.peaks);
    n_samples = length(INFO.smoothed_trace);
    
    Ca_events = NaN(n_peaks, 3);
    for i_peak = 1:n_peaks
        % Get peak position
        spk = INFO.peaks(i_peak);
        % Re-align peak on fluorescence first derivative
        [~, peak_idx] = min(abs(spk - INFO.positive_crossings));
        spk_peak_idx = INFO.positive_crossings(peak_idx);
        
        % Find two negative 0-crossings, one before and one after the peak
        event_onset = max([INFO.negative_crossings(find(INFO.negative_crossings < spk_peak_idx, 1, 'last')), 1]);
        event_end   = min([INFO.negative_crossings(find(INFO.negative_crossings > spk_peak_idx, 1, 'first')) - 1, n_samples]);
        
        % Store results
        Ca_events(i_peak, :) = [event_onset, spk_peak_idx, event_end];
    end
    % Get amplitudes
    Ca_events(:, 4) = INFO.smoothed_trace(Ca_events(:, 2));
    
    % Store results
    INFO.Ca_events{INFO.current_cell} = Ca_events;

    
function stop_ML_execution(~, ~)
    % Close request function to display a question dialog box
    selection = MessageBox('Stop analysis?', 'Close Request Function', 'Yes','No', 'No');
    switch selection
        case 'Yes'
            % Close GUI
            close('force')
        case 'No'
            return
    end


%#ok<*DEFNU,*NASGU,*LOAD>
