function PROJECTIONS = make_trial_projections_Copy(animal_ID, filename_projections, varargin )
if isempty(varargin)
    animal_ID = 'FK_20';
    filename_projections = 'D:\_MATLAB_CaImaging\FK_20_projections.mat';
    reload_data = false;
    projections =  {'mean', 'median', 'max', 'std', 'corr'};
    print_figures = true;
    only_structural_channel =  false;
    data_type = '2p';
else
    
% Parse inputs
    p = inputParser();
    addParameter(p, 'reload_data', false)  % Whether to overwrite existing file
    addParameter(p, 'projections', {'mean', 'median', 'max', 'std', 'corr'})
    addParameter(p, 'print_figures', true)
    addParameter(p, 'only_structural_channel', false)
    addParameter(p, 'data_type', '2p')
    parse(p, varargin{:})
    % Fetch inputs
    reload_data = p.Results.reload_data;
    projections = p.Results.projections;
    print_figures = p.Results.print_figures;
    only_structural_channel = p.Results.only_structural_channel;
    data_type = p.Results.data_type;
end 
% Enable toolbox for printing
if print_figures
    toggle_toolbox('_plotting', 'on')
end

if ~iscell(projections), projections = {projections}; end  % In case user passed a string
if isempty(projections)
    error('make_trial_projections:emptyProjectionType', 'Specifiy at least one projection type')
end
n_projection_type = length(projections);
for p = 1:n_projection_type
    switch projections{p}
        case 'mean', projection_names{p}  = 'mean';
        case 'median', projection_names{p}  = 'median';
        case 'max', projection_names{p}  = 'max';
        case 'std', projection_names{p}  = 'standard_deviation';
        case 'corr', projection_names{p}  = 'correlation';
    end
end

% Get global variables
global GC LOGGER
has_logger = ~isempty(LOGGER);

msg = 'Computing image stacks projections';
if has_logger
    LOGGER.info(msg)
else
    disp(msg)
end

% Read metadata
switch data_type
    case '2p'
        METADATA = SQL_database.read_table_where('trials', {'+','animal_ID','date','stimulus','experimental_condition','has_structural_channel'}, animal_ID, 'animal_ID', 'return_all_sessions',true);
        group_sessions_by = {'date'};
        
    case 'epi'
        METADATA = SQL_database.read_table_where('sessions', {'+','animal_ID'}, animal_ID, 'animal_ID', 'return_all_sessions',true);
        METADATA{:, 'has_structural_channel'} = false;
        group_sessions_by = {'date', 'experiment'};
        % Load preprocessing parameters
        p_filename = get_filename_of('miniscope_movie_parameters', animal_ID);
        miniscope_parameters = load_variable(p_filename, 'p');
        % Get filename of concatenated file
        local_motion_corrected_filename = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_joint.h5']);
end
sessions = uniqueCellRows(METADATA{:, group_sessions_by}, 'rows');
n_sessions = size(sessions, 1);

% Load previous file
if exist(filename_projections, 'file')
    PROJECTIONS = load_variable(filename_projections, 'PROJECTIONS');
else
    PROJECTIONS = struct();
end
% Initialize status of data file
file_in = [];

for i_sess = 1%
    % Extract metadata for this session
    metadata = METADATA(ismemberCellRows(METADATA{:, group_sessions_by}, sessions(i_sess, :)), :);
    switch data_type
        case '2p'
            n_trials = height(metadata);
            
        case 'epi'
    end

    % Reload data into this file, if requested
    session_fieldname = ['session_', strjoin(sessions(i_sess, :), '_')];
    if ~isfield(PROJECTIONS, session_fieldname) || reload_data
        reload_data_this_session = true;
    else
        reload_data_this_session = false;
    end
    % Open file for reading, if not open already
    if reload_data_this_session
        if isempty(file_in)
            switch data_type
                case '2p'
%                     file_in = open_motion_corrected_file(animal_ID, only_structural_channel);
                    file_in = matfile(['L:\motion_correction_data_2p\', animal_ID,'.mat'], 'Writable',false);
                    data_size = size(file_in, 'Y');

                case 'epi'
                    file_in = struct('Y', loadMovie(local_motion_corrected_filename, miniscope_parameters.user.jointExtraction.preprocessed_frame_edges(i_sess, 1), miniscope_parameters.user.jointExtraction.preprocessed_frame_edges(i_sess, 2)));
                    data_size = size(file_in.Y);
                    n_frames = data_size(3);
                    n_trials = ceil(n_frames / GC.epifluorescence_segmentation_n_frames_per_chunk);
                    % Set beginning and end of each trial
                    chunk_start = 1:GC.epifluorescence_segmentation_n_frames_per_chunk:n_frames;
                    chunk_end = [chunk_start(2:end) - 1, n_frames];
                    frames_idx = [chunk_start(:), chunk_end(:)];
            end
            frame_height = data_size(1);
            frame_width = data_size(2);
        end

        % Delete previous runs on this session
        PROJECTIONS.(session_fieldname) = struct();
        % Allocate array
        if only_structural_channel
            empty_array = NaN(frame_height, frame_width);
            log_prefix = 'Structural channel';
        else
            empty_array = NaN(frame_height, frame_width, n_trials);
            log_prefix = 'Functional channel';
        end
        % Initialize output variable
        for p = 1:n_projection_type
            switch projections{p}
                case 'mean'
                    PROJECTIONS.(session_fieldname).mean = empty_array;
                    projection_names{p}  = 'mean';
                case 'median'
                    PROJECTIONS.(session_fieldname).median = empty_array;
                    projection_names{p}  = 'median';
                case 'max'
                    PROJECTIONS.(session_fieldname).max = empty_array;
                    projection_names{p}  = 'max';
                case 'std'
                    PROJECTIONS.(session_fieldname).standard_deviation = empty_array;
                    projection_names{p}  = 'standard_deviation';
                case 'corr'
                    PROJECTIONS.(session_fieldname).correlation = empty_array;
                    projection_names{p}  = 'correlation';
            end
        end

        for itrial = 1:n_trials
            msg = [sessions{i_sess}, ': ', log_prefix, ' sweep ', num2str(itrial), '/', num2str(n_trials)];
            if has_logger
                LOGGER.trace(msg, 'write_to_file',false, 'overwrite_last_message',itrial>1)
            else
                disp(msg)
            end
            % Read frames from data file
            switch data_type
                case '2p'
                    idx = metadata.frames_idx{itrial};
                    idx = regexp(idx, ',', 'split');
                    idx = cellfun(@str2double, idx);
                    % Get the frame rate
                    frame_rate = metadata.frame_rate(itrial);
                
                case 'epi'
                    idx = frames_idx(itrial, :);
                    % Get the frame rate
                    frame_rate = GC.epifluorescence_downsample_to_frame_rate;
            end
            frames = file_in.Y(1:frame_height, 1:frame_width, idx(1):idx(2));

            % Compute projections
            for p = 1:n_projection_type
                switch projections{p}
                    case 'mean'
                        PROJECTIONS.(session_fieldname).mean(:, :, itrial) = nanmean(frames, 3);

                    case 'median'
                        PROJECTIONS.(session_fieldname).median(:, :, itrial) = nanmedian(frames, 3);

                    case 'max'
                        PROJECTIONS.(session_fieldname).max(:, :, itrial) = nanmax(frames, [], 3);

                    case 'std'
                        PROJECTIONS.(session_fieldname).standard_deviation(:, :, itrial) = nanstd(frames, [], 3);

                    case 'corr'
                        % Compute crosscorrelation after smoothing video in time (no NaNs allowed)
                        time_window = ceil(GC.crosscorrelation_time_smoothing_window * frame_rate);
                        frames_smoothed = temporalSmooth(frames, time_window);
                        frames_smoothed(isnan(frames_smoothed)) = 0;
                        PROJECTIONS.(session_fieldname).correlation(:, :, itrial) = CrossCorrImage(frames_smoothed);
                end
            end
        end
    end
    
    if print_figures
        msg = 'Printing projections to PDF file';
        if has_logger
            LOGGER.info(msg)
        else
            disp(msg)
        end
        % Set and make folder
        figure_folder = os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name, 'fields_of_view', animal_ID);
        if ~exist(figure_folder, 'dir')
            mkdir(figure_folder)
        end
        if only_structural_channel
            figure_filename = os.path.join(figure_folder, [strjoin(sessions(i_sess, :), '_'), '_structural_channel.pdf']);
        else
            figure_filename = os.path.join(figure_folder, [strjoin(sessions(i_sess, :), '_'), '.pdf']);
        end
        % Make a page per projection type
        temp_folder = os.path.join(GC.temp_dir, animal_ID, 'projections');
        if exist(temp_folder, 'dir'), rmdir(temp_folder, 's'), end
        mkdir(temp_folder)
%         PDF_pages = {};  % This variable will contain the list of printed PDFs
        for p = 4:n_projection_type
            % Show mean image
            data2show = nanmean(PROJECTIONS.(session_fieldname).(projection_names{p}), 3);
            fig = figure('color','w', 'units','normalized','outerposition',[0 0 1 1], 'Visible','on');
            imagesc(data2show);
            % Adjust axes
            colormap(gray)
            axis square
            set(gca, 'box','on', 'XTick',[], 'YTick',[], 'XColor','w', 'YColor','w')
            title_str = strrep(projection_names{p}, '_', ' ');
            title_str(1) = upper(title_str(1));
            title(title_str, 'FontSize',16, 'FontWeight','Bold', 'Color','w')
            
            % Export figure to .pdf
            filename = 'M:\Thomas\Manuscripts\in preparation\Invivo ACC\Material for Figures\check_points\FOV_FK_20.pdf';
%             filename = os.path.join(temp_folder, [num2str(length(PDF_pages) + 1, '%05i'), '.pdf']);
%             PDF_pages{end+1} = filename;
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters')%,fig)
            close(fig)
        end
        
        % Merge all temporary PDFs in one on the remote server
        if exist(figure_filename, 'file')
            delete(figure_filename)
        end
%         append_pdfs(figure_filename, PDF_pages{:})
%         printFilePath(figure_filename, figure_filename, false)
        % Remove folder containing temporary PDF files from disk
%         rmdir(temp_folder, 's')
    end
    
    % Empty pointer to data file, so that another session can be loaded
    if strcmp(data_type, 'epi')
        file_in = [];
    end
end

% Store data to disk
save(filename_projections, 'PROJECTIONS', '-v7.3')

% Close data file
clear file_in

% Disable toolbox for printing
if print_figures
    toggle_toolbox('_plotting', 'off')
end


%% MLint exceptions
%#ok<*AGROW>
