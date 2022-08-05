function data = read_motion_corrected_data(filename, frames_idx, frame_height, frame_width)
    
    % First input could be handle to matfile
    if ischar(filename)
        % Open file
        f = matfile(filename, 'Writable',false);
    else
        f = filename;
    end
    
    if ~isnumeric(frames_idx)
        % Make sure that frames_idx is a cell array
        if istable(frames_idx)
            frames_idx = table2cell(frames_idx);
        end
        if ~iscell(frames_idx)
            frames_idx = {frames_idx};
        end
        % Convert strings to cells
        frames_idx = cellfun(@(x) regexp(x,',','split'), frames_idx, 'UniformOutput',false);
        % Concatenate cells
        frames_idx = vertcat(frames_idx{:});
        % Convert strings to numbers
        frames_idx = cellfun(@str2double, frames_idx);
    end
    
    % Get all indices of the frames to read
    % Get list of frames belonging to these trials
    frames_list = [];
    for itrial = 1:size(frames_idx,1)
        frames_list = [frames_list, frames_idx(itrial,1):frames_idx(itrial,2)]; %#ok<AGROW>
    end
    
    % Get frame size
    if ~exist('frame_height','var') || isempty(frame_height) || ~exist('frame_width','var') || isempty(frame_width)
        [frame_height, frame_width, ~] = size(f, 'Y');
    end

    % Check whether the list of frames is contiguous
    if any(diff(frames_list) ~= 1)
        % Frames are not contiguous, which will raise the error "Cannot index
        % into 'Y' because ranges for MatFile objects must increase in equally
        % spaced intervals."
        % Allocate output variable
        data = NaN(frame_height, frame_width, length(frames_list));
        start_frame = 1;
        % Read trial by trial
        n_trials = size(frames_idx, 1);
        % Loop through trials and store it in memory
        for itrial = 1:n_trials
            frames2read = frames_idx(itrial,1):frames_idx(itrial,2);
            n_frames2read = length(frames2read);
            data(:,:,start_frame:start_frame+n_frames2read-1) = f.Y(1:frame_height, 1:frame_width, frames2read);
            start_frame = start_frame+n_frames2read;
        end
        
    else
        % Frames are contiguous, so read them in one go
        data = f.Y(1:frame_height, 1:frame_width, frames_list);
    end
