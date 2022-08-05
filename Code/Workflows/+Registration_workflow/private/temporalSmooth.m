function stack_in_raw_out = temporalSmooth(stack_in_raw, window_size, add_last_frame)
    % Get size of data
    [frame_height, frame_width, n_frames] = size(stack_in_raw);

    % Assign default values and check inputs
    if ~exist('add_last_frame','var') || isempty(add_last_frame)
        add_last_frame = true;
    end
    if window_size > n_frames
        window_size = n_frames;
    end

    % Preallocate variable
    stack_in_raw_out = zeros(frame_height, frame_width, n_frames);
    % Loop until the frame before the last
    for iframe = 1:n_frames-1
        % Get index of first and last frame
        start_frame = iframe;
        end_frame = min([start_frame + window_size, n_frames]);
        % Compute mean
        stack_in_raw_out(:,:,iframe) = nanmean(stack_in_raw(:,:,start_frame:end_frame), 3);
    end

    % Add last frame
    if add_last_frame
        stack_in_raw_out(:,:,end) = stack_in_raw(:,:,end);
    else
        stack_in_raw_out(:,:,end) = [];
    end
