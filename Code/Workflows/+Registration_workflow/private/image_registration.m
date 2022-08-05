function image_registration(METADATA, missing_only, skip)
% This functin performs image registration.

% Disable warnings raised in normcorre_batch while the Parallel Pool is active
warning_state = warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary');

% Parse inputs
p = inputParser();
addOptional(p, 'missing_only', false)  % Do all the sessions
addOptional(p, 'skip', false)  % Registration is skipped
parse(p, missing_only, skip)
% Fetch inputs
missing_only = p.Results.missing_only;
skip = p.Results.skip;

% Get LOGGER and general_configs
global LOGGER GC

% Set options for memory-intensive options
corr_batch_size = compute_max_batch_size(100, 85, 10000);
max_buffer_width = compute_max_batch_size(50, 85, 1000);
max_mem_batch_size = compute_max_batch_size(5000, 85, 10000);

% Get name of mouse and whether imaging has a structural channel
animal_ID = unique(METADATA.animal_ID);
animal_ID = animal_ID{1};
has_structural_channel = any(METADATA.has_structural_channel);

% Get info on frame size and rate
frame_height = max(METADATA{:, 'frame_height'});
frame_width  = max(METADATA{:, 'frame_width'});
frame_size   = max([frame_height, frame_width]);

% Perform this analysis by session
sessions = unique(METADATA(:, 'date'), 'rows', 'stable');
n_sessions = height(sessions);

% Make filenames
raw_data_filenames = cell(n_sessions, 1);
motion_corrected_filenames = cell(n_sessions, 1);
if has_structural_channel
    struct_raw_data_filenames = raw_data_filenames;
    struct_motion_corrected_filenames = motion_corrected_filenames;
end
for isess = 1:n_sessions
    raw_data_filenames{isess} = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '__', num2str(isess), '.mat']);
    motion_corrected_filenames{isess} = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '__', num2str(isess), '__motion_corrected_data.mat']);
    if has_structural_channel
        struct_raw_data_filenames{isess} = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '__', num2str(isess), '_structural_channel.mat']);
        struct_motion_corrected_filenames{isess} = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '__', num2str(isess), '__motion_corrected_structural_channel.mat']);
    end
end

% Initialize flag to know whether at least one session has been processed
processed_a_session = false;

% Loop through sessions
for isess = 1:n_sessions
    % If user requested to process only missing sessions, check whether all
    % output files of this loop are present
    if missing_only
        process_this_session = ~exist(raw_data_filenames{isess}, 'file') || ~exist(motion_corrected_filenames{isess}, 'file') || ...
            has_structural_channel && (~exist(struct_raw_data_filenames{isess}, 'file') || ~exist(struct_motion_corrected_filenames{isess}, 'file'));
        if ~process_this_session
            LOGGER.info(['Skipped ', strjoin(sessions{isess,:},' - ')])
            continue
        end
    end

    % Get trials corresponding to this session
    trials_idx = find(ismember(METADATA(:, 'date'), sessions(isess,:), 'rows'));
    tiff_files = table2cell(METADATA(trials_idx, 'tiff_path'));
    % Add root path
    tiff_files = strcat(GC.data_root_path, filesep(), GC.tiff_raw_path, filesep(), tiff_files);

    % Check whether this FOV contains a structural channel
    session_has_structural_channel = max(logical(unique(METADATA{trials_idx, 'has_structural_channel'})));

    % Get frame rate
    frame_rate = max(METADATA{trials_idx, 'frame_rate'});
    % Get total number of frames and trials in this session
    n_frames = sum(METADATA{trials_idx, 'n_frames'});
    n_trials = length(tiff_files);

    % Get beginning and end of each stack
    trials_frames_idx = cumsum(METADATA{trials_idx, 'n_frames'});
    trials_frames_idx = [[1;trials_frames_idx(1:end-1)+1], trials_frames_idx(:)];
    
    % If left by previous runs of the algorithm
    if exist(raw_data_filenames{isess}, 'file'), delete(raw_data_filenames{isess}), end
    if session_has_structural_channel && exist(struct_raw_data_filenames{isess}, 'file'), delete(struct_raw_data_filenames{isess}), end
    % Open file with writing access and without loading it into memory
    all_data_file = matfile(raw_data_filenames{isess}, 'Writable',true);
    if session_has_structural_channel
        all_struct_data_file = matfile(struct_raw_data_filenames{isess}, 'Writable',true);
    end
    
    % Allocate variable data type and then size
    all_data_file.Y = cat(3, single(0), single(0));
    all_data_file.Y(frame_size, frame_size, n_frames) = single(0);
    if session_has_structural_channel
        all_struct_data_file.Y = cat(3, single(NaN), single(NaN));
        all_struct_data_file.Y(frame_size, frame_size, n_frames) = single(NaN);
    end
    LOGGER.info(['Concatenating ', strjoin(sessions{isess,:},' - '), ' (', num2str(n_trials), ' trials, ', num2str(n_frames), ' frames)'])
    % Loop trials for registration
    for itrial = 1:n_trials
        % Log progress
        LOGGER.trace(['Processing trial ', num2str(itrial)], 'overwrite_last_message',itrial>1)
        % Load TIFF and convert it to single because NoRMCorre uses single precision arrays
        try
            TIFF = im2single(load_tiff_stack(tiff_files{itrial}));
        catch ME
            if strcmp(ME.identifier, 'tiffReader:cannotRead')
                TIFF = im2single(ScanImageTiffReader(tiff_files{itrial}).data);
            else
                rethrow(ME)
            end
        end
        % Fix aspect ratio
        TIFF = imresize(TIFF, [frame_size, frame_size], 'nearest', 'Antialiasing',true);
        % Write to memory-mapped file the functional and possibly the structural channel
        trial_has_structural_channel = logical(METADATA{trials_idx(itrial), 'has_structural_channel'});
        if trial_has_structural_channel
            all_data_file.Y(1:frame_size, 1:frame_size, trials_frames_idx(itrial,1):trials_frames_idx(itrial,2)) = TIFF(:, :, 1:2:end);
            all_struct_data_file.Y(1:frame_size, 1:frame_size, trials_frames_idx(itrial,1):trials_frames_idx(itrial,2)) = TIFF(:, :, 2:2:end);
        else
            all_data_file.Y(1:frame_size, 1:frame_size, trials_frames_idx(itrial,1):trials_frames_idx(itrial,2)) = TIFF;
        end
    end

    % Correlate each frame to the mean frame
    LOGGER.trace('Computing correlation to the mean')
    if session_has_structural_channel
        [corr_to_mean, ~, ~] = motion_metrics(all_struct_data_file, [], corr_batch_size);
    else
        [corr_to_mean, ~, ~] = motion_metrics(all_data_file, [], corr_batch_size);
    end
    % Take the most correlated frames to make a template
    n_frames_template = min([n_frames, GC.motion_correction_init_batch]);
    corr_to_mean = [corr_to_mean(:), (1:n_frames)'];
    corr_to_mean = sortrows(corr_to_mean, -1);  % Sort in descending order
    % Take indices of most correlated frames
    template_idx = corr_to_mean(1:n_frames_template, 2);
    template_idx = sort(template_idx);
    % Loop through these frames
    motion_correction_template = zeros(frame_size, frame_size);
    for iframe = 1:n_frames_template
        if session_has_structural_channel
            motion_correction_template = motion_correction_template + all_struct_data_file.Y(1:frame_size, 1:frame_size, template_idx(iframe));
        else
            motion_correction_template = motion_correction_template + all_data_file.Y(1:frame_size, 1:frame_size, template_idx(iframe));
        end
    end
    % Divide by number of frames to get mean
    motion_correction_template = motion_correction_template ./ n_frames_template;
    % Make sure that there are no NaNs
    motion_correction_template(isnan(motion_correction_template)) = 0;
    % Close file
    clear all_data_file all_struct_data_file

    if skip
        LOGGER.info('Skipping motion correction')
    else
        % Make options for NoRMCorre
        NoRMCorre_options = NoRMCorreSetParms('iter',           1, ...
                                              'd1',             frame_size, ...
                                              'd2',             frame_size, ...
                                              'grid_size',      [ceil(frame_size/4), ceil(frame_size/4)], ... % from paper
                                              'overlap_pre',    [ceil(frame_size/16), ceil(frame_size/16)], ...  % from paper
                                              'mot_uf',         4, ...
                                              'max_dev',        ceil(8 ./ (512 ./ [frame_size, frame_size])), ...  % because the default value was 8 for a FOV of 512
                                              'bin_width',      round(frame_rate * GC.motion_correction_bin_width), ...
                                              'max_shift',      ceil(20 ./ (512 ./ [frame_size, frame_size])), ... % because the default value was 20 for a FOV of 512
                                              'us_fac',         100, ...
                                              'init_batch',     round(frame_rate * GC.motion_correction_init_batch), ...
                                              'boundary',       'zero', ...
                                              'use_parallel',   true,...
                                              'phase_flag',     false, ...
                                              'shifts_method',  'FFT', ...
                                              'correct_bidir',  false, ...
                                              'buffer_width',   max_buffer_width, ...
                                              'upd_template',   false, ...  % force the use of the computed template
                                              'memmap',         true, ...
                                              'output_type',    'memmap', ...
                                              'mem_filename',   motion_corrected_filenames{isess}, ...
                                              'mem_batch_size', max_mem_batch_size);
        % Check that number of frames is correct for batch-processing
        if NoRMCorre_options.bin_width > n_frames, NoRMCorre_options.bin_width = n_frames; end
        if NoRMCorre_options.init_batch > n_frames, NoRMCorre_options.init_batch = n_frames; end
        % Get whether these trials need to be fixed for distortions due to
        % bi-directional scanning
        NoRMCorre_options.correct_bidir = any(logical(METADATA{trials_idx, 'needs_fix_for_bidirectional_scanning'}));

        % Store output filename
        if session_has_structural_channel
            % Re-route output
            NoRMCorre_options.mem_filename = struct_motion_corrected_filenames{isess};
        end

        % Delete previous local file to make sure that a new one will be created
        if exist(motion_corrected_filenames{isess}, 'file'), delete(motion_corrected_filenames{isess}), end
        if session_has_structural_channel && exist(struct_motion_corrected_filenames{isess}, 'file'), delete(struct_motion_corrected_filenames{isess}), end
        % Apply non-rigid motion correction
        if session_has_structural_channel
            LOGGER.info('Step 1/2: Performing motion correction on structural channel')
            [~, shifts, ~, ~, col_shift] = normcorre_batch(struct_raw_data_filenames{isess}, NoRMCorre_options, motion_correction_template);
            % Set output to final file
            NoRMCorre_options.mem_filename = motion_corrected_filenames{isess};
            LOGGER.info('Step 2/2: Performing motion correction on functional channel')
            apply_shifts(raw_data_filenames{isess}, shifts, NoRMCorre_options, 0, 0, 0, col_shift);

        else
            LOGGER.info('Performing motion correction')
            normcorre_batch(raw_data_filenames{isess}, NoRMCorre_options, motion_correction_template);
        end
    end
    
    % Toggle flag
    processed_a_session = true;
end

% Get names of input files
if skip
    input_filenames = raw_data_filenames;
else
    input_filenames = motion_corrected_filenames;
end
% Make filename for final file in the local folder
normcorre_local_filename = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '.mat']);
struct_normcorre_local_filename = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_structural_channel.mat']);

% Check whether the files with the concatenated images need be updated
concatenate_functional_channel = processed_a_session || ~exist(normcorre_local_filename,'file');
concatenate_structural_channel = has_structural_channel && (processed_a_session || ~exist(normcorre_local_filename,'file'));

% Delete previous files
if concatenate_functional_channel && exist(normcorre_local_filename,'file'), delete(normcorre_local_filename), end
if concatenate_structural_channel && exist(struct_normcorre_local_filename,'file'), delete(struct_normcorre_local_filename), end

if n_sessions == 1  % Simply copy file
    if concatenate_functional_channel
        copyfile(input_filenames{1}, normcorre_local_filename)
    end
    if has_structural_channel
        if skip
            input_filenames = struct_raw_data_filenames;
        else
            input_filenames = struct_motion_corrected_filenames;
        end
        if concatenate_structural_channel
            copyfile(input_filenames{1}, struct_normcorre_local_filename)
        end
    end

else
    if concatenate_functional_channel
        % Concatenate files
        LOGGER.info('Concatenating sessions into one file. Please be patient')
        % Open file for writing
        normcorre_local_file = matfile(normcorre_local_filename, 'Writable',true);
        % Initialize empty, 3D variable called 'Y'
        normcorre_local_file.Y = cat(3, single(0), single(0));
        normcorre_local_file.Y(frame_size, frame_size, 1) = single(0);
        n_total_frames = 0;
        n_frames_per_batch = 2000;
        for isess = 1:n_sessions
            % Get filename of file to copy
            file_to_copy = input_filenames{isess};
            % Log session filename
            LOGGER.trace(['Copying ', file_to_copy], 'contains_path',true)
            % Open file fore reading and get number of frames of it
            file = matfile(file_to_copy, 'Writable',false);
            n_frames = size(file, 'Y');
            n_frames = n_frames(3);
            % Copy file in batches
            start_frame = n_total_frames + 1;
            n_batches = ceil(n_frames / n_frames_per_batch);
            for ibatch = 1:n_batches        
                % Calculate the first and last frame of the batch
                start_frame_batch = (ibatch-1) * n_frames_per_batch + 1;
                end_frame_batch = min([start_frame_batch + n_frames_per_batch - 1, n_frames]);
                start_frame_target_file = start_frame + start_frame_batch - 1;
                end_frame_target_file = start_frame + start_frame_batch - 1 + (end_frame_batch - start_frame_batch);
                % Copy frames
                normcorre_local_file.Y(1:frame_size, 1:frame_size, start_frame_target_file:end_frame_target_file) = file.Y(1:frame_size, 1:frame_size, start_frame_batch:end_frame_batch);
            end
            % Update total number of frames copied
            n_total_frames = n_total_frames + n_frames;
        end
        % Close files
        clear normcorre_local_file file
    end
    
    if concatenate_structural_channel
        % Open file for writing
        normcorre_local_file = matfile(struct_normcorre_local_filename, 'Writable',true);
        % Initialize empty, 3D variable called 'Y'
        normcorre_local_file.Y = cat(3, single(0), single(0));
        normcorre_local_file.Y(frame_size, frame_size, 1) = single(0);
        n_total_frames = 0;
        n_frames_per_batch = 2000;
        for isess = 1:n_sessions
            % Get filename of file to copy
            file_to_copy = input_filenames{isess};
            % Log session filename
            LOGGER.trace(['Copying ', file_to_copy], 'contains_path',true)
            % Open file fore reading and get number of frames of it
            file = matfile(file_to_copy, 'Writable',false);
            n_frames = size(file, 'Y');
            n_frames = n_frames(3);
            % Copy file in batches
            start_frame = n_total_frames + 1;
            n_batches = ceil(n_frames / n_frames_per_batch);
            for ibatch = 1:n_batches
                % Calculate the first and last frame of the batch
                start_frame_batch = (ibatch-1) * n_frames_per_batch + 1;
                end_frame_batch = min([start_frame_batch + n_frames_per_batch - 1, n_frames]);
                start_frame_target_file = start_frame + start_frame_batch - 1;
                end_frame_target_file = start_frame + start_frame_batch - 1 + (end_frame_batch - start_frame_batch);
                % Copy frames
                normcorre_local_file.Y(1:frame_size, 1:frame_size, start_frame_target_file:end_frame_target_file) = file.Y(1:frame_size, 1:frame_size, start_frame_batch:end_frame_batch);
            end
            % Update total number of frames copied
            n_total_frames = n_total_frames + n_frames;
        end
        % Close files
        clear normcorre_local_file file
    end
end

% Copy concatenated file to server
if concatenate_functional_channel
    LOGGER.info('Copying functional channel to server. Please be patient')
    % Copy file to server
    remote_filename = get_filename_of('motion_corrected', METADATA.animal_ID{1});
    copyfile(normcorre_local_filename, remote_filename)
end
if concatenate_structural_channel
    LOGGER.info('Copying structural channel to server. Please be patient')
    % Copy file to server
    remote_filename = get_filename_of('motion_corrected', METADATA.animal_ID{1}, true);
    copyfile(struct_normcorre_local_filename, remote_filename)
end

LOGGER.trace('Updating database')
% Convert list of indices to comma-separated strings
trials_frames_idx = cumsum(METADATA{:, 'n_frames'});
trials_frames_idx = [[1;trials_frames_idx(1:end-1)+1], trials_frames_idx(:)];
trials_frames_str = mat2cell(trials_frames_idx, ones(height(METADATA),1), 2);
trials_frames_str = cellfun(@(x) strjoin(value2str(x,'%d'),','), trials_frames_str, 'UniformOutput',false);
% Update database with position of each trial in the newly copied file
SQL_database.update('trials', 'frames_idx', trials_frames_str, METADATA.tiff_path, 'tiff_path')

% Re-enable default warnings
warning(warning_state)
