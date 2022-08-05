function DONE = adjust_FOV_from_registration(animal_ID, experiment_type)
%% Preamble
% This functions will take the registered images (either h5 or mat files)
% that were created during the registration procedure. It then adjust for
% FOV transformation (only translation at the moment) and generate a h5
% file that will be used for further analysis (automatic ROI detection,
% using CNMF (-e)).
% Input = - animal_ID
%         - experiment_type : either 'epi' or '2p'
% Output = DONE : a flag that passess to the next script signaling that
% the video was successfully created.

%% Init
%Load the ROI transformation matrices and idx in the ROI info file
clc
global GC
fprintf('%% INITIALIZING FOV ADJUSTMENT %% \n')
debug_mode = 0;
DONE = 0;
warning('off')
%% Load parameters
is_epi = strcmp(experiment_type, 'epi');

% animal_ID = 'MA_9';
% is_epi = 0;
temp_dir = GC.temp_dir;

local_folder = os.path.join(temp_dir, animal_ID);%['D:\_MATLAB_CaImaging\test_CNMFE'];
%     local_folder = ['D:\_MATLAB_CaImaging\test_CNMFE'];
% roi_folder = 'V:\Ca_imaging_pain\3_segmentation'; % this requires that the user first saves the data on the server
roi_folder = local_folder;
% you can load these parameters after prepare data for segmentation
roi_file = os.path.join(roi_folder,[animal_ID, '_ROI_info.mat']);
% os.path.join(temp_dir, animal_ID, [animal_ID, '_ROI_info.mat']);
disp('Loading ROI info')
ROI_info = load(roi_file);
disp('ROI info loaded')

% if the data come form 2p recording you can get it from here
% then it's ok if we don't use the params.mat file
% Read the frames and idx of transformations
% frames_idx = PARAMETERS.frames_idx; % Do later for miniscope (check prepare_data_for_segmentation.m)
if is_epi
    
    param_file =  os.path.join(temp_dir, animal_ID, [animal_ID, '_params.mat']);
    P = load(param_file);
    PARAMETERS = P.PARAMETERS;

    % load p file
    p = load(['V:\Ca_imaging_pain\1_raw_movies\', animal_ID, '\', 'p.mat']);
    frames_idx = PARAMETERS.frames_idx;
    % Load joint h5 file for epi data
    joint_filename = [animal_ID, '_joint.h5'];
    
else
    frames_idx =  get_trial_frames(animal_ID, 'reset_index',false, 'return_all_sessions',true);
    joint_filename =  [animal_ID, '.mat'];
    
end
% Set the path for either h5 or mat file
local_joint_file = os.path.join(GC.temp_dir, animal_ID, joint_filename);

n_chunks = size(frames_idx,1);
n_frames = frames_idx(end,end);
% frame_size = [PARAMETERS.frame_height, PARAMETERS.frame_width];
ROIs_size = size(ROI_info.ROIs);
frame_size = ROIs_size([1,2]);

transf_idx = ROI_info.TRANSFORMATION_IDX;
transf_matrices = ROI_info.TRANSFORMATION_MATRICES;
MAPPING_TABLE = ROI_info.MAPPING_TABLE;

%% Load motion correected file
% in debug mode relate it to another folder otherwise delete these 3 lines:
if ~exist(local_folder, 'dir')
    mkdir(local_folder)
end
dataset_mov = '/1';
% if reading from miniscope data, get all the h5 files
% otherwise get the .mat file
if is_epi
    % In the normal pipeline, there is a _joint.h5 file that contains all
    % the sessions.
    % We need to take this file otherwise (when taking an already analysed
    % video) we have to concatenate all the h5 files
    filepath = 'K:\Ca_imaging_pain\2_motion_corrected_movies\';
    
    if ~exist(local_joint_file, 'file')
        %         disp('Loading previous h5 files and concatenating them, be patient')
        %         keyboard
        %         %%%%%%%%%%%%%%%%%%%% Modify the next lines!! %%%%%%%%%%%%%%%%%%%%%%
        %         files_folder = os.path.join(filepath, animal_ID);
        %         names = dir(files_folder);
        %         names = {names.name};
        %         h5toread = names(endsWith(names, '.h5'));
        %         h5toread = sort(h5toread)';
        %         prep_filename = os.path.join(local_folder, [animal_ID, '_joint.h5']);
        %         dataset_mov = '/1';
        %
        % %         n_frames = P.p.user.jointExtraction.preprocessed_total_n_frames;
        %         % Create h5 file
        %         h5create(prep_filename, '/1', [frame_size, n_frames]);
        %
        %         % Loop through each chunk
        %         n_files = length(h5toread);
        %
        %         for i_file = 1:n_files
        %             disp(['Writing now file nr: ',num2str(i_file)])
        %             this_file = os.path.join(files_folder,h5toread{i_file});
        %             % read the file
        %             start_frame = frames_idx(i_file,1);
        %             end_frame = frames_idx(i_file,2);
        %             n_chunk_frames = end_frame - start_frame + 1;
        %             data = h5read(this_file,dataset_mov,[1,1,frames_idx(i_file,1)],[frame_size(1), frame_size(2),n_frames_this_chunk]);
        % %             h5write(joint_file, dataset_mov, data, [1, 1, start_frame], [frame_size(1), frame_size(2),end_frame])
        %
        %             % Correct for FOV adjustment
        %             this_transf_idx = cell2mat(transf_idx(i_file,2)) + 1; % Since it comes from python, we have to add 1
        %             if  isempty(this_transf_idx) % data do not need to be corrected
        %                 data_to_write = data;
        %                 %             keyboard
        %                 h5write(prep_filename, dataset_mov, data, [1, 1, start_frame], [frame_size(1), frame_size(2),n_chunk_frames])
        %             else % we will correct the data based on the manual alignment
        %                 % Correct for displacement
        %                 mask = single(zeros(size(data)));
        %                 %         this_idx = find(ismember(cell2mat(transf_matrices(:,2)), this_transf_idx));
        %
        %                 frames_to_move = transf_matrices{this_transf_idx, 3}{1};
        %                 frames_to_move = frames_to_move([1,3]); % vertical and horizontal
        %                 %             mask(frames_to_move(1):end -1 , frames_to_move(2):end-1,:) = data(1:end-frames_to_move(1),1:end-frames_to_move(2),:);
        %                 mask(frames_to_move(1)+1:end , frames_to_move(2)+1:end,:) = data(1:end-frames_to_move(1),1:end - frames_to_move(2),:);
        %                 data_to_write = mask;
        %                 % Write on the 5 file
        %                 h5write(prep_filename, dataset_mov, data_to_write, [1, 1, start_frame], [frame_size(1), frame_size(2),n_chunk_frames])
        %
        %             end
        %             chunks_images(:,:,i_chunk) = sum(data_to_write,3);
        %             deleted_pixels(i_chunk,:) = frames_to_move';
        %
        %
        %
        %
        %
        %         end
        is_h5 = 1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        mov_info = h5info(local_joint_file);
        
        movie_file = local_joint_file;
    end
    is_h5 = 1;
else
    % for 2p pipeline
    if debug_mode
        movie_file = os.path.join(filepath, animal_ID,'.mat'); % read it from server
        % or read it from local drive
    else
        movie_file = local_joint_file;
    end
    mov_info = matfile(movie_file);
    is_h5 = 0;
end

%% Create a h5 file in the local folder
prep_filename = os.path.join(local_folder, [animal_ID, '_joint_corrected.h5']);
if exist(prep_filename, 'file')
    DONE = 1;
    return
end
h5create(prep_filename, dataset_mov, [frame_size, n_frames], 'Datatype', 'single');
% deleted_pixels = [];
% chunks_images = [];
% max_chunks = [];
% IDX_trans = [];
% Disp

%% Run FOV adjustment
disp('Reading and writing corrected file')
nbytes = fprintf('processing chunk 0 of %d', n_chunks);

tic
try
    for i_chunk = 1:n_chunks
        %         fprintf('Writing now chunk nr: %2d \n',i_chunk)
        % read a chunk
        this_frame = frames_idx(i_chunk,:);
        start_frame = this_frame(1);
        finish_frame = this_frame(2);
        n_chunk_frames = finish_frame - start_frame+1;
        if is_h5
            data = h5read(movie_file,dataset_mov,[1,1,start_frame],[frame_size,n_chunk_frames]);
        else
            data =  single(mov_info.Y(:, :, this_frame(1):this_frame(2)));
        end
        % Get the coordinates to crop
        this_transf_idx = cell2mat(transf_idx(i_chunk,2)) + 1; % Since it comes from python, we have to add 1
        % Read original file
        if  isempty(this_transf_idx) % data do not need to be corrected
            data_to_write = data;
            h5write(prep_filename, dataset_mov, data, [1, 1, start_frame], [frame_size(1), frame_size(2),n_chunk_frames])
        else % we will correct the data based on the manual alignment
            map_tab = MAPPING_TABLE(i_chunk,:);
            map_reshaped = reshape(~isnan(map_tab), frame_size(1),frame_size(2));
            map_idx_reshaped = reshape(map_tab, frame_size(1),frame_size(2));
            %             map_idx_reshaped = map_idx_reshaped'; % turn it around
            map_idx_reshaped = map_idx_reshaped(:) +1; % vectorize it and start at 1
            map_idx_reshaped_idx = ~isnan(map_idx_reshaped);
            map_idx = find(~isnan(map_idx_reshaped));
            map_idx_find = find(map_idx_reshaped_idx);
            %             map_reshaped = map_reshaped';
            % use this matrix as a mask. Add the pixels from the data to
            % points where ~isnan of map_tab
            %             mask_vector_init = map_reshaped(:); % vectorize to get each pixel
            %             mask_vector_init = logical(mask_vector_init);
            data_to_write = [];
            for i_frame = 1:size(data,3)
                this_data = data(:,:,i_frame)';
                data_vector = this_data(:);
                mask_vector = zeros(size(data_vector));
                %                 where_idx = mask_vector_init;
                % translate the image
                %                 to_idx = find(~isnan(map_tab))';
                % Assign new identity to each pixel
                for ipixel = 1:length(find(map_idx_find)) % nr pixels in image
                    this_pixel = map_idx_find(ipixel);
                    this_pixel_data = map_idx_reshaped(this_pixel)  ;
                    mask_vector(map_idx(ipixel)) = data_vector(this_pixel_data);
                end
                % Reshape and rotate
                data_cropped = reshape(mask_vector, frame_size(2),frame_size(1))';
                % write it back to matrix
                data_to_write(:,:,i_frame) = data_cropped;
            end
            % Write data into the h5 file
            h5write(prep_filename, dataset_mov, data_to_write, [1, 1, start_frame], [frame_size(1), frame_size(2),n_chunk_frames])
        end
        
        % for evaluation
        %         chunks_images(:,:,i_chunk) = mean(data_to_write,3);
        %         deleted_pixels(i_chunk,:) = frames_to_move';
        %         rmax_chunks(:,:,i_chunk) = max(data_to_write,[],3);
        % check the transformation idxs
        %         if isempty(this_transf_idx)
        %             IDX_trans(i_chunk,1) = 0;
        %         else
        %             IDX_trans(i_chunk,1) = this_transf_idx;
        %         end
        while nbytes > 0
            fprintf('\b')
            nbytes = nbytes - 1;
        end
        nbytes = fprintf('processing chunk %d of %d \n', i_chunk, n_chunks);
        
    end
catch ME
    keyboard
    DONE = 0;
    return
end
T = toc;
fprintf('### DONE WITH ADJUSTMENT ### \nTime elapsed: %2f mins \n', T/60 )
DONE = 1;
end

%%
%#ok<*SAGROW>
