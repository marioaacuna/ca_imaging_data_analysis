function DONE = Run_CNMFe(animal_ID, exp_type)

%% clear the workspace and select data 
clc; close all; 
disp('%% RUNNING CNMF-e %%')
global GC LOGGER
%% Set paths
% animal_ID = 'MA_9';
% prep_filename = ['D:\_MATLAB_CaImaging\test_CNMFE\', animal_ID, '.mat'];
prep_filename = os.path.join(GC.temp_dir, animal_ID, [animal_ID,'_joint_corrected.h5']);
% prep_filename = ['D:\_MATLAB_CaImaging\test_CNMFE\', animal_ID, '_210727_SP.h5'];

cnmfe_folder = 'C:\Users\acuna\Documents\CNMF_E';
cnmfe_loaded = 0;
if ~exist('cnmfe_loaded', 'var') || ~cnmfe_loaded
    addpath(fullfile(cnmfe_folder, 'ca_source_extraction'));
    addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'utilities')));
    addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'endoscope')));
    addpath(fullfile(cnmfe_folder, 'GUI'));
    addpath(fullfile(cnmfe_folder, 'GUI', 'gui_callbacks'));
    addpath(fullfile(cnmfe_folder, 'GUI', 'modules'));
    addpath(fullfile(cnmfe_folder, 'OASIS_matlab'));
    addpath(fullfile(cnmfe_folder, 'scripts'));
    cnmfe_loaded = true;
end
oasis_setup()

is_epi = strcmp(exp_type, 'epi');

%% check if file exists already
neuron_path = os.path.join(GC.temp_dir, animal_ID);
local_neuron_filename = os.path.join(neuron_path, 'neuron_temp_results.mat');

overwrite = 0;
exists_temp_file = any(exist(local_neuron_filename, 'file')) ;
if exists_temp_file
    redo_temp = questdlg('Temp file aready exists, do you want to overwrite it (re do it)?');
    if strcmp(redo_temp, 'Yes')
        overwrite =1;
    end
end

if ~exists_temp_file || overwrite
%% Create the Sources2D obj
neuron = Sources2D(); 
%% Check if tun big data or not
prep_filename_info = h5info(prep_filename);
size_prep_file = prep_filename_info.Datasets.Dataspace.Size;
if ~is_epi
    patch_dims = [size_prep_file(1)/2,size_prep_file(2)/2]; % for 2x2 patches
else
    patch_dims = [size_prep_file(1)/4,size_prep_file(2)/4]; % for 4x4 patches
end

dir_prep_file = dir(prep_filename);         
filesize = dir_prep_file.bytes;       
size_of_patch = filesize/(2*(size_prep_file(1)/patch_dims(1)));
%%
if size_of_patch < 5*10^9 % this did not work, due to future mempry problems
    run_big_data = 1;
else
    run_big_data = 1;

end

%%
% nam = prep_filename;          
nam = neuron.select_multiple_files({prep_filename});  %if nam is [], then select data interactively 
%% Get file size
% m = matfile(prep_filename);
% m_info = whos(m);
% s_file = m_info.size;
% patch_dims = [s_file(1)/4,s_file(2)/4]; 
% patch_dims = [128,128]; 
%% parameters  
% -------------------------    COMPUTATION    -------------------------  %
if run_big_data
    pars_envs = struct('memory_size_to_use', 3, ...   % GB, memory space you allow to use in MATLAB
        'memory_size_per_patch', 0.3, ...   % GB, space for loading data within one patch
        'patch_dims', patch_dims, ... % patch size
        'batch_frames', 8000);   % nr of frames per batch
else
    pars_envs = struct('memory_size_to_use', 12, ...   % GB, memory space you allow to use in MATLAB
        'memory_size_per_patch', 12, ...   % GB, space for loading data within one patch
        'patch_dims', patch_dims/2); % patch size
end
% -------------------------      SPATIAL      -------------------------  %
gSig = 1;%0.5;      % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
ssub = 2;           % spatial downsampling factor
with_dendrites = false;   % with dendrites or not 
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';  %#ok<UNRCH>
    updateA_bSiz = 20;
    updateA_dist = neuron.options.dist; 
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse'; %#ok<UNRCH>
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
min_pixel = 80; % min size neuron in px % 120 identified only 19 neurons (with 30, more than 1000)
min_size = 80; % still don't know for what this is


% -------------------------      TEMPORAL     -------------------------  %
if is_epi
    Fs = GC.epifluorescence_downsample_to_frame_rate; % frame rate
    dmin_only = 15;%10;  % merge neurons if their distances are smaller than dmin_only. 
     bd = 30;%1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)

else
%     keyboard % set up later
    Fs = table2array(unique(SQL_database.read_table_where('trials', {'frame_rate'}, animal_ID,'animal_ID')));  % frame rate
    dmin_only = 7;  % merge neurons if their distances are smaller than dmin_only. 
     bd = 25;%1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)

end
tsub = 2;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level 
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true);% optimize the baseline); 
nk = 1;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
                    % when changed, try some integers smaller than total_frame/(Fs*30) 
detrend_method = 'local_min';  % compute the local minimum as an estimation of trend. 

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'ring';%'svd';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
ring_radius = 13;  % when the ring model used, it is the radius of the ring used in the background model. 
                    %otherwise, it's just the width of the overlapping area 

% -------------------------      MERGING      -------------------------  %
show_merge = false;  % if true, manually verify the merging step 
merge_thr = 0.25;%0.5;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'mean';   % method for computing neuron distances {'mean', 'max'}
dmin = 10;%3;       % minimum distances between two neurons. it is used together with merge_thr
% dmin_only = 7;  % merge neurons if their distances are smaller than dmin_only. 

% -------------------------  INITIALIZATION   -------------------------  %
K = 30;%[];             % maximum number of neurons per patch. when K=[], take as many as possible 
answer = MessageBox('Do you think you have few neurons?', 'YES', 'No');
if strcmp(answer, 'Yes')
    min_corr =0.45;     % minimum local correlation for a seeding pixel % 2p = 0.65
    min_pnr = 3;       % minimum peak-to-noise ratio for a seeding pixel
    min_pixel = 30; % min size neuron in px % 120 identified only 19 neurons (with 30, more than 1000)
    min_size = 30; % still don't know for what this is
    gSiz = 30;          % pixel, neuron diameter 

else
    min_corr = 0.65;     % minimum local correlation for a seeding pixel % 2p = 0.65
    min_pnr = 5;       % minimum peak-to-noise ratio for a seeding pixel
    min_pixel = 80; % min size neuron in px % 120 identified only 19 neurons (with 30, more than 1000)
    min_size = 80; % still don't know for what this is
    gSiz = 80;          % pixel, neuron diameter 

end
% min_pixel = 2^2;      % minimum number of nonzero pixels for each neuron
% bd = 25;%1;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames 
save_initialization = false;    % save the initialization procedure as a video. 
use_parallel = true;    % use parallel computation for parallel computing 
show_init = true;   % show initialization results 
choose_params = false; % manually choose parameters 
center_psf = false;  % set the value as true when the background fluctuation is large (usually 1p data) 
                    % set the value as false when the background fluctuation is small (2p)

% ----------------------   MANUAL INTERVENTION  --------------------  %
with_manual_intervention = true; 

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not 
kt = 1;                 % frame intervals 

% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial -------- 
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'tsub', tsub, ...                       % -------- temporal -------- 
    'deconv_options', deconv_options, ...    '
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background -------- 
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization ----- 
    'min_pnr', min_pnr, ...
    'min_size', min_size, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf);
neuron.Fs = Fs; 

%% distribute data and be ready to run source extraction 
if run_big_data
    neuron.getReady_batch(pars_envs);
    choose_params = 0;
else
    neuron.getReady(pars_envs);
    choose_params = 1;
end
%% initialize neurons from the video data within a selected temporal range 
if choose_params
    % change parameters for optimized initialization
    [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
end

if run_big_data
    neuron.initComponents_batch(K, save_initialization, use_parallel);
else
    [center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
    
    % [center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
    % if show_init
    %     figure;
    %     imagesc(Cn, [0, 1]); colormap gray;
    %     hold on;
    %     plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
    % end
end
neuron_init = neuron.copy(); 

%% udpate spatial components for all batches
neuron.update_spatial_batch(use_parallel);

%% udpate temporal components for all bataches
neuron.update_temporal_batch(use_parallel);

%% update background
neuron.update_background_batch(use_parallel);
%% Copy results to a new structure
neuron_bkg = neuron.copy(); % copy results to a new object in case next steps don't work

%% Merge 
% % merge neurons
% n_batches = length(neuron_final.batches);
% for ibatch = 1:n_batches
%     neuron_final.batches{ibatch}.neuron.merge_neurons_dist_corr(1)
%     neuron_final.batches{ibatch}.neuron.merge_close_neighbors(show_merge, dmin_only)
% end

% if you merge neurons per batch before correlation_pnr_batch then it
% doesn't work. 
% It seems that merging in batch mode doesn't work

% neuron_bkg = neuron.copy(); % Copying like this works
%% get the correlation image and PNR image for all neurons
neuron.correlation_pnr_batch();

%% concatenate temporal components
neuron.concatenate_temporal_batch();

%% deal with bad neurons 
% tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
% neuron.delete(tags>0); 
%% make back up of concatenated neurons
neuron_concat = neuron.copy();

% It seems that merging and then deleting is good
% Merge anutomatically
neuron_concat.merge_neurons_dist_corr(show_merge);
neuron_concat.merge_close_neighbors(show_merge, dmin_only);

else
    disp('Loading previous file')
    neuron_concat = load_variable(local_neuron_filename, 'neuron_results');
    if is_epi
        dmin_only = 10;  % merge neurons if their distances are smaller than dmin_only.
    else
        dmin_only = 7;  % merge neurons if their distances are smaller than dmin_only.
    end

end
%% Check neurons
disp('%% Plotting Countours %%')
% Compile images of batches
n_batches = length(neuron_concat.batches);
compiled_batches = [];
for i_batch  = 1: n_batches
    compiled_batches(:,:,i_batch) = neuron_concat.batches{i_batch}.neuron.Cn; 
end

% Plot contours on compiled image
mean_images = median(compiled_batches,3);
figure('pos', [10 10 1200 800]), 
Cn = (mean_images);
hold on
Coor = neuron_concat.Coor;
Aor = neuron_concat.A;
CC = plot_contours(Aor, Cn,[], 1, [], Coor);
selection = [];
selection = MessageBox('Do need to delete cells?', 'Close Request Function', 'Yes','No', 'No');

if strcmp(selection, 'Yes')
    % Then delete
    % check if it's re-loaded from local
    if exists_temp_file && ~strcmp(redo_temp, 'Yes')
        neuron_filename =os.path.join(GC.data_root_path,GC.segmentation_path, [animal_ID, '_Neuron_Source2D.mat']);
        neuron_exist = exist(neuron_filename, 'file');
        show_merge = false;
        if neuron_exist
            neuron = load_variable(os.path.join(GC.data_root_path,GC.segmentation_path, [animal_ID, '_Neuron_Source2D.mat'])); 
            neuron_concat = neuron;
        else
            % copy important parameters into a Sources2D
            neuron = Sources2D;
            neuron.A = neuron_concat.A;
            neuron.C = neuron_concat.C;
            neuron.C_raw = neuron_concat.C_raw;
            neuron.S = neuron_concat.S;
            neuron.batches = neuron_concat.batches;
            neuron.tags = neuron_concat.tags;
            neuron.ids = neuron_concat.ids;
            neuron.Fs = neuron_concat.Fs;
            neuron.P = neuron_concat.P;
            neuron.kernel = neuron_concat.kernel;
            neuron.b = neuron_concat.b;
            neuron.f = neuron_concat.f;
            neuron.W = neuron_concat.W;
            neuron.b0 = neuron_concat.b0;
            neuron.options = neuron_concat.options;
            
            % replace neuron_concat
            neuron_concat = neuron;
            
        end
    end
    neuron_concat.viewNeurons([], neuron_concat.C_raw);
    % Run merge manual in case
%     try
%         cnmfe_manual_merge();
%     catch
%         disp('something went wrong, but dont worry')
%     end
    % run again Automatic merge
    neuron_concat.merge_neurons_dist_corr(show_merge);
    neuron_concat.merge_close_neighbors(show_merge, dmin_only + round(0.3*dmin_only));

end

%% save neuron obj locally
neuron_concat.save_results(local_neuron_filename);
%% save ROI extraction
% we have to save the fluorescence traces trasposed in: ROI_fluorescence in
% the ROI_info_file

ROI_info_filename = get_filename_of('ROI_info', animal_ID);
load(ROI_info_filename) %#ok<LOAD> %load all variables and save them again
ROI_info = CC; % X,Y position in pixel dimension. Convert to d1,d2 coordinates with reshape and draw.

%% DF/F extraction only for 2p data
if strcmp(exp_type, 'epi')
    ROI_fluorescence = neuron_concat.C_raw';
else

    disp('%% DF-Fing the traces %%')

    ROI_DFF = [];
    nbytes = 0;
    data_set = '/1';
    if is_epi
        param_file =  os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_params.mat']);
        P = load(param_file);
        PARAMETERS = P.PARAMETERS;
        %     p = load(['V:\Ca_imaging_pain\1_raw_movies\', animal_ID, '\', 'p.mat']);
        frames_idx = PARAMETERS.frames_idx;
    else
        frames_idx =  get_trial_frames(animal_ID, 'reset_index',false, 'return_all_sessions',true); % for miniscope use the PARAMETERS.mat file
    end
    n_chunks = length(frames_idx);
    for ich = 1:n_chunks
        % Chunk idx
        init_1= frames_idx(ich,:);
        init_1_idx = init_1(1,1);
        fin_1_idx = init_1(1,2);
        h5_info = h5info(prep_filename);
        size_vid = h5_info.Datasets.Dataspace.Size(1:2);
        n_frames  = fin_1_idx - init_1_idx +1;
        Y = h5read(prep_filename,data_set,[1,1,init_1_idx],[size_vid,n_frames]);
        % Y needs to be pixels by time
        Y = reshape(Y, size(Y,1)*size(Y,2), []);
        C = neuron_concat.C_raw(:,init_1_idx:fin_1_idx );
        dff = extract_DF_F(Y,neuron_concat.A,C,[],[],neuron_concat.options);
        %     dff_old =  extract_DF_F_old(Y, neuron_concat.A,C,[],neuron_concat.options);
        ROI_DFF(init_1_idx:fin_1_idx, :) = dff';
        while nbytes > 0
            fprintf('\b')
            nbytes = nbytes - 1;
        end
        nbytes = fprintf('processing chunk %d of %d \n', ich, n_chunks);
    end
    ROI_fluorescence = ROI_DFF;
end
%% UPDATE DATABASE
LOGGER.trace('Updating database')
n_ROIs = size(neuron_concat.C,1);
SQL_database.update('experiments', 'n_ROIs',n_ROIs, animal_ID,'animal_ID')

%% SAVE
LOGGER.trace('Writing ROIs to server ...')
save(ROI_info_filename, 'ROI_info', 'ROI_fluorescence', ...
    'modified_MAPPING', 'MAPPING_TABLE','TRANSFORMATION_MATRICES', 'TRANSFORMATION_IDX', '-v7.3') 
% save Neuron source2D object
final_filename = get_filename_of('source2D', animal_ID);
try
    neuron_concat.save_results(final_filename)
catch ME
end
neurons_filename = get_filename_of('ROI_map_cnmfe', animal_ID);
neuron_concat.save_neurons(neurons_filename)
close all


LOGGER.trace('done', 'append',true)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%#ok<*AGROW>
%#ok<*USENS>
DONE = 1;
end

























