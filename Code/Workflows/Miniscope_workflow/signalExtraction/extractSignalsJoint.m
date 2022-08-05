function extractSignalsJoint(animal_ID, experiments)
%RUNJOINTEXTRACTION loads, registers and concatenates session movies and
%runs joint signal extraction for one subject

global GC LOGGER

% Load parameters
p_filename = get_filename_of('miniscope_movie_parameters', animal_ID);
p = load_variable(p_filename, 'p');

n_sessions = size(experiments, 1);

%% Get number of frames in each movie
% Load helper
local_saveDir = os.path.join(GC.temp_dir, animal_ID);
if ~exist(local_saveDir, 'dir'), mkdir(local_saveDir), end

% Get video height and width
position = p.user.position{1};
video_height = length(position(2):position(2)+position(4)-1);
video_width = length(position(1):position(1)+position(3)-1);

% Initialize helper structure
helper = struct();
helper.joint = struct();
helper.joint.zero_mask = zeros(video_height, video_width, n_sessions);
helper.joint.imref2d = cell(n_sessions, 1);

n_frames_per_session = NaN(n_sessions, 1);
for i_sess = 1:n_sessions    
    this_experiment = experiments(i_sess, :);
    moviePath = get_filename_of('miniscope_preprocessed_movie', animal_ID, this_experiment{1}, this_experiment{2});

    % Get info
    info = h5info(moviePath);
    chunkSize = info.Datasets.Dataspace.Size;
        
    % Register movies with validated registration coordinates
    ov = imref2d(chunkSize([1, 2]));
    n_frames_per_session(i_sess) = chunkSize(3);
    helper.joint.imref2d{i_sess} = ov;
end

% Store info from this loop
helper.preprocessed_total_n_frames = sum(n_frames_per_session);
helper.preprocessed_frame_edges = [[1; cumsum(n_frames_per_session(1:end-1)) + 1], cumsum(n_frames_per_session)];

%% crop, normalize, concatenate and save movies
LOGGER.info([char(datetime('now')) ' Cropping, normalizing, concatenating and saving movies'])

% Compute mean of each frame and store flattened version of data
cropped_width = video_width;
cropped_height = video_height;

% Load transformation matrices
% filename = os.path.join(GC.data_root_path, GC.registered_images, animal_ID, [animal_ID, '_alignedCellMaps.mat']);
% tforms = load_variable(filename, 'tforms');

% Make temporary file
temporary_file = os.path.join(local_saveDir, [animal_ID '_joint.h5']);
% Overwrite temporary files
if exist(temporary_file, 'file')
    delete(temporary_file)
end
h5create(temporary_file, '/1', [cropped_height, cropped_width, helper.preprocessed_total_n_frames], 'Datatype', 'single');

for i_sess = 1:n_sessions
    this_experiment = experiments(i_sess, :);
    LOGGER.info(['Loading movie ', num2str(i_sess), '/', num2str(n_sessions), ': ', this_experiment{1}, ' - ' this_experiment{2}])
    moviePath = get_filename_of('miniscope_preprocessed_movie', animal_ID, this_experiment{1}, this_experiment{2});
    movie = loadMovie(moviePath);

    % Convert movie to 'single' data type
    movie = single(movie);
    
    % Align movie (This is done by the GUI_ROI_segmentation, so there is no need to do it here)
    %movie = imwarp(movie, tforms{i_sess}, 'OutputView',helper.joint.imref2d{i_sess}, 'FillValues',NaN);

    % Normalize movie
    movie_sd = nanstd(movie(:));
    movie_avg = nanmean(movie(:));
    movie = movie - movie_avg;
    movie = movie ./ movie_sd;
    movie(isnan(movie)) = 0;
    n_frames = size(movie, 3);

    h5write(temporary_file, '/1', movie, [1, 1, helper.preprocessed_frame_edges(i_sess, 1)], [size(movie, 1), size(movie, 2), n_frames])
end

% Store helper structure in p
p.user.jointExtraction = helper;
save(p_filename, 'p');

