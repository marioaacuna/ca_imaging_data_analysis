function concatMovies(p, movDir, position, rawMovie_filename)

    global LOGGER
    global GC
    tic;

    % Get names of local folder
    local_saveDir = fileparts(fileparts(rawMovie_filename));
    
    % Get name of subject and session number
    f = strsplit(rawMovie_filename, filesep());
    subject_name = f{end-3};
    session_number = f{end-2};
    % Make filename of helper file
    helper_file = os.path.join(local_saveDir, ['helper_', subject_name, '_', session_number, '.mat']);
    if exist(helper_file, 'file')
        load(helper_file, 'helper');
        do_create_files = false;
         
    else  % Made new one
        do_create_files = true;
        % count number of avi files
        fnames = dir(movDir);
        n_movies = sum(cellfun(@(x) endsWith(x, '.avi'), {fnames.name}'));

        % Get number of frames
        total_nFrames = zeros(n_movies, 1);
        for m = 1:n_movies
            % load movie
            movPath = os.path.join(movDir, ['msCam', num2str(m), '.avi']);
            movieObj = VideoReader(movPath);
            total_nFrames(m) = movieObj.NumberOfFrames;
        end

        % Initialize helper
        helper = struct();
    end
    
    % Initialize or store information in helper
    if ~isfield(helper, 'n_movies'), helper.n_movies = n_movies; end
    if ~isfield(helper, 'total_nFrames'), helper.total_nFrames = total_nFrames; end
    if ~isfield(helper, 'position'), helper.position = position; end
    if ~isfield(helper, 'completed_PCAICA'), helper.completed_PCAICA = false; end
    % Unpack variables
    n_movies = helper.n_movies;
    total_nFrames = helper.total_nFrames;
    position = helper.position;

    cropped_height = length(position(2):position(2)+position(4)-1);
    cropped_width = length(position(1):position(1)+position(3)-1);

    % Add more variables
    if ~isfield(helper, 'analyzed'), helper.analyzed = repmat({''}, n_movies, 1); end
    if ~isfield(helper, 'zero_mask'), helper.zero_mask = false([cropped_height, cropped_width, n_movies]); end
    if ~isfield(helper, 'max_mask'), helper.max_mask = []; end
    if ~isfield(helper, 'first_frame_downsampled'), helper.first_frame_downsampled = 1; end
    if ~isfield(helper, 'last_preprocessed_chunk'), helper.last_preprocessed_chunk = 0; end
    if ~isfield(helper, 'last_preprocessed_chunk_mean_all'), helper.last_preprocessed_chunk_mean_all = 0; end
    if ~isfield(helper, 'cropping_coordinates'), helper.cropping_coordinates = []; end
    if ~isfield(helper, 'last_preprocessed_chunk_make_flat'), helper.last_preprocessed_chunk_make_flat = 0; end
            
    % Set beginning and end of each chunk
    movie_frame_edges = [[1; cumsum(helper.total_nFrames(1:end-1)) + 1], cumsum(helper.total_nFrames)];
    
    if p.user.large_data
        % Temporary files are stored locally
        preprocessed_movie_filename = os.path.join(local_saveDir, 'preprocessedMovie.h5');
        temporary_file = os.path.join(local_saveDir, 'temp.h5');
        temporary_flat_file = os.path.join(local_saveDir, 'temp_flat.h5');
        
        if do_create_files || ~exist(temporary_file, 'file')
            if exist(temporary_file, 'file')
                delete(temporary_file)
            end
            h5create(temporary_file, '/1', [cropped_height, cropped_width, sum(total_nFrames)], 'Datatype', 'single');
        end
    end
    
    % loop through avi files and concatenate them
    for m = 1:n_movies
        % load movie
        movPath = os.path.join(movDir, ['msCam', num2str(m), '.avi']);
        
        if ismember(movPath, helper.analyzed)
            continue
        end
        
        LOGGER.info([char(datetime('now')), ' Concatenating movie ' num2str(m), '/', num2str(n_movies)])
        movieObj = VideoReader(movPath);
        nFrames = movieObj.NumberOfFrames;
        thisStack = squeeze(read(movieObj));
                
        % Crop movie
        thisStack = thisStack(position(2):position(2)+position(4)-1, position(1):position(1)+position(3)-1,:);
        
        % For very large datasets, it is worth to do preprocessing of each movie
        if p.user.large_data
            LOGGER.info('\tFiltering movie: dividing by lowpass')
            movie = filterMovie(p, single(thisStack), 'lowpass', false);
            
            % Check whether to delete the first frame
            first_10_frames = movie(:, :, 1:10);
            n_zeros = squeeze(sum(sum(first_10_frames == 0, 1), 2));
            fraction_zeros = n_zeros ./ (size(thisStack, 1) * size(thisStack, 2));
            frames_to_replace = find(fraction_zeros >= GC.epifluorescence_skip_frames_with_zeros);
            if ~isempty(frames_to_replace)
                n_frames_to_replace = length(frames_to_replace);
                first_good_frame = find(fraction_zeros < GC.epifluorescence_skip_frames_with_zeros, 1, 'first');
                movie(:, :, frames_to_replace) = repmat(movie(:, :, first_good_frame), [1, 1, n_frames_to_replace]);
                LOGGER.warn([num2str(n_frames_to_replace) ' frames have been replaced by frame #', num2str(first_good_frame), ' because had too much noise'])
            end

            LOGGER.info('\tRegistering movie')
            [movie, this_zeroMask] = registerMovie(p, movie);
            % Store mask
            helper.zero_mask(:, :, m) = this_zeroMask;
            
            LOGGER.info('\tComputing dFF')
            movie = dfofMovie(movie);
            
            LOGGER.info('\tWriting movie to disk')
            h5write(temporary_file, '/1', movie, [1, 1, movie_frame_edges(m, 1)], [cropped_height, cropped_width, nFrames]);
        end
        
        % Update helper file
        helper.analyzed{m} = movPath;
        save(helper_file, 'helper', '-v7.3')
    end
    % Clear unused variables
    clear thisStack movie
    
    % Complete preprocessing
    if p.user.large_data
        if do_create_files || ~exist(preprocessed_movie_filename, 'file')
            LOGGER.info([char(datetime('now')), ' Creating file ''preprocessedMovie.h5'''])
            if exist(preprocessed_movie_filename, 'file')
                delete(preprocessed_movie_filename)
            end
            h5create(preprocessed_movie_filename, '/1', [cropped_height, cropped_width, Inf], 'Datatype', 'single', 'ChunkSize',[cropped_height, cropped_width, 100]);
        end
        
        % Make final zeroMask
        final_zeroMask = max(helper.zero_mask, [], 3);
        % Set chunk size and number of chunks
        chunk_size = max(total_nFrames);
        n_chunks = ceil(sum(total_nFrames) / chunk_size);
        % Initialize working status, so later we can display appropriate messages
        worked_on_a_chunk = false;
        % Read chunk
        for i_chunk = 1:n_chunks
            if helper.last_preprocessed_chunk >= i_chunk
                continue
            end
            worked_on_a_chunk = true;
            LOGGER.info(['\tReading chunk ', num2str(i_chunk), '/', num2str(n_chunks)])
            % Read movie
            first_frame = (i_chunk - 1) * chunk_size + 1;
            last_frame = (i_chunk - 1) * chunk_size + chunk_size + p.downsampleTime.factor;
            last_frame = min(last_frame, sum(total_nFrames));
            movie = loadMovie(temporary_file, first_frame, last_frame);

            LOGGER.info('\tDownsampling movie')
            movie = downsampleMovie(p, movie, 'time');

            % Remove padding, because last frame was not used for downsampling
            if ndims(movie) == 3  % if there are more than one frame
                movie(:, :, end) = [];
            end
        
            % Crop
            % LOGGER.info([char(datetime('now')) ' cropping movie'])
            % movie = cropMovie(movie, final_zeroMask); 
            
            % Write to disk
            nFrames = size(movie, 3);
            LOGGER.info(['\tSaving movie in frames ', num2str(helper.first_frame_downsampled), ' to ', num2str(helper.first_frame_downsampled + nFrames)])
            h5write(preprocessed_movie_filename, '/1', movie, [1, 1, helper.first_frame_downsampled], [size(movie, 1), size(movie, 2), nFrames]);
            % Update pointer
            helper.first_frame_downsampled = helper.first_frame_downsampled + nFrames + 1;
            
            % Comput max-mask
            helper.max_mask(:, :, i_chunk) = max(movie == 0, [], 3);
            
            % Store at which chunk we are
            helper.last_preprocessed_chunk = i_chunk;
            save(helper_file, 'helper', '-v7.3')
        end
        
        if worked_on_a_chunk
            % Log time that it took for the preprocessing
            time_end_preprocessing = toc;
            LOGGER.info(['Finished preprocessing in ' num2str(time_end_preprocessing / 60, '%.2f') ' min'])
        end
        
        % Delete unused variables
        clear movie
    end
    
    % Get number of frames after downsampling
    info = h5info(preprocessed_movie_filename);
    n_frames_downsampled = info.Datasets.Dataspace.Size(3);
    
    % Start signal extraction
    if p.user.large_data
        % Compute mean across all frames
        zeroMask = max(helper.max_mask, [], 3);
        zeroMask(zeroMask > 0) = 1;
        % [N, S, W, E] = getCropCoords(zeroMask);
        
        n_chunks = ceil(n_frames_downsampled / chunk_size);
        if helper.last_preprocessed_chunk_mean_all == 0
            helper.inputMean = NaN(n_chunks, 2);
        end

        LOGGER.info('\tComputing grand mean')
        for i_chunk = 1:n_chunks
            if helper.last_preprocessed_chunk_mean_all >= i_chunk
                continue
            end
            % Read movie
            LOGGER.info(['\t\tReading chunk ', num2str(i_chunk), '/', num2str(n_chunks)])
            first_frame = (i_chunk - 1) * chunk_size + 1;
            last_frame = (i_chunk - 1) * chunk_size + chunk_size;
            last_frame = min(last_frame, n_frames_downsampled);
            movie = loadMovie(preprocessed_movie_filename, first_frame, last_frame);

            % Crop movie
            % movie = movie(N:S, W:E, :);
            % Store sum and number of pixels in movie (to ompute mean later)
            helper.inputMean(i_chunk, :) = [sum(movie(:)), numel(movie)];

            % Store where we are
            helper.last_preprocessed_chunk_mean_all = i_chunk;
            save(helper_file, 'helper', '-v7.3')
        end

        LOGGER.info([char(datetime('now')) ' Preparing data for signal extraction'])
        % Compute mean of each frame and store flattened version of data
        % cropped_width = length(W:E);
        % cropped_height = length(N:S);
        cropped_width = size(helper.zero_mask, 2);
        cropped_height = size(helper.zero_mask, 1);
        
        % Create file
        if do_create_files || ~exist(temporary_flat_file, 'file')
            if exist(temporary_flat_file, 'file')
                delete(temporary_flat_file)
            end
            h5create(temporary_flat_file, '/1', [cropped_width * cropped_height, n_frames_downsampled], 'Datatype', 'single');
        end
        n_pixels = cropped_width * cropped_height;

        % Compute mean of all movies
        inputMean = sum(helper.inputMean(:, 1)) / sum(helper.inputMean(:, 2));

        % De-mean and flatten movies
        for i_chunk = 1:n_chunks
            if helper.last_preprocessed_chunk_make_flat >= i_chunk
                continue
            end
            first_frame = (i_chunk - 1) * chunk_size + 1;
            last_frame = (i_chunk - 1) * chunk_size + chunk_size;
            last_frame = min(last_frame, sum(n_frames_downsampled));

            LOGGER.info(['\tReading chunk ', num2str(i_chunk), '/', num2str(n_chunks)])
            movie = loadMovie(preprocessed_movie_filename, first_frame, last_frame);
            nFrames = size(movie, 3);
            
            % Crop movie
            % movie = movie(N:S, W:E, :);

            % Subtract grand mean
            movie = movie - inputMean;
            
            % Reshape movie into [space x time] matrix
            movie = reshape(movie, n_pixels, nFrames);

            % Subtract mean of each frame
            mean_M = mean(movie, 1);
            movie = bsxfun(@minus, movie, mean_M);

            % Write data to disk
            h5write(temporary_flat_file, '/1', movie, [1, first_frame], [n_pixels, nFrames]);

            % Update counter
            helper.last_preprocessed_chunk_make_flat = i_chunk;
            save(helper_file, 'helper', '-v7.3')
        end
        
        clear movie
    end
    
    % Signal extraction
    if p.user.large_data
        if ~helper.completed_PCAICA
            saveDir = os.path.join(local_saveDir, 'extracted');
            if ~exist(saveDir, 'dir')
                mkdir(saveDir)
            end
            resultsPCAICA_filename = os.path.join(saveDir, 'resultsPCAICA.mat');
            cellMap_filename = os.path.join(saveDir, 'cellMap.mat');

            % run PCA
            movie = loadMovie(temporary_flat_file);
            LOGGER.info([char(datetime('now')), ' Signal extraction: PCA'])
            [spatial, temporal, S] = compute_pca(movie, p.PCAICA.nPCs);
            S = diag(S);  % keep only the diagonal of S
            clear movie
            
            % ICA
            LOGGER.info([char(datetime('now')), ' Signal extraction: ICA'])
            % set parameters
            mu = p.PCAICA.mu;
            num_ICs = p.PCAICA.nICs;
            term_tol = p.PCAICA.term_tol;
            max_iter = p.PCAICA.max_iter;

            % run ICA
            ica_mixed = compute_spatiotemporal_ica_input(spatial, temporal, mu);
            ica_W = compute_ica_weights(ica_mixed, num_ICs, term_tol, max_iter)';
            [filters, traces] = compute_ica_pairs(spatial, temporal, S, cropped_height, cropped_width, ica_W);
            traces = permute(traces, [2, 1]);

            % put filters back into right frame and save
            % fFrame = zeros([size(zeroMask), size(filters,3)]);
            % Recalculate coordinates
            % [crop_N, crop_S, crop_W, crop_E] = getCropCoords(zeroMask);
            % fFrame(crop_N:crop_S, crop_W:crop_E, :) = filters;
            % filters = fFrame;
            cellMap = squeeze(max(filters, [], 3));

            % Increaes contrast by removing outliers
            q = quantile(cellMap(:), [.01, .99]);
            cellMap(cellMap < q(1)) = q(1);
            cellMap(cellMap > q(2)) = q(2);

            % save data
            save(resultsPCAICA_filename, 'p', 'traces', 'filters');
            save(cellMap_filename, 'cellMap', 'filters');

            helper.completed_PCAICA = true;
            save(helper_file, 'helper', '-v7.3')
        end
    end
    

%% MLint exceptions
%#ok<*TNMLP,*VIDREAD,*NASGU,*STRNU>
