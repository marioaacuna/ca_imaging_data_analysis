function organizeData(p, sessionDirs, experiments, missing_only, done_preprocessing) 
%ORGANIZEDATA organizes movie files for all sessions of one subject

    global GC LOGGER

    p_filename = get_filename_of('miniscope_movie_parameters', p.subjects{1});
    animal_ID = p.subjects{1};
    
    %% STEP 2: get refFrames from all sessions for choosing ROIs
    
    if isfield(p.user, 'position')
        fprintf('ROI positions already stored for these sessions. Would you like to overwrite, or add more? y/[n]\n')
        while true
            reply = input('','s');
            if ismember(reply, {'y', 'n'}) || isempty(reply)
                break
            else  % delete invalid input
                fprintf(repmat('\b', 1, length(reply)+1))
            end
        end
        redo_positions = strcmp(reply, 'y');
    
    else
        redo_positions = true;
    end
    
    if redo_positions
        refFrames = cell(p.nSessions, 1);
        for i_sess = 1:p.nSessions
            % load first movie in the session directory
            movPath = os.path.join(sessionDirs{i_sess}, 'msCam1.avi');
            movieObj = VideoReader(movPath);
            %movieObj = imread(movPath);
            % take frame number 100 as reference frame
            refFrames{i_sess} = read(movieObj,100);
        end

        % display the referance frames, choose ROIs such that they are
        % approximately the same for each session
        ROIsize = [550, 450];
        
        f = figure('units','normalized','outerposition',[0 0 1 1]);
        n_subplot_cols = 3;
        n_subplot_rows = ceil(p.nSessions / n_subplot_cols);
        for i_sess = 1:p.nSessions
            if done_preprocessing(i_sess)
                subaxis(n_subplot_rows, n_subplot_cols, i_sess, 'ml',.03, 'mr',.01, 'mt',.1, 'mb',.03, 'sv',.05, 'sh',.02);
                imagesc(refFrames{i_sess}); axis image; colormap gray;
                %imagesc(movieObj); axis image; colormap gray;
                title(['Session ' num2str(i_sess)])
                try
                    if ~isempty(p.user.position{i_sess})
                        h{i_sess} = imrect(gca,[p.user.position{i_sess}]);
                    else
                        h{i_sess} = imrect(gca,[10,10,ROIsize(1),ROIsize(2)]);
                    end
                catch
                    h{i_sess} = imrect(gca,[10,10,ROIsize(1),ROIsize(2)]);
                end
            else
                subaxis(n_subplot_rows, n_subplot_cols, i_sess, 'ml',.03, 'mr',.01, 'mt',.1, 'mb',.03, 'sv',.05, 'sh',.02);
                imagesc(refFrames{i_sess}); axis image; colormap gray;
                %imagesc(movieObj); axis image; colormap gray;
                title(['Session ' num2str(i_sess)])
                h{i_sess} = imrect(gca,[10,10,ROIsize(1),ROIsize(2)]);
            end
        end
        sgtitle({'Position the rectangles to select the best image region' ; 'Choose approx. the same area in every session' ; 'Confirm with double click on rectangles (left to right)'})
        
        % Store ROI position
        p.user.position = cell(p.nSessions, 1);
        for i_sess = 1:p.nSessions
            position = wait(h{i_sess});
            p.user.position{i_sess} = round(position);
            % Close figure
            if i_sess == p.nSessions
                close(f)
            end
        end

        % Store p to disk
        save(p_filename, 'p')
    end
    
    %% STEP 3: for every session: load and concat movies, subtract dead pixel values, crop and save

    % Get folder on server where to store data
    root_dir_on_server = os.path.join(GC.registered_images, p.subjects{1});
    if ~exist(root_dir_on_server, 'dir'), mkdir(root_dir_on_server), end
    
    for i_sess = 1:p.nSessions
        % Check whether this session has been already analyzed
        this_experiment = experiments(i_sess, :);
        GC.experiment_name = 'ACC_SNI_anxiety'; % for some reason the experiment name was changed to the date. FIX
        METADATA = SQL_database.read_table_where('sessions', {'done_preprocessing', 'experimental_condition'}, [p.subjects(1), this_experiment{1}], {'animal_ID', 'date'});
        row = find(ismember(METADATA.experiment, this_experiment{2}));
        experimental_condition = unique(METADATA.experimental_condition(row));
        done_preprocessing = METADATA.done_preprocessing(row);
        % Because different stimuli might be associated to the same session,
        % just take the lowest value of all, which means that if by any chance
        % there are both 0s and 1s, it means that it is not complete.
        done_preprocessing = min(logical(done_preprocessing));
        if done_preprocessing && missing_only, continue, end
        
        LOGGER.info([char(datetime('now')) ': Processing ', p.subjects{1}, ' - session ' num2str(i_sess), ' (', this_experiment{1}, ' - ', this_experiment{2}, ')'])

        % Get cropping position of ROI
        position = p.user.position{i_sess};
        movDir = sessionDirs{i_sess};
        % Get path of file
        raw_movie_folder = os.path.join(GC.temp_dir, p.subjects{1}, ['exp', num2str(i_sess)], 'raw');
        if ~exist(raw_movie_folder, 'dir')
            mkdir(raw_movie_folder)
        end
        rawMovie_filename = os.path.join(raw_movie_folder, 'rawMovie.h5');

        % Copy filtering parameters in the right place
        if (p.nSessions > 1) % in case it only has one session, skip it
            p.turboreg.bandpassFreqs = p.user.turboreg.bandpassFreqs{i_sess};
            p.filtering.lowpassFreq = p.user.filtering.lowpassFreq{i_sess};
        end

        % Compute downsampling factor
        timestamp_filename = os.path.join(sessionDirs{i_sess}, 'timestamp.dat');
        timestamp = readtable(timestamp_filename);
        timestamp_interval = diff(timestamp.sysClock);
        frame_rate_this_experiment = 1000 / median(timestamp_interval(2:end));
        p.user.frameRate(i_sess) = frame_rate_this_experiment;
        p.downsampleTime.factor = round(frame_rate_this_experiment / GC.epifluorescence_downsample_to_frame_rate);
        % Store p to disk
        save(p_filename, 'p')
        
        % Run concatenation and analysis
        concatMovies(p, movDir, position, rawMovie_filename);
        preprocessed_movie_filename_on_server = get_filename_of('miniscope_preprocessed_movie', animal_ID, this_experiment{1}, this_experiment{2});
        preprocessed_movie_filename_local = os.path.join(GC.temp_dir, p.subjects{1}, ['exp', num2str(i_sess)], 'preprocessedMovie.h5');
        preprocessed_movie_filename_local_int16 = os.path.join(GC.temp_dir, p.subjects{1}, ['exp', num2str(i_sess)], 'preprocessedMovie_int16.h5');
        % Convert data to a lower precision to save space on disk
        convert_hdf5_to(preprocessed_movie_filename_local, preprocessed_movie_filename_local_int16, 'int16', 'scale',true);

        % Copy data to server
        LOGGER.info([char(datetime('now')) ' Copying preprocessed movie to the server'])
        copyfile(preprocessed_movie_filename_local_int16, preprocessed_movie_filename_on_server)
        
        LOGGER.info([char(datetime('now')) ' Copying cell map to the server'])
        cell_map_filename_on_server = get_filename_of('miniscope_cell_map', animal_ID, this_experiment{1}, this_experiment{2});
        cell_map_filename = os.path.join(GC.temp_dir, p.subjects{1}, ['exp', num2str(i_sess)], 'extracted', 'cellMap.mat');
        resultsPCAICA_filename = os.path.join(GC.temp_dir, p.subjects{1}, ['exp', num2str(i_sess)], 'extracted', 'resultsPCAICA.mat');
        copyfile(cell_map_filename, cell_map_filename_on_server)
        
        % Delete local files
        local_folder_to_delete = os.path.join(GC.temp_dir, p.subjects{1}, ['exp', num2str(i_sess)]);
        if exist(local_folder_to_delete, 'dir'), rmdir(local_folder_to_delete, 's'), end
        
        % Update database
        session_id = SQL_database.read_table_where('sessions', 'session_id', [p.subjects(1), this_experiment{1}, experimental_condition], {'animal_ID', 'date', 'experimental_condition'});
        SQL_database.update('sessions', 'done_preprocessing', {true}, session_id.session_id, 'session_id');
    end
end

%#ok<*NASGU,*TNMLP>
