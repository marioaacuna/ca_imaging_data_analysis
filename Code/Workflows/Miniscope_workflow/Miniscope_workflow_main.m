function Miniscope_workflow_main(animal_ID, missing_only)

% Import global variables with parameters and for logging
global GC LOGGER

toggle_toolbox('_plotting', 'on')

%% Create parameter structure 

% Get info from database
METADATA = SQL_database.read_table_where('sessions', {}, animal_ID,'animal_ID');
experiments = uniqueCellRows(METADATA{:, {'date', 'experiment'}}, 'rows');
n_experiments = size(experiments, 1);
done_preprocessing = unique(SQL_database.read_table_where('sessions', {'experimental_condition','done_preprocessing'}, animal_ID, 'animal_ID', 'split_experimental_condition', false), 'stable');
done_preprocessing = done_preprocessing.done_preprocessing;

% Load existing p structure
p_dir = os.path.join(GC.data_root_path, GC.movie_raw_path, animal_ID);
if ~exist(p_dir, 'dir')
    mkdir(p_dir)
end
p_filename = get_filename_of('miniscope_movie_parameters', animal_ID);
if exist(p_filename, 'file')
    p = load_variable(p_filename, 'p');
    
else  % Make a new structure
    p = params();
    p.nSessions = n_experiments;
    p.experiment = {''};
    p.subjects = {animal_ID};
    p.user.large_data = true;
    p.user.turboreg = struct();
    p.user.frameRate = NaN(n_experiments, 1);
end
p.nSessions = n_experiments;
p.user.experiments = experiments;
p.PCAICA.nPCs = GC.epifluorescence_PCAICA_nPCs;
p.PCAICA.nICs = GC.epifluorescence_PCAICA_nICs;
save(p_filename, 'p');

%% Set filtering parameters per session
if ~isfield(p.user.turboreg, 'bandpassFreqs') || any(cellfun(@(x) isempty(x), p.user.turboreg.bandpassFreqs)) || ...
        ~isfield(p.user.filtering, 'lowpassFreq') || any(cellfun(@(x) isempty(x), p.user.filtering.lowpassFreq)) ||...
        length(p.user.filtering.lowpassFreq) < length (done_preprocessing)
    if ~isfield(p.user.turboreg, 'bandpassFreqs') || any(cellfun(@(x) isempty(x), p.user.turboreg.bandpassFreqs)) || ...
            ~isfield(p.user.filtering, 'lowpassFreq') || any(cellfun(@(x) isempty(x), p.user.filtering.lowpassFreq))
        
        p.user.turboreg.bandpassFreqs = cell(p.nSessions, 1);
        p.user.filtering.lowpassFreq = cell(p.nSessions, 1);
        for i_sess = 1:p.nSessions
            % Read reference frame from this session
            movPath = get_filename_of('miniscope_movie_avi', animal_ID, experiments{i_sess, 1}, experiments{i_sess, 2}, '1');
            movieObj = VideoReader(movPath);
            nFrames = movieObj.NumberOfFrames;
            frame_to_read = 100;
            if frame_to_read > nFrames
                frame_to_read = nFrames;
            end
            
            LOGGER.info(['Used frame #', num2str(frame_to_read), ' from file ''', movPath, ''' as refFrame'], 'contains_path',true)
            refFrame = read(movieObj, frame_to_read);
            
            % test different parameter combinations for lowpass and bandpass filtering on this frame
            filterParams = testFilters(refFrame);
            % set filter parameters
            p.user.turboreg.bandpassFreqs{i_sess} = filterParams(2:3);
            p.user.filtering.lowpassFreq{i_sess} = filterParams(1);
        end
        
        % set the values of the filter parameters to be equal
        % Band pass filter
        filter_values_bp = cell2mat(p.user.turboreg.bandpassFreqs);
        ref_values_bp = round(mean(filter_values_bp));
        % Low pass filter
        filter_values_lp = cell2mat(p.user.filtering.lowpassFreq);
        ref_values_lp = round(mean(filter_values_lp));
        
        % Loop through sessions
        for i_sess = 1:p.nSessions
            p.user.turboreg.bandpassFreqs{i_sess} = ref_values_bp;
            p.user.filtering.lowpassFreq{i_sess} = ref_values_lp;
        end
        
        save(p_filename, 'p');
        
    elseif length(p.user.filtering.lowpassFreq) < length (done_preprocessing)
        for i_sess = 1:p.nSessions
            if ~done_preprocessing(i_sess)
                % Read reference frame from this session
                movPath = get_filename_of('miniscope_movie_avi', animal_ID, experiments{i_sess, 1}, experiments{i_sess, 2}, '1');
                movieObj = VideoReader(movPath);
                nFrames = movieObj.NumberOfFrames;
                frame_to_read = 100;
                if frame_to_read > nFrames
                    frame_to_read = nFrames;
                end
                
                LOGGER.info(['Used frame #', num2str(frame_to_read), ' from file ''', movPath, ''' as refFrame'], 'contains_path',true)
                refFrame = read(movieObj, frame_to_read);
                
                % test different parameter combinations for lowpass and bandpass filtering on this frame
                filterParams = testFilters(refFrame);
                % set filter parameters
                p.user.turboreg.bandpassFreqs{i_sess} = filterParams(2:3);
                p.user.filtering.lowpassFreq{i_sess} = filterParams(1);
            else
                continue
            end
        end
        % set the values of the filter parameters to be equal
        % Band pass filter
        ref_values_bp = p.user.turboreg.bandpassFreqs{1};
        ref_values_lp = p.user.filtering.lowpassFreq{1};
        for i_sess = 1:p.nSessions 
            if ~done_preprocessing(i_sess) 
            p.user.turboreg.bandpassFreqs{i_sess} = ref_values_bp;
            p.user.filtering.lowpassFreq{i_sess} = ref_values_lp;
            else
                continue
            end
        end  
        save(p_filename, 'p');
    end
else
    LOGGER.info('Filters are already set for all experiments')
end


%% Concatenate and extract signals from big data sets
sessionDirs = cell(p.nSessions, 1);
for i_sess = 1:p.nSessions
    movPath = get_filename_of('miniscope_movie_avi', animal_ID, experiments{i_sess, 1}, experiments{i_sess, 2}, '1');
    sessionDirs{i_sess} = fileparts(movPath);
end

organizeData(p, sessionDirs, experiments, missing_only, done_preprocessing) 


toggle_toolbox('_plotting', 'off')


%% MLint exceptions
%#ok<*TNMLP,*VIDREAD>
