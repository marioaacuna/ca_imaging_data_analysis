function ca_trace_extraction(animal_ID, neuron_exist)
    % This function extract calcium traces from images.
    
    % Get logger and general_configs
    global LOGGER GC

    % Make filename for output
    dFF_file = get_filename_of('dFF', animal_ID);
    % Get filenames of input
    ROI_info_file = get_filename_of('ROI_info', animal_ID);
    % Get data type
    data_type = SQL_database.read_table_where('experiments', 'data_type', animal_ID,'animal_ID', 'return_as_table',false);
    
    % Read metadata from database
    switch data_type
        case '2p'
            METADATA = SQL_database.read_table_where('trials', {'+','date','experimental_condition','stimulus'}, animal_ID,'animal_ID');
            % Get list of indices where each trial starts
             if any(strcmp(animal_ID, {'NRN_15a', 'NRN_15b' ,'NRN_16a', 'NRN_16b'}))
                 trials_frames_idx = get_trial_frames(animal_ID, 'reset_index',false, 'return_all_sessions',1);
             else 
                 trials_frames_idx = get_trial_frames(animal_ID, 'reset_index',false, 'return_all_sessions',false);
             end

        case 'epi'
            try
                METADATA = SQL_database.read_epi_trials(animal_ID);
            catch
                METADATA = SQL_database.read_table_where('sessions', {'+'}, animal_ID, 'animal_ID');
                disp('Check, this animal does not have any stimulation (yet)')
            end

            miniscope_parameters = load_variable(get_filename_of('miniscope_movie_parameters', animal_ID), 'p');
            sessions_order = miniscope_parameters.user.experiments;
            group_sessions_by = {'date', 'experiment'};
            sessions = uniqueCellRows(METADATA{:, group_sessions_by}, 'rows');
            preprocessed_frame_edges = miniscope_parameters.user.jointExtraction.preprocessed_frame_edges;
            
            % Make time axis for each session
            frame_rate = GC.epifluorescence_downsample_to_frame_rate;
            time_axis = cell(size(preprocessed_frame_edges, 1), 1);
            for i_sess = 1:length(sessions)% 1:miniscope_parameters.nSessions
                n_frames = preprocessed_frame_edges(i_sess, 2) - preprocessed_frame_edges(i_sess, 1) + 1;
                duration = n_frames / frame_rate;
                time_axis{i_sess} = linspace(0, duration, n_frames);
            end
    end
    n_trials = height(METADATA);
    
    [~, ~, session_idx] = unique(METADATA.session_id, 'stable');
    n_sessions = length(unique(session_idx));
    % Read the Ca2+ indicator used for this animal and get its decay time constant
    Ca_indicator = SQL_database.read_table_where('experiments', 'Calcium_indicator', animal_ID, 'animal_ID', 'return_as_table',false);
    tau_decay = GC.decay_time_constant_calcium_reporter.(Ca_indicator);

    % Load data
    LOGGER.trace('Loading data')
    F = load_variable(ROI_info_file, 'ROI_fluorescence');
    % Get number of ROIs and samples
    n_ROIs = size(F, 2);
    n_samples = size(F, 1);
         
    %% Compute dF/F0 for each ROI
    LOGGER.trace(['Computing ', char(916), 'F/F0'])

    switch data_type
        case '2p'
            % Get frame rate
            frame_rate = METADATA.frame_rate(1);
            if ~neuron_exist
                % Make running-window
                % Window is the maximum of 15s or 40x time decay
                win_tau = 40 * tau_decay / (1000/frame_rate);
                win_time = 15 * frame_rate;
                win = max([win_tau, win_time]);
                win = floor(win/2)*2 + 1;  % Round to nearest odd number
                % Signal padders
                left_pad = abs(ceil(-win/2));
                right_pad = floor(win/2);
            end
        case 'epi'
    end
    
    % Initialize output variables
    dFF = cell(n_ROIs, n_trials);

    for iroi = 1:n_ROIs
        LOGGER.trace(['Processed ', num2str(iroi), ' ROIs (out of ', num2str(n_ROIs), ')'], 'overwrite',iroi>1)

        % Get cell signal
        F_cell = double(F(:, iroi));
        % Replace NaNs with 0s
        F_cell(isnan(F_cell)) = 0;
        
        switch data_type
            case '2p'
                % Initialize
                dFF_cell = NaN(n_samples, 1);
                if ~neuron_exist
                    % Running-window average of the 8th percentile of the ROI's fluorescence, as in Romano et al (DOI:10.1371/journal.pcbi.1005526)
                    % Split data by session
                    for isess = 1:n_sessions
                        frames_this_session = trials_frames_idx(session_idx == isess, :);
                        % Make mask of frames to select
                        mask = false(n_samples, 1);
                        for i_trial = 1:size(frames_this_session, 1)
                            mask(frames_this_session(i_trial,1):frames_this_session(i_trial,2)) = true;
                        end
                        F_session = F(mask, iroi);
                        % Get indices of each window
                        all_windows = repmat(1:length(F_session),win,1).';
                        all_windows = bsxfun(@plus, all_windows, 0:win-1);
                        % Remove lines with out-of-bounds indices
                        length_padded_F = length(F_session) + left_pad + right_pad;
                        all_windows(any(all_windows.' > length_padded_F), :) = [];
                        
                        % Pad signal
                        temp_F = [NaN(left_pad,1); F_session; NaN(right_pad,1)];
                        % Get values in windows
                        temp_F = temp_F(all_windows);
                        % Compute 8th percentile in each window
                        baseline = prctile(temp_F, 8, 2);
                        F_smooth = runline(baseline, win, 1);
                        
                        % Compute (F-F0)/F0
                        dFF_cell(mask) = (F_session - F_smooth) ./ F_smooth;
                        
                        % ### Visualize F ###
                        % figure('color','w'), clf, hold on, plot(F_session), plot(F_smooth,'linewidth',3)
                        % figure('color','w'), clf, hold on, plot(dFF_cell(mask))
                        % ##########################
                    end
                    
                    % Replace non-finite values with 0s [NaN is the result of 0/0, Inf is the result of anything divided by 0]
                    dFF_cell(~isfinite(dFF_cell)) = 0;
                else
                    dFF_cell = F(1:n_samples,iroi);
                end
                % Split data by trial
                for i_trial = 1:n_trials
                    % Get indices of trial start and end
                    trial_start = trials_frames_idx(i_trial,1);
                    trial_end   = trials_frames_idx(i_trial,2);
                    % Get dFF segment for this trial
                    dFF_trial = dFF_cell(trial_start:trial_end);
                    % Store data in memory as a row vector
                    dFF{iroi, i_trial} = dFF_trial(:)';
                end

            case 'epi'
                has_ts =  any(strcmp('timestamps',METADATA.Properties.VariableNames));
                % Split data by trial
                for i_trial = 1:n_trials
                    % Get to which condition this trial belongs
                    info_this_trial = METADATA{i_trial, {'date', 'experiment'}};
                    stimulus_this_trial = METADATA{i_trial, 'stimulus'}; stimulus_this_trial = stimulus_this_trial{1};
                    session_idx = find(ismemberCellRows(sessions_order, info_this_trial));
                    frames_this_session = preprocessed_frame_edges(session_idx, :);
                    n_frames_this_session = frames_this_session(2) - frames_this_session(1) + 1;
                    % Get F of this session
                    F_this_session = F_cell(frames_this_session(1):frames_this_session(2));
                    
                    
                    if has_ts
                        timestamps = METADATA.timestamps(i_trial,:);
                        if all(isnan(timestamps))
                            dFF_trial = F_this_session;
                        else
                            timestamps_sample = NaN(1, 2);
                            for i_ts = 1:2
                                if isnan(timestamps(i_ts)), continue, end
                                [~, timestamps_sample(i_ts)] = min(abs(timestamps(i_ts) - time_axis{session_idx}));
                            end
                            earliest_timestamp = min(timestamps_sample(~isnan(timestamps_sample)));
                            latest_timestamp = max(timestamps_sample(~isnan(timestamps_sample)));
                            baseline_window = ceil(GC.miniscope_trial_duration.(stimulus_this_trial)(1) * frame_rate);
                            evoked_window = ceil(GC.miniscope_trial_duration.(stimulus_this_trial)(2) * frame_rate);
                            n_frames_this_trial = baseline_window + evoked_window + 1;
                            dFF_trial = NaN(n_frames_this_trial, 1);

                            % Make sure that trial selects existing frames
                            frames_idx = earliest_timestamp - baseline_window:latest_timestamp + evoked_window;
                            keep_frames = ~(frames_idx < 1 | frames_idx > n_frames_this_session);
                            frames_idx = frames_idx(keep_frames);
                            % Get dFF segment for this trial
                            dFF_trial(keep_frames) = F_this_session(frames_idx);
                        end
                    else
                        dFF_trial = F_this_session;
                    end
                    
                    % Store data in memory as a row vector
                    dFF{iroi, i_trial} = dFF_trial(:)';
                end
        end
    end
    
    % Write data to disk
    LOGGER.trace('Writing Ca2+ traces to disk ...')
    save(dFF_file, 'dFF', '-v7.3')
    LOGGER.trace('done', 'append',true)
end
%% Mlint exceptions
%#ok<*NASGU,*PFBNS,*PFOUS>
