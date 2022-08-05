%% MAIN
function done = main()

% Allocate output variable
done = true;

% Get logger and general_configs
global GC LOGGER

% Find stimulus type
METADATA = table2cell(unique(SQL_database.read_table_where('sessions', 'stimulus')));
n_stimuli = length(METADATA);

% Make sure main folder exists
stimuli_path = os.path.join(GC.data_root_path, GC.stimuli_path);
if ~exist(stimuli_path, 'dir')
    mkdir(stimuli_path)
end

% Loop through stimuli
for istim = 1:n_stimuli
    % Log progress
    LOGGER.trace(['Making ''', METADATA{istim}, ''''])
    
    % Initialize output variable
    STIMULUS = struct('code',        METADATA{istim}, ...
        'description', '', ...
        'timestamps',  []);
    
    % Make a string with description
    switch STIMULUS.code
        case 'FPS'
            description = 'Forepaw electrical stimulation';
        case 'HPS'
            description = 'Hindpaw electrical stimulation';
        case 'HPS_50%'
            description = 'Hindpaw electrical stimulation (50% intensity)';
        case 'HPS_short'
            description = 'Hindpaw electrical stimulation (shortly after another HPS block)';
        case 'HPSn'
            description = 'Hindpaw electrical stimulation (n: new frame rate, new ScanImage)';
        case 'HPS_2n'
            description = 'Hindpaw electrical stimulation end of session (n: new frame rate, new ScanImage)';
        case 'HPSnl'
            description = 'Hindpaw electrical stimulation (n: new frame rate, l: long duration of recording)';
        case 'SP'
            description = 'Spontaneous activity';
        case 'SPn'
            description = 'Spontaneous activity (n: new frame rate, new ScanImage) ';
        case 'SPnl'
            description = 'Spontaneous activity (n: new frame rate, l: long duration of recording)';
        case 'SP_2'
            description = 'Spontaneous activity (end of session)';
        case 'SP_2n'
            description = 'Spontaneous activity (end of session), new ScanImage';
        case 'acet+'
            description = 'Acetone application to hindpaw';
        case 'acet-'
            description = 'Saline application to hindpaw (negative control for acet+)';
        case 'odor+'
            description = 'Aversive odor (TMT) of fox urine';
        case 'odor-'
            description = 'Normal air flow (negative control for odor+)';
        case 'puff'
            description = 'Aversive air puff';
        case 'puffn'
            description = 'Aversive air puff (n: new frame rate(new ScanImage))';
        case 'puffnl'
            description = 'Aversive air puff (n: new frame rate(new ScanImage), l: long duration of recording)';
        case 'reward+'
            description = 'Nicotine injection i.p.';
        case 'reward-'
            description = 'Saline injection i.p. (negative control for reward+)';
        case 'sound'
            description = 'Neutral loud sound (85 dB)';
        case 'soundn'
            description = 'Neutral loud sound (85 dB) (n: new frame rate(new ScanImage))';
        case 'soundnl'
            description = 'Neutral loud sound (85 dB) (n: new frame rate(new ScanImage), l: long duration of recording)';
        case 'temp_38'
            description = 'Neutral 38 deg ramping temperature';
        case 'temp_42'
            description = 'Possibly nociceptive 42 deg ramping temperature';
        case 'temp_43'
            description = 'Possibly nociceptive 43 deg ramping temperature';
        case 'temp_48'
            description = 'Noxious 48 deg ramping temperature';
        case 'touch'
            description = 'Mechanoceptive stimulus';
        case 'pinprick'
            description = 'Pinprick stimulus';
        case 'EPM'
            description = 'Elevated plus maze';
        case 'cold'
            description = 'Noxious dry ice application';
        case 'heat'
            description = 'Noxious laser heat application';
        case 'heatn'
            description = 'Noxious laser heat application under the 2p,(n: new frame rate(new ScanImage)) ';
        case 'heatnl'
            description = 'Noxious laser heat application under the 2p,(n: new frame rate(new ScanImage), l: long duration of recording) )';
        case 'TMT'
            description = 'Aversive TMT air puff not directly to the face';
        case 'water'
            description = 'Neutral water air puff not directly to the face';

        otherwise
            error('update_timestamps:unknown_stimulus', '%s is an unknown stimulus type', STIMULUS.code)
    end
    STIMULUS.description = description;
    
    % Make time axis, timestamps and stimulus profile
    switch STIMULUS.code
        case {'FPS', 'HPS', 'HPS_50%', 'HPS_short'}
            STIMULUS.timestamps = 10;  % Middle of the trial
            STIMULUS.duration = 0.05;  % 50 ms
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case {'SP', 'SP_2'}
            STIMULUS.timestamps = [];
            STIMULUS.duration = [];
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            % Add info from epi data
            row = height(info) + 1;
            info{row, 'frame_rate'} = GC.epifluorescence_downsample_to_frame_rate;
            info{row, 'n_frames'} = NaN;
            info = sortrows(table2array(info));
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                if ~isnan(n_frames)
                    time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                    
                    % Make stimulus profile
                    stimulus_profile = zeros(n_frames, 1);
                else
                    time_axis = NaN;
                    stimulus_profile = NaN;
                end
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
        case {'SPn', 'SP_2n', 'SPnl'}
            STIMULUS.timestamps = [];
            STIMULUS.duration = [];
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            % Add info from epi data
            %row = height(info) + 1;
            %info{row, 'frame_rate'} = GC.epifluorescence_downsample_to_frame_rate;
            %info{row, 'n_frames'} = NaN;
            info = sortrows(table2array(info));
%             if length(info) > 2, keyboard, end
            
            if any(ismember(info(:), 361))
%                 continue
                info(any(ismember(info, 361),2),:) = [];
            end
            
            if any(ismember(info(:), 180))
%                 continue
               info(any(ismember(info, 361),2),:) = [];
            end
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                if ~isnan(n_frames)
                    time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                    
                    % Make stimulus profile
                    stimulus_profile = zeros(n_frames, 1);
                else
                    time_axis = NaN;
                    stimulus_profile = NaN;
                end
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case {'acet+', 'acet-'}
            STIMULUS.timestamps = 0;
            STIMULUS.duration = 60;
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                try
                    frame_rate = info(i_comb, 1);
                catch, keyboard, end
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Make stimulus profile
                if strcmp(STIMULUS.code, 'acet+')
                    stimulus_profile = ones(n_frames, 1);
                else
                    stimulus_profile = zeros(n_frames, 1);
                end
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case {'odor+', 'odor-'}
            STIMULUS.timestamps = 10;
            STIMULUS.duration = 5;
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case {'puff', 'sound'}
            STIMULUS.timestamps = 10;
            STIMULUS.duration = 0.5;
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            try % for sound there's no epi stimulus yet. 
                info = [info; [GC.epifluorescence_downsample_to_frame_rate, GC.epifluorescence_downsample_to_frame_rate * (GC.miniscope_trial_duration.(STIMULUS.code)(2) + GC.miniscope_trial_duration.(STIMULUS.code)(1)) + 1;]];
            catch
                continue
            end
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                % Get indices of timestamps in time axis
                if n_frames ~= 76
                    timestamps = 10;
                elseif n_frames == 76
                    timestamps = 5;
                end
                [~, timestamp_idx] = min(abs(time_axis - timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
                STIMULUS.stimulus_profile{i_comb, 5} = timestamps;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile', 'timestamps'});
            STIMULUS.stimulus_profile = sortrows(STIMULUS.stimulus_profile, {'frame_rate', 'n_frames'});
        case {'puffn', 'soundn', 'puffnl', 'soundnl', 'heatn', 'heatnl'}
            STIMULUS.timestamps = 10;
            STIMULUS.duration = 0.5;
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
%             if length(info) > 2, keyboard, end
            if any(ismember(info(:), 361))
                info(any(ismember(info, 361),2),:) = [];
            end
            if any(ismember(info(:), 180))
%                 keyboard
%                 continue
                info(any(ismember(info, 180),2),:) = [];
            end
            
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
        case {'reward+', 'reward-'}
            STIMULUS.timestamps = [];
            STIMULUS.duration = [];
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case {'temp_38', 'temp_42', 'temp_43', 'temp_48'}
            speed = 1;  % 1 degree C / s
            total_duration_ramp_plateau = 22;  % s
            baseline_temperature = 32;  % degree C
            max_temperature = strsplit(STIMULUS.code, '_');
            max_temperature = str2double(max_temperature{2});
            stimulus_onset = 10;  % s
            STIMULUS.timestamps = [stimulus_onset, ...  % beginning of ramp up
                stimulus_onset + (max_temperature - baseline_temperature) / speed, ...  % beginning of plateau temperature
                stimulus_onset + total_duration_ramp_plateau, ...  % beginning of ramp down
                stimulus_onset + total_duration_ramp_plateau + (max_temperature - baseline_temperature) / speed, ...  % end of stimulus
                ];
            STIMULUS.duration = STIMULUS.timestamps(4) - STIMULUS.timestamps(1);
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1) + baseline_temperature;
                % Ramp up
                ramp_up = linspace(baseline_temperature, max_temperature, timestamp_idx(2) - timestamp_idx(1) + 1);
                stimulus_profile(timestamp_idx(1):timestamp_idx(2)-1) = ramp_up(1:end-1);
                % Plateau
                stimulus_profile(timestamp_idx(2):timestamp_idx(3)-1) = max_temperature;
                % Ramp down
                ramp_down = linspace(max_temperature, baseline_temperature, timestamp_idx(4) - timestamp_idx(3) + 1);
                stimulus_profile(timestamp_idx(3):timestamp_idx(4)-1) = ramp_down(1:end-1);
                % Normalize between 0 and 1
                stimulus_profile = stimulus_profile - min(stimulus_profile);
                stimulus_profile = stimulus_profile  / max(stimulus_profile);
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case {'touch', 'pinprick'}
            STIMULUS.timestamps = 10;  %10; Middle of the trial
            STIMULUS.duration = 0.5;  % 500 ms
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
            info = [info; [GC.epifluorescence_downsample_to_frame_rate, GC.epifluorescence_downsample_to_frame_rate * (GC.miniscope_trial_duration.(STIMULUS.code)(2) + GC.miniscope_trial_duration.(STIMULUS.code)(1)) + 1;]];
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 5);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                % Get indices of timestamps in time axis
                if n_frames == 170
                    timestamps = 10;
                elseif n_frames == 76
                    timestamps = 5;
                end
                [~, timestamp_idx] = min(abs(time_axis - timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
                STIMULUS.stimulus_profile{i_comb, 5} = timestamps;
            end

            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile', 'timestamps'});
            STIMULUS.stimulus_profile = sortrows(STIMULUS.stimulus_profile, {'frame_rate', 'n_frames'});
            
        case 'EPM'
            STIMULUS.timestamps = [];
            STIMULUS.duration = [];
            % Get number of frames and frame rates for this stimulus
            info = GC.epifluorescence_downsample_to_frame_rate;
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = NaN;
                STIMULUS.stimulus_profile{i_comb, 3} = NaN;
                STIMULUS.stimulus_profile{i_comb, 4} = NaN;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case 'cold'
            STIMULUS.timestamps = GC.miniscope_trial_duration.cold(1);
            STIMULUS.duration = 0;
            % Get number of frames and frame rates for this stimulus
            info = [GC.epifluorescence_downsample_to_frame_rate, ...
                GC.epifluorescence_downsample_to_frame_rate * (GC.miniscope_trial_duration.cold(2) + GC.miniscope_trial_duration.cold(1)) + 1];
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
            
        case 'heat'
            STIMULUS.timestamps = GC.miniscope_trial_duration.heat(1);
            STIMULUS.duration = 0;
            % Get number of frames and frame rates for this stimulus
            info = [GC.epifluorescence_downsample_to_frame_rate, ...
                GC.epifluorescence_downsample_to_frame_rate * (GC.miniscope_trial_duration.heat(2) + GC.miniscope_trial_duration.heat(1)) + 1];
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});
        case {'FPSn', 'HPSn', 'HPS_2n', 'FPSnl', 'HPSnl'}
            STIMULUS.timestamps = 10;  % Middle of the trial
            STIMULUS.duration = 0.05;  % 50 ms
            % Get number of frames and frame rates for this stimulus
            info = unique(SQL_database.read_table_where('trials', {'frame_rate','n_frames'}, STIMULUS.code, 'stimulus'), 'rows');
            info = sortrows(table2array(info));
%             if length(info) > 2, keyboard, end
            if any(ismember(info(:), 361))
                
                info(any(ismember(info, 361),2),:) = [];
            end
            if any(ismember(info(:), 180))
                
                info(any(ismember(info, 180),2),:) = [];
            end
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});

        case {'TMT', 'water'}
            STIMULUS.timestamps = GC.miniscope_trial_duration.TMT(1);
            STIMULUS.duration = 0;
            % Get number of frames and frame rates for this stimulus
            info = [GC.epifluorescence_downsample_to_frame_rate, ...
                GC.epifluorescence_downsample_to_frame_rate * (GC.miniscope_trial_duration.TMT(2) + GC.miniscope_trial_duration.TMT(1)) + 1];
            n_combinations = size(info, 1);
            STIMULUS.stimulus_profile = cell(n_combinations, 4);
            for i_comb = 1:n_combinations
                frame_rate = info(i_comb, 1);
                n_frames = info(i_comb, 2);
                time_axis = linspace(0, n_frames / frame_rate, n_frames)';
                
                % Get indices of timestamps in time axis
                [~, timestamp_idx] = min(abs(time_axis - STIMULUS.timestamps));
                % Make stimulus profile
                stimulus_profile = zeros(n_frames, 1);
                stimulus_profile(timestamp_idx:timestamp_idx + round(STIMULUS.duration * frame_rate)) = 1;
                
                % Store data
                STIMULUS.stimulus_profile{i_comb, 1} = frame_rate;
                STIMULUS.stimulus_profile{i_comb, 2} = n_frames;
                STIMULUS.stimulus_profile{i_comb, 3} = time_axis;
                STIMULUS.stimulus_profile{i_comb, 4} = stimulus_profile;
            end
            STIMULUS.stimulus_profile = cell2table(STIMULUS.stimulus_profile, 'VariableNames',{'frame_rate','n_frames','time_axis','stimulus_profile'});

    end
    
    % Make name of output file
    stimulus_filename = os.path.join(stimuli_path, [STIMULUS.code, '.mat']);
    save(stimulus_filename, 'STIMULUS', '-v7.3')
end

