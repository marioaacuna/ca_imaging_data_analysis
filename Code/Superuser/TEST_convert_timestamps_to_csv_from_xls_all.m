%% Preamble
% This script is for converting many of the xlsx files that came from
% behavior data collected not with our GUI.
%% Set variables
clear, clc
global GC
table_filepath = 'T:\Mario\Behavior\ACC_pain\behavior cams\_stimuli_timestaps_csv_to_database\originals';
csv_filepath = 'T:\Mario\Behavior\ACC_pain\behavior cams\_stimuli_timestaps_csv_to_database\csvs';
experiment = 'pain';
warning('off')
files_in_folder = dir(table_filepath);
files_in_folder = {files_in_folder.name};
files = files_in_folder(endsWith(files_in_folder, '.xlsx'));

%% Loop thorugh files
n_files = length(files);
for iif = 1: n_files
    this_file = files{iif};
    this_name = strsplit(this_file, '.xlsx');
    splits = strsplit(this_name{1}, '_');
    date_recording = splits{1};
    animal_ID = strjoin(splits(2:end), '_');
    video_folder = os.path.join('T:\Mario\Behavior\', 'ACC_pain', GC.behavior_video_subfolder, animal_ID, date_recording, experiment);
    filename = os.path.join(video_folder, 'stimuli_timestamps2.csv');
    original_table = readtable([table_filepath, '\',[date_recording, '_', animal_ID, '.xlsx']]);
    fieldnames_table = fieldnames(original_table);
    disp(animal_ID)
    % set parametres
    % Create exemption if the table is organised in ms and sec only
    if ~any(ismember(fieldnames_table, 'min'))
        time_start = num2cell(original_table.sec);
    else
        time_start = num2cell(original_table.min .* 60 + original_table.sec + original_table.msec .*0.01);
    end
    time_end =  num2cell(NaN(length(time_start),1));
    event = original_table.Event;
    withdraw = original_table.Withdrawal;
    aff = original_table.AffectiveResponse;
    
    % Re set values to match data on database
    withdraw = num2cell(strcmp(withdraw, 'yes'));
    new_aff = aff;
    new_aff(strcmp(aff, 'no')) = {''};
%     new_aff(strcmp(aff, 'other')) = {''};
    % Outputs
    var_names = {'Time_start', 'Time_end', 'Event', 'Withdraw', 'Affective_response'};
    output_array = [time_start, time_end, event, withdraw, new_aff];
    output_table = array2table(output_array,'VariableNames', var_names);
    
    % Check if last value is nan
    has_nan = isnan(cell2mat(output_table.Time_start));
    if sum(has_nan) > 0
        output_table(has_nan,:)=[];
    end
    
    % Write csv into csvs folder
    name_csv = [animal_ID, '_', date_recording, '_pain2.csv'];
    csv_filename = [csv_filepath,'\', name_csv];
    writetable(output_table,csv_filename)
    disp('done csv')
    %% Copy files to original forlder to be next updated using C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Code\Superuser\import_behavior_timestamps_to_database.m
    
    copyfile(csv_filename, filename)
    disp('done copying to original video folder')
    
    %% UPLOAD TO DATABASE
    METADATA = SQL_database.read_table_where('sessions', {'experiment_id', 'session_id', 'experimental_condition', 'stimulus'}, {animal_ID, date_recording}, {'animal_ID', 'date'});
    METADATA = METADATA(ismember(METADATA.experiment, experiment), :);
    
    % Get Data
    data = table2cell(readtable(filename));
    table_data = cell(size(data, 1), 5);
    % Convert timestamps
    time = data(:, 1);
    time_s = cell(size(time));
    for i_ts = 1:length(time_s)
        value = time{i_ts};
        time_minutes = floor(value / 60);
        time_seconds = floor(value - time_minutes * 60);
        time_milliseconds = round((value - time_minutes * 60 - time_seconds) * 1000);
        time_s{i_ts} = [num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];
    end
    table_data(:, 1) = time_s;
    
    if size(data, 2) > 2  % new format
        time = data(:, 2);
        time_s = cell(size(time));
        for i_ts = 1:length(time_s)
            value = time{i_ts};
            if isnan(value)
                time_s{i_ts} = '';
            else
                time_minutes = floor(value / 60);
                time_seconds = floor(value - time_minutes * 60);
                time_milliseconds = round((value - time_minutes * 60 - time_seconds) * 1000);
                time_s{i_ts} = [num2str(time_minutes, '%02i'), ':', num2str(time_seconds, '%02i'), '.', num2str(time_milliseconds, '%03i')];
            end
        end
        table_data(:, 2) = time_s;
        table_data(:, 3:5) = data(:, 3:5);
        
    else  % old format
        disp('Old format')
        table_data(:, 3) = data(:, 2);
        table_data(:, 4) = {false};
    end
    
    % Fill in empty cells
    table_data(cellfun(@(x) isempty(x), table_data)) = {''};
    
    % Make sure that the checkbox contains a logical value
    table_data(:, 4) = num2cell(cellfun(@(x) logical(x), table_data(:, 4)));
    
    table_columns = {'Time_start', 'Time_end', 'Event', 'Withdrawal', 'Affective_response'};
    
    
    %% Convert data
    data = cell2table(table_data, 'VariableNames', table_columns);
    % Convert Time_start to seconds
    time_str = data.Time_start;
    time_s = NaN(size(time_str));
    for i_ts = 1:length(time_s)
        this_time_str = time_str{i_ts};
        try
            output = strsplit(this_time_str, ':');
            time_minutes = str2double(output{1}) * 60;
            time_seconds = str2double(output{2});
            time_s(i_ts) = time_minutes + time_seconds;
        end
    end
    % Put data back
    data.Time_start = time_s;
    
    % Convert Time_end to seconds
    time_str = data.Time_end;
    time_s = NaN(size(time_str));
    for i_ts = 1:length(time_s)
        this_time_str = time_str{i_ts};
        try
            output = strsplit(this_time_str, ':');
            time_minutes = str2double(output{1}) * 60;
            time_seconds = str2double(output{2});
            time_s(i_ts) = time_minutes + time_seconds;
        end
    end
    % Put data back
    data.Time_end = time_s;
    
    
    
    
    %% Import data
    stimuli_detected = unique(data.Event);
    
    for i_stim = 1:length(stimuli_detected)
        this_stimulus = stimuli_detected{i_stim};
        rows = find(ismember(data.Event, this_stimulus));
        
        % --- Timestamps ----------
        timestamps = value2str(data{rows, {'Time_start', 'Time_end'}}, '%.3f');
        % Replace 'NaNs' with ''
        timestamps(ismember(timestamps, 'NaN')) = {''};
        % Add double quotes
        timestamps = cellfun(@(x) ['"', x, '"'], timestamps, 'UniformOutput',false);
        % Convert to a long string
        timestamps_ts = cell(size(timestamps, 1), 1);
        for i_ts = 1:size(timestamps, 1)
            timestamps_ts{i_ts} = ['[', strjoin(timestamps(i_ts, :), ', '), ']'];
        end
        timestamps_ts = ['[', strjoin(timestamps_ts, ', '), ']'];
        
        % --- Withdrawal response ----------
        response = strjoin(value2str(data.Withdrawal(rows), '%i'), '');
        
        % --- Affective response ----------
        all_affective_responses = [{''}, GC.freely_moving_affective_responses(:)'];
        affective_responses = data.Affective_response(rows);
        if iscell(affective_responses)
            affective_responses(ismember(affective_responses, ' ')) = {''};
            affective_responses(ismember(affective_responses, 'NaN')) = {''};
            affective_response = strjoin(value2str(cellfun(@(x) find(ismember(all_affective_responses, x)), affective_responses), '%i'), '');
        else
            affective_responses = repmat({''}, length(rows), 1);
        end
        
        % --- Valid stimuli ----------
        valid = strjoin(value2str(true(length(rows), 1), '%i'), '');
        
        % --- Stimulus type ----------
        if ismember(this_stimulus, GC.freely_moving_stimuli)
            stimulus_type = 'evoked';
        else
            stimulus_type = 'spontaneous';
        end
        
        % Get stimulus id
        session_id = METADATA.session_id(ismember(METADATA.stimulus, this_stimulus));
        experiment_id = METADATA.experiment_id(ismember(METADATA.stimulus, this_stimulus));
        stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id', {experiment_id, session_id}, {'experiment_id', 'session_id'}, 'return_as_table',false);
        if ~isempty(stimulus_id)
            disp(['All done for stimulus ''', this_stimulus, '''', ' and database will be modified'])
%             continue
        else
            disp(['Copying stimulus ''', this_stimulus, ''''])
        end
        
        if isempty(stimulus_id)
            % Add experiment and session id, and get newly assigned stimulus_id
            SQL_database.update('stimulations', {'experiment_id', 'session_id'}, {experiment_id, session_id})
            stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id', {experiment_id, session_id}, {'experiment_id', 'session_id'}, 'return_as_table',false);
        end
        SQL_database.update('stimulations', {'timestamps', 'response', 'valid', 'type', 'affective'}, {timestamps_ts, response, valid, stimulus_type, affective_response}, stimulus_id,'stimulus_id')
    end
    
    disp('Finished')
    
    
    
    
    
    
    
    

end
disp('all done')
warning('on')