clear, clc
global GC

animal_ID = 'MA_31epi'
session_date = '190919'
experiment = 'pain'

video_folder = os.path.join('T:\Mario\Behavior', 'ACC_pain', GC.behavior_video_subfolder, animal_ID, session_date, experiment);
filename = os.path.join(video_folder, 'stimuli_timestamps.csv');
if ~exist(filename, 'file'), disp('File does not exist'), return, end

backup_folder = os.path.join(GC.temp_dir, 'stimuli_timestamps_backup');
if ~exist(backup_folder, 'dir'), mkdir(backup_folder), end
backup_file = os.path.join(backup_folder, [animal_ID, '_', session_date, '_', experiment, '.csv']);
copyfile(filename, backup_file)

% Get session id
METADATA = SQL_database.read_table_where('sessions', {'experiment_id', 'session_id', 'experimental_condition', 'stimulus'}, {animal_ID, session_date}, {'animal_ID', 'date'});
METADATA = METADATA(ismember(METADATA.experiment, experiment), :);


%% Load data
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
        disp(['All done for stimulus ''', this_stimulus, ''''])
        continue
    else
        disp(['Copying stimulus ''', this_stimulus, ''''])
    end
    
    if isempty(stimulus_id)
        % Add experiment and session id, and get newly assigned stimulus_id
        SQL_database.update('stimulations', {'experiment_id', 'session_id'}, {experiment_id, session_id})
        stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id', {experiment_id, session_id}, {'experiment_id', 'session_id'}, 'return_as_table',false);
    end
    SQL_database.update('stimulations', {'timestamps', 'response', 'valid', 'type'}, {timestamps_ts, response, valid, stimulus_type}, stimulus_id,'stimulus_id')
end

disp('Finished')


