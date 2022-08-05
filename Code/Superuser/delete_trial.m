function delete_trial(animal_ID, day_index, session_stimulus, trial_to_delete)
%DELETE_TRIAL Deletes a trial from the database.

% Clear Command Window, and read general_configs and SQL connector
clc
global GC SQL_DATABASE_connection

% Read database
METADATA = SQL_database.read_table_where('trials', {'+','date','stimulus'}, animal_ID,'animal_ID');
% Add day index
[~, ~, METADATA.day_index] = unique(METADATA.date, 'stable');
% Get rows for trials during the chosen session
trial_rows = find(ismember(METADATA.day_index, day_index) & ismember(METADATA.stimulus, session_stimulus));
metadata_rows_to_delete = trial_rows(trial_to_delete);
trial_id_to_delete = METADATA.trial_id(metadata_rows_to_delete);
% Execute SQL query
SQL_cursor = exec(SQL_DATABASE_connection, ['DELETE FROM `nevian2`.`2p_imaging_trials` WHERE `trial_id`=''', num2str(trial_id_to_delete) ,''';']);
close(SQL_cursor);

% Move .tiff file to trash
path_to_file = [GC.data_root_path, GC.tiff_raw_path, METADATA.tiff_path{metadata_rows_to_delete}];
[~, name, ext] = fileparts(path_to_file);
previousState = recycle('on');
movefile(path_to_file, GC.temp_dir)
delete([GC.temp_dir, name, ext])
recycle(previousState)
