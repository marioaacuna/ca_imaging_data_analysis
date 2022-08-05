function delete_session(animal_ID, session_date)
%DELETE_ANIMAL Deletes all entries from database regarding one animal.

% Clear Command Window
clc

% Read trials to delete
METADATA = SQL_database.read_table_where('trials', {'animal_ID', 'trial_id'}, {animal_ID, session_date},{'animal_ID','date'}, 'return_all_sessions',true);
if ~isempty(METADATA)
    disp('Deleting trials')
    SQL_database.delete_table_where('trials', METADATA(:, 'trial_id'))
end

% Read sessions to delete
METADATA = SQL_database.read_table_where('sessions', 'session_id', {animal_ID, session_date},{'animal_ID','date'}, 'return_all_sessions',true);
if ~isempty(METADATA)
    disp('Deleting sessions')
    SQL_database.delete_table_where('sessions', METADATA(:, 'session_id'))
end

disp('done')
