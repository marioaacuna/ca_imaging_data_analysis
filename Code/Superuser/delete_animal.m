function delete_animal(animal_ID)
%DELETE_ANIMAL Deletes all entries from database regarding one animal.

% Clear Command Window
clc

% Read trials to delete
METADATA = SQL_database.read_table_where('trials', {'animal_ID', 'trial_id'}, animal_ID,'animal_ID');
if ~isempty(METADATA)
    disp('Deleting trials')
    SQL_database.delete_table_where('trials', METADATA(:, 'trial_id'))
end

% Read sessions to delete
METADATA = SQL_database.read_table_where('sessions', {'animal_ID', 'session_id'}, animal_ID,'animal_ID');
if ~isempty(METADATA)
    disp('Deleting sessions')
    SQL_database.delete_table_where('sessions', METADATA(:, 'session_id'))
end

% Read experiments to delete
METADATA = SQL_database.read_table_where('experiments', {'animal_ID', 'experiment_id'}, animal_ID,'animal_ID');
if ~isempty(METADATA)
    disp('Deleting experiments')
    SQL_database.delete_table_where('experiments', METADATA(:, 'experiment_id'))
end

disp('done')
