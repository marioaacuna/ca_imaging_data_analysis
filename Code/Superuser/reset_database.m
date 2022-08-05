function reset_database()
%RESET_DATABASE Deletes all entries from database.

% Clear Command Window
clc

% Read general configs and get connector to database
global GC SQL_DATABASE_connection

% Move to folder where functions to interface with SQL database are located
cd([GC.repository_root_path, 'Code', filesep, 'Utilities', filesep, '+SQL_database', filesep, 'private'])

% First, disable checks on foreign keys, which allows us to truncate tables
% while reseeding primary and foreign keys.
execute_SQL_query(SQL_DATABASE_connection, 'SET FOREIGN_KEY_CHECKS=0;');

% Name of tables to delete
table_names = get_table_names_of(SQL_DATABASE_connection);

% Show final message before executing the SQL query
disp('Will reset:')
disp(table_names)
fprintf('\nPress ANY key to continue, Ctrl+C to abort.\nTHERE IS NO WAY BACK!\n')
pause

for t = 1:length(table_names)
    % Delete data from table
    execute_SQL_query(SQL_DATABASE_connection, ['TRUNCATE ', GC.database_name, '.', table_names{t}]);    
end

% Lastly, re-enable checks on foreign keys
execute_SQL_query(SQL_DATABASE_connection, 'SET FOREIGN_KEY_CHECKS=1;');

% Communicate outcome
disp('done')
