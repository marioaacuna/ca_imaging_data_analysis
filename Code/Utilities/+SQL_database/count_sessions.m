function n_sessions = count_sessions(animal_ID)

% Get connector to database
global GC SQL_DATABASE_connection

table_prefix = [GC.database_name, '.', GC.database_table_prefix];

SQL_query = sprintf('SELECT COUNT(DISTINCT(date)) FROM %s_sessions INNER JOIN %s_experiments WHERE (%s_sessions.experiment_id=%s_experiments.experiment_id) AND ((%s_experiments.animal_ID=''%s'' AND %s_sessions.keep=''1''))', ...
    table_prefix, table_prefix, table_prefix, table_prefix, table_prefix, animal_ID, table_prefix);
    
% Execute query
metadata = execute_SQL_query(SQL_DATABASE_connection, SQL_query);
% Extract result
n_sessions = table2array(metadata);
