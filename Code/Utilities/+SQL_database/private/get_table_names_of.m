function table_names = get_table_names_of(DB_conn)
%GET_COLUMN_NAMES_OF Returns the name of the columns in the table.

% Read general configs
global GC

% Make the SQL query
SQL_query = ['SHOW tables FROM ', GC.database_name];
data = execute_SQL_query(DB_conn, SQL_query);
table_names = table2cell(data);

% Keep only tables that are related to calcium imaging experiments
table_names = table_names(cell2mat(cellfun(@(x) ~isempty(regexp(x, ['^', GC.database_table_prefix, '_'],'once')), table_names, 'UniformOutput',false)));
