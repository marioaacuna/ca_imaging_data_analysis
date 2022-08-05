function columns = get_column_names_of(DB_conn, table_name, include_foreign_keys)
%GET_COLUMN_NAMES_OF Returns the name of the columns in the table.

% Read general configs
global GC

if contains(DB_conn.Type, 'ODBC') || ~ispc()
    use_ODBC = true;
else
    use_ODBC = false;
end

% Make the SQL query
SQL_query = ['SHOW columns FROM ', GC.database_name, '.', table_name];
data = execute_SQL_query(DB_conn, SQL_query);
% Get names of columns
if use_ODBC 
    columns = data.Field(:);
else
    columns = data.COLUMN_NAME(:);
end

if exist('include_foreign_keys','var') && ~include_foreign_keys
    if use_ODBC
        columns(ismember(data.Key,'MUL')) = [];
    else
        columns(ismember(data.COLUMN_KEY,'MUL')) = [];
    end
end
