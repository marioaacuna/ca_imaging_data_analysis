function delete_table_where(table_name, where_rows, where_columns)
%DELETE_TABLE Delete metadata from a table of the SQL database.

% Read general_configs and get connector to database
global GC SQL_DATABASE_connection

% Make full path of table in database
table_name = [GC.database_table_prefix, '_', table_name];

% Check that inputs exist and are cells

if ~exist('where_rows', 'var') || isempty(where_rows)
    where_rows = {};
else
    if istable(where_rows)
        % Get names of columns and convert to cell array
        where_columns = where_rows.Properties.VariableNames;
        where_rows = table2cell(where_rows);
    elseif ~iscell(where_rows)
        where_rows = {where_rows};
    end
end
if ~exist('where_columns', 'var') || isempty(where_columns)
    where_columns = {};
else
    if ~iscell(where_columns)
        where_columns = {where_columns};
    end
end

% Make WHERE query
[where_columns, join_where] = make_SQL_query(SQL_DATABASE_connection, where_columns, table_name, false);
% Get unique elements of the uniqe clause
tables2join = unique(join_where(:,1));
join_where  = unique(join_where(:,2));

% Concatenate queries in one statement (and add JOIN clause, if any)
SQL_query_delete = ['DELETE FROM ', GC.database_name, '.', table_name];
if ~isempty(tables2join) && ~isempty(join_where)
    SQL_query_delete = [SQL_query_delete, ' INNER JOIN ', strjoin(tables2join,','), ' WHERE (', strjoin(join_where,' AND '), ')'];
end

% Convert all cells to match to characters
non_char_cells = cellfun(@(x) ~ischar(x), where_rows);
if any(non_char_cells(:))
    idx = find(non_char_cells);
    for ii = 1:length(idx)
        value = where_rows{idx(ii)};
        if isnumeric(value) || islogical(value)
            value = num2str(value);
        end
        where_rows{idx(ii)} = value;
    end
end

% Make where clause
n_rows_to_match = size(where_rows,1);
if n_rows_to_match > 0
    SQL_query_where = cell(n_rows_to_match, 1);
    for irow = 1:n_rows_to_match
        % Prepare data by escaping backslashes in paths, if any
        row_data = cellfun(@(x) regexprep(x, '\\', '\\\'), where_rows(irow,:), 'UniformOutput',false);
        % Get which cells are empty
        is_null = cellfun(@isempty, row_data);
        % Write columns that are not empty
        SQL_query_where{irow} = ['(', strjoin(cellfun(@(x,y) [x, '=''', y, ''''], where_columns(~is_null), row_data(~is_null), 'UniformOutput',false), ' AND ')];
        % Replace empty cells with 'IS NULL'
        if any(is_null)
            SQL_query_where{irow} = [SQL_query_where{irow}, ' AND ', strjoin(cellfun(@(x) [x, ' IS NULL'], where_columns(is_null), 'UniformOutput',false), ' AND ')];
        end
        % Add closing parenthesis
        SQL_query_where{irow} = [SQL_query_where{irow}, ')'];
    end
    % Concatenate conditions
    SQL_query_where = strjoin(SQL_query_where, ' OR ');

    % If we have already a WHERE clause in the SELECT query do not repeat keyword
    if isempty(strfind(SQL_query_delete, 'WHERE')) %#ok<STREMP>
        SQL_query_where = [' WHERE ', SQL_query_where];
    else  % Add parentheses
        SQL_query_where = [' AND (', SQL_query_where, ')'];
    end
    
else
    SQL_query_where = '';
end

% Execute query
execute_SQL_query(SQL_DATABASE_connection, [SQL_query_delete, SQL_query_where]);
