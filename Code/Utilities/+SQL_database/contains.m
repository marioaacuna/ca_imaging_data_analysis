function in_DB = contains(table_name, data_to_match, varargin)
%CONTAINS Checks whether data exist in SQL database.
% Parse inputs
p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'table_name', @ischar)
addRequired(p, 'data_to_match', @(x) istable(x) || iscell(x))
addParameter(p, 'columns', @iscell);
parse(p, table_name, data_to_match, varargin{:});

% Read general_configs and get connector to database
global GC SQL_DATABASE_connection

% Unpack inputs
data_to_match = p.Results.data_to_match;
if ~istable(data_to_match)
    columns = p.Results.columns;
else
    % Get names of columns and convert to cell array
    columns = data_to_match.Properties.VariableNames;
    data_to_match = table2cell(data_to_match);
end
% Get name of table
table_name = [GC.database_table_prefix, '_', p.Results.table_name];

% Initialize output as all false
in_DB = false(size(data_to_match, 1), 1);

% Get the row count from the table
output = table2cell(execute_SQL_query(SQL_DATABASE_connection, ['SELECT COUNT(*) FROM ', GC.database_name, '.', table_name]));
row_count = output{1};

% If table is empty, it will not contain our data
if row_count == 0
    return
end

% Make SELECT query
SQL_query_select = make_SQL_query(SQL_DATABASE_connection, columns, table_name);

% Make the WHERE clause
% Convert all cells to characters
non_char_cells = cellfun(@(x) ~ischar(x), data_to_match);
if any(non_char_cells(:))
    idx = find(non_char_cells);
    for ii = 1:length(idx)
        value = data_to_match{idx(ii)};
        if isnumeric(value) || islogical(value)
            value = num2str(value);
        end
        data_to_match{idx(ii)} = value;
    end
end

% Concatenate all data_to_match
for irow = 1:size(data_to_match,1)
    % Prepare data by escaping backslashes in paths, if any
    row_data = cellfun(@(x) regexprep(x, '\\', '\\\'), data_to_match(irow,:), 'UniformOutput',false);
    
    % Get which cells are empty
    is_null = cellfun(@isempty, row_data);
    % Write columns that are not empty
    SQL_query_where = strjoin(cellfun(@(x,y) [x, '=''', y, ''''], columns(~is_null), row_data(~is_null), 'UniformOutput',false), ' AND ');
    
    % Replace empty cells with 'IS NULL'
    if any(is_null)
        SQL_query_where = [SQL_query_where, ' AND ', strjoin(cellfun(@(x) [x, ' IS NULL'], columns(is_null), 'UniformOutput',false), ' AND ')];
    end
    
    % If we have already a WHERE clause in the SELECT query do not repeat keyword
    if isempty(strfind(SQL_query_select, 'WHERE')) %#ok<STREMP>
        SQL_query_where = ['WHERE (', SQL_query_where, ')'];
    else  % Prepend AND
        SQL_query_where = [' AND ', SQL_query_where];
    end
    
    % Execute query
    output = execute_SQL_query(SQL_DATABASE_connection, [SQL_query_select, ' ', SQL_query_where]);
    if istable(output) && height(output) > 0
        in_DB(irow) = true;
    end
end

%#ok<*AGROW>
