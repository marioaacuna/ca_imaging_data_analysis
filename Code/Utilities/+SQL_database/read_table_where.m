function metadata = read_table_where(table_name, return_columns, where_rows, where_columns, varargin)
%READ_TABLE Reads metadata from a table of the SQL database.

% Get connector to database
global GC SQL_DATABASE_connection

% Fix missing inputs
if ~exist('return_columns', 'var'), return_columns={}; end
if ~exist('where_rows', 'var'), where_rows={}; end
if ~exist('where_columns', 'var'), where_columns={}; end

% Parse inputs
p = inputParser();
addOptional(p, 'return_columns', [], @(x) isempty(x) || ischar(x) || iscell(x))
addOptional(p, 'where_rows', {}, @(x) isempty(x) || ischar(x) || isnumeric(x) || iscell(x) || istable(x))
addOptional(p, 'where_columns', {}, @(x) isempty(x) || ischar(x) || iscell(x))
addParameter(p, 'return_as_table', true, @(x) islogical(x) || (isnumeric(x) && length(x(:))==1 && ismember(x, [0,1])))
addParameter(p, 'split_experimental_condition', true, @(x) islogical(x) || (isnumeric(x) && length(x(:))==1 && ismember(x, [0,1])))
addParameter(p, 'return_all_sessions', false, @(x) islogical(x) || (isnumeric(x) && length(x(:))==1 && ismember(x, [0,1])))
% Concatenate all user inputs
args = {return_columns, where_rows, where_columns};
if ~isempty(varargin)
    args = [args, varargin];
end
% Look for parameter names
param_names = {'return_as_table', 'split_experimental_condition', 'return_all_sessions'};
found_param = false;
for i_arg = 1:length(args)
    try
        found_param = ismember(args{i_arg}, param_names);
    end
    if found_param
        break
    end
end
if found_param
    first_param = i_arg;
    varargin = args(first_param:end);
    return_columns = {};
    where_rows = {};
    where_columns = {};
    for i_arg = 1:first_param-1
        switch i_arg
            case 1
                return_columns = args{i_arg};
            case 2
                where_rows = args{i_arg};
            case 3
                where_columns = args{i_arg};
        end
    end
end

% Parse inputs and assign defaults
parse(p, return_columns, where_rows, where_columns, varargin{:})
% Retrieve inputs
return_columns = p.Results.return_columns;
where_rows = p.Results.where_rows;
where_columns = p.Results.where_columns;
return_as_table = logical(p.Results.return_as_table);
split_experimental_condition = logical(p.Results.split_experimental_condition);
return_all_sessions = logical(p.Results.return_all_sessions);

criteria = [];

% Make full path of table in database
full_table_name = [GC.database_table_prefix, '_', table_name];
% Check that table exists
all_table_names = get_table_names_of(SQL_DATABASE_connection);
if ~ismember(full_table_name, all_table_names)
    error('SQL:unknownTable', 'Unknown table ''%s''', full_table_name)
end

% Get name of columns
columns_to_delete = {};
return_all_columns = false;
additional_columns = {};
if isempty(return_columns)
    return_all_columns = true;
else
    try
        if iscell(return_columns) && strcmp(return_columns{1}, '+')
            return_all_columns = true;
            additional_columns = return_columns(2:end);
        end
    end
end
if return_all_columns
    return_columns = get_column_names_of(SQL_DATABASE_connection, full_table_name);
end
if ~isempty(additional_columns)
    return_columns = [return_columns(:); additional_columns(:)];
end
if ~iscell(return_columns)
    return_columns = {return_columns};
end

% Add other columns used simply for splitting
criteria_to_match = {};
if split_experimental_condition && ...
        ismember('experimental_condition', return_columns) && ...
        ~isempty(GC.experiment_name) && ...
        ~isempty(GC.columns_experimental_condition.(GC.experiment_name))
    add_column_names = GC.columns_experimental_condition.(GC.experiment_name);
    if size(add_column_names, 1) > 1
        criteria = cellfun(@(x) regexp(x, '=', 'split'), add_column_names(:, 1), 'UniformOutput',false);
        criteria_to_match = cellfun(@(x) x{1}, criteria, 'UniformOutput',false);
        columns_experimental_condition = unique(criteria_to_match);
        columns_to_add = setdiff(columns_experimental_condition, return_columns);
        return_columns = [return_columns(:); columns_to_add(:)];
        columns_to_delete = [columns_to_delete; columns_to_add(:)];
    end
end

% Check that inputs exist and are cells
if istable(where_rows)
    % Get names of columns and convert to cell array
    where_columns = where_rows.Properties.VariableNames;
    where_rows = table2cell(where_rows);
elseif ischar(where_rows)
    where_rows = {where_rows};
elseif isnumeric(where_rows)
    if numel(where_rows) == 1
        where_rows = {where_rows};
    else
        where_rows = num2cell(where_rows);
    end
elseif iscell(where_rows)
else
    error('SQL_database:read_table_where', 'where_rows should be a cell array, a table or a string or a number')
end
if ~iscell(where_columns)
    where_columns = {where_columns};
end

% Add check on 'kept' column in not returning all sessions
if ~return_all_sessions && ismember(table_name, {'sessions', 'trials', 'stimulations'})  % Cannot read from higher-order table, such as 'experiments'
    where_columns = [where_columns, 'keep'];
    where_rows(:, end+1) = {1};
end

% Add columns that will be used to sort the ouptut table
original_return_columns = return_columns(:);
% Remove columns that exist only in child tables
switch table_name
    case 'experiments'
        return_columns(ismember(return_columns, 'date') | ismember(return_columns, 'recording_time')) = [];
        sorting_columns = {'animal_ID'};
    case 'sessions'
        return_columns(ismember(return_columns, 'recording_time')) = [];
        sorting_columns = {'animal_ID', 'date'};
    case 'trials'
        sorting_columns = {'animal_ID', 'date', 'recording_time'};
    case 'stimulations'
        sorting_columns = {'animal_ID', 'date'};
end
return_columns = unique([return_columns(:); sorting_columns(:)]);

% Make SELECT and WHERE queries
[select_columns, join_select] = make_SQL_query(SQL_DATABASE_connection, return_columns, full_table_name, false);
[where_columns, join_where] = make_SQL_query(SQL_DATABASE_connection, where_columns, full_table_name, false);
% Get unique elements of the uniqe clause
tables2join = unique([join_select(:,1); join_where(:,1)]);
join_where = unique([join_select(:,2); join_where(:,2)]);

% Concatenate queries in one statement (and add JOIN clause, if any)
SQL_query_select = ['SELECT ', strjoin(select_columns,','), ' FROM ', GC.database_name, '.', full_table_name];
if ~isempty(tables2join) && ~isempty(join_where)
    SQL_query_select = [SQL_query_select, ' INNER JOIN ', strjoin(tables2join,','), ' WHERE (', strjoin(join_where,' AND '), ')'];
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
    if isempty(strfind(SQL_query_select, 'WHERE')) %#ok<STREMP>
        SQL_query_where = [' WHERE ', SQL_query_where];
    else  % Add parentheses
        SQL_query_where = [' AND (', SQL_query_where, ')'];
    end
    
else
    SQL_query_where = '';
end

% Execute query
metadata = execute_SQL_query(SQL_DATABASE_connection, [SQL_query_select, SQL_query_where]);

% Data should be a table
if ~istable(metadata) && (isnumeric(metadata) || (length(metadata) == 1 && strcmp(metadata{1}, 'No Data')))
    metadata = table();
end

if ~isempty(metadata)
    % Sort and delete additional columns
    if ~isempty(sorting_columns)
        columns_to_sort = NaN(length(sorting_columns), 1);
        for icol = 1:length(sorting_columns)
            columns_to_sort(icol) = find(ismember(metadata.Properties.VariableNames, sorting_columns{icol}));
        end
        metadata = natsortrows(metadata, columns_to_sort);

        % Keep only columns that the user did request
        metadata = metadata(:, original_return_columns);
    end
    column_names = metadata.Properties.VariableNames;
    
    % Process special columns, by converting values or splitting them in subcolumns
    if split_experimental_condition && ismember('experimental_condition', column_names)
        metadata = process_experimental_condition(metadata, criteria, criteria_to_match);
    end
    if strcmp(table_name, 'stimulations')
        metadata = process_stimulations(metadata);
    end

    % Delete unused columns
    if ~isempty(columns_to_delete)
        metadata(:, columns_to_delete) = [];
    end
    
    % Convert from table, if user requested it
    if ~return_as_table
        metadata = table2cell(metadata);
        if length(metadata) == 1
            metadata = metadata{1};
        end
    end

else
    if ~return_as_table
        metadata = {};
    end
end
