function update(table_name_without_prefix, columns_to_update, new_values, where_rows, where_columns)
%UPDATE_EXPERIMENTS Update experiments table.

% Read general configs and get connector to database
global GC SQL_DATABASE_connection

% Make full path of table in database
table_name = [GC.database_table_prefix, '_', table_name_without_prefix];

% If not provided, the columns to update are all the columns
if ~exist('columns_to_update', 'var') || isempty(columns_to_update)
    if istable(new_values)
        columns_to_update = new_values.Properties.VariableNames;
    else
        columns_to_update = get_column_names_of(SQL_DATABASE_connection, table_name, true);
    end
else
    % Make sure that input is a cell arrays
    if ~iscell(columns_to_update)
        columns_to_update = {columns_to_update};
    end
end

% Make new values a vertical cell array
if istable(new_values)
    new_values = table2cell(new_values);
end

% Check that inputs exist and are cells
if ~exist('where_rows', 'var') || isempty(where_rows)
    where_rows = {};
else
    if istable(where_rows)
        where_columns = where_rows.Properties.VariableNames;
        where_rows = table2cell(where_rows);
    elseif isnumeric(where_rows)
        where_rows = num2cell(where_rows);
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

% Make WHERE clause
[where_columns, join_where] = make_SQL_query(SQL_DATABASE_connection, where_columns, table_name, false);
% Get unique elements of the uniqe clause
tables2join = unique(join_where(:,1));
join_where = unique(join_where(:,2));

% If we are checking multiple rows but there is only one value to update, it
% means that the same value should be replicated to all rows that meet the
% where clause.
if ~isempty(where_rows)
    if length(where_rows)==1 && strcmpi(where_rows{1}, 'all') && isempty(where_columns)
        n_rows = table2cell(execute_SQL_query(SQL_DATABASE_connection, ['SELECT COUNT(*) FROM ', GC.database_name, '.', table_name]));
        n_rows = n_rows{1};
        update_all_rows = true;
    else
        n_rows = size(where_rows, 1);
        update_all_rows = false;
    end
else
    n_rows = 1;
	update_all_rows = false;
end

% Check here whether we have to update existing rows or simply add new rows
if isempty(where_rows) && isempty(where_columns)
    action = 'insert';
    % Make sure we have to insert data. Check the unique keys for this table.
%     description = execute_SQL_query(SQL_DATABASE_connection, ['SHOW CREATE TABLE ', table_name ,';']);
%     unique_keys = extractBetween(description.CreateTable{1}, 'UNIQUE KEY `', '_UNIQUE');
%     if any(ismember(columns_to_update, unique_keys))
%         % Check whether the data that we are about to insert in these unique
%         % columns do not exist already in the database.
%         data_to_match = new_values(:, ismember(columns_to_update, unique_keys));
%         unique_columns = columns_to_update(ismember(columns_to_update, unique_keys));
%         in_DB = SQL_database.contains(table_name_without_prefix, data_to_match, 'columns',unique_columns);
%         if any(in_DB)
%             action = 'update';
%         end
%         
%         % Insert new data
%         data_to_insert = new_values(~in_DB, :);
%         if ~isempty(data_to_insert)
%             SQL_database.update(table_name_without_prefix, columns_to_update, data_to_insert)
%         end
%         
%         % Continue with updating existing data
%         new_values = new_values(in_DB, :);
%         overlapping_values = ismember(unique_keys, columns_to_update);
%         where_columns = unique_keys(overlapping_values);
%         where_rows = new_values(:, overlapping_values);
%     end
else
    action = 'update';
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
% Convert all values to insert to characters
if ~iscell(new_values)
    if size(new_values, 1) == size(where_rows, 1) && ismatrix(new_values)  % User input an array
        new_values = num2cell(new_values);
    end
end
non_char_cells = cellfun(@(x) ~ischar(x), new_values);
if any(non_char_cells(:))
    idx = find(non_char_cells);
    for ii = 1:length(idx)
        value = new_values{idx(ii)};
        if isnumeric(value) || islogical(value)
            value = num2str(value);
        end
        new_values{idx(ii)} = value;
    end
end

% Replicate value in case user provided only one for all rows to update
if length(new_values) == 1 && ~update_all_rows && n_rows ~= 1
    new_values = repmat(new_values, n_rows, 1);
end
% Get number of values to update
n_values_to_update = size(new_values, 1);

switch action
    case 'insert'
        % Make fixed part of the SQL query
        SQL_query_insert = ['INSERT INTO ', GC.database_name, '.', table_name, ' (', strjoin(columns_to_update, ', '), ')'];

        % Loop through values to update
        for irow = 1:n_values_to_update
            % Prepare data to insert by escaping backslashes in paths, if any
            data_to_insert = cellfun(@(x) regexprep(x, '\\', '\\\'), new_values(irow, :), 'UniformOutput',false);
            % Get which cells are empty
            is_null = cellfun(@isempty, data_to_insert);
            % Write columns that are not empty
            SQL_query_values = data_to_insert;
            SQL_query_values(~is_null) = cellfun(@(x) ['''', x, ''''], data_to_insert(~is_null), 'UniformOutput',false);
            % Replace empty cells with 'NULL'
            if any(is_null)
                SQL_query_values(is_null) = {'NULL'};
            end
            % Concatenate values
            SQL_query_values = ['VALUES (', strjoin(SQL_query_values, ', '), ')'];
            
            % Concatenate statements in one query
            SQL_query = [SQL_query_insert, ' ', SQL_query_values];
            % Execute statement
            execute_SQL_query(SQL_DATABASE_connection, SQL_query);
        end

    case 'update'
        % Make fixed part of the SQL query
        SQL_query_update = ['UPDATE ', GC.database_name, '.', table_name];
        
        % Make a join statement
        if ~isempty(tables2join) && ~isempty(join_where)
            SQL_query_update = [SQL_query_update, ' INNER JOIN ', strjoin(tables2join,','), ' ON (', strjoin(join_where,' AND '), ')'];
        end

        % Loop through values to update
        for irow = 1:n_values_to_update
            % Prepare data to insert by escaping backslashes in paths, if any
            data_to_insert = cellfun(@(x) regexprep(x, '\\', '\\\'), new_values(irow, :), 'UniformOutput',false);
            % Get which cells are empty
            is_null = cellfun(@isempty, data_to_insert);
            % Write columns that are not empty
            SQL_query_set = ['SET ', strjoin(cellfun(@(x,y) [x, '=''', y, ''''], columns_to_update(~is_null), data_to_insert(~is_null), 'UniformOutput',false), ', ')];
            % Add a comma only if there is a mix of null and non-null values
            if any(is_null) && any(~is_null)
                query_conjunction = ', ';
            else
                query_conjunction = '';
            end
            % Replace empty cells with 'NULL'
            if any(is_null)
                SQL_query_set = [SQL_query_set, query_conjunction, strjoin(cellfun(@(x) [x, '=NULL'], columns_to_update(is_null), 'UniformOutput',false), ', ')];
            end
            
            % If updating a subset of rows, make a WHERE statement
            if ~update_all_rows
                % Prepare data to match by escaping backslashes in paths, if any
                data_to_match = cellfun(@(x) regexprep(x, '\\', '\\\'), where_rows(irow,:), 'UniformOutput',false);
                % Get which cells are empty
                is_null = cellfun(@isempty, data_to_match);
                % Write columns that are not empty
                SQL_query_where = ['WHERE (',  strjoin(cellfun(@(x,y) [x, '=''', y, ''''], where_columns(~is_null), data_to_match(~is_null), 'UniformOutput',false), ' AND ')];
                % Replace empty cells with 'IS NULL'
                if any(is_null)
                    SQL_query_where = [SQL_query_where, ' AND ', strjoin(cellfun(@(x) [x, ' IS NULL'], where_columns(is_null), 'UniformOutput',false), ' AND ')];
                end
                % Add closing parenthesis
                SQL_query_where = [SQL_query_where, ')'];
            else
                SQL_query_where = '';
            end
            
            % Concatenate statements in one query
            SQL_query = [SQL_query_update, ' ', SQL_query_set, ' ', SQL_query_where];
            % Execute statement
            execute_SQL_query(SQL_DATABASE_connection, SQL_query);
        end
end

%% Mlint
%#ok<*AGROW>
