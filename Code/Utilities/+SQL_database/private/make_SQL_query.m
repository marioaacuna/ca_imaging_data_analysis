function varargout = make_SQL_query(DB_conn, columns, table_name, return_query)

% Default input
if ~exist('return_query', 'var'), return_query=true; end

% Check that something has to be found
if isempty(columns)
    if return_query
        varargout{1} = '';
    else
        varargout{1} = columns;
        varargout{2} = cell(0,2);
    end
    % Quit function
    return
end

% Read general configs
global GC

% Copy input array
out_columns = columns;

% Pre-allocate output variable
SQL_query_join = cell(0, 2);

% Get the names of the database
db_name = GC.database_name;

% Get the names of the tables
all_table_names = get_table_names_of(DB_conn);

% Get index of target table
table_idx = ismember(all_table_names, table_name);
table_columns = get_column_names_of(DB_conn, all_table_names{table_idx});

% Get columns in parent table
columns_in_table_idx = ismember(columns, table_columns);
columns_in_table = columns(columns_in_table_idx);
% Add name of table to columns
out_columns(columns_in_table_idx) = strcat(table_name, '.', columns_in_table(:));

% Check whether the user wants to select columns that are not in this table
if any(~ismember(columns, table_columns))
    % Get columns in other tables
    columns_not_in_table_idx = ~ismember(columns, table_columns);
    columns_not_in_table = columns(columns_not_in_table_idx);
    
    % Get name of columsn from other tables
    other_table_names = all_table_names(~table_idx);
    n_other_tables = length(other_table_names);
    other_table_columns = cell(n_other_tables, 1);
    for itable = 1:n_other_tables
        other_table_columns{itable} = get_column_names_of(DB_conn, other_table_names{itable}, false);
    end
    % Concatenate columns names
    table_id = cellfun(@(x,y) ones(length(x),1)*y, other_table_columns, num2cell((1:n_other_tables)'), 'UniformOutput',false);
    table_id = num2cell(cat(1, table_id{:}));
    other_table_columns = [cat(1, other_table_columns{:}), table_id];
    % Find columns not in table
    joined_foreign_tables = {};
    for icol = 1:length(columns_not_in_table)
        % Get name of foreign table in lookup table
        foreign_table = other_table_names{other_table_columns{ismember(other_table_columns(:,1), columns_not_in_table{icol}), 2}};
        % Concatenate the name of the foreign table to its column
        columns_not_in_table{icol} = [foreign_table, '.', columns_not_in_table{icol}];
        
        % Make a SQL join statement to look in the foreign table
        if ~ismember(foreign_table, joined_foreign_tables)
            % Get the column name containing the id
            id_name = regexp(foreign_table, '_', 'split');
            id_name = [id_name{end}(1:end-1), '_id'];
            % Add statement to SQL query
            SQL_query_join = [SQL_query_join; [{[db_name, '.', foreign_table]}, {[db_name, '.', table_name, '.', id_name, '=', db_name, '.', foreign_table, '.', id_name]}]];
            
            % Mark the foreign table as already added, so we won't do it twice
            joined_foreign_tables = [joined_foreign_tables; foreign_table];
        end
    end
    % Merge all columns
    out_columns(columns_not_in_table_idx) = columns_not_in_table;

end

% Prepend name of database
out_columns = strcat(db_name, '.', out_columns);

% Make SQL query from columns and join statement
if return_query
    SQL_query = ['SELECT ', strjoin(out_columns, ','), ' FROM ', db_name, '.', table_name];
    if size(SQL_query_join, 1) > 0
        SQL_query = [SQL_query, ' INNER JOIN ', strjoin(SQL_query_join(:,1),','), ' WHERE (', strjoin(SQL_query_join(:,2),' AND '), ')'];
    end
    varargout{1} = SQL_query;

else
    varargout{1} = out_columns;
    varargout{2} = SQL_query_join;
end


%% MLint exception
%#ok<*AGROW>
