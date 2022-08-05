function data = execute_SQL_query(DB_conn, SQL_query)
%EXECUTE_SQL_QUERY Runs an SQL query.

% Run cursor
SQL_cursor = exec(DB_conn, SQL_query);
    
% Check for errors in query
if ~isempty(SQL_cursor.Message)
    if endsWith(SQL_cursor.Message, 'MySQL server has gone away')
        % Re-initialize SQL_DATABASE_connection
        global SQL_DATABASE_connection
        SQL_DATABASE_connection = SQL_database.initialize_db_connection();
        % Run cursor
        SQL_cursor = exec(SQL_DATABASE_connection, SQL_query);        
    
    else  % there was a problem
        % Close cursor but keep error message
        SQL_cursor_Message = SQL_cursor.Message;
        close(SQL_cursor)
        
        % Make error message
        error_message = sprintf('ERROR while executing the following SQL query:\n%s\n\n%s', SQL_query, SQL_cursor_Message);

        % User is asking for an uknown column. Suggest which ones are available in error message
        if contains(SQL_cursor_Message, 'Unknown column')
            % Get table name
            if contains(SQL_query, 'FROM')
                % Get names of table
                start_idx = strfind(SQL_query, 'FROM ');
                if ~isempty(start_idx)
                    str = SQL_query(start_idx+5:end);
                    idx = strfind(str, ' ');
                    if ~isempty(idx)
                        str = str(1:idx(1)-1);
                        % Remove name of database
                        idx = strfind(str, '.');
                        table_name = str(idx(1)+1:end);
                        % Get names of columns
                        known_columns = get_column_names_of(DB_conn, table_name, true);
                        % Append to error message
                        error_message = [error_message, '\n\nKnown columns in ''', table_name, ''' are:\n', strjoin(known_columns, ', '), '\n\n'];
                    end
                end
            end
        end
        % Throw error
        error('execute_SQL_query:SQL_error', error_message)
    end
end

% Fetch data if user requested it
SQL_data = fetch(SQL_cursor);

% Check for errors in query
if ~isempty(SQL_data.Message) && ~strcmp(SQL_data.Message, 'Invalid Result Set')
    error('execute_SQL_query:SQL_error', 'ERROR while executing the following SQL query:\n%s\n\n%s', SQL_query, SQL_cursor.Message)
end

if ~strcmp(SQL_data.Message, 'Invalid Result Set')
    % Retrieve data
    data = SQL_data.Data;
else
    data = [];
end

% Close cursor
close(SQL_cursor)


%% MLint exceptions
%#ok<*TLEV>
