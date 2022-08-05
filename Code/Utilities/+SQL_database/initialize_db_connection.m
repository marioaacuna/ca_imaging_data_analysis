function DB_conn = initialize_db_connection()
%INITIALIZE_DB_CONNECTION Starts connection with database.

% Set new preferences
setdbprefs({'DataReturnFormat', 'NullStringRead', 'NullNumberRead'}, ...
           {'table',            '',               'NaN'})

% Connect to database
if ispc()
    DB_conn = database('neviandb', '', '');
else
    global GC  %#ok<TLEV>

    jdbc_config_table = readtable(GC.jdbc_datasource_config);
    DB_conn = database(jdbc_config_table.database{1}, ...
                       jdbc_config_table.username{1}, ...
                       jdbc_config_table.password{1});%, ...
%                        'Vendor', jdbc_config_table.vendor{1}, ...
%                        'Server', jdbc_config_table.server{1});
end
