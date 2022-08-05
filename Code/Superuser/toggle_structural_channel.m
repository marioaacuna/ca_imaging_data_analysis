function toggle_structural_channel(new_state, animal_ID, sessions_date)
%TOGGLE_STRUCTURAL_CHANNEL Turns on or off the structural_channel property of
%some sessions.

% Clear Command Window
clc

% If no session provided, apply to all
if ~exist('sessions_date', 'var')
    sessions_date = unique(SQL_database.read_table_where('sessions', 'date', animal_ID,'animal_ID', 'return_as_table',false));
end
if ~iscell(sessions_date), sessions_date = {sessions_date}; end
% Log progress
disp(['Will affect the following sessions of mouse ''', animal_ID, ''':'])
disp(sessions_date)

% Convert toggle state to number
switch new_state
    case 'on', new_state = 1;
    case 'off', new_state = 0;
    otherwise, error('toggle_structural_channel:new_state', 'Unknown state %s', new_state)
end
disp('Will toggle the structural channel to:')
disp(new_state)

% Read ID of trials to modify
what_to_look_for = [repmat({animal_ID}, length(sessions_date), 1), sessions_date(:)];
trials_id = SQL_database.read_table_where('trials', 'trial_id', what_to_look_for,{'animal_ID','date'}, 'return_all_sessions',true, 'return_as_table',false);

% Update database
SQL_database.update('trials', 'has_structural_channel', {new_state}, trials_id, 'trial_id')
disp('done')
