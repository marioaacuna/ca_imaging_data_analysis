function toggle_toolbox(toolbox_names, new_state, verbose)
%TOGGLE_TOOLBOX Enables or disables toolboxes.

% Get LOGGER and general_configs
global LOGGER GC

% Set default verbose state to true
if ~exist('verbose', 'var'), verbose = false; end

% Make toolbox_names a cell array
if ischar(toolbox_names)
    toolbox_names = {toolbox_names};
end
n_toolboxes = length(toolbox_names);

% Ignore these warnings
warning('off','MATLAB:rmpath:DirNotFound')
warning('off','MATLAB:dispatcher:nameConflict')  % We take the risk of knowing that this might happen

% Prepare string to print on screen
outcome_message = 'Toolbox';
if n_toolboxes > 1
    outcome_message = [outcome_message, 'es'];
    action_verb1 = 'have been';
else
    action_verb1 = 'has been';
end

for itb = 1:n_toolboxes
    % Make path to toolbox folder
    folder = os.path.join(GC.toolboxes_root_path, toolbox_names{itb});
    % Raise an error if it doesn't exit
    if ~exist(folder,'dir')
        msg = ['The toolbox "', toolbox_names{itb}, '" does not exist'];
        if isempty(LOGGER)
            disp(msg)
        else
            LOGGER.critical(msg)
        end
        error('toggle_toolbox:unknown_toolbox', msg)
    end
	
    % Toggle state
    toolbox_folders = genpath(folder);
    switch new_state
        case 'on', addpath(toolbox_folders);
        case 'off', rmpath(toolbox_folders);
    end
end

if verbose
    % Print message on screen
    switch new_state
        case 'on',  action_verb2 = 'enabled';
        case 'off', action_verb2 = 'disabled';
    end
    % Add quotes to names
    toolbox_names = strcat('''', toolbox_names, '''');
    if n_toolboxes > 1
        outcome_message = [outcome_message, ' ', strjoin(toolbox_names(1:end-1),', '), ' and ', toolbox_names{end}];
    else
        outcome_message = [outcome_message, ' ', toolbox_names{1}];
    end
    outcome_message = [outcome_message, ' ', action_verb1, ' ', action_verb2];
    if isempty(LOGGER)
        disp(outcome_message)
    else
        LOGGER.trace(outcome_message)
    end
end

% Reenable these warnings
warning('on','MATLAB:rmpath:DirNotFound')
warning('on','MATLAB:dispatcher:nameConflict')
