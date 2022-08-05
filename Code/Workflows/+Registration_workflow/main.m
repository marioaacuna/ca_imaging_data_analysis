%% PREAMBLE
% This function launches the Preprocessing workflow, which extract calcium
% traces from image stacks recorded at the two-photon microscope.

%% MAIN
function done = main(INFO)
done = true;  % This is the output of the workflow, which is read by the launcher GUI and used to continue to the next workflow

% Get logger and general_configs
global LOGGER GC

% Unpack INFO to make all variables of this script explicit
experiments = INFO.experiments;
data_types = INFO.sessions_table.data_type(ismember(INFO.sessions_table.animal_ID, INFO.experiments));
% Actions to perform
do_registration = INFO.actions.registration.main;
n_actions_to_perform = sum([do_registration]);
if n_actions_to_perform < 1, return, end

% Options for each action
skip_registration         = INFO.actions.registration.skip;
registration_missing_only = INFO.actions.registration.missing;

% Log beginning of workflow
LOGGER.info('Preprocessing workflow', 'decorate',true)

% Split experiments by animal
all_animal_names = natsort(experiments);

%% IMAGE REGISTRATION
if do_registration
    % Enable toolboxes used by this workflow
    toolboxes_to_use = {'NoRMCorre', 'ScanImageTiffReader', '_plotting'};
    toggle_toolbox(toolboxes_to_use, 'on')

    % Initialize counter
    current_action = 1;

    % Log action
    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Image registration'])

    for iid = 1:length(all_animal_names)
        % Get data for this animal
        animal_ID = all_animal_names{iid};
        this_data_type = data_types{iid};
        if ~strcmp(this_data_type, 'epi')
            % Read metadata from database
            METADATA = SQL_database.read_table_where('trials', {'+','animal_ID','date','stimulus','experimental_condition'}, animal_ID, 'animal_ID', 'return_all_sessions',true, 'split_experimental_condition', false);

            % Check the number of rows
            n_trials = height(METADATA);
            if n_trials < 1
                LOGGER.warn(['No trials in the database for ', animal_ID])
                continue
            end
        end

        % Log dataset to work on
        LOGGER.info(['Processing ', animal_ID])
        
        % Create folder
        folder = os.path.join(GC.temp_dir, animal_ID);
        if ~exist(folder, 'dir')
            mkdir(folder)
        end

        switch this_data_type
            case '2p'
                % Perform motion correction by sub-pixel, piecewise-rigid image registration
                has_structural_channel = max(logical(unique(METADATA{:, 'has_structural_channel'})));
                image_registration(METADATA, registration_missing_only, skip_registration);
        
            case 'epi'
                has_structural_channel = false;
                Miniscope_workflow_main(animal_ID, registration_missing_only)
               GC.experiment_name = INFO.selected_experiment;
        end
    end
    
    % Disable toolboxes used in this section
    toggle_toolbox(toolboxes_to_use, 'off')
    
    % Increment counter
    current_action = current_action + 1;
end


%% MLint exceptions
%#ok<*AGROW,*NASGU,*STRNU>
