%% PREAMBLE
% This function launches the ROI_selection workflow, which runs a specialized
% GUI for extracting regions-of-interest (ROIs) corresponding to cells or axons.

%% MAIN
function done = main(INFO)
done = true;  % This is the output of the workflow, which is read by the launcher GUI and used to continue to the next workflow
% Get logger and general_configs
global LOGGER GC

% Unpack INFO to make all variables of this script explicit
experiments = INFO.experiments;
data_types = INFO.sessions_table.data_type(ismember(INFO.sessions_table.animal_ID, INFO.experiments));

% Log beginning of workflow
LOGGER.info('ROI selection workflow', 'decorate',true)

% Split experiments by animal
all_animal_names = natsort(experiments);
experiment_name = INFO.selected_experiment; 

for iid = 1:length(all_animal_names)
    % Get data for this animal
    animal_ID = all_animal_names{iid};
    this_data_type = data_types{iid};
    
    LOGGER.info(['Preparing data of ', animal_ID, ' (', this_data_type, '-imaging) for segmentation GUI'])
    
    % Prepare data for segmentation, if they don't exist
    Registration_workflow.prepare_data_for_segmentation(animal_ID, this_data_type, experiment_name)
    
    % Launch python
    LOGGER.info('Launching GUI_ROI_segmentation.py')
    disp('YOU MUST SELECT ONE ROI AND SAVE IT')
    filename_params = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_params.mat']);
    % Load pre existing data
    filename_remote = get_filename_of('ROI_info', animal_ID);
    if exist(filename_remote, 'file')
        selection = MessageBox('Do you want to re run segmentation GUI?', 'Close Request Function', 'YES','No', 'No');
        if strcmp(selection, 'YES')
            run_in_python('GUI_ROI_segmentation', filename_params, false)
        end
    else
        run_in_python('GUI_ROI_segmentation', filename_params, false)
    end
    
    % Process GUI's output -----------------------------------------------------
    % Get filenames
    filename_local = os.path.join(GC.temp_dir, animal_ID, [animal_ID, '_ROI_info.mat']);
    
    % Get number of detected ROIs
    ROIs = load_variable(filename_local, 'ROIs');
    n_ROIs = size(ROIs, 3);
    LOGGER.info([num2str(n_ROIs), ' ROIs detected in dataset ''', animal_ID, ''''])
    % Update database
    LOGGER.trace('Updating database')
    SQL_database.update('experiments', 'n_ROIs',n_ROIs, animal_ID,'animal_ID')
    
    % Ask user to copy file to server
    if exist(filename_remote, 'file')
        DEFAULT = 'Yes, overwrite';
    else
        DEFAULT = 'Yes';
    end
    selection = MessageBox('Copy ROI info to server?', 'Close Request Function', DEFAULT,'No', 'No');
    if strcmp(selection, DEFAULT)
        LOGGER.info('Copying ROI info to server ...')
        copyfile(filename_local, filename_remote)
        LOGGER.info('done', 'append',true)
    end
    DEFAULT = 'YES';
    selection = MessageBox('Do you want to continue with CNMF-e analysis?', 'Close Request Function', DEFAULT,'No', 'No');
    if strcmp(selection, DEFAULT)
        LOGGER.info('Adjust FOV based on manual adjustment')
        toggle_toolbox('CNMF_E', 'on')
%         addpath('C:\Users\acuna\Documents\CNMF_E') % Later change it to path in repo
        done_adjustment = adjust_FOV_from_registration(animal_ID,this_data_type);
        
        if done_adjustment
            selection = MessageBox('Do you want to plot the FOV?', 'Close Request Function', DEFAULT,'No', 'No');
            if strcmp(selection, 'YES')
                LOGGER.info('Checking FOV adjustment')
                test_FOV_correction(animal_ID,this_data_type )
            end
        end
        
        LOGGER.info('About to perform CNFM-e')
        done_CNMFe = Run_CNMFe(animal_ID,this_data_type);
        toggle_toolbox('CNMF_E', 'off')
    end

end


%% MLint exceptions
%#ok<*AGROW,*NASGU,*STRNU>
