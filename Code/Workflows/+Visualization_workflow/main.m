%% PREAMBLE
% This function launches the routines to plot deconvolved spikes.

%% MAIN
function done = main(INFO)

% Allocate output variable
done = true;

% Get logger and general_configs
global LOGGER GC

% Unpack INFO to make all variables of this script explicit.
do_response_detection                 = INFO.actions.visualization.response_detection.main;
do_response_detection_only_selected   = INFO.actions.visualization.response_detection.selected;
do_response_detection_summary         = INFO.actions.visualization.response_detection.summary;
do_response_detection_overwrite       = INFO.actions.visualization.response_detection.overwrite;
do_response_decoding                  = INFO.actions.visualization.response_decoding.main;
do_response_decoding_only_selected    = INFO.actions.visualization.response_decoding.selected;
do_response_decoding_summary          = INFO.actions.visualization.response_decoding.summary;
do_response_decoding_overwrite        = INFO.actions.visualization.response_decoding.overwrite;
do_spontaneous_activity               = INFO.actions.visualization.spontaneous_activity.main;
do_spontaneous_activity_only_selected = INFO.actions.visualization.spontaneous_activity.selected;
do_spontaneous_activity_summary       = INFO.actions.visualization.spontaneous_activity.summary;
do_spontaneous_activity_overwrite     = INFO.actions.visualization.spontaneous_activity.overwrite;
n_actions_to_perform = sum([do_spontaneous_activity, do_response_detection*5, do_response_decoding]);
if n_actions_to_perform < 1, return, end

% Make filename of temporary file with instructions to pass to python
temp_filename = [GC.temp_dir, 'plotting_instructions.mat'];
temp_figure_filename = [GC.temp_dir, 'figure.pdf'];
    
% Log beginning of workflow
LOGGER.info('Visualization workflow', 'decorate',true)

% Initialize counter
current_action = 1;

%% MAKE PLOTS
if do_response_detection
    % Get list of animals to analyze
    if do_response_detection_only_selected
        all_animal_IDs = INFO.experiments;
    else
        all_animal_IDs = {};  % all
    end
    
%     % Heatmap of responses
%     LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot cell responses as heatmaps'])
%     try    
%         dataset_response_detection_heatmap(all_animal_IDs, do_response_detection_overwrite)
%     catch ME
%         LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
%     end
%     current_action = current_action + 1;
% 
%     % Response detection
%     LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot response detection heatmaps'])
%     try
%         dataset_response_detection_heatmap(all_animal_IDs, do_response_detection_overwrite)
%     catch ME
%         LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
%     end
%     current_action = current_action + 1;
%     
    % Alluvial plots
    switch INFO.selected_experiment
        case 'ACC_CCI_anesth'
            LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot cells selectivity across sessions'])
            try
                dataset_response_modulation_AUROC_alluvial(all_animal_IDs, do_response_detection_overwrite, do_response_detection_summary, INFO.selected_experiment)
            catch ME
                LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
            end
            current_action = current_action + 1;
        case 'ACC_SNI_anxiety'
            LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot cells selectivity across sessions'])
            try
%                 dataset_response_modulation_AUROC_alluvial_epi(all_animal_IDs, do_response_detection_overwrite, do_response_detection_summary)
%                 dataset_response_modulation_AUROC_alluvial(all_animal_IDs, do_response_detection_overwrite, do_response_detection_summary,INFO.selected_experiment)
                response_detection_plot (all_animal_IDs, do_response_detection_overwrite, INFO.selected_experiment)
            catch ME
                LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
            end
            current_action = current_action + 1;
        case 'CLA_pain'
            LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot cells selectivity across sessions'])
            try
                dataset_response_modulation_AUROC_alluvial(all_animal_IDs, do_response_detection_overwrite, do_response_detection_summary,INFO.selected_experiment)
            catch ME
                LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
            end
            current_action = current_action + 1;

    end
%     % Traces (and AUROC)
%     LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot response traces and AUROC per cell'])
%     try
%         single_cell_response_detection_AUROC(all_animal_IDs, do_response_detection_overwrite)
% 	  catch ME
%         LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
%     end
%     current_action = current_action + 1;
%     
%     % Traces (and AUC)
%     LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot response traces and AUC per cell'])
%     try
%         single_cell_response_detection_AUC(all_animal_IDs, do_response_detection_overwrite)
%     catch ME
%         LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
%     end
%     current_action = current_action + 1;
% 
%     % Scatterplot of AUROC
%     if ~strcmp(GC.experiment_name, 'ACC_pain_LP-211')
%         LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot AUROC to each stimulus'])
%         try
%             dataset_response_detection_AUROC_scatterplot(all_animal_IDs, do_response_detection_overwrite, do_response_detection_summary)
%         catch ME
%             LOGGER.critical(['dataset_response_detection_heatmap\n\n', ME.message])
%         end
%         current_action = current_action + 1;
%     end
end


if do_response_decoding
    % Get list of animals to analyze
    if do_response_decoding_only_selected
        all_animal_IDs = INFO.experiments;
    else
        all_animal_IDs = {};  % all
    end
    
    % Overlap decoding stimulus types
    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot cell responses as heatmaps'])
    decoding_overlap_confusion_matrices(all_animal_IDs, do_response_decoding_overwrite)
    current_action = current_action + 1;
end


if do_spontaneous_activity
    % Get list of animals to analyze
    if do_spontaneous_activity_only_selected
        all_animal_IDs = INFO.experiments;
    else
        all_animal_IDs = {};  % all
    end

    % Cumulative distribution functions of amplitude and frequency
    LOGGER.info(['Step ', num2str(current_action), '/', num2str(n_actions_to_perform), ': Plot spontaneous activity'])
    dataset_spontaneous_activity_distributions(all_animal_IDs, do_spontaneous_activity_overwrite, do_spontaneous_activity_summary)
    current_action = current_action + 1;
end


%% MLint exceptions
%#ok<*AGROW,*NASGU,*STRNU>
