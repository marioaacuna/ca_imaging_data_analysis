function dataset_response_modulation_heatmap(all_animal_IDs, do_overwrite)

% Get general_configs and LOGGER
global GC LOGGER

% Load file containing single-cell response modulation statistics fot this animal
stats_filename = get_filename_of('response_modulation_stats', GC.experiment_name);
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');

% Read general_configs
group_trials_by = GC.group_trials_by_timepoint.(GC.experiment_name);

%% Individual animals
% Get list of animals to analyze
all_animal_names = fieldnames(STATS_response_modulation);
if isempty(all_animal_IDs)
    all_animal_IDs = all_animal_names;
end
all_animal_IDs = natsort(all_animal_IDs);

% List allowed stimuli and compounds for the analysis
allowed_stimuli = GC.response_modulation_heatmap_stimuli.(GC.experiment_name);
allowed_compounds = GC.analysis_allowed_compounds.(GC.experiment_name);

% Loop through animals
for iid = 1:length(all_animal_IDs)
    % Skip if not completely analyzed
    animal_ID = all_animal_IDs{iid};
    if ~ismember(animal_ID, all_animal_names)
        LOGGER.warn([animal_ID, ' not yet analyzed. Skipping'])
        continue
    end
    
    % Create output folder, if it doesn't exist
    output_folder =  [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep(), 'response_modulation', filesep()];
    if ~exist(output_folder, 'dir')
        mkdir(output_folder)
    end
    
    % Skip if output file exists and user requested not to overwrite previous
    % runs of this script
    output_filename = [output_folder, animal_ID , '_response_heatmap.pdf'];
    if exist(output_filename, 'file') && ~do_overwrite
        LOGGER.info(['Plot for ' animal_ID, ' already exists. Skipping'])
        continue
    end
    
    LOGGER.info(['Processing ', animal_ID])
    % Read metadata and spikes for this animal
    spikes_filename = get_filename_of('spikes', animal_ID);
    spikes = load_variable(spikes_filename, 'spikes');
    METADATA = SQL_database.read_table_where('trials', {'+','keep','experimental_condition','stimulus','date'}, animal_ID,'animal_ID', 'return_all_sessions',true);
    % Remove sessions
    trials_to_remove = METADATA{:, 'keep'} == 0;
    METADATA(trials_to_remove, :) = [];
    spikes(:, trials_to_remove) = [];
    % Get data for this animal only
    STATS = STATS_response_modulation.(animal_ID);
    % Skip not allowed stimuli or compounds
    if ismember('compound', STATS.Properties.VariableNames)
        STATS = STATS(ismember(STATS.compound, allowed_compounds), :);
    end
    % Remove spontaneous data
    STATS = STATS(~ismember(STATS.stimulus, 'SP'), :);
    cells_id = unique(STATS.ROI_id);
    n_all_cells = length(cells_id);
    % Convert timepoint info to numbers
    STATS.timepoints = cellfun(@str2double, STATS.(group_trials_by));
    if all(isnan(STATS.timepoints))
        [~, ~, STATS.timepoints] = unique(STATS.(group_trials_by), 'stable');
        [~, ~, METADATA.timepoints] = unique(METADATA.(group_trials_by), 'stable');
        timepoint_names_are_string = true;
    else
        METADATA.timepoints = cellfun(@str2double, METADATA.(group_trials_by));
        timepoint_names_are_string = false;
    end
    % Remove sessions that user discarded from metadata
    STATS(~ismember(STATS.timepoints, METADATA.timepoints), :) = [];
    timepoint_names = unique(STATS.timepoints);
    
    for istim = 1:length(allowed_stimuli)
        stats = STATS(ismember(STATS.stimulus, allowed_stimuli{istim}), :);
        if height(stats) == 0, continue, end
        rows = ismember(table2cell(METADATA(:, 'stimulus')), allowed_stimuli{istim});
        if ismember('compound', METADATA.Properties.VariableNames)
            rows = rows & ismember(table2cell(METADATA(:, 'compound')), allowed_compounds);
        end
        data = spikes(:, rows);
        metadata = METADATA(rows, :);
        
        % Mark selective cells
        stats.is_selective = stats.p <= 0.05;
        stats.selective = sign(stats.AUROC - 0.5);
    
        % Prepare heatmap data
        trial_duration = unique(metadata{:, 'n_frames'});
        % Get timepoints
        sessions = unique(stats.timepoints);
        n_sessions = length(sessions);
        % Initialize variable
        HEATMAP = zeros(n_all_cells, trial_duration, n_sessions);
        AUROC = zeros(n_all_cells, n_sessions) + 0.5;
        IS_SELECTIVE = zeros(n_all_cells, n_sessions);
        
        for isess = 1:n_sessions
            % Get trials corresponding to this timepoint
            trials_idx = ismember(metadata{:, 'timepoints'}, sessions(isess));

            for iroi = 1:n_all_cells
                % Get activity of this cell                
                activity = data(iroi, trials_idx);
                activity = cell2mat(activity(:));
                % Store mean activity in main variable
                HEATMAP(iroi, :, isess) = mean(activity, 1);
                
                % Get AUROC of this cell
                row = ismember(stats.ROI_id, cells_id(iroi)) & ismember(stats.timepoints, sessions(isess));
                if any(row)
                    % Copy AUROC
                    value = stats{row, 'AUROC'};
                    if ~isnan(value)
                        AUROC(iroi, isess) = value;
                    end
                    % Get whether this AUROC is statistically significant
                    value = stats{row, 'p'};
                    if ~isnan(value)
                        IS_SELECTIVE(iroi, isess) = double(value <= 0.05);
                    end
                end
            end
        end
        
        % Sort cells
        switch GC.response_modulation_heatmap_sort
            case 'none'
                % nothing to do
                
            case 'session_1'
                % Sort according to amount of activation during session 1
                sorting_criteria = [sum(HEATMAP(:, :, 1), 2), AUROC(:, 1), IS_SELECTIVE(:, 1)];
                [~, order] = sortrows(sorting_criteria, [-3, -2, 1]);
                % Sort variables
                HEATMAP = HEATMAP(order, :, :);
                AUROC = AUROC(order, :);
                IS_SELECTIVE = IS_SELECTIVE(order, :);
                
            case 'all'
                for isess = 1:n_sessions
                    sorting_criteria = [sum(HEATMAP(:, :, isess), 2), AUROC(:, isess), IS_SELECTIVE(:, isess)];
                    [~, order] = sort(sorting_criteria, 'descend');
                    % Sort variables
                    HEATMAP(:, :, isess) = HEATMAP(order, :, isess);
                    AUROC(:, isess) = AUROC(order, isess);
                    IS_SELECTIVE(:, isess) = IS_SELECTIVE(order, isess);
                end
        end
        
        % Get stimulus timestamps
        frame_rate = unique(metadata{:, 'frame_rate'});
        n_frames = unique(metadata{:, 'n_frames'});
        STIMULUS = Metadata_workflow.load_stimuli(allowed_stimuli{istim}, 'frame_rate',frame_rate, 'n_frames',n_frames);
        
        % Make output variable to pass to python
        PARAMETERS = struct();
        PARAMETERS.data = HEATMAP;
        PARAMETERS.selectivity = AUROC;
        PARAMETERS.is_selective = IS_SELECTIVE;
        if timepoint_names_are_string
            labels = timepoint_names(:)';
        else
            labels = strcat('day', {' '}, value2str(sessions(:)', '%+i', true));
        end
        PARAMETERS.labels = labels;
        PARAMETERS.title = [animal_ID, ' - ', STIMULUS.code];
        PARAMETERS.output_filename = output_filename;
        PARAMETERS.stimulus_timestamps = STIMULUS.timestamps;
        PARAMETERS.time_axis = STIMULUS.time_axis;
        PARAMETERS.font_size = GC.python.font_size_labels;

        % Store data to disk
        temp_data_filename = [GC.temp_dir, 'heatmap_plot_data.mat'];
        save(temp_data_filename, 'PARAMETERS', '-v6')

        % Call python
        run_in_python('response_heatmap', temp_data_filename)
        printFilePath(PARAMETERS.output_filename, PARAMETERS.output_filename, false)
    end   
end

% Delete temporary file
delete(temp_data_filename)


%% MLint exceptions
%#ok<*AGROW,*STRNU>
