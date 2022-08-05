function dataset_response_modulation_AUROC_alluvial(all_animal_IDs, do_overwrite, do_summary)

% Get general_configs and LOGGER
global GC LOGGER

% Load file containing single-cell response modulation statistics fot this animal
% stats_filename = get_filename_of('response_modulation_stats', GC.experiment_name);
% STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');
STATS_response_modulation = load_variable(os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', 'SELECTIVITY_stats__p.mat'), 'SELECTIVITY_stats');

% animal_ID = all_animal_IDs{1};
% detection_filename = get_filename_of('response_detection', animal_ID);
% results = load_variable(detection_filename, 'RESULTS');

group_trials_by = GC.analysis.(GC.experiment_name).group_trials_by_timepoint;
session_name_prefix = GC.analysis.(GC.experiment_name).session_column_prefix;

%% Individual animals
% Get list of animals to analyze
all_animal_names = fieldnames(STATS_response_modulation);
if isempty(all_animal_IDs)
    all_animal_IDs = all_animal_names;
end
all_animal_IDs = natsort(all_animal_IDs);

% List allowed stimuli and compounds for the analysis
allowed_stimuli = GC.analysis.(GC.experiment_name).analysis_allowed_stimuli;
allowed_compounds = GC.analysis.(GC.experiment_name).analysis_allowed_compounds;

% Make variable for all datasets
if do_summary
    keyboard
    n_sessions = zeros(length(all_animal_IDs), 1);
    n_sessions_around_0 = zeros(length(all_animal_IDs), 2);
    for iid = 1:length(all_animal_IDs)
        animal_ID = all_animal_IDs{iid};
        METADATA = SQL_database.read_table_where('sessions', {'+'}, animal_ID, 'animal_ID');
%         STATS = STATS_response_modulation.(animal_ID);
%         STATS.day_from_surgery = cellfun(@str2double, STATS.day_from_surgery);
        sessions = unique(METADATA.day_from_surgery, 'stable');
%         sessions = unique(STATS.day_from_surgery);
        n_sessions(iid) = length(sessions);
        keyboard % FIX next line
        n_sessions_around_0(iid, :) = [length(find(sessions < 0)), length(find(sessions > 0))];
    end
    N_SESSIONS = max(n_sessions);
    n_sessions_around_0 = max(n_sessions_around_0);
    ALL_DATA = cell(0, N_SESSIONS + 1);
end

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
    output_filename = [output_folder, animal_ID , '_class_identity.pdf'];
    if exist(output_filename, 'file') && ~do_overwrite
        LOGGER.info(['Plot for ' animal_ID, ' already exists. Skipping'])
        continue
    end
    
    LOGGER.info(['Processing ', animal_ID])
    % Get data for this animal only
    STATS = STATS_response_modulation.(animal_ID);
    % Skip not allowed stimuli or compounds
    if ismember('compound', STATS.Properties.VariableNames)
        STATS = STATS(ismember(STATS.compound, allowed_compounds), :);
    end
    STATS = STATS(:, ismember(fieldnames(STATS), allowed_stimuli));
    % Remove spontaneous data
%     STATS = STATS(~ismember(fieldnames(STATS), 'SP'), :);
    % Load Keep cells
    try 
        Ca_events_filename = get_filename_of('Ca_events', animal_ID);
        keep_cells = logical(load_variable(Ca_events_filename, 'keep_cell'));
        catch 
            disp('Animal_not fully analyzed')
    end
    STATS = STATS(keep_cells, :);
    n_all_cells = height(STATS);
%     n_all_cells = length(unique(STATS.ROI_id));
%     % Convert timepoint info to numbers
    % Select session to plot
    n_sessions_to_plot = 3;
    METADATA = SQL_database.read_table_where('sessions', {'+'}, animal_ID, 'animal_ID');
%         STATS = STATS_response_modulation.(animal_ID);
%         STATS.day_from_surgery = cellfun(@str2double, STATS.day_from_surgery);
    sessions = unique(METADATA.day_from_surgery, 'stable');
%         sessions = unique(STATS.day_from_surgery);
   
    timepoint_names = sessions;
    timepoint_names_are_string = false;

%     STATS.timepoints = cellfun(@str2double, STATS.(group_trials_by));
%     if all(isnan(STATS.timepoints))
%         [timepoint_names, ~, STATS.timepoints] = unique(STATS.(group_trials_by), 'stable');
%         timepoint_names_are_string = true;
%     else
%         timepoint_names = unique(STATS.timepoints);
%         timepoint_names_are_string = false;
%     end

    % Keep only significant variations
%     STATS = STATS(STATS.p <= 0.05, :); % already selected
%     n_selective = length(unique(STATS.ROI_id));
    stimuli = allowed_stimuli(ismember(allowed_stimuli, fieldnames(STATS)));
    n_selective = NaN(n_all_cells, length(stimuli));
    for i_stim = 1:length(stimuli) 
        this_stim = stimuli{i_stim};
        n_selective(:, i_stim) = sum(~isnan(STATS.(this_stim)(:,1:n_sessions_to_plot)),2);
    end
     n_nonselective = length(find(sum(n_selective, 2) == 0));
     n_selective = n_all_cells - n_nonselective;
    n_not_selective = n_all_cells - n_selective;
    % Skip if there are no selective cells
%     if n_selective == 0
%         LOGGER.error(['There are no selective cells in ''', animal_ID, ''''])
%         continue
%     end

    % Get timepoints
%     sessions = unique(STATS.timepoints);
    n_sessions = length(sessions);
    % Make all combinations of selectivity across days
%     selectivities = setdiff(allowed_stimuli, 'SP');

    % Make empty cell array
%     DATA = cell(n_selective, n_sessions);  % Last column contains the frequencies
    DATA = cell(n_all_cells, n_sessions_to_plot);  % Last column contains the frequencies
    % Loop through ROIs
    %     ROIs = unique(STATS.ROI_id);
     ROIs = 1:height(STATS);
     ROIs = ROIs';
    stim_to_plot = stimuli(ismember(stimuli, {'HPS', 'puff'}));
    for iroi = 1:n_all_cells%n_selective
        % Get ROI
        ROI_id = ROIs(iroi);
        % Get stimulus to which this ROI responds
%         stats = STATS(ismember(STATS.ROI_id, ROI_id), :);
%         [timepoints, ~, timepoint_idx] = unique(stats.day_from_surgery, 'stable');
        timepoints  = sessions;
%         timepoint_idx = 1:length(sessions);
        % Initialize selectivity of this ROI
        selectivity_ROI = repmat({'none'}, 1, n_sessions_to_plot);
        % Loop through days
        for itimepoint = 1: n_sessions_to_plot%length(timepoints)
%             stats_of_this_timepoint = stats(timepoint_idx == itimepoint, :);
            % Get name of preferred stimuli, and whether they increase or suppress 
            % activity in this cell
            preferred_stimulus = '';
            for istim = 1:length(stim_to_plot)
                this_stim = stim_to_plot{istim};
%                 if ismember(selectivities{istim}, stats_of_this_timepoint.stimulus)
                    AUROC = STATS.(this_stim)(ROI_id,itimepoint);
%                     AUROC = stats_of_this_timepoint{ismember(stats_of_this_timepoint.stimulus, selectivities{istim}), 'AUROC'};
                    if ~isnan(AUROC)                    
                           if AUROC > 0.5
                               preference_sign = '+';
                           else
                               preference_sign = '-';
                           end
                    else
                        continue
                    end
                    % Append a space to separate mixed selectivities
                    if ~strcmp(preferred_stimulus, '')
                        preferred_stimulus = [preferred_stimulus, ''];
                    end
                    % Put it all toghether
                    preferred_stimulus = [preferred_stimulus, stimuli{istim}, preference_sign];
%                 end
            end

            % Store selectivity in main variable
            day_index = ismember(sessions, timepoints(itimepoint));
            selectivity_ROI{day_index} = preferred_stimulus;
        end

        % Store selectivity of this ROI in main variable
        DATA(iroi, :) = selectivity_ROI;
    end

    % Transform strings to numbers for correct sorting of the rows
    data = zeros(size(DATA));
    for stim_idx = 1:length(stim_to_plot)%length(allowed_stimuli)
        data(cellfun(@(x) ~isempty(regexp(x, ['^', stim_to_plot{stim_idx}, '[+-]$'], 'once')), DATA)) = stim_idx;

%         data(cellfun(@(x) ~isempty(regexp(x, ['^', allowed_stimuli{stim_idx}, '[+-]$'], 'once')), DATA)) = stim_idx;
    end
    
    % Non selective cells will be at the end
    data(ismember(DATA, 'none')) = Inf;
    % Mark missing data with -1
    data(ismember(DATA, '')) = -1;
    % Move mixed at the end
    mixed_selectivities = unique(DATA(data == 0));
    if length(mixed_selectivities) == 1
        data(data == 0) = length(allowed_stimuli) + 1;
    else
        for imix = 1:length(mixed_selectivities)
            data(ismember(DATA, mixed_selectivities{imix})) = length(allowed_stimuli) + 1 + imix;
        end
    end
    % Remove incomplete rows, sort them and compute frequencies
    bad_rows = any(data' == -1);
    data(bad_rows, :) = [];
    DATA(bad_rows, :) = [];
    [~, sorting_order] = sortrows(data);
    DATA = DATA(sorting_order, :);
    % Get unique combinations of selectivity among ROIs
    DATA = cell2table(DATA, 'VariableNames',strcat(session_name_prefix, value2str(1:n_sessions_to_plot, '%i', true)));
    [data, ~, ic] = unique(DATA, 'rows');
    freq = accumarray(ic, 1);
    % Add a row for cells that never changed selectivity
    data(height(data)+1, :) = repmat({'none'}, 1, size(data,2));
    freq(end+1) = n_not_selective;
    
    % Make output variable
    PARAMETERS = struct();
    PARAMETERS.data = table2cell(data);
    PARAMETERS.freq = freq;
    PARAMETERS.time_step = 3;%find(sessions > 0, 1, 'first');
    if timepoint_names_are_string
        labels = timepoint_names(:)';
    else
        % labels = strcat('day', {' '}, value2str(sessions(:)', '%+i', true));
        labels = strcat('day', {' '}, value2str(1:n_sessions_to_plot(:)'));
    end    
    PARAMETERS.labels = labels;
    PARAMETERS.title = animal_ID;
    PARAMETERS.output_filename = output_filename;
    PARAMETERS.group_colors = {'HPS',[1,0,0]; 'puff',[.05,.48,.75]; 'mixed',[1,.65,0]; 'none',[.7,.7,.7]};

    % Store data to disk
    temp_data_filename = [GC.temp_dir, '\alluvial_plot_data.mat'];
    save(temp_data_filename, 'PARAMETERS', '-v6')

    % Call python
    run_in_python('alluvial_plot', temp_data_filename)
    printFilePath(PARAMETERS.output_filename, PARAMETERS.output_filename, false)
    
    % Delete temporary file
    delete(temp_data_filename)
    
    % Keep data for further processing
    if do_summary
        data_cell = table2cell(data);
        data2save = repmat({''}, size(data_cell, 1), N_SESSIONS);
        % Save data before 0
        target_columns = n_sessions_around_0(1) - sum(sessions < 0) + 1;
        data2save(:, target_columns:n_sessions_around_0(1)) = data_cell(:, sessions < 0);
        % Save data after 0
        target_columns = n_sessions_around_0(2) + sum(sessions > 0) + 1;
        data2save(:, n_sessions_around_0(1) + 1:target_columns) = data_cell(:, sessions > 0);
        % Remove last row and replace it with actual cells that are never selective
        data2save = [data2save(1:end-1, :); repmat({'none'}, n_not_selective, size(data2save,2))];
        % Store cell selectivity together with the experimental group to which they belong
        group_name = SQL_database.read_table_where('experiments', {'experimental_group'}, animal_ID,'animal_ID', 'return_as_table',false);
        ALL_DATA = [ALL_DATA; [repmat({group_name}, size(data2save,1), 1), data2save]];
    end
end


%% Summary of all datasets
if do_summary
    LOGGER.info('Making summary plots of all cells')
    
    [groups, ~, group_idx] = unique(ALL_DATA(:,1));
    for igroup = 1:length(groups)
        % Get data
        data_cell = ALL_DATA(group_idx == igroup, 2:end);
        data = zeros(size(data_cell));
        % Transform strings to numbers for correct sorting of the rows
        for stim_idx = 1:length(allowed_stimuli)
            data(cellfun(@(x) ~isempty(regexp(x, ['^', allowed_stimuli{stim_idx}, '[+-]$'], 'once')), data_cell)) = stim_idx;
        end
        % Non selective cells will be at the end
        data(ismember(data_cell, 'none')) = Inf;
        % Mark missing data with -1
        data(ismember(data_cell, '')) = -1;
        % Move mixed at the end
        mixed_selectivities = unique(data_cell(data == 0));
        if length(mixed_selectivities) == 1
            data(data == 0) = length(allowed_stimuli) + 1;
        else
            for imix = 1:length(mixed_selectivities)
                data(ismember(data_cell, mixed_selectivities{imix})) = length(allowed_stimuli) + 1 + imix;
            end
        end
        % Remove incomplete rows, sort them and compute frequencies
        bad_rows = any(data' == -1);
        data(bad_rows, :) = [];
        data_cell(bad_rows, :) = [];
        [~, sorting_order] = sortrows(data);
        data_cell = data_cell(sorting_order, :);
        [data, ~, ic] = unique(cell2table(data_cell), 'rows', 'stable');
        freq = accumarray(ic, 1);

        % Make output variable
        PARAMETERS = struct();
        PARAMETERS.add_n = false;  % There will be too many rows
        PARAMETERS.split_flows = false;  % There will be too many rows
        PARAMETERS.data = table2cell(data);
        PARAMETERS.freq = freq;
        time_step = n_sessions_around_0(1) + 1;
        labels = [value2str(-time_step+1:-1, '%+i', true), value2str(1:n_sessions-time_step+1, '%+i', true)];
        PARAMETERS.time_step = time_step;
        PARAMETERS.labels = strcat('session ', {' '}, labels);
        PARAMETERS.title = groups{igroup};
        PARAMETERS.output_filename = [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep(), 'response_modulation', filesep(), '_', groups{igroup} , '_class_identity.pdf'];
        PARAMETERS.group_colors = {'HPS',[1,0,0]; 'puff',[.05,.48,.75]; 'mixed',[1,.65,0]; 'none',[.7,.7,.7]};
        % Store data to disk
        temp_data_filename = [GC.temp_dir, 'alluvial_plot_data.mat'];
        save(temp_data_filename, 'PARAMETERS', '-v6')

        % Call python
        run_in_python('alluvial_plot', temp_data_filename)
        printFilePath(PARAMETERS.output_filename, PARAMETERS.output_filename, false)

        % Delete temporary file
        delete(temp_data_filename)
    end
end

%% MLint exceptions
%#ok<*AGROW,*STRNU>
