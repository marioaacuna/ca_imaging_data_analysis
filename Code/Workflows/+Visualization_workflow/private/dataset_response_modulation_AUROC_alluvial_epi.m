function dataset_response_modulation_AUROC_alluvial_epi(all_animal_IDs, do_overwrite, do_summary)

% Get general_configs and LOGGER
global GC LOGGER
GC.experiment_name = 'ACC_SNI_anxiety';
% Load file containing single-cell response modulation statistics fot this animal
% stats_filename = get_filename_of('response_detection_epi', GC.ex);
% STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');
% group_trials_by = GC.group_trials_by_timepoint.(GC.experiment_name);
% session_name_prefix = GC.session_name_prefix.(GC.experiment_name);

%% Individual animals
% Get list of animals to analyze
all_animal_names =(all_animal_IDs);
if isempty(all_animal_IDs)
    all_animal_IDs = all_animal_names;
end
all_animal_IDs = natsort(all_animal_IDs);

% % List allowed stimuli and compounds for the analysis
% allowed_stimuli = GC.analysis_allowed_stimuli.(GC.experiment_name);
% allowed_compounds = GC.analysis_allowed_compounds.(GC.experiment_name);

% Make variable for all datasets
if do_summary
    n_sessions = zeros(length(all_animal_IDs), 1);
    n_sessions_around_0 = zeros(length(all_animal_IDs), 2);
    for iid = 1:length(all_animal_IDs)
        animal_ID = all_animal_IDs{iid};
        detection_filename = get_filename_of('response_detection_epi', animal_ID);
        detection_results = load_variable(detection_filename, 'RESULTS');
        session_names = fieldnames(detection_results);
        n_sessions (iid) = length(session_names);
        C = cell(n_sessions(iid),1);
        M = zeros(n_sessions(iid),1);
        P = zeros(n_sessions(iid),1);
        for i_sess = 1: n_sessions(iid)
            C {i_sess} = strsplit(char(session_names(i_sess)), '__');
            M (i_sess) = any(startsWith(C{i_sess}, 'minus'), 2);
            P (i_sess) = any(startsWith(C{i_sess}, 'plus'), 2);
        end   
        n_sessions_around_0(iid, :) = [sum(M), sum(P)];
    end
    N_SESSIONS = max(n_sessions);
    n_sessions_around_0 = max(n_sessions_around_0);
    ALL_DATA = cell(0, N_SESSIONS + 1);
end
%%
% Loop through animals
for iid = 1:length(all_animal_IDs)
    % Skip if not completely analyzed
    animal_ID = all_animal_IDs{iid};
    if ~ismember(animal_ID, all_animal_names)
        LOGGER.warn([animal_ID, ' not yet analyzed. Skipping'])
        continue
    end
    
    % Create output folder, if it doesn't exist
    output_folder =  os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name,[ 'response_modulation', filesep()]);
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
    METADATA = SQL_database.read_epi_trials(animal_ID);
%     % Get data for this animal only
%     STATS = STATS_response_modulation.(animal_ID);
%     % Skip not allowed stimuli or compounds
%     if ismember('compound', STATS.Properties.VariableNames)
%         STATS = STATS(ismember(STATS.compound, allowed_compounds), :);
%     end
%     STATS = STATS(ismember(STATS.stimulus, allowed_stimuli), :);
%     % Remove spontaneous data
%     STATS = STATS(~ismember(STATS.stimulus, 'SP'), :);
%     n_all_cells = length(unique(STATS.ROI_id));
%     % Convert timepoint info to numbers
%     STATS.timepoints = cellfun(@str2double, STATS.(group_trials_by));
%     if all(isnan(STATS.timepoints))
%         [timepoint_names, ~, STATS.timepoints] = unique(STATS.(group_trials_by), 'stable');
%         timepoint_names_are_string = true;
%     else
%         timepoint_names = unique(STATS.timepoints);
%         timepoint_names_are_string = false;
%     end
%     % Keep only significant variations
%     STATS = STATS(STATS.p <= 0.05, :);
%     n_selective = length(unique(STATS.ROI_id));
%     n_not_selective = n_all_cells - n_selective;
    sessions_metadata = unique(METADATA.day_from_surgery(ismember(METADATA.type, {'evoked'}), :), 'stable');
    results_filename = get_filename_of('response_detection_epi', animal_ID);
%     results_filename = get_filename_of('response_detection_p_0_01', animal_ID, GC.experiment_name);
    RESULTS = load_variable (results_filename, 'RESULTS');
    sessions = fieldnames(RESULTS);
    n_sessions = length(sessions);
    stimuli = {'cold', 'heat'}; 
    STATS = cell (length(stimuli), n_sessions);
    selectivity = cell (length(stimuli), 1);
    n_cells = SQL_database.read_table_where('experiments', {'n_ROIs'}, animal_ID,'animal_ID', 'return_as_table', false);
    for i_stim = 1: length(stimuli)
        STATS {i_stim}  = zeros (n_cells, n_sessions);
        for i_sess = 1:n_sessions
            selectivity_this_sess = RESULTS.(sessions{i_sess}).(stimuli{i_stim}).selectivity.excited;
            STATS {i_stim, i_sess} = selectivity_this_sess ;
        end
        selectivity{i_stim} = horzcat(STATS {i_stim, :});
    end
    
    % Make empty cell array
    ROIs = num2cell(1:1:n_cells)';
    data_cold = selectivity{1};
    data_heat = selectivity{2};
    data = sum(cat(3,data_heat, data_cold), 3);
    data_cold = num2cell (data_cold);
    data_cold ((cellfun(@(x) isempty(setxor(1,x)), data_cold))) = {'cold'};
    data_cold ((cellfun(@(x) isempty(setxor(0,x)), data_cold))) = {'none'};
    data_heat = num2cell (data_heat);
    data_heat ((cellfun(@(x) isempty(setxor(1,x)), data_heat))) = {'heat'};
    data_heat ((cellfun(@(x) isempty(setxor(0,x)), data_heat))) = {'none'};
    data_cold = cat(2, ROIs, data_cold);
    data_heat = cat(2, ROIs, data_heat);
    data_all = cat (1, data_cold,data_heat);
    data_all = sortrows(data_all);
    col_names = {'ROIs', sessions{:}};
    data_all = cell2table (data_all, 'VariableNames', col_names);

    DATA = cell(n_cells, n_sessions);  % Last column contains the frequencies
   

    % Loop through ROIs
    ROIs = unique(data_all.ROIs);
    for iroi = 1 : n_cells
        % Get ROI
        ROI_id = ROIs(iroi);
        % Get stimulus to which this ROI responds
        stats = data_all(ismemberCellRows(data_all.ROIs, ROI_id),:);%STATS(ismember(STATS.ROI_id, ROI_id), :);
%         [timepoints, ~, timepoint_idx] = unique(stats.timepoints, 'stable');

        % Initialize selectivity of this ROI
        selectivity_ROI = repmat({'none'}, 1, n_sessions);
        timepoints = n_sessions;
        % Loop through days
        for itimepoint = 1:timepoints
            stats_of_this_timepoint = stats.(col_names{itimepoint + 1});
            % Get name of preferred stimuli, and whether they increase or suppress 
            % activity in this cell
            preferred_stimulus = '';
            
%             is_selective = ~ismember(stats_of_this_timepoint(:) ,'none');
%             char(stats_of_this_timepoint(is_selective))
            
            for i_row = 1:length(stimuli)
                 if ~ismember(stats_of_this_timepoint(i_row) ,'none') 
                     preferred_this_stimulus = [stats_of_this_timepoint{i_row}, '+'];
%                     if ~ismember(stats_of_this_timepoint(:) ,'none') 
%                         preferred_this_stimulus = [stimuli{i_row}, '+'];
                 else
                     preferred_this_stimulus = '';
                 end

                    % Append a space to separate mixed selectivities
                    if ~strcmp(preferred_stimulus, '')
                        preferred_stimulus = [preferred_stimulus, ''];
                    end
                    % Put it all toghether
                    preferred_stimulus = [preferred_this_stimulus, preferred_stimulus];
                
            end

            % Store selectivity in main variable
%             day_index = ismember(sessions, sessions(itimepoint))';
            selectivity_ROI{itimepoint} = preferred_stimulus;
        end

        % Store selectivity of this ROI in main variable
        DATA(iroi, :) = selectivity_ROI;
    end
    DATA(ismember(DATA, {''})) = {'none'};
% 
%     % Transform strings to numbers for correct sorting of the rows
%     data = zeros(size(DATA));
%     for stim_idx = 1:length(allowed_stimuli)
%         data(cellfun(@(x) ~isempty(regexp(x, ['^', allowed_stimuli{stim_idx}, '[+-]$'], 'once')), DATA)) = stim_idx;
%     end
%     % Non selective cells will be at the end
%     data(ismember(DATA, 'none')) = Inf;
%     % Mark missing data with -1
%     data(ismember(DATA, '')) = -1;
%     % Move mixed at the end
    mixed_selectivities = unique(DATA(data == 2));
%     if length(mixed_selectivities) == 1
%         data(data == 2) = length(stimuli) + 1;
%     else
%         for imix = 1:length(mixed_selectivities)
%             data(ismember(DATA, mixed_selectivities{imix})) = length(allowed_stimuli) + 1 + imix;
%         end
%     end
%     % Remove incomplete rows, sort them and compute frequencies
%     bad_rows = any(data' == -1);
%     data(bad_rows, :) = [];
%     DATA(bad_rows, :) = [];
    [~, sorting_order] = sortrows(data);
    DATA = DATA (sorting_order, :);
    % Get unique combinations of selectivity among ROIs
    DATA = cell2table(DATA, 'VariableNames',strcat(sessions));
%     data = DATA{1};
    
    [data_x, ~, ic] = unique(DATA, 'rows');
    freq = accumarray(ic, 1);
    % Add a row for cells that never changed selectivity
    data_x(height(data_x)+1, :) = repmat({'none'}, 1, size(data_x,2));
    freq(end+1) = 0; 
    
    % Make output variable
    PARAMETERS = struct();
    PARAMETERS.data = table2cell(data_x);%table2cell(data);
    PARAMETERS.freq = freq;
    PARAMETERS.time_step = 0;%find(P, 1, 'first');
    timepoint_names_are_string = 0;
    if timepoint_names_are_string
        labels = {sessions_metadata(:)'};
    else
%         labels = strcat('day', {' '}, value2str([sessions_metadata{:}], '%+i', true));
        labels = strcat('day', {' '}, char(sessions_metadata{:}));
    end    
    PARAMETERS.labels = labels;
    PARAMETERS.title = animal_ID;
    PARAMETERS.output_filename = output_filename;
    PARAMETERS.group_colors = {'cold+',[0,0,1]; 'heat+',[1,0,0]; 'mixed',[1,.65,0]; 'none',[.7,.7,.7]};

    % Store data to disk
    temp_data_filename = [GC.temp_dir, '\', 'alluvial_plot_data.mat'];
    save(temp_data_filename, 'PARAMETERS', '-v6')

    % Call python
    run_in_python('alluvial_plot', temp_data_filename)
    printFilePath(PARAMETERS.output_filename, PARAMETERS.output_filename, false)
    
    % Delete temporary file
    delete(temp_data_filename)
    
    % Keep data for further processing
    if do_summary
        data_cell = table2cell(data_x);
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
