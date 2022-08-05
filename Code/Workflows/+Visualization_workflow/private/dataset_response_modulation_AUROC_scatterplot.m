function dataset_response_modulation_AUROC_scatterplot(all_animal_IDs, do_overwrite, do_summary)

% Get general_configs and LOGGER
global GC LOGGER

% Load file containing single-cell response modulation statistics fot this animal
stats_filename = get_filename_of('response_modulation_stats', GC.experiment_name);
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');

%% Individual animals
% Get list of animals to analyze
all_animal_names = fieldnames(STATS_response_modulation);
if isempty(all_animal_IDs)
    all_animal_IDs = all_animal_names;
end
all_animal_IDs = natsort(all_animal_IDs);

% List allowed stimuli and compounds for the analysis
allowed_stimuli = {'HPS', 'puff'};
n_stimuli = length(allowed_stimuli);
allowed_compounds = GC.analysis_allowed_compounds.(GC.experiment_name);

% List colors for different selectivities
group_colors = {'HPS',[1,0,0]; 'puff',[.05,.48,.75]; 'mixed',[1,.65,0]; 'none',[.7,.7,.7]};

% Make variable for all datasets
if do_summary
    n_sessions = zeros(length(all_animal_IDs), 1);
    n_sessions_around_0 = zeros(length(all_animal_IDs), 2);
    for iid = 1:length(all_animal_IDs)
        animal_ID = all_animal_IDs{iid};
        STATS = STATS_response_modulation.(animal_ID);
        STATS.day_from_surgery = cellfun(@str2double, STATS.day_from_surgery);
        sessions = unique(STATS.day_from_surgery);
        n_sessions(iid) = length(sessions);
        n_sessions_around_0(iid, :) = [length(find(sessions < 0)), length(find(sessions > 0))];
    end
    N_SESSIONS = max(n_sessions);
    n_sessions_around_0 = max(n_sessions_around_0);
    ALL_DATA = cell(0, N_SESSIONS + 1);
    ALL_GROUPS = cell(0, N_SESSIONS + 1);
end

% Loop through animals
for iid = 1:length(all_animal_IDs)
    % Skip if not completely analyzed
    animal_ID = all_animal_IDs{iid};
    if ~ismember(animal_ID, all_animal_names)
        LOGGER.warn([animal_ID, ' not yet analyzed. Skipping'])
        continue
    end
    
    % Skip if output file exists and user requested not to overwrite previous
    % runs of this script
    output_filename = [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep(), 'response_modulation', filesep(), animal_ID , '_stimulus_selectivity_scatterplot.pdf'];
    if exist(output_filename, 'file') && ~do_overwrite
        LOGGER.info(['Plot for ' animal_ID, ' already exists. Skipping'])
        continue
    end

    LOGGER.info(['Processing ', animal_ID])
    % Get data for this animal only
    STATS = STATS_response_modulation.(animal_ID);
    % Skip not allowed stimuli or compounds
    STATS = STATS(ismember(STATS.compound, allowed_compounds), :);
    STATS = STATS(ismember(STATS.stimulus, allowed_stimuli), :);
    
    % Get timepoints
    sessions = sort(cellfun(@str2double, unique(STATS.day_from_surgery)));
    n_sessions = length(sessions);

    % Loop through sessions
    ROI_ids = unique(STATS.ROI_id);
    DATA = repmat({NaN(length(ROI_ids), 2)}, n_sessions, 1);
    GROUPS = repmat({zeros(length(ROI_ids), 1)}, n_sessions, 1);
    for isess = 1:n_sessions
        data = STATS(ismember(STATS.day_from_surgery, num2str(sessions(isess), '%+i')), :);
        P = NaN(length(ROI_ids), 2);
        for row = 1:height(data)
            stimulus = data.stimulus{row};
            if ismember(stimulus, allowed_stimuli)
                DATA{isess}(ismember(ROI_ids, data.ROI_id(row)), ismember(allowed_stimuli, stimulus)) = data.AUROC(row);
                P(ismember(ROI_ids, data.ROI_id(row)), ismember(allowed_stimuli, stimulus)) = data.p(row);
            end
        end
        % Set colors depending on selectivity
        selectivity = P <= 0.05;
        for istim = 1:n_stimuli
            rows = selectivity(:, istim);
            GROUPS{isess}(rows, :) = istim;
        end
        double_selective = all(selectivity');
        GROUPS{isess}(double_selective, :) = n_stimuli + 1;
        GROUPS{isess}(GROUPS{isess} == 0) = n_stimuli + 2;
        
        % Sort points
        [GROUPS{isess}, sorting_order] = sort(GROUPS{isess});
        DATA{isess} = DATA{isess}(sorting_order, :);
    end   
    % Remove NaNs (where there was a problem with the baseline)
    for isess = 1:n_sessions
        bad_rows = any(isnan(DATA{isess})');
        DATA{isess}(bad_rows, :) = [];
        GROUPS{isess}(bad_rows, :) = [];
    end
    % Remove timepoints where comparison cannot be performed
    bad_rows = cellfun(@(x) size(x,1)==0, DATA);
    DATA(bad_rows) = [];
    GROUPS(bad_rows) = [];
    sessions(bad_rows) = [];
    
    % Make plot
    PARAMETERS = struct();
    PARAMETERS.data = DATA;
    PARAMETERS.groups = GROUPS;
    PARAMETERS.group_colors = group_colors;
    PARAMETERS.stimuli = allowed_stimuli;
    PARAMETERS.labels = strcat('day', {' '}, value2str(sessions, '%+i', true));
    PARAMETERS.output_filename = output_filename;
    PARAMETERS.title = animal_ID;
    PARAMETERS.add_n = true;
    PARAMETERS.add_average = false;
    % Store data to disk
    temp_data_filename = [GC.temp_dir, 'scatterplot_data.mat'];
    save(temp_data_filename, 'PARAMETERS', '-v6')
    
    % Call python
    run_in_python('response_modulation_scatterplot_AUROC', temp_data_filename)
    printFilePath(PARAMETERS.output_filename, PARAMETERS.output_filename, false)
    
    % Delete temporary file
    delete(temp_data_filename)

    % Keep data in memory for further processing
    if do_summary
        experimental_group = SQL_database.read_table_where('experiments', 'experimental_group', animal_ID,'animal_ID', 'return_as_table',false);
        ALL_DATA(end+1, [1, find(~bad_rows(:)')+1]) = [{experimental_group}, DATA(:)'];
        ALL_GROUPS(end+1, [1, find(~bad_rows(:)')+1]) = [{experimental_group}, GROUPS(:)'];
    end
end


%% Summary of all datasets
if do_summary
    LOGGER.info('Making summary plots of all cells')
    
    STATS = struct();
    [experimental_groups, ~, group_idx] = unique(ALL_DATA(:,1));
    for igroup = 1:length(experimental_groups)
        % Concatenate data from this group
        data_cell = ALL_DATA(group_idx == igroup, 2:end);
        group_cell = ALL_GROUPS(group_idx == igroup, 2:end);
        data = cell(1, N_SESSIONS);
        groups = cell(1, N_SESSIONS);
        for isess = 1:N_SESSIONS
            data{isess} = cat(1, data_cell{:, isess});
            groups{isess} = cat(1, group_cell{:, isess});
        end
                
        % Make plot
        PARAMETERS = struct();
        PARAMETERS.data = data;
        PARAMETERS.groups = groups;
        PARAMETERS.group_colors = group_colors;
        PARAMETERS.group_colors{end, 2} = [PARAMETERS.group_colors{end, 2}, 0];
        time_step = n_sessions_around_0(1) + 1;
        labels = [value2str(-time_step+1:-1, '%+i', true), value2str(1:N_SESSIONS-time_step+1, '%+i', true)];
        PARAMETERS.labels = strcat('session ', {' '}, labels);
        PARAMETERS.title = experimental_groups{igroup};
        PARAMETERS.stimuli = allowed_stimuli;
        PARAMETERS.output_filename = [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep(), 'response_modulation', filesep(), '_', experimental_groups{igroup}, '_stimulus_selectivity_scatterplot.pdf'];
        % Add other parameters that work better on summary
        PARAMETERS.add_n = true;
        PARAMETERS.add_average = false;
        PARAMETERS.scatterplot_size = 20;
        % Store data to disk
        temp_data_filename = [GC.temp_dir, 'scatterplot_data.mat'];
        save(temp_data_filename, 'PARAMETERS', '-v6')

        % Call python
        run_in_python('response_modulation_scatterplot_AUROC', temp_data_filename)
        printFilePath(PARAMETERS.output_filename, PARAMETERS.output_filename, false)

        % Delete temporary file
        delete(temp_data_filename)

        continue
        % Do statistics
        STATS.(experimental_groups{igroup}).vs_mixed = zeros(N_SESSIONS, 3, n_stimuli);
        STATS.(experimental_groups{igroup}).vs_other_stimuli = zeros(N_SESSIONS, 3);
        for isess = 1:N_SESSIONS
            % Compare single selectivity to mixed selectivity
            selective_for_one_stimulus = cell(n_stimuli, 1);
            for istim = 1:n_stimuli
                selective_for_one_stimulus{istim} = ismember(groups{isess}, istim);
            end
            mixed_selectivity = ismember(groups{isess}, n_stimuli+1);
            mixed_selectivity_AUROC = data{isess}(mixed_selectivity, :);
            for istim = 1:n_stimuli
                s1 = data{isess}(selective_for_one_stimulus{istim}, istim);
                s2 = mixed_selectivity_AUROC(:, istim);
                s1 = s1(s1 >= 0.5);
                s2 = s2(s2 >= 0.5);
                p = bootstrap_ttest([s1; s2], [zeros(length(s1),1); ones(length(s2),1)], false, [], false);
                STATS.(experimental_groups{igroup}).vs_mixed(isess, :, istim) = [median(s1), median(s2), p];
            end
            
            % Compare stimuli against each other
            s1 = data{isess}(selective_for_one_stimulus{1}, 1);
            s2 = data{isess}(selective_for_one_stimulus{2}, 2);
            s1 = s1(s1 >= 0.5);
            s2 = s2(s2 >= 0.5);
            p = bootstrap_ttest([s1; s2], [zeros(length(s1),1); ones(length(s2),1)], false, [], false);
            STATS.(experimental_groups{igroup}).vs_other_stimuli(isess, :) = [median(s1), median(s2), p];
        end
        
        % Compare days per stimulus
        P = NaN(0, 4);
        for istim = 1:n_stimuli
            s1 = [];
            for isess = 1:n_sessions_around_0(1)
                s1 = [s1; data{isess}(ismember(groups{isess}, [istim, n_stimuli+1]), istim)];
            end
            s1 = s1(s1 >= 0.5);
            median_s1 = median(s1);
            for isess = 1:n_sessions_around_0(2)
                session_idx = isess + n_sessions_around_0(1);
                s2 = data{session_idx}(ismember(groups{session_idx}, [istim, n_stimuli+1]), istim);
                s2 = s2(s2 >= 0.5);
                p = bootstrap_ttest([s1; s2], [zeros(length(s1),1); ones(length(s2),1)], false, [], false);
                P(end+1, :) = [istim, median_s1, median(s2), p];
            end
            rows = P(:,1) == istim;
            [~, P(rows,4)] = adjust_Pvalues(P(rows,4));
        end
        STATS.(experimental_groups{igroup}).days_per_stimulus = array2table(P, 'VariableNames', {'stimulus', 'baseline', 'after', 'p'});
    end
    
    % Compare groups
    data = cell(length(experimental_groups), N_SESSIONS);
    groups = cell(length(experimental_groups), N_SESSIONS);
    for igroup = 1:length(experimental_groups)
        % Concatenate data from this group
        data_cell = ALL_DATA(group_idx == igroup, 2:end);
        group_cell = ALL_GROUPS(group_idx == igroup, 2:end);
        for isess = 1:N_SESSIONS
            data{igroup, isess} = cat(1, data_cell{:, isess});
            groups{igroup, isess} = cat(1, group_cell{:, isess});
        end
    end
    
    P = NaN(0, 4);
    for istim = 1:n_stimuli
        % Baseline
        s1 = [];
        s2 = [];
        for isess = 1:n_sessions_around_0(1)
            s1 = [s1; data{1, isess}(ismember(groups{1, isess}, [istim, n_stimuli+1]), istim)];
            s2 = [s2; data{2, isess}(ismember(groups{2, isess}, [istim, n_stimuli+1]), istim)];
        end
        s1 = s1(s1 >= 0.5);
        s2 = s2(s2 >= 0.5);
        p = bootstrap_ttest([s1; s2], [zeros(length(s1),1); ones(length(s2),1)], false, [], false);
        P(end+1, :) = [istim, median(s1), median(s2), p];
        % Days after surgery
        for isess = 1:n_sessions_around_0(2)
            session_idx = isess + n_sessions_around_0(1);
            s1 = data{1, session_idx}(ismember(groups{1, session_idx}, [istim, n_stimuli+1]), istim);            
            s2 = data{2, session_idx}(ismember(groups{2, session_idx}, [istim, n_stimuli+1]), istim);
            s1 = s1(s1 >= 0.5);
            s2 = s2(s2 >= 0.5);
            p = bootstrap_ttest([s1; s2], [zeros(length(s1),1); ones(length(s2),1)], false, [], false);
            P(end+1, :) = [istim, median(s1), median(s2), p];
        end
        rows = P(:,1) == istim;
        [~, P(rows,4)] = adjust_Pvalues(P(rows,4));
    end
    STATS.CCI_vs_sham = array2table(P, 'VariableNames', {'stimulus', 'CCI', 'sham', 'p'});
end

%% MLint exceptions
%#ok<*AGROW,*STRNU>
