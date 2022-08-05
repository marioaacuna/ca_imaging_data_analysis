function dataset_spontaneous_activity_distributions(all_animal_IDs, do_overwrite, do_summary)

% Get general_configs and LOGGER
global GC LOGGER

% Load file containing single-cell response modulation statistics fot this animal
stats_filename = get_filename_of('spontaneous_activity_stats', GC.experiment_name);
DISTRIBUTIONS = load_variable(stats_filename, 'DISTRIBUTIONS');

%% Individual animals
% Get list of animals to analyze
all_animal_names = fieldnames(DISTRIBUTIONS);
if isempty(all_animal_IDs)
    all_animal_IDs = all_animal_names;
end
all_animal_IDs = natsort(all_animal_IDs);

% Make names of columns to keep in the DISTRIBUTION table
attribute_names = GC.spontaneous_activity_attributes;
n_attributes = length(attribute_names);
all_columns = {};
for iid = 1:length(all_animal_IDs)
    if isfield(DISTRIBUTIONS, all_animal_IDs{iid})
        all_columns = [all_columns, DISTRIBUTIONS.(all_animal_IDs{iid}).Properties.VariableNames];
    end
end
all_columns = unique(all_columns);
timepoints = all_columns(~cellfun(@isempty, regexp(all_columns, [attribute_names{1}, '_dist_'])));
timepoints = cellfun(@(x) strsplit(x, [attribute_names{1}, '_dist_']), timepoints, 'UniformOutput',false);
timepoints = cellfun(@(x) ['_dist_', x{2}], timepoints, 'UniformOutput',false);
distribution_columns = {};
for iab = 1:n_attributes
    for ii = 1:length(timepoints)
        distribution_columns{end+1} = [attribute_names{iab}, timepoints{ii}];
    end
end
columns = ['ROI_id', distribution_columns, 'stable_baseline', 'change_after'];
% Get number of timepoints after
n_timepoints_after = length(timepoints);
n_timepoints_after = n_timepoints_after - sum(~cellfun(@isempty, regexp(timepoints, '_before')));
attributes = struct();
if ismember('amplitude', attribute_names)
    attributes.amplitude = GC.spontaneous_amplitude_bins;
end
if ismember('interval', attribute_names)
    attributes.interval = GC.spontaneous_interval_bins;
end

%% Loop through animals
% Make instructions for python
GRAPH_PARAMETERS = struct();
% Write colors
GRAPH_PARAMETERS.group_colors = GC.experimental_groups_colors.(GC.experiment_name);
GRAPH_PARAMETERS.font_size    = GC.python.font_size_labels;
% Plot flags
if strcmp(GC.experiment_name, 'ACC_CCI_anesth')
    GRAPH_PARAMETERS.discard_unstable_baseline = false;
    GRAPH_PARAMETERS.discard_unchanged_after = false;
elseif strcmp(GC.experiment_name, 'ACC_pain_LP-211')
    GRAPH_PARAMETERS.discard_unstable_baseline = true;
    GRAPH_PARAMETERS.discard_unchanged_after = false;
end
GRAPH_PARAMETERS.overlap_timepoints = true;
% Write file
temp_filename = [GC.temp_dir, 'plotting_instructions_spontaneous_activity.mat'];
save(temp_filename, 'GRAPH_PARAMETERS', '-v6');

DATA_all = cell(length(all_animal_IDs), 1);
for iid = 1:length(all_animal_IDs)
    if ~isfield(DISTRIBUTIONS, all_animal_IDs{iid})
        LOGGER.trace([all_animal_IDs{iid}, ' not yet analyzed. Skipping'])
        continue
    end
    LOGGER.trace(['Processing ', all_animal_IDs{iid}])

    columns_this_dataset = DISTRIBUTIONS.(all_animal_IDs{iid}).Properties.VariableNames;
    columns_not_present_in_dataset = setdiff(columns, columns_this_dataset);
    columns_present_in_dataset = columns(ismember(columns, columns_this_dataset));
    % Get info and average distributions
    data = DISTRIBUTIONS.(all_animal_IDs{iid})(:, columns_present_in_dataset);
    for icol = 1:length(columns_not_present_in_dataset)
        this_column_name = columns_not_present_in_dataset{icol};
        this_attribute_column = attribute_names{~cellfun(@isempty, regexp(this_column_name, attribute_names))};
        data(:, this_column_name) = repmat({NaN(0, length(attributes.(this_attribute_column)))}, height(data), 1);
    end
    for jj = 1:length(distribution_columns)
        data{:, distribution_columns{jj}} = cellfun(@(x) mean(x,1).', data{:, distribution_columns{jj}}, 'UniformOutput',false);
    end
    if size(data.change_after, 2) < n_timepoints_after
        data.change_after = [data.change_after, zeros(size(data.change_after,1), n_timepoints_after - size(data.change_after,2))];
    end
    % Re-order columns
    data = data(:, columns);
    % Append animal name and its experimental group
    data(:, 'animal_ID') = {all_animal_IDs{iid}};
    if strcmp(GC.experiment_name, 'ACC_CCI_anesth')
        group_name = SQL_database.read_table_where('experiments', 'experimental_group', all_animal_IDs{iid},'animal_ID', 'return_as_table',false);
    elseif strcmp(GC.experiment_name, 'ACC_pain_LP-211')
        experimental_group = SQL_database.read_table_where('trials', {'experimental_group','experimental_condition'}, all_animal_IDs{iid},'animal_ID');
        experimental_group = experimental_group(:, {'experimental_group','compound'});
        group_name = unique(experimental_group, 'rows', 'stable');
        group_name = table2cell(group_name);
        group_name = strjoin(group_name, ', ');
    end
    data(:, 'experimental_group') = {group_name};
    % Append data to main table
    DATA_all{iid} = data;
end

% Make folder for figures, if it doesn't exist
output_folder = [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep(), 'spontaneous_activity'];
if ~exist(output_folder, 'dir')
    mkdir(output_folder)
end

% Make a plot per dataset
all_animal_IDs_table = cell2table(all_animal_IDs, 'VariableNames',{'animal_ID'});
for iid = 1:length(all_animal_IDs)
    animal_ID = all_animal_IDs_table{iid, 'animal_ID'};
    animal_ID = animal_ID{1};

    output_filename = [output_folder, filesep(), animal_ID, '_CDF.pdf'];
    % Skip if not allowed to overwrite
    if exist(output_filename, 'file') && ~do_overwrite
        LOGGER.info(['Plot for ' animal_ID, ' already exists. Skipping'])
        continue
    end
    
    LOGGER.info(['Plotting ', animal_ID])
    % Slice data for this animal only
    DATA = DATA_all{iid};
    % Get timepoints
    these_columns = DISTRIBUTIONS.(animal_ID).Properties.VariableNames;
    these_timepoints = these_columns(~cellfun(@isempty, regexp(these_columns, [attribute_names{1}, '_dist_'])));
    these_timepoints = cellfun(@(x) strsplit(x, [attribute_names{1}, '_dist_']), these_timepoints, 'UniformOutput',false);
    these_timepoints = cellfun(@(x) ['_dist_', x{2}], these_timepoints, 'UniformOutput',false);
    % Discard timepoints not present in this dataset
    timepoints_to_discard = setdiff(timepoints, these_timepoints);
    if ~isempty(timepoints_to_discard)
        cols_to_discard = [];
        for icol = 1:length(timepoints_to_discard)
            for iat = 1:length(attribute_names)
                column_to_discard = [attribute_names{iat}, timepoints_to_discard{icol}];
                DATA(:, column_to_discard) = [];
            end
            cols_to_discard(end+1) = str2double(timepoints_to_discard{icol}(end));
        end
        DATA.change_after(:, cols_to_discard) = [];
    end
    % Add a grouping index
    [groups, ~, DATA.grouping_index] = unique(DATA.experimental_group);
    % Convert data to a structure that can be opened in python
    DATA = table2struct_of_arrays(DATA);
    DATA.allowed_groups = groups;
    % Add other columns depending on features analyzed
    if ismember('amplitude', attribute_names)
        DATA.attributes.amplitude = GC.spontaneous_amplitude_bins;
    end
    if ismember('interval', attribute_names)
        DATA.attributes.interval = GC.spontaneous_interval_bins;
    end
    % Store data to disk
    temp_data_filename = [GC.temp_dir, 'spontaneous_activity.mat'];
    save(temp_data_filename, 'DATA', '-v6')

    % Run python script
    run_in_python('spontaneous_activity_distributions', temp_filename, temp_data_filename, output_filename)
    printFilePath(output_filename, output_filename, false)
end

%% Summary
if do_summary
    % Make instructions for python
    GRAPH_PARAMETERS = struct();
    % Write colors
    GRAPH_PARAMETERS.group_colors = GC.experimental_groups_colors.(GC.experiment_name);
    GRAPH_PARAMETERS.font_size    = GC.python.font_size_labels;
    % Plot flags
    if strcmp(GC.experiment_name, 'ACC_CCI_anesth')
        GRAPH_PARAMETERS.discard_unstable_baseline = false;
        GRAPH_PARAMETERS.discard_unchanged_after = false;
    elseif strcmp(GC.experiment_name, 'ACC_pain_LP-211')
        GRAPH_PARAMETERS.discard_unstable_baseline = true;
        GRAPH_PARAMETERS.discard_unchanged_after = false;
    end
    GRAPH_PARAMETERS.overlap_timepoints = false;
    % Write file
    temp_filename = [GC.temp_dir, 'plotting_instructions_spontaneous_activity.mat'];
    save(temp_filename, 'GRAPH_PARAMETERS', '-v6');

    % Convert data to a structure that can be opened in python
    if strcmp(GC.experiment_name, 'ACC_CCI_anesth')
        groups = GC.spontaneous_groups_to_plot;
    elseif strcmp(GC.experiment_name, 'ACC_pain_LP-211')
        groups = unique(DATA_all.experimental_group);
    end
    % Get columns for table
    all_columns = {};
    for iid = 1:length(all_animal_IDs)
        if isfield(DISTRIBUTIONS, all_animal_IDs{iid})
            METADATA = SQL_database.read_table_where('experiments', 'experimental_group', all_animal_IDs{iid}, 'animal_ID', 'return_as_table',false);
            if ~ismember(METADATA, groups), continue, end
            all_columns = [all_columns, DISTRIBUTIONS.(all_animal_IDs{iid}).Properties.VariableNames];
        end
    end
    all_columns = unique(all_columns);
    timepoints = all_columns(~cellfun(@isempty, regexp(all_columns, [attribute_names{1}, '_dist_'])));
    timepoints = cellfun(@(x) strsplit(x, [attribute_names{1}, '_dist_']), timepoints, 'UniformOutput',false);
    timepoints = cellfun(@(x) ['_dist_', x{2}], timepoints, 'UniformOutput',false);
    distribution_columns = {};
    for iab = 1:n_attributes
        for ii = 1:length(timepoints)
            distribution_columns{end+1} = [attribute_names{iab}, timepoints{ii}];
        end
    end
    columns = ['ROI_id', 'experimental_group', distribution_columns, 'stable_baseline', 'change_after'];

    % Concatenate data
    DATA_all_redux = table();
    for iid = 1:length(all_animal_IDs)
        METADATA = SQL_database.read_table_where('experiments', 'experimental_group', all_animal_IDs{iid}, 'animal_ID', 'return_as_table',false);
        if ~ismember(METADATA, groups), continue, end
        DATA_all_redux = [DATA_all_redux; DATA_all{iid}];
    end
    DATA_all_redux = DATA_all_redux(:, columns);
    % Add a grouping index
    [groups, ~, DATA_all_redux.grouping_index] = unique(DATA_all_redux.experimental_group);
    % Convert data to a structure that can be opened in python
    DATA = table2struct_of_arrays(DATA_all_redux);
    DATA.allowed_groups = groups;
    if ismember('amplitude', attribute_names)
        DATA.attributes.amplitude = GC.spontaneous_amplitude_bins;
    end
    if ismember('interval', attribute_names)
        DATA.attributes.interval = GC.spontaneous_interval_bins;
    end
    % Store data to disk
    temp_data_filename = [GC.temp_dir, 'spontaneous_activity.mat'];
    save(temp_data_filename, 'DATA', '-v6')

    % Run python
    output_filename = [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep(), 'spontaneous_activity_CDF.pdf'];
    run_in_python('spontaneous_activity_distributions', temp_filename, temp_data_filename, output_filename)
    printFilePath(output_filename, output_filename, false)
end

% Delete temporary files
if exist(temp_filename, 'file')
    delete(temp_filename)
end
if exist(temp_data_filename, 'file')
    delete(temp_data_filename)
end


%% MLint exceptions
%#ok<*AGROW,*STRNU>
