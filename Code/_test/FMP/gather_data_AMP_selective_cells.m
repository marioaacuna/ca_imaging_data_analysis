clear, clc, close all

global GC
GC.experiment_name = 'ACC_CCI_anesth';

toggle_toolbox('_plotting', 'on')
FP = figure_properties();

output_folder = [GC.data_root_path, GC.plots_path, 'Figures_paper_FMP', filesep()];

amplitude_type = 'peak';  % Take the difference between peak fluorescence in an event minus the fluorescence at the beginning of the event
% amplitude_type = 'auc';  % Take the sum of the fluorescence in the event
% amplitude_type = 'z';  % z-score of whole trace
% amplitude_type = 'mad';  % Like z-score but it divides median by the MAD
% amplitude_type = 'spikes';  % Deconvolved event amplitude

average_type = 'mean';

stimuli_to_determine_selectivity = {'touch'};  % If more than one stimulus, only cells that have been presented with all stimuli will be analyzed
name = '';

debug_mode = 1;

%% Read file containing the data of interest
disp(['Loading event statistics: SP - ', amplitude_type, ' - ', average_type])
stats_filename = [GC.repository_root_path, '\Figures_paper_FMP\_data\EVENT_stats_SP_', amplitude_type, '.mat'];
EVENT_STATS = load_variable(stats_filename, 'EVENT_STATS');

disp('Loading selectivity')
stats_filename = [GC.repository_root_path, '\Figures_paper_FMP\_data\SELECTIVITY_stats__p.mat'];
SELECTIVITY_stats = load_variable(stats_filename, 'SELECTIVITY_stats');

% Set parameters
experimental_groups_to_analyze = {'sham', 'CCI'};
timepoints_to_analyze = {'S1', 'S2'};

filename_figure = [output_folder, 'fig3_scatter__spontaneous_amplitude_selective_for_', strjoin(stimuli_to_determine_selectivity, '_')];

%% Get data
n_stimuli = length(stimuli_to_determine_selectivity);
disp(['Stimuli to analyze for selectivity: ', strjoin(stimuli_to_determine_selectivity, ', ')])
n_experimental_groups = length(experimental_groups_to_analyze);
n_timepoints = length(timepoints_to_analyze);

filename_figures = strcat(filename_figure, strcat('__', timepoints_to_analyze, '_surgery'), '.pdf');

promoter = {};
promoter(end+1) = {'hSyn'};
promoter(end+1) = {'CaMKII'};
disp(['Analyzing experiments with promoters: ', strjoin(promoter, ' and ')])

animals_to_analyze = animal_list(name);
animal_names = fieldnames(EVENT_STATS);
animal_names = intersect(animal_names, animals_to_analyze);
animals_to_delete = {};
for i_mouse = 1:length(animal_names)
    % Get statistics
    stats = EVENT_STATS.(animal_names{i_mouse});
    if ~isempty(stats), continue, end
    animals_to_delete{end+1} = animal_names{i_mouse};
end

animal_names(ismember(animal_names, animals_to_delete)) = [];

METADATA = SQL_database.read_table_where('experiments', {'animal_ID', 'experimental_group', 'promoter'}, animal_names,'animal_ID');
METADATA = METADATA(ismember(METADATA.animal_ID, animal_names), :);

animal_names = METADATA.animal_ID(ismember(METADATA.promoter, promoter));
METADATA = METADATA(ismember(METADATA.promoter, promoter), :);
n_animals = length(animal_names);

group_trials_by = GC.analysis.(GC.experiment_name).group_trials_by_timepoint;
allowed_compounds = GC.analysis.(GC.experiment_name).analysis_allowed_compounds;%{ismember(GC.analysis.(GC.experiment_name).analysis_allowed_compounds(:,1), GC.experiment_name), 2};

amplitude_column_names = {'animal_ID', 'group', 'cell_ID', 'selectivity_type_S3','selectivity_type_S4', 'S3','S4'};
AMPLITUDE = cell(1, n_timepoints);
for i_mouse = 1:n_animals
    if ~ismember(animal_names{i_mouse}, fieldnames(SELECTIVITY_stats)) || ~ismember(animal_names{i_mouse}, fieldnames(EVENT_STATS)), continue, end
    % Get the experimental group for this animal
    experimental_group = SQL_database.read_table_where('experiments', 'experimental_group', animal_names{i_mouse},'animal_ID', 'return_as_table',false);
    if ~ismember(experimental_group, experimental_groups_to_analyze), continue, end
    
    % Read metadata
    metadata = SQL_database.read_table_where('sessions', {'animal_ID', 'stimulus', 'date', 'experimental_condition'}, animal_names{i_mouse},'animal_ID');
    metadata = metadata(ismember(metadata.compound, allowed_compounds), :);
    if isempty(metadata), continue, end

    % Get trials of interest
    metadata.day_from_surgery = str2double(metadata.day_from_surgery);
    [timepoints, ~, ~] = unique(metadata{:, group_trials_by}, 'stable');
    dates = cell(length(timepoints), 1);
    for i_tp = 1:length(timepoints)
        temp = unique(metadata{ismember(metadata.day_from_surgery, timepoints(i_tp)), 'date'});
        dates{i_tp} = temp{1};
    end
    
    disp([animal_names{i_mouse}, ' - ', experimental_group])
    
    for i_time = 1:n_timepoints
        this_timepoint = timepoints_to_analyze{i_time};
        
        
        sessions_of_interest = find(timepoints > 0 , 2, 'first');
      
       
        good_sessions = [];
        for i_sess = 1:length(sessions_of_interest)
            stimuli_this_session = metadata{ismember(metadata.day_from_surgery, timepoints(sessions_of_interest(i_sess))), 'stimulus'};
            if all(ismember(stimuli_to_determine_selectivity, stimuli_this_session))
                good_sessions(end+1) = sessions_of_interest(i_sess);
            end
        end
       
        sessions_of_interest = good_sessions;
        if isempty(sessions_of_interest), continue, end
        n_sessions = length(sessions_of_interest);
        
        % Get trials of interest
        metadata = SQL_database.read_table_where('trials', {'stimulus', 'experimental_condition', 'date'}, animal_names{i_mouse},'animal_ID');
        metadata = metadata(ismember(metadata.compound, allowed_compounds), :);
        metadata.day_from_surgery = str2double(metadata.day_from_surgery);

        % Get statistics
        selectivity_stats = SELECTIVITY_stats.(animal_names{i_mouse});
        amplitude_stats = EVENT_STATS.(animal_names{i_mouse});
        
        try
            Ca_events_filename = get_filename_of('Ca_events', animal_names{i_mouse});
            keep_cells = load_variable(Ca_events_filename, 'keep_cell');
            keep_cells = logical(keep_cells);
            %         data= (data(keep_cells,:));
        catch
            disp('manual inspection of ROI is missing')
            %         data= (data);
        end
        
        mean_event_amplitude = cell2mat(cellfun(@(x) x.', amplitude_stats{keep_cells, [average_type, '_event_amplitude']}, 'UniformOutput',false));
        mean_event_amplitude(mean_event_amplitude == 0) = NaN;
        n_cells = sum(keep_cells);% height(selectivity_stats);
        
        % Get amplitude according to selectivity of cell
        is_cell_selective = NaN(n_cells, length(timepoints_to_analyze));
        for i_sess = 1:length(sessions_of_interest)
            selectivity_this_day = selectivity_stats{keep_cells, stimuli_to_determine_selectivity}(:, sessions_of_interest(i_sess));
            values_selectivity = zeros(size(selectivity_this_day, 1), 1);
            values_selectivity(selectivity_this_day > .5) = 1;  % responding
            %             values_selectivity(selectivity_this_day < .5) = -1;  % inhibited
            is_cell_selective(:, i_sess) = values_selectivity;
        end
        %         is_cell_selective = abs(is_cell_selective) > 0;
        
        %         session_pairs = nchoosek(1:length(sessions_of_interest), 2);
        %         n_pairs = size(session_pairs, 1);
        %         amplitude = NaN(n_cells * n_pairs, 3);  % 3 columns: selectivity type, amplitude session 1, amplitude session 2
        amplitude = NaN(n_cells ,2);
        session_1 = sessions_of_interest(1);
       
        if n_sessions == 2
            session_2 = sessions_of_interest(2);
            amplitude = mean_event_amplitude(:, [session_1, session_2]);
        elseif n_sessions == 1
            
            amplitude(:,n_sessions) = mean_event_amplitude(:, [session_1]);
        else
            keyboard
        end
        % Store data
        amplitude_metadata = [repmat([{animal_names{i_mouse}, experimental_group}], n_cells , 1), num2cell(repmat((1:n_cells)', 1))];
        AMPLITUDE{i_time} = [AMPLITUDE{i_time}; [amplitude_metadata, num2cell(is_cell_selective), num2cell(amplitude)]];
    end
end

% Convert cell arrays to table
for i_time = 1:n_timepoints
    AMPLITUDE{i_time} = cell2table(AMPLITUDE{i_time}, 'VariableNames',amplitude_column_names);
end
%% Gather data
% select cells only selective in sessions 3 and 4


data_S3 = AMPLITUDE{1};
data_S4 = AMPLITUDE{2};

data_sham_S3 = data_S3.S3(ismember(data_S3.selectivity_type_S3, 1) & ismember(data_S3.group, 'sham'));
data_CCI_S3 =  data_S3.S3(ismember(data_S3.selectivity_type_S3, 1) & ismember(data_S3.group, 'CCI'));

data_sham_S4 = data_S4.S4(ismember(data_S4.selectivity_type_S4, 1) & ismember(data_S3.group, 'sham'));
data_CCI_S4 = data_S4.S4(ismember(data_S4.selectivity_type_S4, 1) & ismember(data_S3.group, 'CCI'));


proportion_CCI_S3 = nansum(data_S3.selectivity_type_S3(ismember(data_S3.group, 'CCI'))) ./ height(data_S3);
proportion_sham_S3 = nansum(data_S3.selectivity_type_S3(ismember(data_S3.group, 'sham'))) ./ height(data_S3);


proportion_CCI_S4 = nansum(data_S4.selectivity_type_S4(ismember(data_S3.group, 'CCI'))) ./ height(data_S4);
proportion_sham_S4 = nansum(data_S4.selectivity_type_S4(ismember(data_S3.group, 'sham'))) ./ height(data_S4);


[h,p] = wilcoxon_ranksum (data_sham_S4, data_CCI_S4)