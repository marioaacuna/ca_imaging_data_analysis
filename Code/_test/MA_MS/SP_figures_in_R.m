%%
clear, clc

global GC
disp('Loading event statistics')

GC.experiment_type = 'ACC_SNI_anxiety';
stats_filename = get_filename_of('spontaneous_activity_stats_epi', GC.experiment_type);
EVENT_STATS = load_variable(stats_filename, 'EVENT_STATS_spontaneous_activi_sessy');
animals = fieldnames(EVENT_STATS);
n_animals = length (animals);
combine_sham_and_CCI = false;

experimental_groups = SQL_database.read_table_where ('experiments', {'animal_ID', 'experimental_group'}, animals, 'animal_ID', 'return_as_table', false);
experimental_conditions = experimental_groups(:, 2);

experimental_groups = unique(experimental_conditions);
%
[~, ~, idx] = unique(experimental_groups(:, 1));
experimental_groups(:, 2) = num2cell(idx);

average_type = 'median';
sessions_to_analyse = 4;
output_folder = GC.temp_dir;% ['M:\Thomas\Manuscripts\in preparation\Invivo ACC\Material for Figures\check_points\'];

disp('Gathering data')

% gather data
AMPLITUDE = cell(n_animals, 1);
FREQUENCY = cell(n_animals, 1);
TIMEPOINTS = cell(n_animals, 3);
TIMEPOINTS(:, 1) = animals;
TIMEPOINTS(:, 2) = experimental_conditions;

%%
for iid = 1: n_animals
    
    animal_ID = animals{iid};
    n_cells = size(EVENT_STATS.(animal_ID),1);
    
    event_stats = EVENT_STATS.(animal_ID);
    if height(event_stats) == 0, continue, end
    METADATA = SQL_database.read_epi_trials (animal_ID);
    dates = unique(METADATA.date, 'stable');

%     dates_EPM = METADATA.date(ismember (METADATA.stimulus, 'EPM'));
    dates_to_take =dates(5:end);
    
%     day_from_surgery = METADATA.day_from_surgery(ismember (METADATA.stimulus, 'EPM'));
    sessions_to_keep = find (ismember(dates, dates_to_take));
    TIMEPOINTS{iid, 3} = dates(sessions_to_keep);
    
    mean_event_amplitude = cell2mat(cellfun(@(x) x.', event_stats{:, [average_type, '_event_ampli_sessude']}, 'UniformOutput',false));
    mean_event_amplitude(mean_event_amplitude == 0) = NaN;
    AMPLITUDE{iid} = mean_event_amplitude(:, sessions_to_keep);

    mean_event_freq =  cell2mat(cellfun(@(x) x.', event_stats{:, ['mean', '_event_rate']}, 'UniformOutput',false));
    FREQUENCY{iid} = mean_event_freq(:, sessions_to_keep);
   
end
disp ('done loading')
%% Plot amplitude and freq
n_columns = sessions_to_analyse;
disp('Re-arranging data for plotting')

% Mark groups
DATA_amplitude = NaN(0, n_columns + 1);
DATA_frequency = NaN(0, n_columns + 1);
for i_mouse = 1:n_animals
    timepoints = TIMEPOINTS{i_mouse, 3};
    if isempty(timepoints), continue, end
    
    amplitudes = NaN(size(AMPLITUDE{i_mouse}, 1), n_columns);
    frequencies = amplitudes;
    %     sessions_before_surgery = find(timepoints < 0);
    
    amplitudes = AMPLITUDE{i_mouse};
    frequencies = FREQUENCY{i_mouse};
    
    
%     sessions_after_surgery = find(timepoints > 0);
%     if length(sessions_after_surgery) == 2
%         amplitudes(:, 3:4) = AMPLITUDE{i_mouse}(:, sessions_after_surgery);
%         frequencies(:, 3:4) = FREQUENCY{i_mouse}(:, sessions_after_surgery);
%     elseif length(sessions_after_surgery) == 1
%         amplitudes(:, 3) = AMPLITUDE{i_mouse}(:, sessions_after_surgery);
%         frequencies(:, 3) = FREQUENCY{i_mouse}(:, sessions_after_surgery);
   
    
    % Store data
    group_idx = experimental_groups{ismember(experimental_groups(:, 1), TIMEPOINTS(i_mouse, 2)), 2};
    DATA_amplitude = [DATA_amplitude; [ones(size(amplitudes, 1), 1)  * group_idx, amplitudes]];
    DATA_frequency = [DATA_frequency; [ones(size(frequencies, 1), 1) * group_idx, frequencies]];
end


column_names = {'day_1', 'day_2', 'day_1_post','day_2_post'};


% group = ones (n_cells,1);
DATA_amplitude_table = array2table(DATA_amplitude, 'VariableNames',['group', column_names]);
DATA_frequency_table = array2table(DATA_frequency, 'VariableNames',['group', column_names]);
DATA_amplitude_table.group = experimental_groups(DATA_amplitude_table{:, 'group'}, 1);
DATA_frequency_table.group = experimental_groups(DATA_frequency_table{:, 'group'}, 1);

% if combine_sham_and_CCI
% 	DATA_amplitude_table{:, 'group'} = {'S1'} ;
%     DATA_frequency_table{:, 'group'} = {'S1'} ;
% end


amplitude_range = table2array(DATA_amplitude_table(:, 2:end));
amplitude_range = amplitude_range(:);
amplitude_range = [nanmin(amplitude_range), nanmax(amplitude_range)];
amplitude_range = [floor(amplitude_range(1) / 10) * 10, ceil(amplitude_range(2) / 10) * 5];
amplitude_range_str = ['"', strjoin(value2str(amplitude_range, '%.2f'), ', '), '"'];
frequency_range = table2array(DATA_frequency_table(:, 2:end));
frequency_range = frequency_range(:);
frequency_range = [nanmin(frequency_range), nanmax(frequency_range)];
frequency_range_str = ['"', strjoin(value2str(frequency_range, '%.5f'), ', '), '"'];
% Cap amplitude range
amplitude_range = [0, 7];
amplitude_range_str = ['"', strjoin(value2str(amplitude_range, '%.1f'), ', '), '"'];
frequency_range = [0, 0.2];
frequency_range_str = ['"', strjoin(value2str(frequency_range, '%.1f'), ', '), '"'];



filename_amplitude = os.path.join(output_folder, 'amplitude.csv');
filename_frequency = os.path.join(output_folder, 'frequency.csv');
writetable(DATA_amplitude_table, filename_amplitude);
writetable(DATA_frequency_table, filename_frequency);
% Run analysis in R

% output_folder = [GC.data_root_path, GC.plots_path, 'miniscope_experiments'];
output_filename = os.path.join(output_folder, 'amplitude_all_MS.pdf');
run_in_R(os.path.join('Figures_paper_FMP', 'amplitude_frequency_distributions_all_days'), ['--file_input ', filename_amplitude], ['--file_output ', output_filename], '--data_type amplitude', ['--data_limits ', amplitude_range_str], '--log_scale FALSE'); printFilePath(output_filename, [], 0)
output_filename = os.path.join(output_folder, 'frequency_MS.pdf');
run_in_R(os.path.join('Figures_paper_FMP', 'amplitude_frequency_distributions_all_days'), ['--file_input ', filename_frequency], ['--file_output ', output_filename], '--data_type frequency', ['--data_limits ', frequency_range_str]); printFilePath(output_filename, [], 0)

