function gather_data_AMP_times_FREQ()

clear, clc
global GC
load_repository();
GC.experiment_name = 'ACC_CCI_anesth';
% output_folder = os.path.join(GC.data_root_path, GC.plots_path, 'Figures_paper_FMP');
output_folder = GC.temp_dir;% ['M:\Thomas\Manuscripts\in preparation\Invivo ACC\Material for Figures\check_points\'];
combine_sham_and_CCI = false;
discard_naive_from_session_comparisons = true;
% use_all_events = true;
amplitude_type = 'peak';



% Read file containing the data of interest
disp('Loading event statistics')
% if use_all_events
%     stats_filename = os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', ['EVENT_stats_', amplitude_type, '.mat']);
%     EVENT_STATS = load_variable(stats_filename, 'EVENT_STATS');
%     column = 'mean_event_amplitude';
% else
%     stats_filename = get_filename_of('spontaneous_activity_stats', GC.experiment_name);
%     EVENT_STATS = load_variable(stats_filename, 'EVENT_STATS_spontaneous_activity');
% end

data_type = 'old';
% 
if strcmp(data_type, 'new')
    column = 'sum_event_amplitude';
    stats_filename = os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', 'EVENT_stats_SP_auc_new.mat');    
else
    column = 'mean_event_amplitude';
    stats_filename = os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', 'EVENT_stats_SP_AUC_old.mat');
end
EVENT_STATS = load_variable(stats_filename, 'EVENT_STATS');

animal_names = animal_list();
animals_event_stats = fieldnames(EVENT_STATS);
n_animals = length(animal_names);

% % Load Selectivity
% stimuli_to_determine_selectivity = {'temp_48'};
% input_filename = os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', 'SELECTIVITY_stats__p.mat');
% SELECTIVITY_stats = load_variable(input_filename, 'SELECTIVITY_stats');


experimental_groups = SQL_database.read_table_where('experiments', {'animal_ID', 'experimental_group'}, animal_names,'animal_ID', 'return_as_table',false);
experimental_groups = natsortrows(experimental_groups, 1);
% animal_names = experimental_groups(:, 1);
experimental_conditions = experimental_groups(:, 2);
experimental_groups = unique(experimental_conditions);
% Add conde of each group
[~, ~, idx] = unique(experimental_groups(:, 1));
experimental_groups(:, 2) = num2cell(idx);

group_trials_by = GC.analysis.ACC_CCI_anesth.group_trials_by_timepoint;
allowed_compounds =GC.analysis.ACC_CCI_anesth.analysis_allowed_compounds;

disp('Gathering data')

AMPLITUDE = cell(n_animals, 1);
FREQUENCY = cell(n_animals, 1);
TIMEPOINTS = cell(n_animals, 3);
TIMEPOINTS(:, 1) = animal_names;
TIMEPOINTS(:, 2) = experimental_conditions;

for iid = 1:n_animals
    % Get statistics
    animal_ID = (animal_names{iid});
    stats = EVENT_STATS.(animal_names{iid});
    if height(stats) == 0, continue, end

    % Read metadata
    METADATA = SQL_database.read_table_where('trials', {'experimental_condition'}, animal_names{iid},'animal_ID');
    METADATA = METADATA(ismember(METADATA.compound, allowed_compounds), :);
    if isempty(METADATA), continue, end
%     analyzed_stimuli = SELECTIVITY_stats.(animal_ID).Properties.VariableNames(:);
%     analyzed_stimuli(ismember(analyzed_stimuli, 'cell')) = [];
%     analyzed_stimuli(ismember(analyzed_stimuli, 'temp_42')) = {'temp_43'};
%     if ~all(ismember(stimuli_overlap_converted, analyzed_stimuli)), continue, end
%     sessions{:, 'animal_ID'} = {animal_ID};

    [timepoints, ~, ~] = unique(METADATA{:, group_trials_by}, 'stable');
    timepoints = str2double(timepoints);
%     sessions_to_keep_idx = find(timepoints < 0, 2, 'last');%[find(timepoints < 0, 2, 'last'); find(timepoints > 0, 2, 'first')];
    sessions_to_keep_idx = [find(timepoints < 0, 2, 'last'); find(timepoints > 0, 2, 'first')];
    if isempty(sessions_to_keep_idx), continue, end
%     sessions_to_keep_idx = sessions_to_keep_idx(1);
    TIMEPOINTS{iid, 3} = timepoints(sessions_to_keep_idx);

    % Mean event amplitude
    
    Ca_events_filename = get_filename_of('Ca_events', animal_ID);
    keep_cells = load_variable(Ca_events_filename, 'keep_cell');
    keep_cells = logical(keep_cells);
    
    mean_event_amplitude = cell2mat(cellfun(@(x) x(:).', stats{keep_cells, column}, 'UniformOutput',false));
    AMPLITUDE{iid} = mean_event_amplitude(:, sessions_to_keep_idx);
    
    % Mean event frequency
    mean_event_frequency = cell2mat(cellfun(@(x) x(:).', stats{keep_cells, 'mean_event_rate'}, 'UniformOutput',false));
    FREQUENCY{iid} = mean_event_frequency(:, sessions_to_keep_idx);
end

disp ('done loading')
%% Plot amplitude and frequency on day 1
% Day 1 is defined as the first session before surgery
disp('Re-arranging data for plotting')

n_columns = 4; % for one day 1

% Mark groups
DATA_amplitude = NaN(0, n_columns + 1);
DATA_frequency = NaN(0, n_columns + 1);
for i_mouse = 1:n_animals
    timepoints = TIMEPOINTS{i_mouse, 3};
    if isempty(timepoints), continue, end
    
    amplitudes = NaN(size(AMPLITUDE{i_mouse}, 1), n_columns);
    frequencies = amplitudes;
    sessions_before_surgery = find(timepoints < 0);
    if length(sessions_before_surgery) == 2
        amplitudes(:, 1:2) = AMPLITUDE{i_mouse}(:, sessions_before_surgery);
        frequencies(:, 1:2) = FREQUENCY{i_mouse}(:, sessions_before_surgery);
    elseif length(sessions_before_surgery) == 1
        amplitudes(:, 1) = AMPLITUDE{i_mouse}(:, sessions_before_surgery(1));
        frequencies(:, 1) = FREQUENCY{i_mouse}(:, sessions_before_surgery(1));
    end
    
    sessions_after_surgery = find(timepoints > 0);
    if length(sessions_after_surgery) == 2
        amplitudes(:, 3:4) = AMPLITUDE{i_mouse}(:, sessions_after_surgery);
        frequencies(:, 3:4) = FREQUENCY{i_mouse}(:, sessions_after_surgery);
    elseif length(sessions_after_surgery) == 1
        amplitudes(:, 3) = AMPLITUDE{i_mouse}(:, sessions_after_surgery);
        frequencies(:, 3) = FREQUENCY{i_mouse}(:, sessions_after_surgery);
    end
    
    % Store data
    group_idx = experimental_groups{ismember(experimental_groups(:, 1), TIMEPOINTS(i_mouse, 2)), 2};
    DATA_amplitude = [DATA_amplitude; [ones(size(amplitudes, 1), 1)  * group_idx, amplitudes]];
    DATA_frequency = [DATA_frequency; [ones(size(frequencies, 1), 1) * group_idx, frequencies]];
end

% DATA_amplitude(isnan(DATA_amplitude)) = 0;
% DATA_frequency(isnan(DATA_frequency)) = 0;
DATA_amplitude(DATA_amplitude == 0) = NaN;
DATA_frequency(DATA_frequency == 0) = NaN;

disp ('done re-arranging')

%%
% Replace group index with name
column_names = {'day_1', 'day_2', 'day_1_post', 'day_2_post'};
DATA_amplitude_table = array2table(DATA_amplitude, 'VariableNames',['group', column_names]);
DATA_frequency_table = array2table(DATA_frequency, 'VariableNames',['group', column_names]);
DATA_amplitude_table.group = experimental_groups(DATA_amplitude_table{:, 'group'}, 1);
DATA_frequency_table.group = experimental_groups(DATA_frequency_table{:, 'group'}, 1);

if discard_naive_from_session_comparisons
    DATA_amplitude_table(ismember(DATA_amplitude_table.group, 'naive'), :) = [];
    DATA_frequency_table(ismember(DATA_frequency_table.group, 'naive'), :) = [];
end

if combine_sham_and_CCI
	DATA_amplitude_table{:, 'group'} = {'S1'} ;
    DATA_frequency_table{:, 'group'} = {'S1'} ;
end
% DATA_amplitude_table(DATA_amplitude_table.day_1 == 0, :) = [];
% DATA_frequency_table(DATA_frequency_table.day_1 == 0, :) = [];

data_ampl_sham = table2array(DATA_amplitude_table(ismember(DATA_amplitude_table.group, {'sham'}), 2:end));
data_ampl_cci = table2array( DATA_amplitude_table(ismember(DATA_amplitude_table.group, {'CCI'}), 2:end));

data_freq_sham = table2array(DATA_frequency_table(ismember(DATA_amplitude_table.group, {'sham'}), 2:end));
data_freq_cci= table2array( DATA_frequency_table(ismember(DATA_amplitude_table.group, {'CCI'}), 2:end));

data_sham = data_ampl_sham.* data_freq_sham;
data_cci = data_ampl_cci .* data_freq_cci;


%%
% 
% %%
% % Compute data range
% amplitude_range = table2array(DATA_amplitude_table(:, 2:end));
% amplitude_range = amplitude_range(:);
% amplitude_range = [nanmin(amplitude_range), nanmax(amplitude_range)];
% amplitude_range = [floor(amplitude_range(1) / 10) * 10, ceil(amplitude_range(2) / 10) * 10];
% amplitude_range_str = ['"', strjoin(value2str(amplitude_range, '%.1f'), ', '), '"'];
% frequency_range = table2array(DATA_frequency_table(:, 2:end));
% frequency_range = frequency_range(:);
% frequency_range = [nanmin(frequency_range), nanmax(frequency_range)];
% frequency_range_str = ['"', strjoin(value2str(frequency_range, '%.1f'), ', '), '"'];
% 
% % Cap amplitude range
% amplitude_range = [0, 10];
% amplitude_range_str = ['"', strjoin(value2str(amplitude_range, '%.1f'), ', '), '"'];
% frequency_range = [0, 0.2];
% frequency_range_str = ['"', strjoin(value2str(frequency_range, '%.1f'), ', '), '"'];
% 
% 
% % Replace group identity with "all"
% % DATA_amplitude_table{:, 'group'} = {'all'};
% % DATA_frequency_table{:, 'group'} = {'all'};
% 
% % Write data to disk so that it can be passed to the R script plotting
% filename_amplitude = os.path.join(GC.temp_dir, 'amplitude.csv');
% filename_frequency = os.path.join(GC.temp_dir, 'frequency.csv');
% writetable(DATA_amplitude_table, filename_amplitude);
% writetable(DATA_frequency_table, filename_frequency);
% % Run analysis in R
% if strcmp(data_type, 'new')
%     output_filename = os.path.join(output_folder, 'amplitude_session_1_new.pdf');
% else
%     output_filename = os.path.join(output_folder, 'amplitude_session_1.pdf');
% end
% run_in_R(os.path.join('Figures_paper_FMP', 'amplitude_frequency_distributions_all_days'), ['--file_input ', filename_amplitude], ['--file_output ', output_filename], '--data_type amplitude', ['--data_limits ', amplitude_range_str], '--log_scale FALSE'); printFilePath(output_filename, [], 0)
% if strcmp(data_type, 'new')
%     output_filename = os.path.join(output_folder, 'frequency_session_1_new.pdf');
% else
%     output_filename = os.path.join(output_folder, 'frequency_session_1.pdf');
% end
% run_in_R(os.path.join('Figures_paper_FMP', 'amplitude_frequency_distributions_all_days'), ['--file_input ', filename_frequency], ['--file_output ', output_filename], '--data_type frequency', ['--data_limits ', frequency_range_str]); printFilePath(output_filename, [], 0)
% 
% 
% %% Do statistics
% % clc
% % figure(1);
% % set(gcf, 'color','w')
% % clf
% % 
% % A = DATA_amplitude_table(:, {'group', 'day_1', 'day_2'});
% % A(any(isnan(A{:, {'day_1', 'day_2'}}), 2), :) = [];
% % sham = A{ismember(A{:, 'group'}, 'sham'), {'day_1', 'day_2'}};
% % CCI = A{ismember(A{:, 'group'}, 'CCI'), {'day_1', 'day_2'}};
% % shamd = sham(:, 2) - sham(:, 1);
% % CCId = CCI(:, 2) - CCI(:, 1);
% % p_day12 = ranksum(CCId, shamd)
% % median_CCId_12 = median(CCId)
% % median_shamd_12 = median(shamd)
% % 
% % subplot(1, 2, 1), cla, hold on
% % bw = .001;
% % x = [CCId(:); shamd(:)];
% % x = linspace(min(x), max(x), 1000);
% % ksdensity(CCId, x, 'Bandwidth',bw)
% % ksdensity(shamd, x, 'Bandwidth',bw)
% % legend('CCI', 'sham', 'AutoUpdate','off')
% % h = plot([0, 0], [0, 1000], 'color','k', 'linestyle','--', 'YLimInclude','off'); uistack(h, 'bottom')
% % xlim([-.5, .5])
% % title('Amplitude difference day 2 - day 1')
% % 
% % 
% % A = DATA_amplitude_table(:, {'group', 'day_1', 'day_2', 'day_1_post'});
% % no_day2 = isnan(A.day_2); has_day1 = ~isnan(A.day_1);
% % A(has_day1 & no_day2, 'day_2') = A(has_day1 & no_day2, 'day_1');
% % A(:, 'day_1') = [];
% % A(any(isnan(A{:, {'day_2', 'day_1_post'}}), 2), :) = [];
% % sham = A{ismember(A{:, 'group'}, 'sham'), {'day_2', 'day_1_post'}};
% % CCI = A{ismember(A{:, 'group'}, 'CCI'), {'day_2', 'day_1_post'}};
% % shamd = sham(:, 2) - sham(:, 1);
% % CCId = CCI(:, 2) - CCI(:, 1);
% % p_day23 = ranksum(CCId, shamd)
% % median_CCId_23 = median(CCId)
% % median_shamd_23 = median(shamd)
% % 
% % subplot(1, 2, 2), cla, hold on
% % x = [CCId(:); shamd(:)];
% % x = linspace(min(x), max(x), 1000);
% % ksdensity(CCId, x, 'Bandwidth',bw)
% % ksdensity(shamd, x, 'Bandwidth',bw)
% % h = plot([0, 0], [0, 1000], 'color','k', 'linestyle','--', 'YLimInclude','off'); uistack(h, 'bottom')
% % xlim([-.5, .5])
% % title('Amplitude difference day 1 post-surgery - last day before surgery')
% 
% 
% % Delete temporary files
% delete(filename_amplitude)
% delete(filename_frequency)
% % if exist('Rplots.pdf', 'file'), delete('Rplots.pdf'), end


%% MLint exceptions
%#ok<*AGROW>
