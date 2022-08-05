clc
clear
global GC
GC.experiment_name = 'ACC_CCI_anesth';
group_trials_by =GC.analysis.(GC.experiment_name).group_trials_by_timepoint;%GC.group_trials_by_timepoint.(GC.experiment_name);
toggle_toolbox('Fractal_dimensions','on')
take_raw_data = true;
session_to_take = 4;
saving_data = 1;

data_type = 'dFF';
% amplitude_type = 'peak';  % Take the difference between peak fluorescence in an event minus the fluorescence at the beginning of the event
% amplitude_type = 'auc';  % Take the sum of the fluorescence in the event
% amplitude_type = 'z';  % z-score of whole trace
% amplitude_type = 'mad';  % Like z-score but it divides median by the MAD
% amplitude_type = 'spikes';  % Deconvolved event amplitude

stimulus_of_interest = 'SP';  % '' means from all stimuli

% Get names of allowed compounds
allowed_compounds ={''};% GC.analysis_allowed_compounds{ismember(GC.analysis_allowed_compounds(:,1), GC.experiment_name), 2};


animal_names = SQL_database.read_table_where('experiments', {'animal_ID', 'experimental_group'}, animal_list(),'animal_ID');
n_animals = height(animal_names);

TRACES = struct();
% filename = ['EVENT_stats_', stimulus_of_interest, '_', amplitude_type, '.mat'];
% EVENT_STATS = load_variable(filename, 'EVENT_STATS');
these_traces_sham = NaN(0, 5100);
these_traces_CCI = these_traces_sham;

%%
for i_mouse = 1:n_animals
    animal_ID = animal_names{i_mouse, 'animal_ID'};
    animal_ID = animal_ID{1};
    disp(animal_ID)
        % Read timepoints from database
    METADATA = SQL_database.read_table_where('trials', {'stimulus','frame_rate','n_frames','experimental_condition','date','frames_idx', 'experimental_group'}, animal_ID,'animal_ID');
    this_group = cell2mat(unique(METADATA.experimental_group, 'stable'));
    all_trials_idx = 1:height(METADATA);
    % Keep only trials with allowed compounds
    % select only one session
    trials_idx = ismember(METADATA.compound, allowed_compounds) & ismember(METADATA.stimulus, stimulus_of_interest);
    frames_idx = cellfun(@(x) strsplit(x, ','), METADATA.frames_idx, 'UniformOutput',false);
    frames_idx = vertcat(frames_idx{:});
    frames_idx = str2double(frames_idx);
    frames_idx = frames_idx(trials_idx, :);
    %         METADATA = METADATA(trials_idx, :);
    all_trials_idx = all_trials_idx(trials_idx);
    bad_trials_idx = find(~trials_idx);
    
    % Split data by timepoint
    [sessions, ~, session_idx] = unique(METADATA{:, group_trials_by}, 'stable');
    sessions_d = str2double(sessions);
    n_sessions = length(sessions);
    sessions_to_keep_idx = [find(sessions_d < 0, 2, 'last'); find(sessions_d > 0, 2, 'first')];
    try
        sessions_to_keep_idx = sessions_to_keep_idx(session_to_take);
    catch
        disp(['no session, animal ', animal_ID])
    continue 
    end

    % load traces 
    traces_filename = get_filename_of('spikes', animal_ID);
    if ~exist(traces_filename, 'file'), continue, end
    Ca_events = load_variable(traces_filename, 'Ca_events');
    if strcmp(data_type, 'spikes')
        dFF = load_variable(traces_filename, 'spikes');
    else
        dFF_filename = get_filename_of('dFF', animal_ID);
        dFF = load_variable(dFF_filename, 'dFF');
    end
    Ca_events_filename = get_filename_of('Ca_events', animal_ID);
    keep_cells = load_variable(Ca_events_filename, 'keep_cell');
    keep_cells = logical(keep_cells);
    Ca_events =  Ca_events(keep_cells,:);
    dFF = dFF(keep_cells,:);
    n_ROIs = size(dFF, 1);
    n_trials = size(dFF, 2);
    % Mark trials of dFF
    dFF_trials = cell(n_trials, 1);
    for itrial = 1:n_trials
        dFF_trials{itrial} = ones(length(dFF{1, itrial}), 1) * itrial;
    end
    dFF_trials = cell2mat(dFF_trials);

    % Allocate traces for this animal
    trials_this_sess = ismember(METADATA.stimulus, stimulus_of_interest) & ismember(METADATA.day_from_surgery, sessions(sessions_to_keep_idx));
    dFF_these_cells = dFF(:, trials_this_sess);
    dFF_these_cells = cell2mat(dFF_these_cells);
    if length(dFF_these_cells) > 5100
        dFF_these_cells = dFF_these_cells (:, 1:5100);
    end
 

%     dFF_these_cells = dFF_these_cells - median(dFF_these_cells);
    
    if strcmp (this_group, 'sham')
        these_traces_sham = [these_traces_sham; dFF_these_cells ];
        TRACES.(this_group) = these_traces_sham;
    else
        these_traces_CCI = [these_traces_CCI; dFF_these_cells ];
        TRACES.(this_group) = these_traces_CCI;
    end
    
end



%% calculate cum muean and cum std
dFF_sham = TRACES.sham;
dFF_CCI = TRACES.CCI; 

% % Surrogates
% n_surrogates = size(y,1);
% surrogate_length = size(y,2);
% y_rnd = Analysis_workflow.generate_surrogate_spike_trains(y, 'n_surrogates',n_surrogates, 'surrogate_length',surrogate_length, 'plot_distributions',false); % shuffled surrogates

[m_sham, s_sham] = f_cummeanstd(dFF_sham);
[m_cci, s_cci] = f_cummeanstd(dFF_CCI);
%% plot
mean_sham = mean(m_sham(:));
mean_sham1 = m_sham - mean_sham;

mean_CCI = mean(m_cci(:));
mean_CCI1 = m_cci - mean_CCI;


figure ('color', 'w')
errorbar(mean(m_sham), sem(m_sham), '-o', 'MarkerFaceColor', 'b')
hold on
errorbar(mean(m_cci), sem(m_cci), '-o',  'MarkerFaceColor', 'r')
xlim([0,100])


figure ('color', 'w')
plot(log10(mean(s_sham)), 'b')
hold on
plot(log10(mean(s_cci)), 'r')

disp('done')
%%
%#ok<*AGROW>