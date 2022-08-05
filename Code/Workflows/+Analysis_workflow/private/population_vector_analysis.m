function results = population_vector_analysis(animal_ID,spikes, METADATA, session_date, stimuli_to_analyze, varargin)

% Get GC and LOGGER
global GC LOGGER

% Parse additional user inputs
p = inputParser(); 
addOptional(p, 'only_evoked',false)
parse(p, varargin{:});
only_evoked = p.Results.only_evoked;

% Count number of cells and stimuli to analyze
n_stimuli_to_analyze = length(stimuli_to_analyze);
n_cells = size(spikes, 1);

% Get other parameters from general_configs
time_bin_activity_n_frames = 1;%GC.assemblies_time_bin_activity_n_frames;

% Gather data from sessions and stimuli of interest
data = [];
bin_table_columns = {'bin', 'trial', 'stimulus', 'period'};
all_bin_edges = cell(0, length(bin_table_columns));
n_bins_per_trial = zeros(n_stimuli_to_analyze, 1);
n_frames_per_condition = zeros(n_stimuli_to_analyze, 1);
timestamps_bin = NaN(n_stimuli_to_analyze, 1);
all_timestamps = cell(n_stimuli_to_analyze, 1);
for i_stim = 1:n_stimuli_to_analyze
    % Get trials corresponding to this session
    trials_this_session = ismember(METADATA.date, session_date);
    trials_stimulus = ismember(METADATA.stimulus, stimuli_to_analyze{i_stim});
    trials_idx = trials_this_session & trials_stimulus; 
    if strcmp(stimuli_to_analyze{i_stim}, 'temp_43') && any(ismember(METADATA.stimulus, 'temp_42'))
       this_stim = 'temp_42';
       trials_stimulus = ismember(METADATA.stimulus, this_stim);
       trials_idx = trials_this_session & trials_stimulus; 
    end
    % Check whether to stop the analysis
    stop_analysis = ~any(trials_idx);
    if stop_analysis, break, end
    % Get spikes on these trials
    these_spikes = spikes(:, trials_idx);

    % Retrieve information about the stimulus
    frame_rate = unique(METADATA.frame_rate(trials_idx));
    n_frames_per_trial = unique(METADATA.n_frames(trials_idx));
    STIMULUS = Metadata_workflow.load_stimuli(stimuli_to_analyze{i_stim}, 'frame_rate',frame_rate, 'n_frames',n_frames_per_trial);
    if ~isempty(STIMULUS.timestamps)
        if ~startsWith(stimuli_to_analyze{i_stim}, 'temp_')
            timestamp = STIMULUS.timestamps(1);
            timestamp_frame = round(timestamp * frame_rate);
            all_timestamps{i_stim} = STIMULUS.timestamps;
        else
            timestamps = STIMULUS.timestamps;
            timestamp_frame = round(timestamps * frame_rate);
            timestamp_frame = timestamp_frame(1);
            all_timestamps{i_stim} = STIMULUS.timestamps;
            
        end
    else
        timestamp_frame = [];
        if only_evoked, continue, end
    end

    % Bin data
    [this_data, ~, bin_edges] = make_tensor(these_spikes, time_bin_activity_n_frames, timestamp_frame);
    LOGGER.info(['\tResponses to ', stimuli_to_analyze{i_stim}, ' binned every ', num2str(1 / frame_rate * time_bin_activity_n_frames * 1000, '%.f'), 'ms'])
    if ~isempty(timestamp_frame)
        % Get bin where the timestamp falls
        [~, timestamp_bin] = min(abs(bin_edges(:, 1) - timestamp_frame));
        % Remove data of epoch that is not of interest
        if only_evoked
            this_data = this_data(:, timestamp_bin:end, :);
        end
    else
        timestamp_bin = NaN;
    end

    % Concatenate data
    this_data = reshape(this_data, n_cells, []);
    data = [data, this_data];
    % Store info
    n_frames_per_condition(i_stim, 1) = size(this_data, 2);
    timestamps_bin(i_stim) = timestamp_bin;
    n_bins_per_trial(i_stim) = size(bin_edges, 1);
    % Store info on bins
    n_bins = size(bin_edges, 1);
    n_trials = sum(trials_idx);
    table_bins = num2cell([repmat((1:n_bins)', n_trials, 1), reshape(repmat((1:n_trials), n_bins, 1), [], 1)]);
    table_bins(:, 3) = stimuli_to_analyze(i_stim);
    if  ~isempty(timestamp_frame)
        if length(STIMULUS.timestamps) == 1
            counti = 0;
            for it = 1: n_trials
                counti = counti * n_bins;
                baseline = counti + 1 : counti + timestamps_bin(i_stim);
                evoked = counti + timestamps_bin(i_stim) + 1 : counti + n_bins_per_trial(i_stim);
                table_bins(baseline, 4) = {1};
                table_bins(evoked, 4)   = {2};
                counti = it ;
            end
        else
            counti = 0;
            for it = 1: n_trials
                counti = counti * n_bins;
                baseline_ts = round(all_timestamps{i_stim}(1) * frame_rate);
                evoked_ts = round(all_timestamps{i_stim}(2) * frame_rate);
                baseline = counti + 1 : counti + baseline_ts;
                evoked = counti +evoked_ts + 1 : counti + baseline_ts;
                table_bins(baseline, 4) = {1};
                table_bins(evoked, 4)   = {2};
                counti = it ;
            end
            
        end
    else
        table_bins(:, 4) = {0};
    end
    all_bin_edges = [all_bin_edges; table_bins];
    
end
% Quit analysis
if stop_analysis
    results = [];
    LOGGER.warn('Not all stimuli are present. Skipped')
    return
end

% Convert to table
all_bin_edges = cell2table(all_bin_edges, 'VariableNames',bin_table_columns);

% Get beginning and end of each condition in final data array
n_frames_per_condition = [[1; n_frames_per_condition(1:end-1) + 1], cumsum(n_frames_per_condition)];


%% SVD ensemble analysis
% Make parameters for analysis
param = struct();
param.pks = eps;
param.plot_sequences = 0;
param.plot_ROC = 0;
param.plot_activation_maps = 0;
param.pks = 3;
% param.jcut = 0.04;
% If automatic detection use below
param.ticut = [];
param.jcut = [];


%% keep track of frames for each stimulus
param.all_bin_edges =all_bin_edges;

%%
if param.plot_activation_maps
    roi_filename = get_filename_of('ROI_info', animal_ID);
    ROI_info = load_variable(roi_filename, 'ROI_info');
    Ca_events_filename = get_filename_of('Ca_events', animal_ID);
    keep_cell = logical(load_variable(Ca_events_filename, 'keep_cell'));
    ROI_info_keep = ROI_info(keep_cell,:);
    coords = [];
    for i_cell = 1 :length (ROI_info_keep)
        coords(i_cell,:) = [mean(ROI_info_keep{i_cell, 5}(:,1)),mean(ROI_info_keep{i_cell, 5}(:,2))];
    end
    param.coords = coords;
end
% Binarize data
data_1 = double(data > 0);
% Perform analysis
[core_svd, state_pks, param] = findSVDensemble(data_1, param, all_bin_edges);
% Store results in output variable
results = struct();
results.core_svd = core_svd;
results.state_pks = state_pks;
results.param = param;
% Add other info
results.info = struct();
results.info.n_frames_per_condition = n_frames_per_condition;
results.info.n_bins_per_trial = n_bins_per_trial;
results.info.stimuli_to_analyze = stimuli_to_analyze;
results.info.timestamps_bin = timestamps_bin;
results.info.all_timestamps = all_timestamps;
results.info.bin_edges = all_bin_edges;


%% MLint exceptions
%#ok<*AGROW>
