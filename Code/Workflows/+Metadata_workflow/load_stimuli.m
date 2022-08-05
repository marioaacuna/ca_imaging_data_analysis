function STIMULUS = load_stimuli(stimulus_name, varargin)

% Read general_configs
global GC

% Initialize empty variable
STIMULUS = struct();

% Parse inputs
p = inputParser();
addParameter(p, 'frame_rate', []);
addParameter(p, 'n_frames', []);
addParameter(p, 'pre_stimulus_window', []);
addParameter(p, 'post_stimulus_window', []);
addParameter(p, 'ignore_if_missing', false);
parse(p, varargin{:})
frame_rate = p.Results.frame_rate;
n_frames = p.Results.n_frames;
pre_stimulus_window = p.Results.pre_stimulus_window;
post_stimulus_window = p.Results.post_stimulus_window;
ignore_if_missing = p.Results.ignore_if_missing;

if isempty(frame_rate) || isempty(n_frames)
    error('MetadataWorkflow:load_stimuli', 'Missing inputs')
end

% Get path to filename of stimuli
stimuli_path = os.path.join(GC.data_root_path, GC.stimuli_path);

% Get name of stimulus and the file containing info on it
stimulus_filename = os.path.join(stimuli_path, [stimulus_name, '.mat']);
if ~exist(stimulus_filename, 'file') && ignore_if_missing
    return
end
% Copy stimulus into memory
STIMULUS = load_variable(stimulus_filename, 'STIMULUS');
% Unpack stimulus profile
row = table2array(STIMULUS.stimulus_profile(:, 'frame_rate')) == frame_rate & table2array(STIMULUS.stimulus_profile(:, 'n_frames')) == n_frames;
if ~any(row)
    error('MetadataWorkflow:load_stimuli', 'Combination of n_frames and frame_rate does not exist for %s', stimulus_name)
end
ta = STIMULUS.stimulus_profile{row, 'time_axis'};
sp = STIMULUS.stimulus_profile{row, 'stimulus_profile'};
STIMULUS.time_axis = ta{1};
STIMULUS.stimulus_profile = sp{1};

% Pre-stimulus window
if ~isempty(pre_stimulus_window)
    if isempty(STIMULUS.timestamps)  % Spontaneous activity has no stimulus onset
        % Get entire period
        baseline_idx = 1:length(STIMULUS.time_axis);
        % Store values
        STIMULUS.baseline_window = baseline_idx;
    else
        if length(STIMULUS.timestamps) == 1  % "Instantaneous" stimulus
            [~,ts] = min(abs(STIMULUS.time_axis - STIMULUS.timestamps(1)));
        else  % Stimulus has many timestamps, but consider only the first one (i.e., the stimulus onset)
            ts = zeros(length(STIMULUS.timestamps),1);
            for ii = 1:length(STIMULUS.timestamps)
                [~,ts(ii)] = min(abs(STIMULUS.time_axis - STIMULUS.timestamps(ii)));
            end
        end
        % Get baseline indices
        baseline_idx = ts(1)-ceil(pre_stimulus_window*frame_rate):ts(1)-1;
        % Remove out-of-bound indices
        baseline_idx(baseline_idx<1) = [];
        baseline_idx(baseline_idx>length(STIMULUS.time_axis)) = [];
        % Store values
        STIMULUS.baseline_window = baseline_idx;
    end
end

% Post-stimulus window
if ~isempty(post_stimulus_window)
    if isempty(STIMULUS.timestamps)  % Spontaneous activity has no "response" interval
        % Simply keep an empty value
        STIMULUS.response_window = [];
    else
        if length(STIMULUS.timestamps) == 1  % "Instantaneous" stimulus
            [~,ts] = min(abs(STIMULUS.time_axis - STIMULUS.timestamps(1)));
        else  % Stimulus has many timestamps, but consider only the first one (i.e., the stimulus onset)
            ts = zeros(length(STIMULUS.timestamps),1);
            for ii = 1:length(STIMULUS.timestamps)
                [~,ts(ii)] = min(abs(STIMULUS.time_axis - STIMULUS.timestamps(ii)));
            end
        end
        % Get the indices of the response window for this stimulus
        response_idx = ts(1):ts(end)+ceil(post_stimulus_window * frame_rate)-1;
        % Remove out-of-bound indices
        response_idx(response_idx<1) = [];
        response_idx(response_idx>length(STIMULUS.time_axis)) = [];
        % Store values
        STIMULUS.response_window = response_idx;
    end
end
