function varargout = generate_surrogate_spike_trains(training_data, varargin)

% Parse inputs
p = inputParser();
% Number of surrogate spike trains
addParameter(p, 'n_surrogates', 10000)
% Number of frames in the final surrogates
addParameter(p, 'surrogate_length', 100)
% Provide range for amplitudes
addParameter(p, 'amplitude_range', [])
% Whether to compute the characteristics of the surrogate spike trains
addParameter(p, 'return_distributions', false)
% Whether to plot the distributions
addParameter(p, 'plot_distributions', false)
parse(p, varargin{:})
% Unpack parsed inputs
n_surrogates = p.Results.n_surrogates;
surrogate_length = p.Results.surrogate_length;
amplitude_range = p.Results.amplitude_range;
return_distributions = p.Results.return_distributions;
plot_distributions = p.Results.plot_distributions;

% Make sure the amplitude range is made of 2 different numbers
if ~isempty(amplitude_range) && amplitude_range(1)==amplitude_range(2)
    amplitude_range(1) = 0;
    if all(amplitude_range==0)
        amplitude_range = [];
    end
end

% Allocate output variables
distributions = struct();
surrogate_spike_trains = zeros(n_surrogates, surrogate_length);

%% Compute feature distributions of train datasets
% Make a cell
if ~iscell(training_data)
    training_data = {training_data};
end

% Get the number of training trials
training_data = training_data(:);
n_trials = size(training_data,1);

% Get non-empty bins
events_location = cellfun(@(x) find(x==0), training_data, 'UniformOutput',false);
% Compute distance between bins
inter_event_intervals = cell(n_trials,1);
parfor itrial = 1:n_trials
    % Convert indices to a table
    event_table = idx2range(events_location{itrial});
    % If only one "no-activity" segment, it means that there is no activity
    % whatsoever. Skip trial
    if size(event_table,1) == 1, continue, end
    % Keep duration of inter-trial intervals
    inter_event_intervals{itrial} = event_table(:,3);
end
inter_event_intervals = cell2mat(inter_event_intervals);
% If all trials were empty, there were no spikes at all. Abort and return empty
if isempty(inter_event_intervals)
    varargout(1:nargout) = {[]};
    return
end
% Compute and store distribution
bins = 1:max(inter_event_intervals);
[CDF, bins, PDF] = make_CDF(inter_event_intervals, bins, true);
distributions.interval.x = bins;
distributions.interval.y = CDF;
if return_distributions || plot_distributions
    distributions.interval.PDF = PDF;
end

% Get amplitude of each firing event
activity_bins = cellfun(@(x) x(x>0).', training_data, 'UniformOutput',false);
activity_bins = cell2mat(activity_bins);
% Compute and store distribution
if isempty(amplitude_range)
    bins = linspace(min(activity_bins), max(activity_bins), 100);
else
    bins = linspace(amplitude_range(1), amplitude_range(2), 100);
end
[CDF, bins, PDF] = make_CDF(activity_bins, bins, true);
distributions.amplitude.x = bins;
distributions.amplitude.y = CDF;
if return_distributions || plot_distributions
    distributions.amplitude.PDF = PDF;
end

%% Generate surrogate data
% Make anonymous functions to pass to parfor-loop
fcn_draw_intervals  = @(x) round(compute_ICDF(surrogate_length, distributions.interval.y, distributions.interval.x, true));
fcn_draw_amplitudes = @(x) compute_ICDF(x, distributions.amplitude.y, distributions.amplitude.x, true);

% Loop through surrogates
parfor isur = 1:n_surrogates
    % Generate many, random intervals
    intervals = fcn_draw_intervals();
    if isempty(intervals), continue, end  % No event in this trial
    % Check the span of all these intervals
    event_location = cumsum(intervals);
    % Keep only up to n intervals that can fit in the surrogate trial
    idx = find(event_location <= surrogate_length, 1, 'last');
    if isempty(idx), continue, end  % The cell would have been active too sparsely to be seen in this trial
    % Discard other intervals
    event_location = event_location(1:idx);

    % Generate random amplitudes
    amplitudes = fcn_draw_amplitudes(length(event_location));

    % Place events in the surrogate trial
    new_train = zeros(1, surrogate_length);
    new_train(event_location) = amplitudes;
    surrogate_spike_trains(isur, :) = new_train;
end

% Return output
if nargin >= 1
    varargout{1} = surrogate_spike_trains;
end

%% Compute the characteristics of the surrogate spike trains
if return_distributions || plot_distributions
    % Convert the surrogate spike trains to a cell array
    surrogate_spike_trains_cell = mat2cell(surrogate_spike_trains, ones(n_surrogates,1), surrogate_length);
    
    % Get non-empty bins
    events_location = cellfun(@(x) find(x==0), surrogate_spike_trains_cell, 'UniformOutput',false);
    % Compute distance between bins
    inter_event_intervals = cell(n_surrogates,1);
    parfor itrial = 1:n_surrogates
        event_table = idx2range(events_location{itrial});
        inter_event_intervals{itrial} = event_table(:,3);
    end
    inter_event_intervals = cell2mat(inter_event_intervals);
    % Compute and store distribution
    bins = 1:max(inter_event_intervals);
    [CDF, bins, PDF] = make_CDF(inter_event_intervals, bins, true);
    distributions.surrogate_interval.x = bins;
    distributions.surrogate_interval.y = CDF;
    distributions.surrogate_interval.PDF = PDF;

    % Get amplitude of each firing event
    activity_bins = cellfun(@(x) x(x>0).', surrogate_spike_trains_cell, 'UniformOutput',false);
    activity_bins = cell2mat(activity_bins);
    % Compute and store distribution
    if isempty(amplitude_range)
        bins = linspace(min(activity_bins), max(activity_bins), 100);
    else
        bins = linspace(amplitude_range(1), amplitude_range(2), 100);
    end
    [CDF, bins, PDF] = make_CDF(activity_bins, bins, true);
    distributions.surrogate_amplitude.x = bins;
    distributions.surrogate_amplitude.y = CDF;
    distributions.surrogate_amplitude.PDF = PDF;
end

% Return output
if return_distributions && nargin >= 2
    varargout{2} = distributions;
end

%% Compare distributions of firing statistics
if plot_distributions
    n_samples = 100000;  % Number of draws in case using asymptotic PDF
    
    figure('color','w')
    clf
    % Amplitude
    subplot(1,2,1)
    cla, hold on
    y = compute_ICDF(n_samples, distributions.amplitude.y, distributions.amplitude.x, true);
    hist_y = hist(y, distributions.amplitude.x);
    x = distributions.amplitude.x;
    y = hist_y/length(y);
    plot(x, y, 'k', 'linewidth',2)
    % Surrogate amplitude
    if use_function
        y = compute_ICDF(n_samples, distributions.surrogate_amplitude.y, distributions.surrogate_amplitude.x, true);
        hist_y = hist(y, distributions.surrogate_amplitude.x);
        x = distributions.surrogate_amplitude.x;
        y = hist_y/length(y);
    else
        x = distributions.surrogate_amplitude.PDF(:,2);
        y = distributions.surrogate_amplitude.PDF(:,1);
    end
    plot(x, y, 'r', 'linewidth',2)
    % Make string for title
    if use_function
        str = ' (asymptotic)';
    else
        str = ' (empirical)';
    end
    xlim([0, .5]), box off, title(['Amplitude',str]), set(gca,'TickDir','out'), ylabel('Probability')
    
    % Inter-event interval
    subplot(1,2,2)
    cla, hold on
    y = compute_ICDF(n_samples, distributions.interval.y, distributions.interval.x, true);
    hist_y = hist(y, distributions.interval.x);
    x = distributions.interval.x;
    y = hist_y/length(y);
    plot(x, y, 'k', 'linewidth',2)
    % Surrogate inter-event interval
    if use_function
        y = compute_ICDF(n_samples, distributions.surrogate_interval.y, distributions.surrogate_interval.x, true);
        hist_y = hist(y, distributions.surrogate_interval.x);
        x = distributions.surrogate_interval.x;
        y = hist_y/length(y);
    else
        x = distributions.surrogate_interval.PDF(:,2);
        y = distributions.surrogate_interval.PDF(:,1);
    end
    plot(x, y, 'r', 'linewidth',2)
    xlim([0,100]), box off, title(['Inter-event interval (frames)', str]), set(gca,'TickDir','out')

    % Add legend
    legend('original','surrogate')
end

%% Mlin exceptions
%#ok<*NASGU,*UNRCH>
