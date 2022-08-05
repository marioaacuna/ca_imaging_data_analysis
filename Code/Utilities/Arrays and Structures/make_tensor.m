function [X, T, bin_edges] = make_tensor(data, time_bin_activity_n_frames, timestamp, average_function)

if ~exist('average_function', 'var'), average_function = 'mean'; end

switch average_function
    case 'mean', fn = @nanmean;
    case 'median', fn = @nanmedian;
    case 'sum', fn = @nansum;
    case 'max', fn = @max;
end
        
[N, K] = size(data);
all_T = cellfun(@length, data);
T = unique(all_T(:));
if length(T) > 1
    T = min(T);
end

% Find edges of each bin
if exist('timestamp','var') && ~isempty(timestamp)
    before_ts = timestamp:-time_bin_activity_n_frames:1;
    after_ts  = timestamp:time_bin_activity_n_frames:T;
    bin_edges = unique([before_ts(:); after_ts(:)]);
else
    bin_edges = (1:time_bin_activity_n_frames:T)';
end
bin_edges = [bin_edges, bin_edges+time_bin_activity_n_frames-1];
bin_edges(any(bin_edges' > T | bin_edges' < 1), :) = [];
n_bins = size(bin_edges, 1);

% Create tensor
X = zeros(N, n_bins, K);
for iroi = 1:N
    for itrial = 1:K
        this_trial = data{iroi, itrial};
        % Bin activity
        this_trial_binned = NaN(n_bins, 1);
        for ibin = 1:n_bins
            this_trial_binned(ibin) = fn(this_trial(bin_edges(ibin,1):bin_edges(ibin,2)));
        end
        X(iroi, 1:length(this_trial_binned), itrial) = this_trial_binned;
    end
end
