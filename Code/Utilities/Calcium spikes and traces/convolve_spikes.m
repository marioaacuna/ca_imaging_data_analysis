function tr = convolve_spikes(spikes, frame_rate, tau_decay, varargin)

if isempty (varargin)
    tau_growth = 0.1;
else
    tau_growth = varargin{1};
end

% Allocate variable to hold output traces
tr = cell(size(spikes));

% Make time axis
trial_duration = max(reshape(cell2mat(cellfun(@length, spikes, 'UniformOutput',false)),[],1));
t = 0 : 1/frame_rate : trial_duration/frame_rate - 1/frame_rate;
% tau_growth = 0.1;
tau_growth_frames = round(frame_rate * tau_growth);
U = (1-exp(-t./tau_growth)) .* exp(-t./(tau_decay/1000));
U = (U-min(U)) ./ (max(U)-min(U)); % normalize in range [0 1]

% Loop through trials
for irow = 1:size(spikes, 1)
    for icol = 1:size(spikes, 2)
        % Get response and convolve it with unitary response
        this_r = spikes{irow, icol};
        this_r_convolved = conv(this_r, U, 'full');
        % Remove extra datapoints and store result
        tr{irow,icol} = this_r_convolved(1+tau_growth_frames:length(this_r)+tau_growth_frames);
    end
end
