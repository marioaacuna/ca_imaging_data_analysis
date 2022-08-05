function [spikes, Ca_events] = deconvolve_spikes_cell(animal_ID, F, frame_rate, verbose)

% Get general_configs
global GC

if ~exist('verbose', 'var'), verbose = true; end

% Algorithm parameters
SNR_threshold = GC.calcium_transient_detection_SNR;
smoothing_width = floor(1000 / 1000 * frame_rate / 2) * 2 + 1;
adjust_peak_location_win = 5;
artefact_window = 5;

% Get number of trials and samples
F = F(:);
n_samples = length(F);

% Smooth fluorescence trace and compute its first-order derivative
[~, g] = sgolay(1, smoothing_width);
% Smooth the data
diff_filter = 1;
F_smoothed = conv(F, factorial(diff_filter-1)/(-(1/frame_rate))^(diff_filter-1) * g(:,diff_filter), 'same');
% Apply 1st order derivative to smoothed data
diff_filter = 2;
F_derivative = conv(F_smoothed, factorial(diff_filter-1)/(-(1/frame_rate))^(diff_filter-1) * g(:,diff_filter), 'same');
% Find 0-crossings in 1st derivative (i.e., sign of product of consecutive
% samples is negative)
zx = find(sign(F_derivative(1:end-1).*F_derivative(2:end)) < 0);
% Remove spurious points
zx(zx<1) = [];
zx(zx>=length(F_derivative)) = [];
% Get the sign of points around 0-crossings
yx = [F_derivative(zx) F_derivative(zx+1)];
% Keep transitions from rising to falling
pos_zx = zx(yx(:,1)>=0 & yx(:,2)<0);
% Keep transitions from falling to rising
neg_zx = zx(yx(:,1)<0 & yx(:,2)>=0);

% Remove samples containing artefacts
artefact_idx = [];
% Get spurious_peaks
spurious_peaks = idx2range(find(abs(F_derivative)>100 | abs(F)>100 | ~isfinite(F_smoothed) | ~isfinite(F_derivative)));
for iart = 1:size(spurious_peaks,1)
    % Find closest 0-crossing on left and right side of each peak
    event_onset = max([zx(find(zx < spurious_peaks(iart,1), 1, 'last')) - artefact_window, 1]);
    event_end   = min([zx(find(zx > spurious_peaks(iart,2), 1, 'first')) - 1 + artefact_window, n_samples]);
    artefact_idx = [artefact_idx, event_onset:event_end];
end
good_samples = setdiff(1:n_samples, artefact_idx);
F_denoised = F(good_samples);

% Skip ROI if no data available
if all(F_denoised==0) || isempty(F_denoised)
    % Fill empty output variable
    spikes = cell(1, n_trials);
    for itrial = 1:n_trials
        spikes{itrial} = zeros(1, length(frames_idx(itrial,1):frames_idx(itrial,2)));
    end
    Ca_events = zeros(0,3);
    return
end

% Estimate baseline noise
F_noise = F_denoised;
F_noise(isnan(F_noise)) = [];
bins = linspace(min(F_noise), max(F_noise), 1000);
% Compute KS-density
[smoothDist, x] = ksdensity(F_noise, bins);
[~, indPeak] = max(smoothDist);
% De-mean
mean_F_noise = x(indPeak);
F_denoised = F_denoised - mean_F_noise;

%% Find peak of each calcium event
if verbose
    disp('Searching for calcium events peak with CWT')
end

% Make filenames of input and output files for R script
filename_data_input  = [GC.temp_dir, 'peak_detection_data_input_', animal_ID, '.csv'];
filename_data_output = [GC.temp_dir, 'peak_detection_data_output_', animal_ID, '.csv'];

% Write data to local disk to be passed to R
csvwrite(filename_data_input, F_denoised);
% Call R
run_in_R('CWTpeak_finder', ['--file_input ', filename_data_input], ['--file_output ', filename_data_output], ['--frame_rate ', num2str(frame_rate, '%.1f')], ['--SNR ', num2str(SNR_threshold)]);

fileID = fopen(filename_data_output, 'r');
lines = textscan(fileID, '%s', 'delimiter','\n');
fclose(fileID);
CWT_peaks = regexp(lines{1}, ',', 'split');
CWT_peaks = cell2mat(cellfun(@str2double, CWT_peaks{1}, 'UniformOutput',false));
CWT_peaks = CWT_peaks(:);
% Allocate 4 columns in the variables CWT_peaks
CWT_peaks(:,2:5) = 0;
% Convert time to frames
CWT_peaks(:,3) = round(CWT_peaks(:,1) * frame_rate);
% Get position of good indices
good_indices = setdiff(1:n_samples, artefact_idx);
% CWT peaks refer to the good indices, so realign values
CWT_peaks(:,3) = good_indices(CWT_peaks(:,3));
clear good_indices

%% Infer the most likely discretized spike train underlying the fluorescence trace with OASIS
if verbose
    disp('Deconvolving spikes with OASIS')
end

[~, s, opt] = deconvolveCa(F, 'method','constrained_foopsi', 'optimize_smin',1, 'optimize_pars',1, 'optimize_b',1);

% Visualize deconvolution
if false
    figure('color','w')
    clf, hold on,
    plot(F,'k'),
    fill([1:length(F), fliplr(1:length(F))], [s(:)', zeros(1,length(s))],'b','EdgeColor','none')
end

% Store data in a row vector
spikes = s(:).';

%% Compare approaches
if verbose
    disp('Compare approaches')
end

% Maks noise in original trace
F_masked_noise = F;
F_masked_noise(artefact_idx) = NaN;

% From each CWT peak find the event in the fluorescence trace
n_events = size(CWT_peaks,1);
for iev = 1:n_events
    % Get peak position
    spk = CWT_peaks(iev,3);
    % Re-align peak on fluorescence first derivative 
    [~, peak_idx] = min(abs(spk - pos_zx));
    spk_peak_idx = pos_zx(peak_idx);

    % Find two negative 0-crossings, one before and one after the peak
    event_onset = max([neg_zx(find(neg_zx < spk_peak_idx, 1, 'last')), 1]);
    event_end   = min([neg_zx(find(neg_zx > spk_peak_idx, 1, 'first')) - 1, n_samples]);

    % Re-calculate peak position in raw trace
    [~, spk_peak_idx] = max(F_masked_noise(event_onset:event_end));
    spk_peak_idx = spk_peak_idx + event_onset - 1;
    
    % Store results
    CWT_peaks(iev, 2:4) = [event_onset, spk_peak_idx, event_end];

    % Get whether there is a spiking event or not during the event
    CWT_peaks(iev, 5) = any(spikes(event_onset:event_end) ~= 0);
end

% Estimate baseline noise
F_noise = F_masked_noise;
for iev = 1:n_events
    F_noise(CWT_peaks(iev,2):CWT_peaks(iev,4)) = NaN;
end
F_noise(isnan(F_noise)) = [];
bins = linspace(min(F_noise), max(F_noise), 1000);
% Compute KS-density
[smoothDist, x] = ksdensity(F_noise, bins);
[~, indPeak] = max(smoothDist);
% De-mean
mean_F_noise = x(indPeak);
F_noise = F_noise - mean_F_noise;
% Compute KS-density
[smoothDist, x] = ksdensity(F_noise, bins);
[~, indPeak] = max(smoothDist);
xFit = x(1:indPeak);
dataToFit = smoothDist(1:indPeak);
% Fit Gaussian curve to negative portion of the distribution. Normalize data
% between 0 and 1 to skip fitting the scaling parameter.
[xData, yData] = prepareCurveData(xFit, dataToFit);
yData = yData - min(yData);
yData = yData / max(yData);
ft = fittype('exp( -(x)^2 / (2*sigma^2) )');
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = 0;
opts.Upper = 100;
opts.StartPoint = 0.5;
[fitresult, ~] = fit(xData, yData, ft, opts);
sigma = fitresult.sigma;
% Set threhsold using the sigma of the fitted curve
noise_threhsold = sigma * 3;
% Add info on peak amplitude and whether it crosses the threshold
CWT_peaks(:, 6) = F_masked_noise(CWT_peaks(:, 3));
CWT_peaks(:, 7) = CWT_peaks(:, 6) >= noise_threhsold;

% Remove spikes that do not fall under a Calcium transient or not cross the
% noise threshold
calcium_index = zeros(n_samples,1);
for iev = 1:n_events
    if CWT_peaks(iev,5) == 1 || CWT_peaks(iev,7) == 1
        calcium_index(CWT_peaks(iev,2):CWT_peaks(iev,4)) = 1;
    end
end
spikes = spikes .* calcium_index(:).';

%% Prepare outputs
% Keep event timestamps of each peak
Ca_events = CWT_peaks(:, 2:4);
Ca_events(CWT_peaks(:,5)==0 | CWT_peaks(:,7)==0, :) = [];

% Delete last output and log files
if exist(filename_data_input, 'file'), delete(filename_data_input), end
if exist(filename_data_output, 'file'), delete(filename_data_output), end

%% MLint exceptions
%#ok<*TRYNC,*NASGU,*ASGLU,*UNRCH,*AGROW>
