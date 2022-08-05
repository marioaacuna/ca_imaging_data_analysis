function [spikes, Ca_events] = deconvolve_spiking_events(animal_ID, F, frame_rate, frames_idx, SPIKES)

% Get logger and general_configs
global LOGGER GC

% Algorithm parameters
SNR_threshold = GC.calcium_transient_detection_SNR;
smoothing_width = floor(1000 / 1000 * frame_rate / 2) * 2 + 1;
adjust_peak_location_win = 5;
artefact_window = 5;

% Get number of trials and samples
n_trials = size(frames_idx, 1);
n_samples = max(frames_idx(:));

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
F_denoised_s = smooth(F_denoised,6, 'sgolay');


% %% Find peak of each calcium event
% LOGGER.trace('Searching for calcium events peak with CWT')
% % Old implementation
% % Make filenames of input and output files for R script
filename_data_input  = os.path.join(GC.temp_dir, ['peak_detection_data_input_', animal_ID, '.csv']);
filename_data_output = os.path.join(GC.temp_dir, ['peak_detection_data_output_', animal_ID, '.csv']);
% Loop through trials
CWT_peaks = [];
Ca_indicator = 'GCaMP6f';
tau_decay = GC.decay_time_constant_calcium_reporter.(Ca_indicator);
win_tau = 40 * tau_decay / (1000/frame_rate);
win_time = 15 * frame_rate;
win = max([win_tau, win_time]);
win = floor(win/2)*2 + 1;  % Round to nearest odd number
% win = 1027;%length(F_denoised_trial)/2;
%% Parameters
% set parameters
baseline = prctile(F_denoised_s, 8, 2);
F_smooth = runline(baseline, win, 1);
F_temp = (F_denoised_s - F_smooth);
rand_nr = F_temp(randperm(length(F_denoised_s), 1000)); % generate random signal
FP_peaks = [];

% F_noise = rssq(rand_nr(:));
% thresh = 2* std(rand_nr);
% 
is_neuron = sum(SPIKES,2)~=1;

if ~is_neuron
    % % R Version
    LOGGER.trace('Doing CWT')
    csvwrite(filename_data_input, F_denoised);
    run_in_R('CWTpeak_finder', ['--file_input ', filename_data_input], ['--file_output ', filename_data_output], ['--frame_rate ', num2str(frame_rate)], ['--SNR ', num2str(SNR_threshold)]);
    fileID = fopen(filename_data_output, 'r');
    lines = textscan(fileID, '%s', 'delimiter','\n');
    fclose(fileID);
    CWT_peaks = regexp(lines{1}, ',', 'split');
    CWT_peaks = cell2mat(cellfun(@str2double, CWT_peaks{1}, 'UniformOutput',false));
    CWT_peaks = CWT_peaks(:);
else
    %% new Version
    LOGGER.trace('Taking spikes from cnmf-e data')
    
    % warning('off') % in case no peaks are found, just skip
    % for itrial = 1:n_trials
    %     F_denoised_trial = F_denoised_s(frames_idx(itrial,1):frames_idx(itrial,2));
    %     % Write data to local disk to be passed to R
    %     win_i = 129;
    %     baselinei = prctile(F_denoised_trial, 8, 2);
    %     F_smooth_i = runline(baselinei, win_i, 1);
    %     F_temp_i = (F_denoised_trial - F_smooth_i);
    % %     [~,FP_peaks_temp] = findpeaks(F_denoised_trial, 'MinPeakHeight', abs((1+SNR_threshold)*median(F_temp(:))), 'MinPeakWidth', ceil(0.28*frame_rate));
    %     [~,FP_peaks_temp] = findpeaks(F_denoised_trial, 'MinPeakHeight', 3*std(F_denoised_trial), 'MinPeakWidth', ceil(0.6*frame_rate));
    %     % Convert to real frame idx
    %     if ~isempty(FP_peaks_temp)
    %         f_idx_fp = frames_idx(itrial,1) + FP_peaks_temp;
    %         FP_peaks = [FP_peaks; f_idx_fp];
    %     else
    %         continue
    %     end
    % end
    
    
    %% get the deconvolved spikes done by CNMFE protocol
    spikes_cnfme = SPIKES;
    SP_peaks = find(spikes_cnfme>0)';
    
end



%% Infer the most likely discretized spike train underlying the fluorescence trace with OASIS
LOGGER.trace('Deconvolving spikes with OASIS')

% Pre-allocate output variable
spikes = cell(1, n_trials);

for itrial = 1:n_trials
    % Pre-allocate output variables
    s = zeros(1, length(frames_idx(itrial,1):frames_idx(itrial,2)));

    % Get segment for this trial removing bad frames
    frames = setdiff(frames_idx(itrial,1):frames_idx(itrial,2), artefact_idx);
    dFF = F(frames);
    if length(dFF) <= 8  % OASIS needs more than 8 samples to run 
        % Store data in a row vector
        spikes{1, itrial} = s;
        continue
    end
    
    % Remove NaNs
    good_indices = isfinite(dFF);
    dFF_for_OASIS = dFF(good_indices);
    % Run OASIS
    if ~isempty(dFF_for_OASIS)
        [~, s_out, opt] = deconvolveCa(dFF_for_OASIS, 'method','constrained_foopsi', 'optimize_smin',1, 'optimize_pars',1, 'optimize_b',1);
        % Reposition values in main array
        s(good_indices) = s_out;

        % Visualize deconvolution
        if false
            figure('color','w')
            clf, hold on,
            plot(dFF,'k'),
            fill([1:length(dFF), fliplr(1:length(dFF))], [s(:)', zeros(1,length(s))],'b','EdgeColor','none')
        end
    end

    % Store data in a row vector
    spikes{1, itrial} = s(:).';
end

% Concatenate spikes of all trials
S = cell2mat(spikes).';


%% Check correlation of spikes to traces for both implementations
% get the correlation to the F trace
% CORR = check_correlation_spikes(CWT_peaks, FP_peaks,  S, F, F_denoised, artefact_idx, pos_zx, neg_zx, frame_rate);

% CWT_peaks = FP_peaks;
%%
% u = unique(CWT_peaks);
% CWT_peaks = FP_peaks;
if is_neuron
    CWT_peaks = SP_peaks;
    % Allocate 4 columns in the variables CWT_peaks
    CWT_peaks(:,2:5) = 0;
    % Convert time to frames
    % CWT_peaks(:,3) =round(CWT_peaks(:,1) * frame_rate); % data coming from CWT
    CWT_peaks(:,3) = CWT_peaks(:,1); %round(CWT_peaks(:,1) * frame_rate); % Data
    % Get position of good indices
else
    % Allocate 4 columns in the variables CWT_peaks
    CWT_peaks(:,2:5) = 0;
    % Convert time to frames
    % CWT_peaks(:,3) =round(CWT_peaks(:,1) * frame_rate); % data coming from CWT
    CWT_peaks(:,3) = round(CWT_peaks(:,1) * frame_rate); % Data
    % Get position of good indices
end


good_indices = setdiff(1:n_samples, artefact_idx);
% CWT peaks refer to the good indices, so realign values
CWT_peaks(:,3) = good_indices(CWT_peaks(:,3));
clear good_indices

%% Compare approaches
LOGGER.trace('Compare approaches')

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
    CWT_peaks(iev, 5) = any(S(event_onset:event_end) ~= 0);
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
S_denoised = S .* calcium_index;

%% Prepare outputs
spikes = cell(1, n_trials);
for itrial = 1:n_trials
    % Get denoised spiking activity and store it as row vector
    spikes{1,itrial} = reshape(S_denoised(frames_idx(itrial,1):frames_idx(itrial,2)), 1,[]);
end

% Keep event timestamps of each peak
Ca_events = CWT_peaks(:, 2:4);
Ca_events(CWT_peaks(:,5)==0 | CWT_peaks(:,7)==0, :) = [];

% plot calcium events
if false
    figure
   for iev = 1:length(Ca_events)
       this_x = Ca_events(iev, 1):Ca_events(iev,3);
       plot(F(this_x))
       hold on
       drawnow
       pause(0.1)
       cla
   end
    
end

% Delete last output and log files
if exist(filename_data_input, 'file'), delete(filename_data_input), end
if exist(filename_data_output, 'file'), delete(filename_data_output), end

%% MLint exceptions
%#ok<*TRYNC,*NASGU,*ASGLU,*UNRCH,*AGROW>
