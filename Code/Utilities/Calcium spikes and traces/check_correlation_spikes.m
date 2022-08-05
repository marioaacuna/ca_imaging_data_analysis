function [cor2trace] = check_correlation_spikes(CWT,FP, S, F, F_denoised, artefact_idx, pos_zx, neg_zx,frame_rate)

fnxs = {'CWT', 'FP'};
F_masked_noise  = F;
F_masked_noise(artefact_idx) = NaN;
n_samples = length(F_masked_noise);

% initiate output
cor2trace = NaN(length(fnxs),1);
for i_f = 1:length(fnxs)
    
    this_f = fnxs{i_f};
    switch this_f
        case 'CWT'
            CWT_peaks = CWT;
            
            % Allocate 4 columns in the variables CWT_peaks
            CWT_peaks(:,2:5) = 0;
            % Convert time to frames
            CWT_peaks(:,3) =round(CWT_peaks(:,1) * frame_rate); % data coming from CWT
            % Get position of good indices
            good_indices = setdiff(1:n_samples, artefact_idx);
            % CWT peaks refer to the good indices, so realign values
            CWT_peaks(:,3) = good_indices(CWT_peaks(:,3));
            clear good_indices

            
        case 'FP'
            CWT_peaks = FP;
            CWT_peaks(:,2:5) = 0;
            CWT_peaks(:,3) = FP; %round(CWT_peaks(:,1) * frame_rate); % Data
            % Get position of good indices
            good_indices = setdiff(1:n_samples, artefact_idx);
            % CWT peaks refer to the good indices, so realign values
            CWT_peaks(:,3) = good_indices(CWT_peaks(:,3));
            clear good_indices
            
            
    end
    
    
    
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
    
    
    %% Get correlation
    CORR= corrcoef(S_denoised, F_denoised);
    cor2trace(i_f) = CORR(2);
    
    
    
end
cor2trace = array2table(cor2trace', 'VariableNames', fnxs);

end