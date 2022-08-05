function z_max = fit_zmax_to_data(rPC_proj, fallback_value)

% Make parametes grid
Us = linspace(0.01, 2, 100);

% Compute solution for all parameters
all_local_minima = NaN(length(Us), 1);
for ip = 1:length(Us)
    % Compute density histogram
    [densityNorm, x] = ksdensity(rPC_proj, 'Bandwidth',Us(ip)); 
    % Invert upside down
    y = -densityNorm;
    % Find first peak
    [peak_values, peak] = findpeaks(y, x, 'NPeaks',1);
    % Skip if no local minimum
    if isempty(peak)
        continue
    end
    % In case the left shoulder was considered (it normally isn't)
    peak = peak(peak_values<0);
    % Keep only the first peak
    all_local_minima(ip) = peak(1);
end

% Remove NaNs
all_local_minima(isnan(all_local_minima)) = [];
% If nothing is left, use default value
if isempty(all_local_minima)
    z_max = fallback_value;

else
    % Compute optimal z_max
    z_max = mean(all_local_minima);
end
