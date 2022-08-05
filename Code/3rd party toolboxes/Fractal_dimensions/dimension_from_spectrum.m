function [D] = dimension_from_spectrum(welch)

y_data_log = log10(welch);
x_data_log = log10(1:size(welch, 1))';

maxSlope = cell(0,1);
windowWidth = 500; % Whatever...
maxk = 1;
i_it = 0;
for k = 1 : windowWidth % length(x_data_log) - windowWidth
    it = 1 + i_it;
    % Get a subset of the data that is just within the window.
    windowedx = x_data_log(k:k+windowWidth-1);
    windowedy = y_data_log(k:k+windowWidth-1);
    % Fit a line to that chunk of data.
    coefficients = polyfit(windowedx, windowedy, 1);
    % Get the slope.
    slope = coefficients(1);
    % See if this slope is steeper than any prior slope.
    if slope < 0
        % It is steeper.  Keep track of this one.
        maxSlope{it} = slope;
        maxk = k;
        i_it = 0 + it;
    end
end
[P, idx ]= min(cell2mat(maxSlope(:,1:20)));
D = 2 + (1/4 * P); % otherwise use 1/1



%% 
if nargout < 1
    figure('color', 'w')
    k = idx;
    it = 1 + i_it;
    % Get a subset of the data that is just within the window.
    windowedx = x_data_log(k:k+windowWidth-1);
    windowedy = y_data_log(k:k+windowWidth-1);
    % Fit a line to that chunk of data.
    coefficients = polyfit(windowedx, windowedy, 1);
    % Get the slope.
    slope = coefficients(1);
    % See if this slope is steeper than any prior slope.
    if slope < 0
        % It is steeper.  Keep track of this one.
        maxSlope{it} = slope;
        maxk = k;
        i_it = 0 + it;
    end
    
    % P = polyfit(x_data_log, y_data_log, 1);
    % f1 = polyval(P,x_data_log);
    f1 = polyval(coefficients,windowedx);
    plot (x_data_log, y_data_log)
    hold on
    plot(windowedx, f1)
    H = 2 - D;
    text(2.5,0, ['H = ', num2str(H)])
    xlabel('Log10(Freq)')
    ylabel('Log10(Power)')
end
