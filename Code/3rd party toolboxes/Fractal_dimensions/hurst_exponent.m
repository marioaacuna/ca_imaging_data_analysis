% The Hurst exponent
%--------------------------------------------------------------------------
% The first 20 lines of code are a small test driver.
% You can delete or comment out this part when you are done validating the 
% function to your satisfaction.
%
% Bill Davidson, quellen@yahoo.com
% 13 Nov 2005

% function []=hurst_exponent()
% disp('testing Hurst calculation');
% 
% n=100;
% data=rand(1,n);
% plot(data);
% 
% hurst=estimate_hurst_exponent(data);
% 
% [s,err]=sprintf('Hurst exponent = %.2f',hurst);disp(s);

%--------------------------------------------------------------------------
% % Hurst introduces the concept of the “Hurst exponent” H, which may
% be understood as quantifying the character of the randomness exhibited in a timeseries
% structure via an autocorrelation measurement. Specifically, a Hurst exponent
% of H ¼ 0:5 describes a process that is purely random, such that the value of the trace
% at time ti is entirely independent of the value at time tj, i 6¼ j. By contrast, Hurst
% exponents in the range 0:5<H< 1 represent traces exhibiting positive autocorrelations,
% while Hurst exponents in the range 0 <H<0:5 represent traces exhibiting
% negative autocorrelations.
%--------------------------------------------------------------------------
% This function does dispersional analysis on a data series, then does a 
% Matlab polyfit to a log-log plot to estimate the Hurst exponent of the 
% series.
%
% This algorithm is far faster than a full-blown implementation of Hurst's
% algorithm.  I got the idea from a 2000 PhD dissertation by Hendrik J 
% Blok, and I make no guarantees whatsoever about the rigor of this approach
% or the accuracy of results.  Use it at your own risk.
%
% Bill Davidson
% 21 Oct 2003

function [h] = hurst_exponent(data0, steepest)   % data set

data=data0;         % make a local copy

[npoints,M]=size(data0);

yvals=zeros(npoints,1);
xvals=zeros(npoints,1);
data2=zeros(npoints,1);

index=0;
binsize=1;

while npoints >= 8
%     y = std(data);
    index = index+1;
    xvals(index) = binsize;
%     yvals(index)= binsize*y;
    % Get range
%     data = detrend(data);
%     figure, plot(data)
%     keyboard
%     R = max(data) - min(data);
%      yvals(index) = R/y;
    yvals(index) = RScalc(data0,binsize);
    npoints = fix(npoints/2);
    binsize =  binsize * 2; 
    for ipoints = 1:npoints % average adjacent points in pairs
        data2(ipoints) = (data(2 * ipoints) + data((2*ipoints) - 1)) * 0.5;
%           data2(ipoints) = data(ipoints);
    end
    data = data2(1:npoints);
end

xvals = xvals(1:index);
yvals = yvals(1:index);

logx = -log10(xvals);
logy = log10(yvals);


% Do steepest line
if steepest
    maxSlope = cell(0,1);
    windowWidth = 3; % Whatever...
    maxk = 1;
    i_it = 0;
    for k = 1 : length(logx) - windowWidth
        it = 1 + i_it;
        % Get a subset of the data that is just within the window.
        windowedx = logx(k:k+windowWidth-1);
        windowedy = logy(k:k+windowWidth-1);
        % Fit a line to that chunk of data.
        coefficients = polyfit(windowedx, windowedy, 1);
        % Get the slope.
        slope = coefficients(1);
        % See if this slope is steeper than any prior slope.
        if slope > 0
            % It is steeper.  Keep track of this one.
            maxSlope{it} = slope;
            maxk = k;
            i_it = 0 + it;
        end
    end
    
    % Select only the max slope
    [P, idx ]= max(cell2mat(maxSlope)); % (use idx to for plotting the steepest line)
    h=P;
else
    p2=polyfit(logx,logy,1);
    h=p2(1); % Hurst exponent is the slope of the linear fit of log-log plot
end
if 1.0001 < h && h < 1.05
    h = 1;
end
%plot
if nargout < 1
    k = idx;
     it = 1 + i_it;
    % Get a subset of the data that is just within the window.
    windowedx = logx(k:k+windowWidth-1);
    windowedy = logy(k:k+windowWidth-1);
    % Fit a line to that chunk of data.
    coefficients = polyfit(windowedx, windowedy, 1);
    % Get the slope.
    slope = coefficients(1);
    % See if this slope is steeper than any prior slope.
    if slope > 0
        % It is steeper.  Keep track of this one.
        maxSlope{it} = slope;
        maxk = k;
        i_it = 0 + it;
    end
    
    figure('color', 'w'), plot(logx, logy, 'ob-')
    xlabel('log(tau)')
    ylabel('log(R/S)')
    hold on
    f1 = polyval(coefficients, logx);
    plot(logx, f1, 'color', 'r')
end

function rs = RScalc(Z,n)
% Calculate (R/S)_n for given n
m = length(Z)/n;
Y = reshape(Z,n,m);
E = mean(Y);
S = std(Y);
for i=1:m
    Y(:,i) = Y(:,i) - E(i);
end
Y = cumsum(Y);
% Find the ranges of cummulative series
MM = max(Y) - min(Y);
% Rescale the ranges by the standard deviations
CS = MM./S;
rs = mean(CS);
return;
