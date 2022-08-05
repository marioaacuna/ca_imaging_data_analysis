function [slope, intercept, p] = TheilSen_Regression(x,y)

x = x(:);
y = y(:);
N = size(x,1);

% Comb = combnk(1:N,2);
% deltay = diff(y(Comb),1,2);
% deltax = diff(x(Comb),1,2);
% 
% theil = deltay ./ deltax;
% slope = nanmedian(theil);
% 
% intercept = nanmedian(y-(slope.*x));

C = nan(N);
for i=1:N
    % Accumulate slopes
    C(i,i:end) = (y(i) - y(i:end)) ./ (x(i) - x(i:end));
end
slope = nanmedian(C(:));  % calculate slope estimate
intercept = nanmedian(y - slope .* x);  % calculate intercept


% Compute p-values
if nargout > 2
    % Source: https://stackoverflow.com/a/42677750
    warning('off', 'MATLAB:singularMatrix')
    
    params = [slope, intercept];
    yHat = x * slope + intercept;
    MSE = sum((y - yHat) .^ 2) / (size(x, 1) - size(x, 2));
    var_b = MSE * diag(inv(dot(x', x)));
    sd_b = sqrt(var_b);
    ts_b = params / sd_b;
    p = 2 .* (1 - tcdf(abs(ts_b), (size(x, 1) - 1)));
    
    warning('on', 'MATLAB:singularMatrix')
end
