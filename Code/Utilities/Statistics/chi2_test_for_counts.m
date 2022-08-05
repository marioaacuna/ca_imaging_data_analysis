function [p, stats] = chi2_test_for_counts(counts1, counts2, varargin)

% Parse inputs
p = inputParser();
addParameter(p, 'data_is_binned', false)  % data has not been binned
parse(p, varargin{:})
data_is_binned = p.Results.data_is_binned;

% Make data into row vectors
counts1 = counts1(:).';
counts2 = counts2(:).';

if data_is_binned
    % If data is binned, data is histograms. Therefore, they should have the
    % same length
    if length(counts1) ~= length(counts2)
        error('chi2_test_for_counts:UnequalLength', 'The two populations should have the same number of bins')
    end

    bins = 1:length(counts1);
    
else  % Data is not binned
    % Compute bins
    all_data = [counts1, counts2];
    bins = min(all_data):max(all_data);
    
    % Bin data
    counts1 = hist(counts1, bins);
    counts2 = hist(counts2, bins); 
end

observed_distribution = [counts1; counts2];
row_total = sum(observed_distribution, 2);
column_total = sum(observed_distribution, 1);
grand_total = sum(observed_distribution(:));
expected_distribution = (row_total * column_total) ./ grand_total;

% Compute chi-square statistics, and degrees of freedom
cstat = (observed_distribution - expected_distribution).^2 ./ expected_distribution;
cstat(~isfinite(cstat)) = 0;
df = size(cstat, 2) - 1;
cstat = sum(cstat(:));
% Gather additional information
stats = struct('chi2stat',cstat, 'df',df, 'Ctrs',bins, 'O',observed_distribution, 'E',expected_distribution);

% Compute p-value
p = gammainc(cstat/2, df/2, 'upper');
