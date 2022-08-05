function result = sem(varargin)

% Parse inputs
p = inputParser();
addRequired(p, 'sample1', @(x) ~isempty(x))
addOptional(p, 'dim', 1)
addOptional(p, 'sample2', [])
parse(p, varargin{:})
sample1 = p.Results.sample1;
dim = p.Results.dim;
sample2 = p.Results.sample2;

n1 = sum(isfinite(sample1), dim);
sd1 = nanstd(sample1, 0, dim);

if isempty(sample2)  % 1-sample case
    result = sd1 ./ sqrt(n1);
    
else
    n2 = sum(isfinite(sample2), dim);
    sd2 = nanstd(sample2, 0, dim);
    sd_pooled = sqrt(((n1-1) .* sd1.^2 + (n2-1) .* sd2.^2) ./ (n1 + n2 - 2));
    
    % https://www.statisticshowto.datasciencecentral.com/find-pooled-sample-standard-error/
    result = sd_pooled .* sqrt(1 ./ n1 + 1 ./ n2);
end

