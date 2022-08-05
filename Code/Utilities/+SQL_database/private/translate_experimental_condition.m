function conds = translate_experimental_condition(experimental_conditions, col_idx)
% TRANSLATE_EXPERIMENTAL_CONDITION This function transforms a semicolon-separated
% array of strings into individual parts, corresponding to different
% experimental variables.

% Set default values
if ~exist('col_idx','var'), col_idx=[]; end

% Make sure input is a cell array
if istable(experimental_conditions)
    experimental_conditions = table2cell(experimental_conditions);
end

% Split at the semicolon
conds_split = cellfun(@(x) regexp(x, ';', 'split'), experimental_conditions, 'UniformOutput',false);
% Check that we obtained the same number of columns
n_columns = unique(cellfun(@length, conds_split));
if length(n_columns) == 1
    conds = vertcat(conds_split{:});
else
    conds = repmat({''}, size(conds_split,1), max(n_columns));
    for irow = 1:size(conds,1)
        conds(irow,1:length(conds_split{irow})) = conds_split{irow};
    end
end

% If user requested, return only selected columns
if ~isempty(col_idx)
    % If col_idx is out of range, let MATLAB raise an error!
    conds = conds(:, col_idx);
end
