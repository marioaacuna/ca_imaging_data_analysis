function [p, x2] = chi2_contingency_table(O, E)

% Compute degrees of freedom
[n_rows, n_cols] = size(O);
if n_rows>1 && n_cols>1
    dof = (n_rows - 1) * (n_cols - 1);
else
    if n_rows == 1
        dof = n_cols-1;
    elseif n_cols == 1
        dof = n_rows-1;
    end
end

if ~exist('E','var')
    % Accumulate total values
    T_rows = sum(O,2);
    T_cols = sum(O,1);
    T = sum(O(:));

    % Calculate expected values
    E = zeros(size(O));
    for irow = 1:size(E,1)
        for icol = 1:size(E,2)
            E(irow, icol) = T_rows(irow) * T_cols(icol) / T;
        end
    end
end

% Calculate chi square statistic
O = O(:); E = E(:);
x2_scores = (O-E).^2 ./ E;
x2_scores(~isfinite(x2_scores)) = [];
x2 = nansum(x2_scores);

% Compute p-value
p = 1 - chi2cdf(x2, dof);
