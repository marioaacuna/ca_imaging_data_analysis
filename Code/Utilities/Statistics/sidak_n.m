function p = sidak_n(p, n)

% Convert to array
if iscell(p)
    input_is_cell = true;
    p = cell2mat(p);
else
    input_is_cell = false;
end

% Apply Sidak correction for n tests
p = 1 - (1 - p) .^ n;
p(p > 1) = 1;
p (~isfinite(p)) = 1;

% Reconvert ouput
if input_is_cell
    p = num2cell(p);
end
