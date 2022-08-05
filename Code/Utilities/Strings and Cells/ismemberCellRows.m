function results = ismemberCellRows(cell_array, cells_to_search, varargin)
%ismemberCellRows Function to get whether a certain cell array is member of 
%another one. It can handle columns with mixed data types.
% 
% Run ismemberCellRows without inputs for a demo.

% Accept varargin to maintain compatibility with ismember, but ignore them
varargin = {};

% Run demo
if nargin == 0, run_demo(), return, end


%% FUNCTION
% Get number of variables
n_variables = size(cell_array, 2);
n_samples = size(cell_array, 1);

% Initialize output variable with all true, because in the following loop we'll
% test for all conditions to be true (using the AND operator), so if any one
% condition will turn out to be false, that sample will marker as false.
results = true(n_samples, 1);

for i_var = 1:n_variables
    % Test ismember
    results = results & ismember(cell_array(:, i_var), cells_to_search(:, i_var));
end


%% DEMO
function run_demo()
    % Make and show input
%     cell_array_in = {'1', '1'; '1', '1'; '1', 1; '1', 1; '2', NaN; '2', NaN};
%     fprintf('\nThis is a demonstration of the function <strong>uniqueCellRows</strong>\n\n')
%     fprintf('Given this cell array with mixed data types in the second column:\n')
%     disp(cell_array_in)
%     
%     % Compute and show result
%     output = uniqueCellRows(cell_array_in);
%     fprintf('\nThe following are the unique rows of it:\n');
%     disp(output)

