function cell_array_out = uniqueCellRows(cell_array_in, varargin)
%uniqueCellRows Function to get the unique rows in a cell array. It can handle
%columns with mixed data types.
% 
% The output is always in 'stable' format, which means that the output will
% always list unique rows in the same order as they appeared.
% 
% Run uniqueCellRows without inputs for a demo.

% Parse the parameter 'stable'/'sorted'.
sorting_options = {'stable', 'sorted'};
default_sorting = 'stable';

if any(ismember(sorting_options, varargin))
    idx = find(ismember(sorting_options, varargin), 1, 'last');
    sorting = sorting_options{idx};
else
    sorting = default_sorting;
end

% Run demo
if nargin == 0, run_demo(), return, end


%% FUNCTION
% Get number of variables
n_variables = size(cell_array_in, 2);

% Initialize internal variables
variable_map = cell(1, n_variables);
variable_type = cell(1, n_variables);
cell_array_in_converted = NaN(size(cell_array_in));

for i_var = 1:n_variables
    % Get unique values and their indices
    try
        % Use builtin function
        [unique_values, ~, indices] = unique(cell_array_in(:, i_var), 'stable');
        % Mark variable type
        this_variable_type = 'string';
    
    catch ME
        if strcmpi(ME.identifier, 'MATLAB:UNIQUE:InputClass')  % it occurs when cells contain numeric values or content is of mixed data types            
            try
                % Transform cells to matrix
                [unique_values, ~, indices] = unique(cell2mat(cell_array_in(:, i_var)), 'stable');
                % Mark variable type
                this_variable_type = 'numeric';
                
            catch subME
                if strcmpi(subME.identifier, 'MATLAB:cell2mat:MixedDataTypes')  % it occurs when cells contain mixed data types
                    % Get data type of each value
                    values_in = cell_array_in(:, i_var);
                    data_types = cellfun(@class, values_in, 'UniformOutput',false);
                    % Convert everything to a string and then to character array
                    values_in = cellfun(@string, values_in, 'UniformOutput',false);
                    % Replace missing values with "NaN"
                    values_in(cellfun(@(x) ismissing(x), values_in)) = {"NaN"};
                    values_in = cellfun(@char, values_in, 'UniformOutput',false);
                    % Concatenate with data type and then select unique
                    % combinations of value and type
                    [~, ~, values_in_indices] = unique(values_in, 'stable');
                    [~, ~, data_types_indices] = unique(data_types, 'stable');
                    % The unique pairs of indices correspond the unique values
                    % in the mixed array
                    [~, ~, indices] = unique([values_in_indices, data_types_indices], 'stable', 'rows');
                    % Take position of each unique value
                    [~, first_index, ~] = unique(indices, 'stable');                    
                    unique_values = cell_array_in(first_index, i_var);
                    % Mark variable type
                    this_variable_type = 'mixed';
            
                else
                    rethrow(subME)
                end
            end
        else
            rethrow(ME)
        end
    end
    
    % Store information
    variable_map{i_var} = unique_values;
    cell_array_in_converted(:, i_var) = indices;
    variable_type{i_var} = this_variable_type;
end

% Get unique rows of indices array
unique_array_out = unique(cell_array_in_converted, 'rows');

% Use the index to retrieve the actual value
cell_array_out = cell(size(unique_array_out));
for i_var = 1:n_variables
    switch variable_type{i_var}
        case {'string', 'mixed'}
            cell_array_out(:, i_var) = variable_map{i_var}(unique_array_out(:, i_var));
            
        case 'numeric'
            cell_array_out(:, i_var) = num2cell(variable_map{i_var}(unique_array_out(:, i_var)));
    end
end

% Sort output, if requested
if strcmp(sorting, 'sorted')
    cell_array_out = sortrows(cell_array_out);
end


%% DEMO
function run_demo()
    % Make and show input
    cell_array_in = {'1', '1'; '1', '1'; '1', 1; '1', 1; '2', NaN; '2', NaN};
    fprintf('\nThis is a demonstration of the function <strong>uniqueCellRows</strong>\n\n')
    fprintf('Given this cell array with mixed data types in the second column:\n')
    disp(cell_array_in)
    
    % Compute and show result
    output = uniqueCellRows(cell_array_in);
    fprintf('\nThe following are the unique rows of it:\n');
    disp(output)

