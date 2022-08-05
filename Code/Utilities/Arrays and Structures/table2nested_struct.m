function struct_out = table2nested_struct(table_in, variables_order, value)

if nargin == 0, demo, end

all_columns = table_in.Properties.VariableNames;
if ~exist('variables_order', 'var') || isempty(variables_order) 
    % Use the order of appearance in the table
    variables_order = all_columns;
    
else
    keyboard
end

% Get the value to introduce at the deepest level
if ~exist('value', 'var') || isempty(value) 
    value = [];
end

% Arrange columns
table_in = table_in(:, variables_order);
n_rows = height(table_in);

struct_out = struct();
for i_row = 1:n_rows
    this_row = table_in{i_row, :};
    struct_out = setfield(struct_out, this_row{:}, value);
end



function demo()
    table_in = cell(9, 3);
    table_in(:, 1) = {'a','a','a','a','a','b','b','b','b'};
    table_in(:, 2) = {'x','x','y','y','z','x','y','y','y'};
    table_in(:, 3) = {'i','j','j','k','i','j','i','j','k'};
    table_in = cell2table(table_in, 'VariableNames',{'A','B','C'});
