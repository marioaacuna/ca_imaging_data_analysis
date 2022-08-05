function data_out = table2struct_of_arrays(data_in)

% Initialize output variable
data_out = struct();
% Get names of columns
column_names = data_in.Properties.VariableNames;
for iname = 1:length(column_names)
    % Extract values
    values = data_in.(column_names{iname});
    if ~iscell(values)
        values = num2cell(values);
    end
    % Concatenate numeric values in a matrix
    if isnumeric(values{1})
        if all(size(values{1})==1)
            data_out.(column_names{iname}) = cell2mat(values);
        else
            data_out.(column_names{iname}) = cell2mat(values')';
        end
        
    else  % Copy cells
        data_out.(column_names{iname}) = values;
    end
end
