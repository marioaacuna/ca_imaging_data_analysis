function write_data_to_excel(data_to_write, filename, sheet_name)

% Clear sheet
try
    [~, ~, Raw] = xlsread(filename, sheet_name);
    [Raw{:, :}] = deal(NaN);
    xlswrite(filename, Raw, sheet_name);
end

% Convert tables to cell arrays
if istable(data_to_write)
    data_to_write = [data_to_write.Properties.VariableNames; table2cell(data_to_write)];
end

% Write new data
xlswrite(filename, data_to_write, sheet_name);
