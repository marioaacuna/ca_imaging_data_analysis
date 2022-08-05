function data = flatten_stimulations_table(metadata)

% Make sure input is correct
if ~istable(metadata)
    error('Metadata_workflow:flatten_stimulations_table', 'input must be a table')
end

column_names = metadata.Properties.VariableNames;
has_timestamps = ismember('timestamps', column_names);
has_response = ismember('response', column_names);
has_valid = ismember('valid', column_names);
has_affective = ismember('affective', column_names);
if ~any([has_timestamps, has_response, has_valid, has_affective])
	error('Metadata_workflow:flatten_stimulations_table', 'input table must contain at least one of the following columns: timestamps, response, valid, affective')
end
other_columns_idx = ~ismember(column_names, {'timestamps', 'response', 'valid', 'affective'});
timestamps_idx = ismember(column_names, 'timestamps');
response_idx = ismember(column_names, 'response');
valid_idx = ismember(column_names, 'valid');
affective_idx = ismember(column_names, 'affective');

% Get total number of stimuli
if has_timestamps
    use_column = 'timestamps';
elseif has_response
    use_column = 'response';
elseif has_valid
    use_column = 'valid';
elseif has_affective
    use_column = 'affective';
end
n_events_per_stimulus = cellfun(@(x) size(x, 1), metadata{:, use_column});
n_stimuli = size(n_events_per_stimulus, 1);

% Find edges of each row in final table
row_idx = cumsum(n_events_per_stimulus);
row_idx = [[1; row_idx(1:end-1) + 1], row_idx];

%%
data = cell(n_stimuli, length(column_names));
for i_stim = 1:n_stimuli
    this_event = table2cell(metadata(i_stim, :));
    % Simply replicate all columns apart from those with timestamps
    this_data = cell(n_events_per_stimulus(i_stim), length(column_names));
    this_data(:, other_columns_idx) = repmat(this_event(other_columns_idx), n_events_per_stimulus(i_stim), 1);
    % Unpack columns
    if has_timestamps
        this_data(:, timestamps_idx) = mat2cell(this_event{timestamps_idx}, ones(n_events_per_stimulus(i_stim), 1), 2);
    end    
    if has_response
        this_data(:, response_idx) = num2cell(this_event{response_idx});
    end
    if has_valid
        this_data(:, valid_idx) = num2cell(this_event{valid_idx});
    end
    if has_affective
        this_data(:, affective_idx) = this_event{affective_idx};
    end
    % Concatenate results
    data(row_idx(i_stim, 1):row_idx(i_stim, 2), :) = this_data;
end

% Convert to table
data = cell2table(data, 'VariableNames',column_names);

