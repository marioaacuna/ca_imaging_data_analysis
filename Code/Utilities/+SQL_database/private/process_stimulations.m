function metadata = process_stimulations(metadata)

if isempty(metadata), return, end

global GC

% Get names of all columns
column_names = metadata.Properties.VariableNames;

% Convert timestamps to numbers
if ismember('timestamps', column_names)
    for i_row = 1:size(metadata, 1)
        timestamps_ts = metadata{i_row, 'timestamps'};
        timestamps_ts = timestamps_ts{1};
        if ~isempty(timestamps_ts)
            % Decode string
            x = jsondecode(timestamps_ts);
            x = [x{:}]';
            % Convert to numbers
            x = str2double(x);
        else
            x = zeros(0, 2);
        end
        metadata{i_row, 'timestamps'} = {x};
    end
end

% Convert response to logical
if ismember('response', column_names)
    for i_row = 1:size(metadata, 1)
        response = metadata{i_row, 'response'};
        response = response{1};
        if ~isempty(response)
            % Convert to logical
            x = logical(str2double(cellstr(response(:))));
        else
            x = false(0);
        end
        metadata{i_row, 'response'} = {x};
    end
end

% Convert valid to logical
if ismember('valid', column_names)
    for i_row = 1:size(metadata, 1)
        value = metadata{i_row, 'valid'};
        value = value{1};
        if ~isempty(value)
            % Convert to logical
            x = logical(str2double(cellstr(value(:))));
        else
            x = false(0);
        end
        metadata{i_row, 'valid'} = {x};
    end
end

% Convert affective to strings
if ismember('affective', column_names)
    all_affective_responses = [{''}, GC.freely_moving_affective_responses(:)'];
    for i_row = 1:size(metadata, 1)
        affective_response = metadata{i_row, 'affective'};
        affective_response = affective_response{1};
        if ~isempty(affective_response)
            x = str2double(cellstr(affective_response(:)));
            % Convert to logical
            x = all_affective_responses(x)';
        else
            x = repmat({''}, size(metadata.timestamps{i_row}, 1), 1);
        end
        % Convert to number
        metadata{i_row, 'affective'} = {x};
    end
end

