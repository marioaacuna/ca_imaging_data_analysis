function trial_frames = get_trial_frames(animal_ID, varargin)

% Parse inputs
p = inputParser();
addParameter(p, 'reset_index', false, @(x) islogical(x) || (isnumeric(x) && length(x(:))==1 && ismember(x, [0,1])))
addParameter(p, 'return_all_sessions', false, @(x) islogical(x) || (isnumeric(x) && length(x(:))==1 && ismember(x, [0,1])))
% Parse inputs and assign defaults
parse(p, varargin{:})
% Retrieve inputs
reset_index = p.Results.reset_index;
return_all_sessions = logical(p.Results.return_all_sessions);

% Read metadata
METADATA = SQL_database.read_table_where('trials', {'n_frames', 'frames_idx', 'date', 'recording_time'}, animal_ID,'animal_ID', 'return_all_sessions',return_all_sessions);

% Convert indices from string to array
if reset_index
    trial_frames = cumsum(METADATA{:, 'n_frames'});
    trial_frames = [[1; trial_frames(1:end-1)+1], trial_frames(:)];
else
    % Convert strings to cells
    trial_frames = cellfun(@(x) regexp(x,',','split'), METADATA.frames_idx, 'UniformOutput',false);
    % Concatenate cells and convert strings to numbers
    trial_frames = cellfun(@(x) [str2double(x{1}), str2double(x{2})], trial_frames, 'UniformOutput',false);
    % Convert to array
    trial_frames = cell2mat(trial_frames);
end
