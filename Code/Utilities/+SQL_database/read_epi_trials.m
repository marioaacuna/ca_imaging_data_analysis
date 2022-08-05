function METADATA = read_epi_trials(animal_ID, varargin)

p = inputParser();
addParameter(p, 'date', '');
addParameter(p, 'stimulus', '');
addParameter(p, 'type', '');
addParameter(p, 'valid', true);
parse(p, varargin{:});
date = p.Results.date;
stimulus = p.Results.stimulus;
type = p.Results.type;
valid = p.Results.valid;

% Make table of trials
try
    METADATA = SQL_database.read_table_where('stimulations', {'session_id', 'date', 'stimulus', 'timestamps', 'response', 'type', 'valid', 'affective', 'experimental_condition'}, animal_ID,'animal_ID');
    METADATA = Metadata_workflow.flatten_stimulations_table(METADATA);
catch
    METADATA = SQL_database.read_table_where('sessions', {'+'}, animal_ID, 'animal_ID');
    return
end
conditions_with_timestamps = uniqueCellRows(METADATA{:, {'date', 'stimulus'}});

% Get conditions without timestamps
METADATA_sessions = SQL_database.read_table_where('sessions', {'date', 'stimulus'}, animal_ID,'animal_ID');
METADATA_done_prepro = logical(table2array(SQL_database.read_table_where('sessions',{'done_preprocessing'}, animal_ID,'animal_ID')));
all_conditions = uniqueCellRows(METADATA_sessions{METADATA_done_prepro, {'date', 'stimulus'}});
missing_conditions = all_conditions(~ismemberCellRows(all_conditions, conditions_with_timestamps), :);
% Add these conditions
if ~isempty(missing_conditions) % this is in case there is no SP in the session (no timestamped stim)
    missing_conditions(:, 3) = {animal_ID};
    METADATA_missing_conds = SQL_database.read_table_where('sessions', {'session_id', 'date', 'stimulus', 'experimental_condition'}, missing_conditions,{'date', 'stimulus', 'animal_ID'});
    n_missing_conds = height(METADATA_missing_conds);
    METADATA_missing_conds.timestamps = NaN(n_missing_conds, 2);
    METADATA_missing_conds.response = NaN(n_missing_conds, 1);
    METADATA_missing_conds{:, 'type'} = {'other'};
    METADATA_missing_conds.type(ismember(METADATA_missing_conds.stimulus, 'SP')) = {'spontaneous'};
    METADATA_missing_conds.valid = true(n_missing_conds, 1);
    METADATA_missing_conds{:, 'affective'} = {''};
    % Concatenate tables
    METADATA = [METADATA; METADATA_missing_conds];
    
end
METADATA = sortrows(METADATA, {'date', 'timestamps'});
% Filter table
rows_to_keep = true(height(METADATA), 1);
if ~isempty(date)
    rows_to_keep = rows_to_keep & ismember(METADATA.date, date);
end
if ~isempty(stimulus)
    rows_to_keep = rows_to_keep & ismember(METADATA.stimulus, stimulus);
end
if ~isempty(type)
    rows_to_keep = rows_to_keep & ismember(METADATA.type, type);
end
if ~isempty(valid)
    rows_to_keep = rows_to_keep & (METADATA.valid == valid);
end
METADATA = METADATA(rows_to_keep, :);
