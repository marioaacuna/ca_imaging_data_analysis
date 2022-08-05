animals = {'MA_34epi';...
    'MA_35epi';...
    'MA_36epi';...
    'MA_37epi';...
    'MA_38epi';...
    'MA_39epi';...
    'MA_40epi'};
METADATA_ALL = SQL_database.read_table_where('sessions', {'experiment_id', 'session_id', 'experimental_condition', 'stimulus'}, animals, {'animal_ID'});
METADATA_ALL = METADATA_ALL(ismember(METADATA_ALL.experiment, experiment), :);

session_id = METADATA_ALL.session_id(ismember(METADATA_ALL.stimulus, 'puff'));
experiment_id = METADATA_ALL.experiment_id(ismember(METADATA_ALL.stimulus, 'puff'));
stimulus_id_all = SQL_database.read_table_where('stimulations', 'stimulus_id', [experiment_id, session_id], {'experiment_id', 'session_id'}, 'return_as_table',false);


for iid = 1:length(animals)
    animal_ID = animals{iid};
    NEW_METADATA = SQL_database.read_epi_trials(animal_ID);
%     dates_of_interest = unique(NEW_METADATA.date(ismember(NEW_METADATA.session_id,session_id)));
    sessions_of_interest = unique(NEW_METADATA.session_id(ismember(NEW_METADATA.session_id, session_id) & ~isnan(NEW_METADATA.timestamps)));
    stimulus_id = SQL_database.read_table_where('stimulations', 'stimulus_id',  sessions_of_interest, {'session_id'}, 'return_as_table',false);

    SQL_database.update('stimulations', {'type'}, {'evoked'}, stimulus_id,'stimulus_id')
end
