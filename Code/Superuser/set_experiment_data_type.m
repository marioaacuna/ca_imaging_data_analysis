function set_experiment_data_type(data_type, animal_ID)
%ADD_TWO_PHOTON_MINISCOPE_KIND Inserts data type into database for a given
%experiment.

% Clear Command Window
clc

% Log progress
disp(['Will assign ''', data_type, ''' to mouse ''' animal_ID, ''''])

% Update database
SQL_database.update('experiments', 'data_type', {data_type}, animal_ID, 'animal_ID')
disp('done')

