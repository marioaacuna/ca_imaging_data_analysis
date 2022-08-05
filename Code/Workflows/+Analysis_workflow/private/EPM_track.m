function done = EPM_track(animal_ID, missing_only)
% This function will run the python script to plot heatmaps ans trayectory
% of animal on the EPM. It also calculates proportion of activated cells in 
% different arms
done = true;

% Get logger and general_configs
global LOGGER GC
GC.experiment_name = 'ACC_SNI_anxiety';

% read Metadata
METADATA = SQL_database.read_epi_trials(animal_ID);
EPM_idx = find(ismember(METADATA.stimulus,{'EPM'}));
EPM_dates = cell(METADATA.date(EPM_idx));
root_path = os.path.join(GC.data_root_path, '6_data', GC.experiment_name, 'EPM_results', animal_ID);
LOGGER.trace(['Computing EPM tracking of ', mat2str(animal_ID)])

for i_sess = 1: length (EPM_dates)
    filename_this_sess = os.path.join(root_path, ['EPM_position_and_F_resampled_',EPM_dates{i_sess},'_', animal_ID, '.csv']);
    if exist(filename_this_sess) && missing_only
        disp ([animal_ID, 'session ', EPM_dates{i_sess}, ' already analysed'])
        continue
    else
        excel_file = ['T:\Mario\Behavior\EPM\raw_files\', EPM_dates{i_sess},'_', animal_ID,'.xlsx'];
        if exist(excel_file)
            run_in_python('EPM_track.py', cell2mat(EPM_dates(i_sess)), animal_ID, EPM_idx (i_sess))
        end
    end
end

