%% compute how many neuorons are activated on the open arm per session
clear, clc

global GC

do_saving = 1;
output_filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'proportions_EPM_responders.mat');

GC.experiment_name =  'ACC_SNI_anxiety';
% load parametrers
animal_list = unique(SQL_database.read_table_where ('sessions', {'animal_ID'}, 'EPM', 'stimulus', 'return_as_table', false), 'stable');
animals_to_use = animal_list(~ismember(animal_list, {'MA_28epi', 'MA_31epi'})); % MA28epi, died, MA_31, not good cells, MA_32 and MA_33, not ready yet (20.1.20)
METADATA_all_animals = SQL_database.read_table_where('sessions', {'+', 'animal_ID'}, animal_list, 'animal_ID');
%%
PERCENT_OPEN = struct ();

for iid = 1:length (animals_to_use)
    animal_ID = animals_to_use{iid};
    METADATA = SQL_database.read_epi_trials(animal_ID);
    dates_EPM = METADATA.date((ismember(METADATA.experiment(:), 'EPM')));
    
    sessions_EPM = METADATA.day_from_surgery((ismember(METADATA.experiment(:), 'EPM')));
    column_names = cell (1,length(sessions_EPM));
    for i_sess = 1 : length(sessions_EPM)
        column_names{1, i_sess} = make_variable_name(['session_',sessions_EPM{i_sess}]);
    end
    n_sessions = length (dates_EPM);
    n_cells = SQL_database.read_table_where('experiments',{'n_ROIs'}, animal_ID, 'animal_ID', 'return_as_table', false);
    percent_cells_open_arm = cell(1, n_sessions);
    for i_sess = 1: n_sessions
        date = dates_EPM{i_sess};
        
        % load traces
        trials_idx = find(ismemberCellRows(METADATA.date, {date}) & ismemberCellRows(METADATA.type(:), {'evoked'})) ;
%         data = all_traces_this_animal(:, trials_idx);
        
        % load session
        day_from_surgery = cell2mat(unique (METADATA.day_from_surgery(trials_idx,: )));
        selected_session =  make_variable_name(['session_',day_from_surgery ]);
        
        % Load EPM data
        filename_EPM_data = os.path.join(GC.data_root_path, GC.aggregated_data_path,GC.experiment_name, ['EPM_results\',animal_ID,'\epm_idx_cells_open_',date,'_', animal_ID,'.json']);
        if ~exist(filename_EPM_data), disp (['EPM data animal ', animal_ID, ', session ', date, ' has not been analysed yet']),continue, end %#ok<EXIST>
        txt = jsondecode(fileread(filename_EPM_data));
        cells = txt.idxCellsRespondingInOpen;
        cells_idx = fieldnames(cells);
        n_EPM_cells = length (cells_idx);
        
        percent_cells_open_arm {1,i_sess} = n_EPM_cells ./ n_cells;
    end
    percent_cells_open_arm_table = cell2table(cell(percent_cells_open_arm), 'VariableNames', column_names);
    PERCENT_OPEN.(animal_ID) = percent_cells_open_arm_table;  
end

if do_saving
    save (output_filename, 'PERCENT_OPEN',  '-v7.3')
    disp(['saved in ', output_filename])
end
disp ('done!')
