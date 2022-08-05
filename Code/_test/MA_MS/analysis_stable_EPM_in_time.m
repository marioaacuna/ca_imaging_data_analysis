function analysis_stable_EPM_in_time ()

%%
clear, clc

global GC LOGGER

saving = 1;
% type_analysis = '';

GC.experiment_name =  'ACC_SNI_anxiety';
% load parametrers
animal_list = unique(SQL_database.read_table_where ('sessions', {'animal_ID'}, 'EPM', 'stimulus', 'return_as_table', false), 'stable');
animals_to_use = animal_list(~ismember(animal_list, {'MA_28epi', 'MA_31epi'})); % MA28epi, died, MA_31, not good cells, MA_32 and MA_33, not ready yet (20.1.20)
% METADATA_all_animals = SQL_database.read_table_where('sessions', {'+', 'animal_ID'}, animal_list, 'animal_ID');
EPM_stable = struct();

% output
root_path = [GC.data_root_path ,'\',GC.aggregated_data_path,'\', GC.experiment_name, '\', 'percent_stable_cells.mat\'];
if ~exist (root_path, 'dir')
    mkdir(root_path)
end



for iid = 1:length(animals_to_use)
    %Load EPM data
    animal_ID = animals_to_use{iid};
    LOGGER.info(['computing proportion of stable EPM cells of ', animal_ID])
%     n_cells = SQL_database.read_table_where('experiments',{'n_ROIs'}, animal_ID, 'animal_ID', 'return_as_table', false);
    METADATA = SQL_database.read_epi_trials(animal_ID);
    dates_EPM = METADATA.date((ismember(METADATA.experiment(:), 'EPM')));
    % look through files to see if it exists
    for i_date = length(dates_EPM)
        filename_EPM_data = os.path.join(GC.data_root_path, GC.aggregated_data_path,GC.experiment_name, ['EPM_results\',animal_ID,'\epm_idx_cells_open_',date,'_', animal_ID,'.json']);
        if ~exist(filename_EPM_data), dates_EPM(i_date) = []; end %#ok<EXIST>
    end
    
    n_sessions = length (dates_EPM);
    session_pairs =  nchoosek(1:length(dates_EPM), 2);
    n_pair = size(session_pairs,1);
    if n_pair > 1, keyboard, end
    
    percent_stable  = cell (1,n_pair);
    for i_pair = 1%: n_pair_sess
        session1 = dates_EPM(session_pairs(i_pair, 1));
        session2 = dates_EPM(session_pairs(i_pair, 2));
        cells_all_sess = cell(1, n_sessions);
        EPM_all_cells = 0;
        for i_sess = 1:n_sessions
            % Load EPM data
            filename_EPM_data = os.path.join(GC.data_root_path, GC.aggregated_data_path,GC.experiment_name, ['EPM_results\',animal_ID,'\epm_idx_cells_open_',dates_EPM{i_sess},'_', animal_ID,'.json']);
            if ~exist(filename_EPM_data), disp (['EPM data animal ', animal_ID, ', session ', dates_EPM(i_sess), ' has not been analysed yet']),continue, end %#ok<EXIST>
            txt = jsondecode(fileread(filename_EPM_data));
            cells = txt.idxCellsRespondingInOpen;
            cells_idx = fieldnames(cells);
            n_EPM_cells = length (cells_idx);
            cells_id_EPM = NaN (n_EPM_cells,1);
            for i_cell = 1:n_EPM_cells
                cells_id_EPM (i_cell) = str2num(cell2mat( strsplit((cells.(cell2mat(cells_idx(i_cell)))), 'cell_') ));
            end
            cells_all_sess{i_sess} = cells_id_EPM;
            EPM_all_cells = sum([EPM_all_cells,n_EPM_cells]);
        end
        percent_stable {1,i_pair} =(length(intersect(cells_all_sess{:})) ./ EPM_all_cells ).*100;
    end
    EPM_stable.(animal_ID) = percent_stable;
    
end

if saving
    output_filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'proportions_EPM_stable.mat');
    save (output_filename, 'EPM_stable',  '-v7.3')
    disp(['saved in ', output_filename])
    
end

%%
clc
all_points = NaN (length (animals_to_use),1);
figure ('Color', 'w', 'Pos',[103, 94, 200, 300])
for iid = 1 : length(animals_to_use)
    this_point = cell2mat(EPM_stable.(animals_to_use {iid}));
    plot(this_point, 'o', 'MarkerFaceColor', 'r')
    hold on
    all_points (iid)  = this_point;
end
plot (mean(all_points, 1), 'o', 'MarkerFaceColor', 'k')
ylim ([0, 50])
xlabel ('Stable', 'FontSize', 11)
errorbar (mean(all_points, 1),sem(all_points, 1) , 'b')
xticks ([])
xlim ([.9, 1.1])
ylabel ('% Stable cells accross sessions', 'FontSize', 11)


    
    
