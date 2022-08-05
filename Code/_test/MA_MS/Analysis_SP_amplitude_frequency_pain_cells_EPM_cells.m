%%
global GC


do_export_file = 1;

disp('Loading event statistics')

GC.experiment_type = 'ACC_SNI_anxiety';
stats_filename = get_filename_of('spontaneous_activity_stats_epi', GC.experiment_type);
EVENT_STATS = load_variable(stats_filename, 'EVENT_STATS_spontaneous_activi_sessy');
animals = fieldnames(EVENT_STATS);
n_animals = length (animals);
experimental_groups = {'pain';'EPM'};
[~, ~, idx] = unique(experimental_groups(:,1), 'stable');
experimental_groups(:, 2) = num2cell(idx);

average_type = 'mean';
sessions_to_analyse = 2;
output_folder = GC.temp_dir;% ['M:\Thomas\Manuscripts\in preparation\Invivo ACC\Material for Figures\check_points\'];


% gather data
AMPLITUDE_EPM = cell(n_animals, 2);
AMPLITUDE_pain = cell(n_animals, 2);
FREQUENCY_EPM = cell(n_animals, 2);
FREQUENCY_pain = cell(n_animals, 2);
TIMEPOINTS = cell(n_animals, 3);
TIMEPOINTS(:, 1) = animals;
TIMEPOINTS(:, 2) = {'all'};
%%

for iid = 1: n_animals
    
    animal_ID = animals{iid};
    
    % Load SP data
    event_stats = EVENT_STATS.(animal_ID);
    if height(event_stats) == 0, continue, end
    METADATA = SQL_database.read_epi_trials (animal_ID);
%   dates = unique(METADATA.date, 'stable');
    dates_EPM = METADATA.date(ismember (METADATA.stimulus, 'EPM'));
    day_from_surgery = METADATA.day_from_surgery(ismember (METADATA.stimulus, 'EPM'));
    sessions_to_keep = find(cell2mat(dates_EPM), sessions_to_analyse, 'first');
    dates_EPM  = dates_EPM (sessions_to_keep);
    TIMEPOINTS{iid, 3} = day_from_surgery(sessions_to_keep);

    % Load  pain data
    
    pain_result_filename = get_filename_of('response_detection_epi', animal_ID);
    RESULTS = load_variable(pain_result_filename, 'RESULTS');
    session_fieldnames = fieldnames(RESULTS);
    selected_session =  session_fieldnames(sessions_to_keep);
    n_sess = length (selected_session);
%     session_idx = find(ismember(session_fieldnames, sessions_to_keep));
    trials_idx = find(ismember(METADATA.date, dates_EPM) & ismemberCellRows(METADATA.type(:), {'evoked'}));%& ismember(METADATA.day_from_surgery, day_from_surgery(sessions_to_keep)));
%     trials_idx = find(ismemberCellRows(METADATA.day_from_surgery, day_from_surgery(sessions_to_keep)) & ismemberCellRows(METADATA.type(:), {'evoked'})) ;
    stimuli = unique(METADATA.stimulus(trials_idx,:));
    n_stim = length(stimuli);
    excited_pain = cell(n_sess, n_stim);
    for i_sess = 1 : n_sess
        for i_stim = 1 : n_stim
            excited_pain{i_sess, i_stim} = RESULTS.(selected_session{i_sess}).(stimuli{i_stim}).selectivity.excited;
        end
    end
    
    % Load EPM data
    cells_id_EPM = cell (n_sess,1);
    for i_sess = 1 : n_sess
        
        filename_EPM_data = os.path.join(GC.data_root_path, GC.aggregated_data_path,GC.experiment_name, 'EPM_results',animal_ID,['epm_idx_cells_open_',dates_EPM{i_sess},'_',animal_ID,'.json']);
        txt = jsondecode(fileread(filename_EPM_data));
        cells = txt.idxCellsRespondingInOpen;
        cells_idx = fieldnames(cells);
        n_EPM_cells = length (cells_idx);
        cells_id_EPM_this_sess = cell (n_EPM_cells,1);
        for i_cell = 1:n_EPM_cells
            cells_id_EPM_this_sess {i_cell, 1} = str2num(cell2mat( strsplit((cells.(cell2mat(cells_idx(i_cell)))), 'cell_') ));
        end
        cells_id_EPM {i_sess} = cells_id_EPM_this_sess;
               
    end
    
    n_cells = size(EVENT_STATS.(animal_ID),1);
    
%     group_idx = experimental_groups{ismember(experimental_groups(:, 1), TIMEPOINTS(iid, 2)), 2};
    mean_event_amplitude = cell2mat(cellfun(@(x) x.', event_stats{:, [average_type, '_event_ampli_sessude']}, 'UniformOutput',false));
    mean_event_amplitude(mean_event_amplitude == 0) = NaN;
    mean_event_freq =  cell2mat(cellfun(@(x) x.', event_stats{:, ['mean', '_event_rate']}, 'UniformOutput',false));

    for i_sess = 1: n_sess
        cells_id_pain = find(sum ([excited_pain{i_sess, :}], 2));
        AMPLITUDE_EPM{iid, i_sess} = mean_event_amplitude(cell2mat(cells_id_EPM{i_sess}));
        AMPLITUDE_pain{iid, i_sess} = mean_event_amplitude(cells_id_pain);
        FREQUENCY_EPM{iid, i_sess} = mean_event_freq(cell2mat(cells_id_EPM{i_sess}));
        FREQUENCY_pain{iid, i_sess} = mean_event_freq(cells_id_pain);
    end
end
disp ('done loading')
%% Re arrange data

DATA_amplitude_EPM = cell (1,n_sess);
DATA_amplitude_pain = cell (1,n_sess);
DATA_frequency_EPM = cell (1,n_sess);
DATA_frequency_pain = cell (1,n_sess);
for i_sess = 1:n_sess
    DATA_amplitude_EPM{1,i_sess} = cell2mat(AMPLITUDE_EPM(:,i_sess));
    DATA_amplitude_pain{1,i_sess}  = cell2mat(AMPLITUDE_pain(:,i_sess));
    DATA_frequency_EPM{1,i_sess}  = cell2mat(FREQUENCY_EPM(:,i_sess));
    DATA_frequency_pain{1,i_sess}  = cell2mat(FREQUENCY_pain(:,i_sess));
    
end

if do_export_file
    disp('put the folder')
    keyboard
end