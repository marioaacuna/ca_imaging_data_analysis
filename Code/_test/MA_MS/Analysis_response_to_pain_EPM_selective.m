%% INFO
% this script analyses the AUC of the evoked response period of neurons
% that are actictivated (excited) by the stim as well as neuros that have higher
% activity in the open arm of the EPM. 
%%
clear, clc

global GC LOGGER

saving = 0;
type_analysis = '';

GC.experiment_name =  'ACC_SNI_anxiety';
% load parametrers
animal_list = unique(SQL_database.read_table_where ('sessions', {'animal_ID'}, 'EPM', 'stimulus', 'return_as_table', false), 'stable');
animals_to_use = animal_list(~ismember(animal_list, {'MA_28epi', 'MA_31epi'})); % MA28epi, died, MA_31, not good cells, MA_32 and MA_33, not ready yet (20.1.20)
METADATA_all_animals = SQL_database.read_table_where('sessions', {'+', 'animal_ID'}, animal_list, 'animal_ID');
%%
AUC_ANALYSIS = struct ();
for iid = 1:length (animals_to_use)
    animal_ID = animals_to_use{iid};
    LOGGER.info(['computing EPM + stimuli responders of ', animal_ID])
    n_cells = SQL_database.read_table_where('experiments',{'n_ROIs'}, animal_ID, 'animal_ID', 'return_as_table', false);
    METADATA = SQL_database.read_epi_trials(animal_ID);
    frame_rate = GC.epifluorescence_downsample_to_frame_rate;
    dates = unique(SQL_database.read_table_where('sessions', {'date'}, animal_ID, 'animal_ID'));
    filename_traces = get_filename_of('dFF', animal_ID);
    all_traces_this_animal = load_variable (filename_traces, 'dFF');
    dates_EPM = METADATA.date((ismember(METADATA.experiment(:), 'EPM')));
    n_sessions = length (dates_EPM);
    % Load RESULTS
    pain_result_filename = get_filename_of('response_detection_epi', animal_ID);
    if ~exist(pain_result_filename), disp (['dFF of ', animal_ID, ',  has not been analysed yet']),continue, end %#ok<EXIST>

    RESULTS = load_variable(pain_result_filename, 'RESULTS');
    session_fieldnames = fieldnames(RESULTS);
    % output
    root_path = [GC.data_root_path ,'\',GC.aggregated_data_path,'\', GC.experiment_name, '\', 'analysis_pain_EPM\'];
    if ~exist (root_path, 'dir')
        mkdir(root_path)
    end
    
    analysis_filename = get_filename_of ('analysis_AUC');
   
    
    %%
    for i_sess = 1: n_sessions
        date = dates_EPM{i_sess};
        
        % load traces
        trials_idx = find(ismemberCellRows(METADATA.date, {date}) & ismemberCellRows(METADATA.type(:), {'evoked'})) ;
        data = all_traces_this_animal(:, trials_idx);
        
        % load session
        day_from_surgery = cell2mat(unique (METADATA.day_from_surgery(trials_idx,: )));
        selected_session =  make_variable_name(['session_',day_from_surgery ]);
        
        % Load  pain data
        session_idx = find(ismember(session_fieldnames, selected_session));
        stimuli = unique(METADATA.stimulus(trials_idx,:));
        n_stim = length(stimuli);
        selectivity_pain = cell(1, n_stim);
        
        for i_stim = 1 : n_stim
            selectivity_pain{:, i_stim} = RESULTS.(selected_session).(stimuli{i_stim}).selectivity.excited;
        end
        
        % load EPM data
        filename_EPM_data = os.path.join(GC.data_root_path, GC.aggregated_data_path,GC.experiment_name, ['EPM_results\',animal_ID,'\epm_idx_cells_open_',date,'_', animal_ID,'.json']);
        if ~exist(filename_EPM_data), disp (['EPM data animal ', animal_ID, ', session ', date, ' has not been analysed yet']),continue, end %#ok<EXIST>
        txt = jsondecode(fileread(filename_EPM_data));
        cells = txt.idxCellsRespondingInOpen;
        cells_idx = fieldnames(cells);
        n_EPM_cells = length (cells_idx);
        cells_id_EPM = cell (n_EPM_cells,1);
        for i_cell = 1:n_EPM_cells
            cells_id_EPM {i_cell, :}= str2num(cell2mat( strsplit((cells.(cell2mat(cells_idx(i_cell)))), 'cell_') ));
            
        end
        
        % get the intersect cells
        idx_cells_intersect = cell (n_stim,1);
        for i_stim = 1 : n_stim
            idx_cells_intersect{i_stim} = intersect(find(selectivity_pain{i_stim}), cell2mat(cells_id_EPM(:,1)));
        end
        
        % analysis
        
        for i_stim = 1: n_stim
            cells_idx_this_stim = selectivity_pain{i_stim};
            cells_intersect_this_stim = cell2mat(idx_cells_intersect(i_stim));
            name_stim = stimuli{i_stim};
            n_int_cells = length(cells_intersect_this_stim);
            trials_idx = ismember(METADATA.stimulus, stimuli(i_stim)) & ismember(METADATA.date, date);
            METADATA_this_stim = METADATA(trials_idx,:);
            data_this_stim = all_traces_this_animal(:, trials_idx);
            
            trial_duration = GC.miniscope_trial_duration.(name_stim);
            baseline_timestamp = 1: (frame_rate * trial_duration(1));
            evoked = (1 + (frame_rate * trial_duration (1)))  : ((trial_duration(2) + trial_duration(1)) * frame_rate);
            
            n_frames = cellfun(@length, data_this_stim);
            n_frames = unique(n_frames(:));
            time_axis = linspace(0,15,n_frames);
            data_to_analyse = NaN(n_cells, n_frames);
            for i_cell = 1:n_cells
                data_this_cell = cell2mat(data_this_stim(i_cell, :).');
                data_this_cell_r = zeros (size(data_this_cell));
                
                for i_trial = 1 : size(data_this_cell,1)
                    data_this_cell_r (i_trial,:) =  data_this_cell(i_trial,:) - median(data_this_cell(i_trial,:));
                end
                
                mean_baseline = nanmean(data_this_cell_r(:, baseline_timestamp), 2);
                std_baseline = std(data_this_cell_r(:, baseline_timestamp, :), 0, 2);
                data_this_cell_z = (data_this_cell_r - mean_baseline) ./ std_baseline;
                
                data_this_cell_z = nanmean(data_this_cell_z, 1);
                data_to_analyse(i_cell, :) = data_this_cell_z;
                %             data_to_analyse (i_cell, :) = data_this_cell_dFF; % data to analyse is dFF
                
            end
            
            % analise AUC
            AUC = table ();
            
            % responsive cells
            data_dFF = data_to_analyse(cells_idx_this_stim == 1, :); % reponsive cells
            n_these_cells = size (data_dFF,1);
            data_smooth = NaN(n_these_cells, n_frames);
            data_avg = NaN(n_cells, 1);
            if ~isempty(data_dFF)
                for i_cell = 1:n_these_cells
                    data_smooth(i_cell,:) = smooth(data_dFF(i_cell,:), 1);
                    data_avg (i_cell) = sum(data_smooth(i_cell, evoked),2);%max(data_smooth(i_cell, evoked));
                end
            end
            AUC.resp = data_avg;
            AUC_ANALYSIS.(animal_ID).(selected_session).(name_stim) = AUC;
            
            
            cell_int_logic = zeros(n_cells,1);
            cell_int_logic(cells_intersect_this_stim) = 1;
            
            % cells intersect EPM & this stim
            data_dFF = data_to_analyse(cell_int_logic == 1, :);
            n_these_cells = size (data_dFF,1);
            data_smooth = NaN(n_these_cells, n_frames);
            data_avg = NaN(n_cells, 1);
            if ~isempty(data_dFF)
                for i_cell = 1:n_these_cells
                    data_smooth(i_cell,:) = smooth(data_dFF(i_cell,:), 1);
                    data_avg (i_cell) = sum(data_smooth(i_cell, evoked),2);%max(data_smooth(i_cell, evoked));
                end
            end
            AUC.intersect = data_avg;
            AUC_ANALYSIS.(animal_ID).(selected_session).(name_stim) = AUC;
            
            % non_resp cells
            data_dFF = data_to_analyse(cells_idx_this_stim == 0, :);
            n_these_cells = size (data_dFF,1);
            data_smooth = NaN(n_these_cells, n_frames);
            data_avg = NaN(n_cells, 1);
            if ~isempty(data_dFF)
                for i_cell = 1:n_these_cells
                    data_smooth(i_cell,:) = smooth(data_dFF(i_cell,:), 1);
                    data_avg (i_cell) = sum(data_smooth(i_cell, evoked),2);%max(data_smooth(i_cell, evoked));
                end
            end
            AUC.no_resp = data_avg;
            AUC_ANALYSIS.(animal_ID).(selected_session).(name_stim) = AUC;
            
            % EPM responding cells
            cell_EPM_logic = zeros(n_cells,1);
            cell_EPM_logic(cell2mat(cells_id_EPM(:,1))) = 1;
            
            data_dFF = data_to_analyse(cell_EPM_logic == 1, :); 
            
            n_these_cells = size (data_dFF,1);
            data_smooth = NaN(n_these_cells, n_frames);
            data_avg = NaN(n_cells, 1);
            if ~isempty(data_dFF)
                for i_cell = 1:n_these_cells
                    data_smooth(i_cell,:) = smooth(data_dFF(i_cell,:), 1);
                    data_avg (i_cell) = sum(data_smooth(i_cell, evoked),2);%max(data_smooth(i_cell, evoked));
                end
            end
            
            AUC.open_arm_cells = data_avg;
            AUC_ANALYSIS.(animal_ID).(selected_session).(name_stim) = AUC;
        end
    end
    
  
    if saving
        
        LOGGER.info (['saving AUC analysis for EPM and pain data of ', animal_ID])
        save(analysis_filename, 'AUC_ANALYSIS', '-v7.3')
    end
    
    disp ('done!')
end
disp ('done with all')