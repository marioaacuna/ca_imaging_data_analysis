clear
global GC
clc, close all hidden
GC.experiment_name = 'ACC_CCI_anesth';
experiment_type = GC.experiment_name;
GC.response_decoding_stats_filename_all = 'response_decoding_stats_ACC.mat';
debug_mode = 1;
do_plotting = 1;
animals = animal_list();
animal_names = SQL_database.read_table_where('experiments', {'animal_ID', 'experimental_group'}, animals, 'animal_ID', 'return_as_table',false);
analysis_type_str = 'all';
analysis_type = ''; % leave empy for all
metrics = {'PR', 'ACC', 'F1', 'MCC', 'AUC'}; % select any of the options of the cm output. For AUPR run PR.
animals_to_discard = {'MA_14', 'MA_15'};
data_dir_root = ['V:\Ca_imaging_pain\6_data\ACC_CCI_anesth\response_decoding'];

for i_m = 1:length (metrics)
    
    metric = metrics{i_m};
%     stats_filename = get_filename_of(['response_decoding_stats_ACC_', analysis_type_str],experiment_type);
    stats_filename = os.path.join([GC.repository_root_path,'\', 'Figures_paper_FMP'],'_data', 'response_decoding_SCORES_stats', ['response_decoding_stats_', metric,'.mat']);
    if exist(stats_filename, 'file') == 2 && do_plotting
        disp (['metric ', metric, ' already analyzed']),
        [STATS_DECODING] = fig3_stimulus_decoding(stats_filename, metric, debug_mode);
%         continue,
    elseif  exist(stats_filename, 'file') == 1
        STATS_response_decoding_ACC = load_variable(stats_filename, STATS_response_decoding_ACC);
    else
        %%
        STATS_response_decoding_ACC = struct();
        for iid = 1 : length(animals)
            % Run score for conf matrx
            animal_ID = animals{iid};
            if ismember(animal_ID, animals_to_discard), continue, end
            disp(animal_ID)
            run_in_python('score_confusion_matrix', ['root_folder=',data_dir_root], ['analysis_type=', analysis_type], ['animal_list=', animal_ID], ['metric=', metric])
            % Run score for GENERALIZABILITY
            data_dir = ['V:\Ca_imaging_pain\6_data\ACC_CCI_anesth\response_decoding', '\',animal_ID,'\'];
            
            output_filename_results          = os.path.join(data_dir, [animal_ID, '_crossday_results_ACC.json']);
            output_filename_confusion_matrix = os.path.join(data_dir, [animal_ID, '_crossday_confusion_matrix_ACC.json']);
            % Run in Python
            run_in_python('decode_activations_different_day_ACC', ['folder=', data_dir], ['output_filename_results=', output_filename_results], ['output_filename_confusion_matrix=', output_filename_confusion_matrix], ['metric=', metric])
            
            % get ACCURACY
            % Read results
            txt = jsondecode(fileread(output_filename_results));
            % Convert output to arrays
            fn = fieldnames(txt);
            for i_field = 1:length(fn)
                values = txt.(fn{i_field});
                if isempty(values), continue, end
                values = cellfun(@(x) x', values, 'UniformOutput',false);
                values = vertcat(values{:});
                values = str2double(values);
                txt.(fn{i_field}) = values;
            end
            % Move results into output structure
            CLASSIFIER_GENERALIZABILITY = txt;
            STATS_response_decoding_ACC.(animal_ID) = CLASSIFIER_GENERALIZABILITY;
        end
        save(stats_filename, 'STATS_response_decoding_ACC', '-v7.3')
        
        %% Run the figure script
        [STATS_DECODING] = fig3_stimulus_decoding(stats_filename, metric, debug_mode);
        
        %%
        save(stats_filename, 'STATS_response_decoding_ACC','STATS_DECODING', '-v7.3')
        disp([metric, ' saved'])
        
    end
end
disp('saved and finished')