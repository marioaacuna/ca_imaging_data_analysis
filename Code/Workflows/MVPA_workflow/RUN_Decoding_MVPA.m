function [CLASSIFIER_GENERALIZABILITY, concat_cm, days, info] = RUN_Decoding_MVPA(animal_ID, spikes, METADATA, STIMULI, allowed_stimuli, CATEGORIES,  experiment_type, do_random_labels,analysis_type_str ,this_cell_subset_info, this_cell_subset, classifier)

    %% Run MVPA for 2p calcium imaging data
    % Downloaded from https://github.com/treder/MVPA-Light
    % citation :  https://www.frontiersin.org/articles/10.3389/fnins.2020.00289/full
    % This script runs Multivariate pattern analysis (MVPA) and crossvalidation
    % for k folds
    convert_data_to_gaussian = 1;
    if ~ispc
        %     load('/Users/marioacuna/Desktop/Test_classifier_DATA_PARAMETERS.mat')
        addpath(genpath('/Users/marioacuna/Documents/MATLAB/MVPA-Light-master'))
    end
%     do_random = do_random_labels;
    skip_same_day_decoding = 0;
    skip_diff_day_decoding = 1;

    % Get LOGGER and general_configs
    global LOGGER GC CFG_MVPA
      
    % Get options from general_configs
    time_bin_activity = GC.response_decoding_window_width;
    n_permutations_significance = GC.decoding_n_permutations_significance;
    baseline_window = GC.response_decoding_baseline_window;

    % Get the number of ROIs
    n_ROIs = size(spikes, 1);
    % Get frame rate
    frame_rate = max(METADATA.frame_rate);
    time_bin_activity_n_frames = round(time_bin_activity * frame_rate);
    baseline_window_n_bins = round(baseline_window / time_bin_activity);

    % Get the name of each condition
    stimuli_name = fieldnames(STIMULI);
    columns_experimental_condition = GC.columns_experimental_condition.(GC.experiment_name);
    if size(columns_experimental_condition, 1) > 1
        criteria = cellfun(@(x) strsplit(x, '='), columns_experimental_condition(:, 1), 'UniformOutput',false);
        criteria = vertcat(criteria{:});
        these_columns_experimental_condition = {};
        for i_rule = 1:size(criteria, 1)
            criterion_to_check = criteria{i_rule, 1};
            value = SQL_database.read_table_where('experiments', criterion_to_check, {animal_ID, experiment_type},{'animal_ID', 'experiment_name'}, 'return_as_table',false);
            if ischar(value)
                did_match = strcmp(value, criteria{i_rule, 2});
            else
                did_match = value == criteria{i_rule, 2};
            end
            if did_match
                these_columns_experimental_condition = columns_experimental_condition{i_rule, 2};
                break
            end
        end
        columns_experimental_condition = these_columns_experimental_condition;
    end
    if strcmp(value, 'epi')
        columns_experimental_condition = {'day_from_surgery', 'compound'};
    end
    conds_trial = METADATA(:, ['stimulus', columns_experimental_condition, 'keep']);
    conds_trial.day_from_surgery = str2double(conds_trial.day_from_surgery);
    % Keep only a subset of trials corresponding to allowed stimuli and compounds
    allowed_compounds = GC.analysis.(experiment_type).analysis_allowed_compounds;
    allowed_trials = ismember(conds_trial.compound, allowed_compounds) & conds_trial.keep;

    % Get the unique combinations of conditions
    stimuli_name_original = {};
    for i_stim = 1:length(stimuli_name), stimuli_name_original{i_stim} = make_variable_name(stimuli_name{i_stim}, true, true); end
    conds = unique(conds_trial(allowed_trials & ismember(conds_trial.stimulus, stimuli_name_original),:), 'rows', 'stable');
    conds.keep = [];
    conds.compound = [];
    [day_from_surgery, ~, days_index] = unique(conds.day_from_surgery);

    % Select only stimuli present in all sessions
    stimuli_across = unique(conds.stimulus, 'stable');
    flag_stim = [];
    for i_stim = 1 : length(stimuli_across)
        this_stim = stimuli_across(i_stim);
        flag_stim(i_stim) = sum(ismember(conds.stimulus, this_stim)) == length(day_from_surgery);
    end
    flag_stim= logical(flag_stim);
    stimuli_across = stimuli_across(flag_stim);

    % Check number of trials for all conditions (this must be the same to pass
    % to classifier
    % count trials per stimuli per condition
    count_trials = struct();
    count_this_trials = [];

    % create an array of n_stim x n_conds
    for i_stim = 1 : length(stimuli_across)
        fname = make_variable_name(stimuli_across{i_stim}, true);
        if strcmp(fname, 'temp_42'), fname = 'temp_43';end
        for i_cond = 1 : length(day_from_surgery)
            count_this_trials(i_stim, i_cond) = sum(allowed_trials & ...
                ismember(conds_trial{:, 'stimulus'}, stimuli_across(i_stim)) &...
                ismember(conds_trial{:, 'day_from_surgery'}, day_from_surgery(i_cond)));
        end
        count_trials.(fname)= min(count_this_trials(i_stim,:));
    end

    % % detect the min nr trials
    % min_trials_per_stim = [];
    % for i_stim = 1: length(stimuli_across)
    %     this_stim =
    %     min_trials_per_stim(i_stim,:) = min(count_trials(i_stim,:));
    % end




    % Keep only last 2 sessions before surgery and the last 2 after
    % sessions_before_surgery = find(day_from_surgery <= 0, 2, 'last');
    % sessions_after_surgery  = find(day_from_surgery >  0, 2, 'first');
    % sessions_to_analyze = [sessions_before_surgery; sessions_after_surgery];
    sessions_to_analyze = 1:length(day_from_surgery);
    n_conditions = length(sessions_to_analyze);

    temp_dir = [GC.temp_dir, '\', animal_ID, filesep];
    if ~exist(temp_dir, 'dir'), mkdir(temp_dir), end

%     output_filename_results          = os.path.join(data_dir, [animal_ID, '_crossday_results.json']);
%     output_filename_confusion_matrix = os.path.join(data_dir, [animal_ID, '_crossday_confusion_matrix.json']);

    %% Decoding analysis
    % Make a structure containing the configuration parameters to pass to the
    % classifier
    cfg = CFG_MVPA ;
    cfg.classifier = classifier;
    %% Check this
%     cfg.classifier =  'naive_bayes';
    cfg.metric     = {'accuracy' 'confusion'};
    cfg.cv         = 'kfold';
    cfg.k          = 5;
    cfg.hyperparameter=[];
    cfg.feedback = 1;
    % cfg.sample_dimension = 1;
    cfg.preprocess = {};
    cfg.preprocess_param = {};
    cfg.append = 0;
    %%
    X = [];
    L = X;
    concat_cm = [];
    info = struct();

    % Loop through conditions and analyze the data
    % It won't perform generalizability, since the analysis will be done
    % per session
    
    
    % Check subsets
    switch analysis_type_str
        case 'all'
%             keyboard % check if we are using the correct name str
            sessions_to_analyze = 1:length(day_from_surgery);
            n_conditions = length(sessions_to_analyze);
            this_cell_subset = repmat(this_cell_subset, n_conditions,1);
        case {'H_noci', 'H_noci_random'}
            
            n_subs = length(this_cell_subset_info);
            day_from_surgery = NaN(n_subs,1);
            for i_inf = 1:n_subs
                 day_from_surgery(i_inf,1) = this_cell_subset_info{i_inf,1}{1};
                
            end
            
            sessions_to_analyze = 1:length(day_from_surgery);
            n_conditions = length(sessions_to_analyze);
            
    end

    
    for ii = 1:n_conditions
        if skip_same_day_decoding, continue, end

        i_cond = sessions_to_analyze(ii);
        LOGGER.info(['Gathering condition ', num2str(i_cond)])

        % Get stimuli in this condition
        this_day = day_from_surgery(ii);
        stimuli_this_condition = natsort(conds{conds.day_from_surgery == this_day, 'stimulus'});
%         stimuli_this_condition = stimuli_across;
        n_stimuli =  length(stimuli_this_condition);

        % Initialize structure containing data and parameters for analysis
        PARAMETERS = struct();
        PARAMETERS.data = struct();
        evoked_activity_max_latency = 6;%GC.response_decoding_evoked_window;
        evoked_activity_max_latency_bins = ceil(evoked_activity_max_latency / time_bin_activity);
        PARAMETERS.data_shape = [];%struct();
        data_shape_noci = [];
        data_shape_ave = [];
        data_shape_inn = [];

        for i_stim = 1:n_stimuli
            % Get data
            trial_idx = allowed_trials & ...
                ismember(conds_trial{:, 'stimulus'}, stimuli_this_condition{i_stim}) & ...
                ismember(conds_trial{:, 'day_from_surgery'}, day_from_surgery(i_cond));
            
            % get the spikes from subset (if it's not 'all' then it takes the real subset, otherwise all)
            cells_to_take = this_cell_subset{ii};
            data = spikes(cells_to_take, trial_idx);
            n_ROIs = size(data,1);
            % Bin data to speed-up computations
            fname = make_variable_name(stimuli_this_condition{i_stim}, true);
            % Skip SP and SP_2
            if ismember(fname, {'SP', 'SP_2', 'temp_38', 'odor_plus'}), continue, end

            if isempty(STIMULI.(fname).timestamps)
                timestamps = [];
                stimulus_onset = [];
            else
                timestamps = round(STIMULI.(fname).timestamps * frame_rate);
                if startsWith(fname, 'temp_')
                    stimulus_onset = timestamps(2);
                    % Read all the frames after onset
%                     stimulus_onset = timestamps(1);
                    time_bin_activity =  0.712; % 50;
                    time_bin_activity_n_frames = round(time_bin_activity * frame_rate);

                else
                    time_bin_activity = 0.712; % CFG_MVPA.time_bin_activity; % GC.response_decoding_window_width;
                    time_bin_activity_n_frames = round(time_bin_activity * frame_rate);
                    stimulus_onset = timestamps(1);
                end
            end
            % Bin and slice data
            [data_binned, ~, bin_edges] = make_tensor(data, time_bin_activity_n_frames, stimulus_onset);
            n_bins = size(bin_edges, 1);
            if ~isempty(timestamps)
                % Get the bin where the stimulus onset falls in
                stimulus_onset_bin = find(bin_edges(:, 1) == stimulus_onset);
                % Convert timestamps to bin indices
                timestamps_bin = NaN(length(timestamps), 1);
                for it = 1:length(timestamps)
                    timestamps_bin(it) = find(timestamps(it) >= bin_edges(:,1) & timestamps(it) <= bin_edges(:,2));
                end

                if startsWith(fname, 'temp_')
                    % Get first and last bin to cut
                    first_bin = timestamps_bin(2);
                    last_bin = timestamps_bin(2) + CFG_MVPA.last_bin; % try and check
                else
                    % Get first and last bin to cut
                    first_bin = timestamps_bin(1);
                    last_bin = timestamps_bin(1) +  CFG_MVPA.last_bin; % timestamps_bin(1) + 1;
                end

                % Make sure that these indices are not out of bounds
                first_bin = max([first_bin, 1]);
                last_bin = min([last_bin, n_bins]);
                % Slice data
                data_binned = data_binned(:, first_bin:last_bin, :);
            end

            % Convert temp_42 to temp_43, and SP_2 to SP
            if strcmp(fname, 'temp_42'), fname = 'temp_43'; end
            if strcmp(fname, 'SP_2'), fname = 'SP'; end

            if ~isempty(CATEGORIES)
                categories = fieldnames(CATEGORIES);
                categories(ismember(categories, {'Properties','Row','Variables'})) = [];

                for i_cat = 1:width(CATEGORIES)
                    this_cat = categories{i_cat};
                    switch this_cat
                        case 'noci'
                            is_noci = any(ismember(CATEGORIES.(this_cat), fname));
                        case 'ave'
                            is_ave = any(ismember(CATEGORIES.(this_cat), fname));
                        case 'inn'
                            is_inn = any(ismember(CATEGORIES.(this_cat), fname));
                    end
                end

                % Check if it's noci
                if is_noci
                    this_cat_str = 'noci';
                    if ~isfield(PARAMETERS.data, this_cat_str)
                        PARAMETERS.data.(this_cat_str) = reshape(data_binned, n_ROIs, []).';
                    else
                        PARAMETERS.data.(this_cat_str) = [PARAMETERS.data.(this_cat_str); reshape(data_binned, n_ROIs, []).'];
                    end
                    % Add index of which stimulus belongs where
                    [~, n_bins, ~] = size(data_binned);
                    n_trials = count_trials.(fname);
                    data_shape_noci = [data_shape_noci ;[n_bins * n_trials]];


                end
                % Check if it's ave
                if is_ave
                    this_cat_str = 'ave';
                    if ~isfield(PARAMETERS.data, this_cat_str)
                        PARAMETERS.data.(this_cat_str) = reshape(data_binned, n_ROIs, []).';
                    else
                        PARAMETERS.data.(this_cat_str) = [PARAMETERS.data.(this_cat_str); reshape(data_binned, n_ROIs, []).'];
                    end
                    % Add index of which stimulus belongs where
                    [~, n_bins, ~] = size(data_binned);
                    n_trials = count_trials.(fname);
                    data_shape_ave = [data_shape_ave ;[n_bins * n_trials]];
                end

                % Check if it's ave
                if is_inn
                    this_cat_str = 'inn';
                    if ~isfield(PARAMETERS.data, this_cat_str)
                        PARAMETERS.data.(this_cat_str) = reshape(data_binned, n_ROIs, []).';
                    else
                        PARAMETERS.data.(this_cat_str) = [PARAMETERS.data.(this_cat_str); reshape(data_binned, n_ROIs, []).'];
                    end
                    % Add index of which stimulus belongs where
                    [~, n_bins, ~] = size(data_binned);
                    n_trials = count_trials.(fname);

                    data_shape_inn = [data_shape_inn ;[n_bins * n_trials]];
                end
            else

                if ~isfield(PARAMETERS.data, fname)
                    PARAMETERS.data.(fname) = reshape(data_binned, n_ROIs, []).';
                else
                    PARAMETERS.data.(fname) = [PARAMETERS.data.(fname); reshape(data_binned, n_ROIs, []).'];
                end
                % Add index of which stimulus belongs where
                [~, n_bins, n_trials] = size(data_binned);
%                 keyboard % fix data_shape to =[];
                PARAMETERS.data_shape = [PARAMETERS.data_shape; [n_bins, n_trials]];
            end
        end
    

        % Reshape data for categories
        if ~isempty(CATEGORIES)
            if sum(data_shape_noci)<=30; row = 1;else, row = 2;end
            PARAMETERS.data_shape.noci = [row, sum(data_shape_noci)/row] ;
            if sum(data_shape_ave)<=30; row = 1;else, row = 2;end
            PARAMETERS.data_shape.ave = [row, sum(data_shape_ave)/row] ;
            if sum(data_shape_inn)<=30; row = 1;else, row = 2;end
            PARAMETERS.data_shape.inn = [row, sum(data_shape_inn)/row] ;
%             PARAMETERS.data_shape(any(PARAMETERS.data_shape == 0,2),:) = [];
        end


        % Set the filename of the output, which will be read back into MATLAB
%         PARAMETERS.output_filenames = struct();
%         PARAMETERS.output_filenames.results    = os.path.join(data_dir, [[animal_ID, '_naive_bayes'], '_cond', num2str(i_cond), '_decoding_results.json']);
%         PARAMETERS.output_filenames.classifier = os.path.join(data_dir, [animal_ID, '_cond', num2str(i_cond), '_decoding_classifier.p']);


        % reshape Data
        X_data = [];
        l_data = X_data;
        these_stim = fieldnames(PARAMETERS.data);
        if ~skip_diff_day_decoding
            keyboard
            for i_stim  = 1:length(these_stim)
                fname = these_stim{i_stim};
                %         n_trials_this_stim = count_trials.(fname);
                these_numbers = PARAMETERS.data_shape.(fname);
                n_rows = these_numbers(1) * these_numbers(2);
                this_data = PARAMETERS.data.(fname);
                this_data = this_data(1:n_rows,:);
                X_data = [X_data; this_data];
                this_label = i_stim*ones(size(this_data, 1),1);
                l_data = [l_data;this_label];
            end
        else
            for i_stim  = 1:length(these_stim)
                fname = these_stim{i_stim};
%                 stim_idx = find(ismember(allowed_stimuli, fname));
                %         n_trials_this_stim = count_trials.(fname);
                these_numbers = PARAMETERS.data_shape(i_stim,:);
%                 n_rows = these_numbers(1) * these_numbers(2);
                this_data = PARAMETERS.data.(fname);
                X_data = [X_data; this_data];
                this_label = i_stim*ones(size(this_data, 1),1);
                l_data = [l_data;this_label];
            end
        end
        
       
        % Run Decoder for single sessions
        if skip_diff_day_decoding
            
            
            if ~do_random_labels
                MVPA_labels = l_data;
            else
                a = l_data;
                MVPA_labels = a(randperm(length(a)));
            end
            
            
            % Set prior
            classes = unique(MVPA_labels, 'stable');
            prior = [];
            for i_class = 1 : length(classes)
                prior (i_class) = sum(ismember(MVPA_labels, i_class))/length(MVPA_labels); %#ok<AGROW>
            end
            % Run Decoder
            cfg.hyperparameter.prior = prior;
            cfg.sample_dimension = 1;
            cfg.feature_dimension = 2;
            cfg.generalization_dimension = [];
            cfg.append = false;
            cfg.dimension_names = {'trials', 'cells'};
            % Run Classifier
            if convert_data_to_gaussian
                X_data = log10(X_data); % create a gaussian distribution of data
                X_data(X_data == -inf) = 0;
            end
            try
                [perf, ~, ~] = mv_classify(cfg,X_data,MVPA_labels);
            catch
                keyboard
            end
            
            PERF = perf{1,1};
            CM = perf{2,1};
               
            %% Allocate data to all stimuli when there is less stim than the allowed stimuli
            CM_all = NaN(length(allowed_stimuli));
            n_elements = numel(CM);
            initial_idx = find(ismember(allowed_stimuli, these_stim));
            for ie = 1:n_elements
               n_rows = size(CM,1);
               n_cols = size(CM,2);
               
               for ir  = 1:n_rows
                    this_idx_r = initial_idx(ir);
                  for ic = 1:n_cols 
                      this_idx_c = initial_idx(ic);
                   CM_all(this_idx_r,this_idx_c ) = CM(ir,ic);
                  end
               end
            end
            
            %%
%             keyboard
%             CLASSIFIER_GENERALIZABILITY = PERF;
%             confusion_matrices_table =CM;
%             classes = these_stim;
%             info = cell(n_conditions,1);
%             row_col_names = cell(1,n_conditions);
%             for i_cond = 1:n_conditions
%                 cond_name = ['cond_', num2str(i_cond)];
%                 row_col_names{1, i_cond} = cond_name;
%                 info{i_cond} = {animal_ID;cond_name};
%             end


        end
        concat_cm(:,:, i_cond) = CM_all; %#ok<AGROW>
        this_day = day_from_surgery(i_cond);
        this_info.stimuli = these_stim;
        this_info.day = this_day;
        info.(['cond', num2str(i_cond)]) = this_info;
        
    end
    days = day_from_surgery;
    CLASSIFIER_GENERALIZABILITY = [];  % To be fixed later
    
%     
%     X(:,:,i_cond) = X_data;
%     L(:,:,i_cond)= l_data;
%     keep_on = find(l_data>1); % if lables is only 1
%     if keep_on && ~skip_diff_day_decoding
% 
%         % Set hyperparameters
%         classes = unique(l_data, 'stable');
%         prior = [];
%         for i_class = 1 : length(classes)
%             prior (i_class) = sum(ismember(l_data, i_class))/length(l_data);
% 
%         end
%         % prior = [sum(ismember(l_data, 1)) / length(l_data),...
%         %     sum(ismember(l_data, 2)) / length(l_data),...
%         %     sum(ismember(l_data, 3)) / length(l_data)];
%         cfg.hyperparameter.prior = prior;
%         cfg.sample_dimension = 1;
%         cfg.feature_dimension = 2;
%         cfg.generalization_dimension = 3;
%         cfg.append = false;
%         cfg.dimension_names = {'trials', 'cells', 'sessions'};
%         % Run Classifier
%         try
%             [perf, result, testlabel] = mv_classify(cfg,X,l_data);
%         catch
%             keyboard
%         end
%         % Plot Results
%         % mv_plot_result(result);
% 
%         % Load Confusion matrices and performance
%         PERF = perf{1,1};
%         CM = perf{2,1};
% 
%         % reshape Confusion matrices
%         confusion_matrices = cell(n_conditions);
%         for i_train_sess = 1 : n_conditions
%             for i_test_sess = 1: n_conditions
%                 this_sess = CM(i_train_sess, :, :, i_test_sess);
%                 this_sess_reshaped = reshape(this_sess, size(this_sess, 2), size(this_sess,3));
%                 confusion_matrices{i_train_sess, i_test_sess} = this_sess_reshaped;
%             end
%         end
% 
%         % set condition names
%         info = cell(n_conditions,1);
%         row_col_names = cell(1,n_conditions);
%         for i_cond = 1:n_conditions
%             cond_name = ['cond_', num2str(i_cond)];
%             row_col_names{1, i_cond} = cond_name;
%             info{i_cond} = {animal_ID;cond_name};
%         end
%         confusion_matrices_table = cell2table(confusion_matrices, 'RowNames', row_col_names, 'VariableNames', row_col_names);
%         CLASSIFIER_GENERALIZABILITY = PERF;
%         % stats_cf_test =confusion_matrix_stats(confusion_matrices_table.cond_1{1});
%         classes = these_stim;
% 
%         % % Save the outout as Json file
%         % X = struct('performance', PERF, 'confusion_matrix', confusion_matrices_table, 'classes', {classes}, 'animal_ID', {animal_ID});
%         % S = jsonencode(X);
%         % json_filename =  '/Users/marioacuna/Desktop/json_trial.json';
%         % fid = fopen(json_filename, 'w');
%         % if fid == -1, error('Cannot create JSON file'); end
%         % fwrite(fid, S, 'char');
%         % fclose(fid);
%         %
% 
%         % %% test extraction of CM
%         % % To read JSON
%         % txt = jsondecode(fileread(json_filename));
%         % % get values of CM from same day decoding (diagonals)
%         % classes = txt.classes;
%         % % Cm accross time and concatenate
%         % n_conditions = size(txt.confusion_matrix,1);
%         % CM = [];
%         % for i_cond = 1:n_conditions
%         %     cond_str = ['cond_', num2str(i_cond)];
%         %     this_CM = array2table(txt.confusion_matrix.(cond_str));
%         %     CM(:,:,i_cond) = [];
%         %
%         %
%         % end
% 
%         % %% Test extraction from table
%         % % Loop first through animals to have something similar to what we have for
%         % % decoding confusion matrix for the logistic regression
%         %
%         % % Set the stimuli to analyse
%         % STIMULI = {'HPS', 'FPS', 'puff', 'touch', 'temp_48'};
%         % i_count = 0;
%         % n_animals = 1;
%         % n_sess = 5; % just to test ()
%         % same_day_CM = NaN(length(STIMULI), length(STIMULI),n_animals * n_sess);
%         % % for i_animal  = 1 : n_animals
%         % this_animal_cm= confusion_matrices_table;
%         % this_animal_classes = classes;
%         % n_conditions = size(this_animal_cm,1);
%         % for i_cond = 1:n_conditions
%         %     i_count = i_count + 1;
%         %     cond_str = ['cond_', num2str(i_cond)];
%         %     this_CM = this_animal_cm.(cond_str){i_cond};
%         %     stim_ids = find(ismember(STIMULI, this_animal_classes));
%         %     same_day_CM(stim_ids,stim_ids, i_count) = this_CM;
%         % end
%         %
%         %
%         % % Plot test
%         % data_to_draw = nanmedian(same_day_CM,3);
%         %
%         % figure, heatmap(data_to_draw)
%     else
%         disp('Found only one class')
%         CLASSIFIER_GENERALIZABILITY = [];
%         confusion_matrices_table=[];
%         classes = [];
%         info = [];
%     end
end
    
    
    
    
    
    
