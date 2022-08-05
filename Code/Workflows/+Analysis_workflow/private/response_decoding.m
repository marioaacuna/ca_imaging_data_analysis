function [CLASSIFIER_GENERALIZABILITY, CONFUSION_MATRICES] = response_decoding(animal_ID, spikes, METADATA, STIMULI, data_dir, experiment_type)

skip_same_day_decoding = 0;
skip_diff_day_decoding = 0;

% Get LOGGER and general_configs
global LOGGER GC

% Get options from general_configs
n_permutations_significance = GC.decoding_n_permutations_significance;
baseline_window = GC.response_decoding_baseline_window;

% Get the number of ROIs
n_ROIs = size(spikes, 1);
% Get frame rate
frame_rate = max(METADATA.frame_rate);
% baseline_window_n_bins = round(baseline_window / time_bin_activity);

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
% Keep only last 2 sessions before surgery and the last 2 after
% sessions_before_surgery = find(day_from_surgery <= 0, 2, 'last');
% sessions_after_surgery  = find(day_from_surgery >  0, 2, 'first');
% sessions_to_analyze = [sessions_before_surgery; sessions_after_surgery];
sessions_to_analyze = 1:length(day_from_surgery);
n_conditions = length(sessions_to_analyze);

temp_dir = [GC.temp_dir, '\', animal_ID, filesep]; 
if ~exist(temp_dir, 'dir'), mkdir(temp_dir), end

output_filename_results          = os.path.join(data_dir, [animal_ID, '_crossday_results.json']);
output_filename_confusion_matrix = os.path.join(data_dir, [animal_ID, '_crossday_confusion_matrix.json']);

%% Decoding analysis
for ii = 1:n_conditions
    if skip_same_day_decoding, continue, end
    
    i_cond = sessions_to_analyze(ii);
    LOGGER.info(['Analyzing condition ', num2str(i_cond)])
    
    % Get stimuli in this condition
    stimuli_this_condition = natsort(conds{days_index == i_cond, 'stimulus'});
    n_stimuli =  length(stimuli_this_condition);

    % Initialize structure containing data and parameters for analysis
    PARAMETERS = struct();
    PARAMETERS.decoder_implementation = 'logistic regression';%'logistic regression'; % 'random forest'
    PARAMETERS.data = struct();
%     evoked_activity_max_latency = 6;%GC.response_decoding_evoked_window;
%     evoked_activity_max_latency_bins = ceil(evoked_activity_max_latency / time_bin_activity);
    PARAMETERS.data_scaling = 'minmax';
    PARAMETERS.CV_scheme = 'k-fold';
    PARAMETERS.fraction_training = 1 - 1 /  GC.decoding_CV_k_folds;
    PARAMETERS.max_n_splits_training = GC.decoding_CV_k_folds;
    PARAMETERS.fraction_training_for_tuning = 0;%0.3;
    PARAMETERS.max_n_splits_tuning = GC.decoding_CV_k_folds;
    PARAMETERS.fraction_training_for_calibration = 0.8;
    PARAMETERS.max_n_splits_calibration = GC.decoding_CV_k_folds;
    PARAMETERS.v_fraction = 0.1;
    PARAMETERS.shuffle_split_repeats = 3;
    PARAMETERS.data_shape = [];

    for i_stim = 1:n_stimuli
        % Get data
        trial_idx = allowed_trials & ...
            ismember(conds_trial{:, 'stimulus'}, stimuli_this_condition{i_stim}) & ...
            ismember(conds_trial{:, 'day_from_surgery'}, day_from_surgery(i_cond));
        data = spikes(:, trial_idx);
        
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
%                 stimulus_onset = timestamps(2);
                % Read all the frames after onset
                stimulus_onset = timestamps(1);
                time_bin_activity = 50;
                time_bin_activity_n_frames = round(time_bin_activity * frame_rate);
 
            else
                if strcmp (experiment_type, 'ACC_SNI_anxiety')
                    time_bin_activity = 0.4; % 3
                else
                    time_bin_activity = GC.response_decoding_window_width;
                end
                time_bin_activity_n_frames = round(time_bin_activity * frame_rate);
                stimulus_onset = timestamps(1);
            end
        end
        % Bin and slice data
        [data_binned, ~, bin_edges] = make_tensor(data, time_bin_activity_n_frames, stimulus_onset, 'mean');
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
                last_bin = timestamps_bin(2) + 1;
            else
                % Get first and last bin to cut
                first_bin = timestamps_bin(1);
                if strcmp (experiment_type, 'ACC_SNI_anxiety')
                    last_bin = timestamps_bin(1) + 3;
                else
                    last_bin = timestamps_bin(1) + 1;
                end
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
        
        if ~isfield(PARAMETERS.data, fname)
            PARAMETERS.data.(fname) = reshape(data_binned, n_ROIs, []).';
        else
            PARAMETERS.data.(fname) = [PARAMETERS.data.(fname); reshape(data_binned, n_ROIs, []).'];
        end
        % Add index of which stimulus belongs where
        [~, n_bins, n_trials] = size(data_binned);
        PARAMETERS.data_shape = [PARAMETERS.data_shape; [n_bins, n_trials]];
    end
    %% load previous resuts to pass to the MVPA
%     MPVA_datafile = ['U:\Analysis_Matlab_tests\MVPA_DATA_',animal_ID ,'3s_1b.mat'];
%     if exist(MPVA_datafile, 'file')
%         MVPA_DATA = load_variable(MPVA_datafile,'MVPA_DATA');
%         fnames= fieldnames(PARAMETERS.data);
%         for iif = 1: length(fnames)
%             this_name = fnames{iif};
%             this_data = PARAMETERS.data.(this_name);
%             % Concatenate data
%             data_mvpa = [MVPA_DATA.data.(this_name); this_data];
%             MVPA_DATA.data.(this_name) = data_mvpa;
%             % concatenate days from surgery
%             days_previous = MVPA_DATA.day_from_surgery.(this_name);
%             days_now = repmat(day_from_surgery(i_cond),size(this_data,1),1);
%             MVPA_DATA.day_from_surgery.(this_name) = [days_previous; days_now];
%         end
%         this_pre_ds = [PARAMETERS.data_shape, zeros(size(PARAMETERS.data_shape,1),1) + i_cond];
%         data_shape_mvpa = [MVPA_DATA.data_shape; this_pre_ds];
% %         MVPA_DATA.data = data_mvpa;
%         MVPA_DATA.data_shape = data_shape_mvpa;
%     else
%         MVPA_DATA = struct();
%         MVPA_DATA.data = PARAMETERS.data;
%         this_data_shape = [PARAMETERS.data_shape, repmat(day_from_surgery(i_cond), size(PARAMETERS.data_shape,1),1)];
%         MVPA_DATA.data_shape = this_data_shape;
%         % create a variable for concatenating the days from surgery
%         fnanes = fieldnames(MVPA_DATA.data);
%         for iif = 1: length(fnanes)
%             this_fname = fnanes{iif};
%             this_data =  MVPA_DATA.data.(this_fname);
%             MVPA_DATA.day_from_surgery.(this_fname) = repmat(day_from_surgery(i_cond), size(this_data,1),1);
%         end
%     end
%     save(MPVA_datafile, 'MVPA_DATA')
%     disp('Saved MVPA')
    %%
    % Set which significance tests to perform
    PARAMETERS.significance_test_n_permutations = n_permutations_significance;
    % Set the filename of the output, which will be read back into MATLAB
    PARAMETERS.output_filenames = struct();
    PARAMETERS.output_filenames.results    = os.path.join(data_dir, [animal_ID, '_cond', num2str(i_cond), '_decoding_results.json']);
    PARAMETERS.output_filenames.classifier = os.path.join(data_dir, [animal_ID, '_cond', num2str(i_cond), '_decoding_classifier.p']);
   
    % Store data to disk
    temp_data_filename = [temp_dir, animal_ID, '_cond', num2str(i_cond), '_decoding_data.mat'];
    save(temp_data_filename, 'PARAMETERS', '-v6')
    
    % Run decoding analysis in python
    try
        run_in_python('decode_activations_same_day', ['input_filename=', temp_data_filename], 'verbose=True')
    catch %#ok<CTCH>
        keyboard
        clc
    end
    % Delete accessory file
    delete(temp_data_filename)
end

% Run across-day predictions to compare generalizability of learned features
if ~skip_diff_day_decoding
    run_in_python('decode_activations_different_day', ['folder=', data_dir], ['output_filename_results=', output_filename_results], ['output_filename_confusion_matrix=', output_filename_confusion_matrix], 'verbose=True')
end

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

% Load confusion matrices
CONFUSION_MATRICES = struct();
for ii = 1:n_conditions
    i_cond = sessions_to_analyze(ii);
    filename = os.path.join(data_dir, [animal_ID, '_cond', num2str(i_cond), '_decoding_results.json']);
    txt = jsondecode(fileread(filename));
    confusion_matrix = txt.performance.confusion_matrix;
    confusion_matrix = cellfun(@(x) x', confusion_matrix, 'UniformOutput',false);
    confusion_matrix = vertcat(confusion_matrix{:});
    % Replace temp_42 with temp_43
    confusion_matrix(1, ismember(confusion_matrix, 'temp_42')) = {'temp_43'};
    confusion_matrix(ismember(confusion_matrix, 'temp_42'), :) = {'temp_43'};
    % Convert to table
    confusion_matrix = array2table(str2double(confusion_matrix(2:end, 2:end)), 'VariableNames',confusion_matrix(1, 2:end), 'RowNames',confusion_matrix(2:end, 1));
    % Store values
    CONFUSION_MATRICES.(['cond_', num2str(i_cond)]) = confusion_matrix;
end

if n_conditions > 1
    % Get the confusion matrices for across-day performance
    confusion_matrix = jsondecode(fileread(output_filename_confusion_matrix));
    if ~isempty(confusion_matrix)
        confusion_matrix = cellfun(@(x) x', confusion_matrix, 'UniformOutput',false);
        confusion_matrix = vertcat(confusion_matrix{:});
        confusion_matrix(:, 1:2) = num2cell(str2double(confusion_matrix(:, 1:2)));
        for i_row = 1:size(confusion_matrix, 1)
            values = confusion_matrix{i_row, 3};
            values = cellfun(@(x) x', values, 'UniformOutput',false);
            values = vertcat(values{:});
            values = array2table(str2double(values(2:end, 2:end)), 'VariableNames',values(1, 2:end), 'RowNames',values(2:end, 1));
            confusion_matrix{i_row, 3} = values;
        end
    end
    CONFUSION_MATRICES.crossday = confusion_matrix;
else
    CONFUSION_MATRICES.crossday = [];
end

if exist(temp_dir, 'dir'), try rmdir(temp_dir), end, end

%% MLint exceptions
%#ok<*PFBNS,*UNRCH,*PFOUS,*AGROW,*STRNU,*NUSED,*NASGU,*NODEF>
