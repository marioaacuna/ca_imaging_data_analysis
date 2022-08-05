function STATS = response_modulation_ROC_stats(spikes, METADATA, STIMULI)

% Get LOGGER and general_configs
global LOGGER GC

% Get options from general_configs
n_surrogates = GC.response_modulation_n_surrogates;
surrogate_length = GC.response_modulation_surrogate_length;

% Get the number of ROIs
n_ROIs = size(spikes, 1);
% Get frame rate
frame_rate = max(METADATA.frame_rate);

% For each stimulus, get baseline and response activity
stimuli_name = fieldnames(STIMULI);
% Get the name of each condition
columns_experimental_condition = GC.columns_experimental_condition{ismember(GC.columns_experimental_condition(:,1),GC.experiment_name),2};
conds_trial = METADATA(:, ['stimulus', columns_experimental_condition, 'keep']);
% Keep only a subset of trials corresponding to the allowed stimuli
allowed_trials = ismember(conds_trial.stimulus, stimuli_name) & conds_trial.keep;

% Get the unique combinations of conditions
conds = unique(conds_trial(allowed_trials,:), 'rows', 'stable');
conds.keep = [];
n_conds = height(conds);

%% ANALYZE ROIs
% Allocate table
table_columns = ['ROI_id','stimulus',columns_experimental_condition,'AUROC','sensitivity','specificity','activity','p'];
RESULTS = repmat({''}, 0, length(table_columns));

% Loop through ROIS
for iroi = 1:n_ROIs
    % Log current ROI
    LOGGER.info(['ROI ', num2str(iroi), '/', num2str(n_ROIs)])
    % Get all data for this ROI
    all_data = spikes(iroi,:);
    
    % Loop through all conditions and get activity in baseline and response windows
    for icond = 1:n_conds
        % Get name of stimulus and condition
        stimulus_name = conds.stimulus{icond};
        cond_name = strjoin(table2cell(conds(icond,:)), ';');
        % Skip 'spontaneous'
        if strcmp(stimulus_name,'SP'), continue, end

        % Get trials and data for this condition
        trials_idx = allowed_trials;
        for icol = 1:size(conds,2)
            trials_idx = trials_idx & ismember(conds_trial(:,icol), conds(icond,icol));
        end
        trials_idx = trials_idx & conds_trial.keep;
        cond_data = all_data(trials_idx);
        if isempty(cond_data), continue, end  % No allowed trials for this condition
        % Stack trials on top of each other
        cond_data = vertcat(cond_data{:});
        n_trials = size(cond_data, 1);

        % Log beginning of algorithm
        LOGGER.trace(['Computing ', num2str(n_surrogates), ', ', num2str(round(surrogate_length/frame_rate)) 's-long surrogates for ', cond_name])
        % Find spontaneous activity for that day
        switch GC.experiment_name
            case 'ACC_CCI_anesth'
                training_data = all_data(ismember(table2cell(conds_trial(:,'stimulus')),'SP') & ismember(conds_trial(:,'day_from_surgery'),conds(icond,'day_from_surgery')));
            case 'Anaesth_check'
                training_data = all_data(ismember(table2cell(conds_trial(:,'stimulus')),'SP') & ismember(conds_trial(:,'anesthesia_level'),conds(icond,'anesthesia_level')));
            case 'ACC_pain_LP-211'
                training_data = all_data(ismember(table2cell(conds_trial(:,'stimulus')),'SP') & ismember(conds_trial(:,'compound_phase'),conds(icond,'compound_phase')));
        end
        training_data = horzcat(training_data{:});
        % Generate surrogate spike trains that preserve cell's firing statistics
        surrogate_spike_trains = Analysis_workflow.generate_surrogate_spike_trains(training_data, 'n_surrogates',n_surrogates, 'surrogate_length',surrogate_length, 'plot_distributions',false);
        % If the algorithm returned nothing, it means that it was impossible to
        % determine firing properties of the cell. Skip condition
        if isempty(surrogate_spike_trains) || all(surrogate_spike_trains(:)==0)
            LOGGER.info('No spontaneous activity to generate surrogate trials. Skip condition')
            continue
        end

        % Get windows of interest
        baseline_window = STIMULI.(stimulus_name).baseline_window;
        response_window = STIMULI.(stimulus_name).response_window;
        % Compute length of the segment of interest
        all_windows_size = length(baseline_window) + length(response_window);
        % Get activity in baseline and response windows
        baseline_activity = sum(cond_data(:,baseline_window), 2) ./ length(baseline_window);
        evoked_activity   = sum(cond_data(:,response_window), 2) ./ length(response_window);
        % Perform ROC analysis on original data
        ROC_stats = Analysis_workflow.ROC_analysis([baseline_activity, evoked_activity]);

        % The baseline window now starts from 1. Realign windows
        offset = baseline_window(1) - 1;
        baseline_window = baseline_window - offset;
        response_window = response_window - offset;

        % Pick a random segment of each surrogate trial
        latest_starting_frame = surrogate_length - all_windows_size;
        % Get a random start in each trial
        start_idx = randi(latest_starting_frame, [n_surrogates,1]);
        % Add trial length
        col_idx = bsxfun(@plus, start_idx, 0:all_windows_size-1);
        row_idx = repmat(1:n_surrogates, all_windows_size,1)';
        % Get the corresponding subset of surrogate data
        surrogate_redux = surrogate_spike_trains(sub2ind([n_surrogates,surrogate_length], row_idx, col_idx));

        % Get data in corresponding "baseline" and "response" windows and
        % normalize by window duration
        surrogate_baseline_activity = sum(surrogate_redux(:,baseline_window), 2) ./ length(baseline_window);
        surrogate_evoked_activity   = sum(surrogate_redux(:,response_window), 2) ./ length(response_window);
        % Concatenate data so it is easy to manipulate
        surrogate_activity = [surrogate_baseline_activity(:), surrogate_evoked_activity(:)];
        % Split the surrogate data in small datasets that have the same number
        % of trials as in the original data
        all_trials_idx = Shuffle(1:n_surrogates);
        n_batches = floor(n_surrogates / n_trials);
        n_surrogates_to_use = n_trials*n_batches;
        all_trials_idx = reshape(all_trials_idx(1:n_surrogates_to_use), n_batches, []);
        % Allocate output variable
        AUROC_surrogate = NaN(n_batches,1);
        parfor ibatch = 1:n_batches
            % Perform ROC analysis on surrogate data
            ROC_stats_surrogate = Analysis_workflow.ROC_analysis(surrogate_activity(all_trials_idx(ibatch,:), :));
            % Keep only AUROC
            AUROC_surrogate(ibatch) = ROC_stats_surrogate.AUROC;
        end

        % Compute two-tailed p-value
        p_left  = (1 + sum(AUROC_surrogate <= ROC_stats.AUROC)) / (n_batches+1);
        p_right = (1 + sum(AUROC_surrogate >= ROC_stats.AUROC)) / (n_batches+1);
        p = 2 * min([p_left, p_right]);
        % Fix any problem in the calculation of p
        if isnan(p) || p>1
            p = 1;
        end
        LOGGER.trace(['AUROC: ', value2str(ROC_stats.AUROC,'%.2f'), ' (p=', num2str(p,'%.5f'), ')'])

        % Concatenate results
        data2add = [{iroi},table2cell(conds(icond,:)),{ROC_stats.AUROC},{ROC_stats.TPR},{ROC_stats.SPC},{[baseline_activity, evoked_activity]},{p}];
        RESULTS(end+1, :) = data2add;
    end
end

% Convert all results to a table
STATS = cell2table(RESULTS, 'VariableName',table_columns);

%% MLint exceptions
%#ok<*PFBNS,*UNRCH,*PFOUS,*AGROW>
