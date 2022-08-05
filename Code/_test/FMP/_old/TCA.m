%% Tensor Component Analysis
% Use TCA to factorize the stimulus response and find the peak of maximum
% activity following stimulus oneset. This will also provide scores (i.e., factors)
% corresponding to the contribution to such response, which will be used as
% features for SVM classification.

if false

% Initialize output variable
TCA_factors = struct();
for stim = 1:n_stimuli
    stimulus = stimuli_name{stim};
    TCA_factors.(stimulus) = struct();
    
    % Perform TCA on stimulus-evoked activity
    rows = ismember(conds_trial.stimulus, stimulus) & allowed_trials;
    condition_trials = str2double(conds_trial.day_from_surgery(rows));
    parameters = struct('tensor_rank',tensor_rank, 'time_bin_activity_n_frames',time_bin_activity_n_frames);
    parameters.stimulus_structure       = STIMULI.(stimulus);
    parameters.condition_trials         = condition_trials;
    parameters.stimulus_onset           = round(STIMULI.(stimulus).timestamps(1) * frame_rate);
    parameters.response_window_n_frames = response_window_n_frames;
    parameters.force_keep_top_n_factors = 0;
    parameters.is_spontaneous           = false;
    LOGGER.info(['Performing TCA of responses to ', stimulus])
    TCA_factors.(stimulus).evoked = perform_TCA(spikes(:, rows), parameters); 
    
    % Perform TCA of spontaneous activity
    rows = ismember(conds_trial.stimulus, 'SP') & ismember(conds_trial.compound, allowed_compounds) & conds_trial.keep;
    condition_trials = str2double(conds_trial.day_from_surgery(rows));
    parameters = struct('tensor_rank',tensor_rank, 'time_bin_activity_n_frames',time_bin_activity_n_frames);
    parameters.stimulus_structure       = STIMULI.SP;
    parameters.condition_trials         = condition_trials;
    parameters.stimulus_onset           = [];
    parameters.response_window_n_frames = [];
    parameters.force_keep_top_n_factors = tensor_rank;
    parameters.is_spontaneous           = true;
    parameters.n_bins                   = size(TCA_factors.(stimulus).evoked.TCs.time, 1);
    parameters.n_trials                 = size(TCA_factors.(stimulus).evoked.TCs.trial, 1);
    LOGGER.info('Performing TCA of spontaneous activity')
    TCA_factors.(stimulus).spontaneous = perform_TCA(spikes(:, rows), parameters);
end

save('C:\Users\deluna\Downloads\test_TCA.mat', 'TCA_factors')
end
TCA_factors = load_variable('C:\Users\deluna\Downloads\test_TCA.mat', 'TCA_factors');
end


%% Tensor Component Analysis

% function output = perform_TCA(spikes_data, parameters)
%     % Unpack inputs
%     stimulus_structure         = parameters.stimulus_structure;
%     time_bin_activity_n_frames = parameters.time_bin_activity_n_frames;
%     condition_trials           = parameters.condition_trials;
%     stimulus_onset             = parameters.stimulus_onset;
%     response_window_n_frames   = parameters.response_window_n_frames;
%     tensor_rank                = parameters.tensor_rank;
%     force_keep_top_n_factors   = parameters.force_keep_top_n_factors;
%     is_spontaneous             = parameters.is_spontaneous;
%     
%     if is_spontaneous
%         n_bins = parameters.n_bins;
%         n_trials = parameters.n_trials;
%         [X, ~] = make_tensor_spontaneous_data(spikes_data, time_bin_activity_n_frames, condition_trials, n_bins, n_trials, false, false);
%     else
%         % Transform data to tensor format
%         [X, ~, bin_edges] = make_tensor(spikes_data, time_bin_activity_n_frames, stimulus_onset, condition_trials, false);
%         [~, n_bins, n_trials] = size(X);
%     end
%     X = fmt(X);
% 
%     % Get response window
%     if is_spontaneous
%         response_window_bins = 1:n_bins;  % use all bins
%         n_response_window_bins = length(response_window_bins);
%         baseline_bins = [];
%         n_baseline_bins = length(baseline_bins);
%         
%     else
%         response_window_frames = stimulus_structure.response_window(1) + response_window_n_frames;
%         response_window_bins = find(response_window_frames(1)>=bin_edges(:,1), 1,'last'):find(response_window_frames(2)>=bin_edges(:,2), 1,'last');
%         n_response_window_bins = length(response_window_bins);
%         % Get bins corresponding to baseline
%         last_baseline_bin = round(stimulus_structure.baseline_window(end) / time_bin_activity_n_frames);
%         if ismember(last_baseline_bin, response_window_bins)
%             last_baseline_bin = last_baseline_bin - 1;
%         end
%         baseline_bins = last_baseline_bin-n_response_window_bins+1:last_baseline_bin;
%         baseline_bins(baseline_bins < 1) = [];
%         n_baseline_bins = length(baseline_bins);
%     end
% 
%     % Build model specifications
%     model = struct();
%     U0 = cpd_gevd(X, tensor_rank);
%     model.variables.neuron = U0{1};
%     model.variables.time   = U0{2};
%     model.variables.trial  = U0{3};
%     model.factors.neuron = {'neuron', @struct_nonneg};
%     model.factors.time   = {'time', @struct_nonneg};
%     model.factors.trial  = {'trial', @struct_nonneg};
%     model.factorizations.tensor.data = X;
%     model.factorizations.tensor.cpd = {'neuron', 'time', 'trial'};
%     options = struct();
%     options.Algorithm = @nls_lm;
%     options.AlgorithmOptions.LineSearch = @cpd_els;  % Exact line search
%     options.MaxIter = 1000;
%     options.TolFun  = 1e-12;  % Function tolerance stop criterion
%     options.TolX    = 1e-12;  % Step size tolerance stop criterion
%     % Perform CPD
%     TCA_model = sdf_nls(model, options);
%     
%     % Compute area under the curve of each factor within the response window
%     loadings = vecnorm(TCA_model.factors.time(response_window_bins, :), 2, 1);
%     [~, sorting_order_factors] = sort(loadings(:), 'descend');
% 
%     if force_keep_top_n_factors == 0
%         % Project the data matrix (X) on combinations of factors that contain the
%         % top-n factors. For each reconstruction (X_hat), compute the amount of
%         % activity that is reconstructed in the response and in the baseline
%         % windows. Finally, keep only the factors that make the contribution to the
%         % response window increase.
%         X_response = X(:, response_window_bins, :);
%         X_baseline = X(:, baseline_bins, :);
%         reconstruction_performance = zeros(tensor_rank, 2);
%         
%         for irank = 1:tensor_rank
%             % Project data onto top-n factors
%             factors_to_use = sorting_order_factors(1:irank);
%             X_hat = zeros(size(X));
%             for k = 1:n_trials
%                 X_hat(:, :, k) = TCA_model.factors.neuron(:, factors_to_use) * diag(TCA_model.factors.trial(k, factors_to_use)) * TCA_model.factors.time(:, factors_to_use).';
%             end
%             
%             % Split evoked from baseline activity
%             X_hat_response_window = X_hat(:, response_window_bins, :);
%             X_hat_baseline_window = X_hat(:, baseline_bins, :);
% 
%             % Compute residual of reconstruction, and ratio of X accounted for by X_hat
%             x = (X_response - X_hat_response_window);
%             added_to_response = 1 - sum(x(:)) / sum(X_response(:));
%             x = (X_baseline - X_hat_baseline_window);
%             added_to_baseline = 1 - sum(x(:)) / sum(X_baseline(:));
%             % Store data
%             reconstruction_performance(irank, :) = [added_to_response, added_to_baseline];
%         end
%         % The relative amount of each window has to be normalized by the number of
%         % bins, and we keep only the difference of the two.
%         distance = reconstruction_performance(:,1)./n_response_window_bins - reconstruction_performance(:,2)./n_baseline_bins;
%         % Find last factor that keeps increasing steadily the contribution to the
%         % evoked activity.
%         best_rank = find(diff(distance)>0 == 0, 1, 'first');
%         if isempty(best_rank)
%             best_rank = tensor_rank;
%         end
%         factors_to_keep = sorting_order_factors(1:best_rank);
%     
%     else  % Use user input to determine output rank
%         factors_to_keep = sorting_order_factors(1:force_keep_top_n_factors);
%     end 
%     
%     % Make output structure
%     output = struct();
%     output.TCs = TCA_model.factors;
%     output.factors_to_keep = factors_to_keep;
% end
