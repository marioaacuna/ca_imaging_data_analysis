function spontaneous_activity_stats(all_animal_names, stats_filename)

% Get LOGGER and general_configs
global LOGGER GC

% Get info to compute distributions
event_attributes_name = GC.spontaneous_activity_attributes;
n_attributes = length(event_attributes_name);
% Make flags to know whhich attributes to analyze
analyze_interval = ismember('interval', event_attributes_name);
analyze_amplitude = ismember('amplitude', event_attributes_name);
% Get parameters for analysis
PARAMETERS = struct();
if analyze_amplitude
    PARAMETERS.amplitude.bins = GC.spontaneous_amplitude_bins;
end
if analyze_interval
    PARAMETERS.interval.bins = GC.spontaneous_interval_bins;
end
LOGGER.info(['Computing distributions of the following attributes: ', strjoin(event_attributes_name, ', ')])

% Load previous runs of this analysis
if exist(stats_filename, 'file')
    f = load(stats_filename);
    STATS_spontaneous_activity = f.STATS_spontaneous_activity;
    DISTRIBUTIONS = f.DISTRIBUTIONS;
    EVENT_STATS_spontaneous_activity = f.STATS_spontaneous_activity;
    clear f

else
    % Initialize output variables
    STATS_spontaneous_activity = struct();
    DISTRIBUTIONS = struct();
    EVENT_STATS_spontaneous_activity = struct();
end

% Get names of allowed compounts
allowed_compounds = GC.analysis_allowed_compounds.(GC.experiment_name);

% Loop through datasets
for iid = 1:length(all_animal_names)
    animal_ID = all_animal_names{iid};
    
    % Skip if spikes file does not exist
    spikes_filename = get_filename_of('spikes', animal_ID);
    if ~exist(spikes_filename,'file')
        LOGGER.info([animal_ID, ' not fully processed. Skipping'])
        continue
    end
    
    LOGGER.info(['Processing ' animal_ID])
    % Load files
    spikes = load_variable(spikes_filename, 'spikes');
    n_ROIs = size(spikes,1);

    % Read metadata
    METADATA = SQL_database.read_table_where('trials', {'stimulus','frame_rate','experimental_condition','date','recording_time'}, animal_ID,'animal_ID');
    % Read type of surgery for this animal
    experimental_group = SQL_database.read_table_where('experiments', 'experimental_group', animal_ID,'animal_ID', 'return_as_table',false);
    
    % Get trials marked as SP (that is, spontaneous)
    SP_idx = ismember(METADATA.stimulus, 'SP');
    % Keep only trials with allowed compounds
    SP_idx = SP_idx & ismember(METADATA.compound, allowed_compounds);
    trials_idx = find(SP_idx);
    % Remove unused trials
    METADATA = METADATA(trials_idx, :);
    spikes = spikes(:, trials_idx);
    
    % Split data by timepoint
    group_trials_by = GC.group_trials_by_timepoint.(GC.experiment_name);
    [timepoints, ~, timepoint_idx] = unique(METADATA{:, group_trials_by}, 'stable');
    n_timepoints = length(timepoints);
    switch GC.experiment_name
        case {'ACC_CCI_anesth','CLA_pain'}
            timepoints = str2double(timepoints);
        case 'ACC_pain_LP-211'
            active_timepoints = find(~cellfun(@isempty, regexp(timepoints, 'active_')));
            active_timepoints_str = timepoints(active_timepoints);
            timepoints_old = timepoints;
            timepoints = zeros(n_timepoints, 1);
            timepoints(active_timepoints) = cellfun(@(x) str2double(x(length('active_')+1:end)), active_timepoints_str);
            timepoints(ismember(timepoints_old, 'baseline')) = -1;
        otherwise
            error('!!! Unknown experiment')
    end

    % Allocate empty output variable
    SPONTANEOUS_ACTIVITY = table();
    SPONTANEOUS_ACTIVITY_basic_stats = table();
    for iroi = 1:n_ROIs
        % Log progress
        LOGGER.trace(['Analyzing ROI ', num2str(iroi), ' (out of ', num2str(n_ROIs), ')'], 'overwrite_last_message',iroi>1)
        % Allocate variables for this ROI
        mean_activity = zeros(1, n_timepoints);
        basic_stats_mean = NaN(n_timepoints, 1);
        basic_stats_rate = NaN(n_timepoints, 1);
        basic_stats_median = NaN(n_timepoints, 5);
        event_attributes = struct();
        for iab = 1:n_attributes
            event_attributes.(event_attributes_name{iab}) = cell(n_timepoints,1);
        end
        
        % Loop through days
        for it = 1:n_timepoints
            % Get and concatenate trials
            data = spikes(iroi, timepoint_idx==it);
            data = horzcat(data{:});
            % Compute sum of all data on this day
            mean_activity(it) = mean(data);
            
            % Find location of gaps between events
            gap_location = find(data == 0);
            gap_table = idx2range(gap_location);
            if analyze_interval
                event_attributes.interval{it} = gap_table(:,3);
            end
            
            % Get amplitude of events
            if analyze_amplitude
                activity = data(data > 0);
                if ~isempty(activity)
                    event_attributes.amplitude{it} = activity;
                else
                    event_attributes.amplitude{it} = 0;
                end
            end
            
            % Compute basic statistics of firing: mean amplitude event and
            % frequency
            activity = data(data > 0);
            data_duration = length(data) / unique(METADATA{:, 'frame_rate'});
            basic_stats_mean(it) = mean(activity);
            basic_stats_rate(it) = length(activity) / data_duration;
            basic_stats_median(it, :) = quantile(activity, [0, .25, .5, .75, 1]);
        end
        % When the cell was not active, replace NaN with 0
        basic_stats_mean(isnan(basic_stats_mean)) = 0;
        basic_stats_median(isnan(basic_stats_median)) = 0;
        
        % Compute probability distributions functions (PDF) across all days
        distributions = struct();
        for iab = 1:n_attributes
            attr_name = event_attributes_name{iab};
            distributions.(attr_name) = zeros(length(PARAMETERS.(attr_name).bins), n_timepoints);
        end
        for it = 1:n_timepoints
            for iab = 1:n_attributes
                % Get parameters
                attr_name = event_attributes_name{iab};
                bins = PARAMETERS.(attr_name).bins;
                values = event_attributes.(attr_name){it};
                % Compute PDF
                [~, ~, PDF]  = make_CDF(values, bins, true);
                distributions.(attr_name)(:,it) = cumsum(PDF(:,1));
            end
        end
        
        % Perform statistics between 'before' and 'after' surgery
        P_table = struct();
        for iab = 1:n_attributes
            P_table.(event_attributes_name{iab}) = [];
        end
        if n_timepoints > 1
            pairs = nchoosek(1:n_timepoints, 2);
            for iab = 1:n_attributes
                attr_name = event_attributes_name{iab};
                % Perform a Kolmogorov-Smirnov tests on all pairs of sessions
                P = ones(size(pairs,1), 4);
                for ipair = 1:size(pairs,1)
                    [~, p] = kstest2(event_attributes.(attr_name){pairs(ipair, 1)}, event_attributes.(attr_name){pairs(ipair, 2)});
                    P(ipair, 1:3) = [pairs(ipair,:), p];
                end
                % Correct p-values for multiple comparisons
                [~, P(:,4)] = adjust_Pvalues(P(:,3));

                % Count number of significant tests in baseline, pre-post, and post comparisons
                timepoints_before = find(timepoints<0);
                timepoints_after = find(timepoints>0);
                % Before
                p = P(all(ismember(P(:,1:2), timepoints_before)'), 4);
                if isempty(p)
                    P_table.(attr_name)(1) = 0;
                else
                    P_table.(attr_name)(1) = sum(p<=0.05) / length(p) * 100;
                end
                % After
                p = P(all(ismember(P(:,1:2), timepoints_after)'), 4);
                if isempty(p)
                    P_table.(attr_name)(end+1) = 0;
                else
                    P_table.(attr_name)(end+1) = sum(p<=0.05) / length(p) * 100;
                end
                % Pre-post
                for ii = 1:length(timepoints_after)
                    p = P((ismember(P(:,1),timepoints_before) & ismember(P(:,2), timepoints_after(ii))) | (ismember(P(:,2),timepoints_before) & ismember(P(:,1),timepoints_after(ii))),4);
                    if isempty(p)
                        P_table.(attr_name)(end+1) = 0;
                    else
                        P_table.(attr_name)(end+1) = sum(p<=0.05) / length(p) * 100;
                    end
                end
            end
            
        else
            for iab = 1:n_attributes
                P_table.(event_attributes_name{iab}) = 0;
            end
        end
        
        % Store data
        Sbs = table();
        Sbs.ROI_id = iroi;
        Sbs.mean_event_amplitude = {basic_stats_mean};
        Sbs.mean_event_rate = {basic_stats_rate};
        Sbs.median_event_amplitude = {basic_stats_median};
                
        S = table();
        S.ROI_id = iroi;
        S.mean_activity = {mean_activity(:)'};
        if analyze_amplitude
            S.amplitude_tests = {P_table.amplitude(:)'};
        end
        if analyze_interval
            S.interval_tests = {P_table.interval(:)'};
        end

        % Store distributions separately
        timepoints_before = find(timepoints < 0);
        timepoints_after = find(timepoints > 0);
        for iab = 1:n_attributes
            attr_name = event_attributes_name{iab};
            % Compute mean distributions
            S.([attr_name, '_dist_before']) = {distributions.(attr_name)(:, timepoints_before).'};
            for ii = 1:length(timepoints_after)
                S.([attr_name, '_dist_after_', num2str(ii)]) = {distributions.(attr_name)(:, timepoints_after(ii)).'};
            end
        end
        % Concatenate with output variable
        SPONTANEOUS_ACTIVITY = [SPONTANEOUS_ACTIVITY; S];
        SPONTANEOUS_ACTIVITY_basic_stats = [SPONTANEOUS_ACTIVITY_basic_stats; Sbs];
    end
    
    % Assess whether attributes change between timepoints
    LOGGER.trace('Computing change between time-points')
    % Gather info on which cells have a stable baseline or do change
    % distribution following surgery
    attribute_changing = table();
    attribute_changing.stable_baseline = true(n_ROIs, 1);
    attribute_changing.change_after = false(n_ROIs, length(timepoints_after));
    for iroi = 1:n_ROIs
        % Collect results from all tests
        p = [];
        for iab = 1:n_attributes
            p = [p; SPONTANEOUS_ACTIVITY.([event_attributes_name{iab}, '_tests']){iroi};];
        end
        p = p.';
        % Change between baseline sessions
        if any(p(1,:) > 0)
            attribute_changing.stable_baseline(iroi) = false;
        end
        % Change between sessions after
        p_others = p(3:end,:);
        if ~isempty(p_others)
            attribute_changing.change_after(iroi, :) = sum(p_others > 0, 2);
        end
    end

    attribute_change = {};
    stability_outcomes = [true, false];
    for it = 1:length(timepoints_after)
        timepoint_after_surgery = it;
        day_after_surgery = timepoints_after(timepoint_after_surgery);

        % Find out the direction of change in each cell and compute some sort of
        % discrimination performance
        for ii = 1:length(stability_outcomes)
            is_stable = stability_outcomes(ii);
            data = []; groups = [];
            % TP (stable baseline & change) or FP (unstable baseline & change)
            baseline_stability = attribute_changing.stable_baseline == is_stable;
            rows = baseline_stability & attribute_changing.change_after(timepoint_after_surgery);
            if any(rows)
                mean_activity = cell2mat(SPONTANEOUS_ACTIVITY.mean_activity(rows));
                activity = [mean(mean_activity(:, timepoints_before), 2), mean(mean_activity(:, day_after_surgery), 2)];
                activity = activity(:,2) - activity(:,1);
                type_of_change = attribute_changing.change_after(rows, timepoint_after_surgery);
                y = activity>0;  data=[data;activity(y)]; groups=[groups;[ones(sum(y),1)*1, type_of_change(y)]];
                y = activity<=0; data=[data;activity(y)]; groups=[groups;[ones(sum(y),1)*2, type_of_change(y)]];
            end
            % TN (stable baseline & no change) or FN (unstable baseline & no change)
            rows = baseline_stability & ~attribute_changing.change_after(timepoint_after_surgery);
            if any(rows)
                mean_activity = cell2mat(SPONTANEOUS_ACTIVITY.mean_activity(rows));
                activity = [mean(mean_activity(:, timepoints_before), 2), mean(mean_activity(:, day_after_surgery), 2)];
                activity = activity(:,2) - activity(:,1);
                type_of_change = attribute_changing.change_after(rows, timepoint_after_surgery);
                y = activity>0;  data=[data;activity(y)]; groups=[groups;[ones(sum(y),1)*3, type_of_change(y)]];
                y = activity<=0; data=[data;activity(y)]; groups=[groups;[ones(sum(y),1)*4, type_of_change(y)]];
            end
            % Store data
            if ~isempty(data)
                did_change_after = ismember(groups(:,1),[1,2]);
                attribute_change(end+1,:) = {it, is_stable, hist(groups(data>=0, 2), 0:3), hist(groups(data<0, 2), 0:3), [data(did_change_after), groups(did_change_after, 2)]};
            end
        end
    end
    
    % For each ditribution mark whether the baseline is stable and whether there
    % was a change compared to after surgery
    SPONTANEOUS_ACTIVITY.stable_baseline = attribute_changing.stable_baseline;
    SPONTANEOUS_ACTIVITY.change_after = attribute_changing.change_after;
    
    % Store data
    if ~isempty(attribute_change)
        d = [repmat({experimental_group}, size(attribute_change,1),1), attribute_change];
        d = cell2table(d, 'VariableNames',{'experimental_group', 'timepoint', 'stable_baseline', 'n_above_0', 'n_below_0', 'activity_change'});
        STATS_spontaneous_activity.(animal_ID) = d;
    else
        STATS_spontaneous_activity.(animal_ID) = NaN;
    end
    DISTRIBUTIONS.(animal_ID) = SPONTANEOUS_ACTIVITY;
    EVENT_STATS_spontaneous_activity.(animal_ID) = SPONTANEOUS_ACTIVITY_basic_stats;
end

% Write file to disk
LOGGER.info('Storing data to disk')
save(stats_filename, 'STATS_spontaneous_activity','DISTRIBUTIONS','EVENT_STATS_spontaneous_activity','PARAMETERS', '-v7.3')
LOGGER.info('(ok)', 'append',true)

%% MLint exceptions
%#ok<*AGROW>
