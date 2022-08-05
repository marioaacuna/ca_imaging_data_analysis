function spontaneous_activity_stats_epi(all_animal_names, stats_filename)

% Get LOGGER and general_configs
global LOGGER GC

% Get info to compute distributions
event_attributes_name =GC.spontaneous_activity_attributes;
n_attributes = length(event_attributes_name);
% Make flags to know whhich attributes to analyze
analyze_interval = ismember('interval', event_attributes_name);
analyze_amplitude = ismember('amplitude', event_attributes_name);
% Get parameters for analysis
PARAMETERS = struct();
LOGGER.info(['Computing distributions of the following attributes: ', strjoin(event_attributes_name, ', ')])

% Load previous runs of this analysis
if exist(stats_filename, 'file')
    f = load(stats_filename);
%     STATS_spontaneous_activi_sessy = f.STATS_spontaneous_activi_sessy;
%     DISTRIBUTIONS = f.DISTRIBUTIONS;
    EVENT_STATS_spontaneous_activi_sessy = f.EVENT_STATS_spontaneous_activi_sessy;
    clear f

else
    % Ini_sessialize output variables
%     STATS_spontaneous_activi_sessy = struct();
%     DISTRIBUTIONS = struct();
    EVENT_STATS_spontaneous_activi_sessy = struct();
end

% Get names of allowed compounts
allowed_compounds =GC.analysis.(GC.experiment_name).analysis_allowed_compounds;


% Loop through datasets
for iid = 1:length(all_animal_names)
    animal_ID = all_animal_names{iid};
    experimental_group = SQL_database.read_table_where('experiments', 'experimental_group', animal_ID,'animal_ID', 'return_as_table',false);
    
    % load traces per animal
    spikes_filename = get_filename_of('spikes', animal_ID);
    traces_filename = get_filename_of('dFF', animal_ID);
    % Skip if spikes file does not exist
    if ~exist(spikes_filename,'file')
        LOGGER.info([animal_ID, ' not fully processed. Skipping'])
        continue
    end
     
    LOGGER.info(['Processing ' animal_ID])
    % Load files
    traces = load_variable(traces_filename, 'dFF');
    spikes = load_variable(spikes_filename, 'spikes');
    n_ROIs = size(spikes,1);
    % read metadata
    METADATA = SQL_database.read_epi_trials (animal_ID);
    SP_idx = find (ismemberCellRows(METADATA.stimulus(:), {'SP'}));
    
    % Remove unused trials
    METADATA = METADATA(SP_idx, :);
    n_sessions =  length (SP_idx);
    SP_spikes = cell (1, n_sessions);
    SP_traces = cell (1, n_sessions);
    for i_sess = 1 : length (SP_idx)    
        SP_spikes{1,i_sess} = spikes(:,SP_idx(i_sess));  
        SP_traces{1,i_sess} = traces(:,SP_idx(i_sess));  
    end

    type_analysis = 'traces';
    % Get trials marked as SP (that is, spontaneous)
    % Keep only trials wi_sessh allowed compounds
%     SP_idx = SP_idx & ismember(METADATA.compound, allowed_compounds);
       
    % Allocate empty output variable
    SPONTANEOUS_ACTIVi_sessY = table();
    SPONTANEOUS_ACTIVi_sessY_basic_stats = table();  
    
    for i_roi = 1:n_ROIs
        % Log progress
%         LOGGER.trace(['Analyzing ROI ', num2str(iroi), ' (out of ', num2str(n_ROIs), ')'], 'overwri_sesse_last_message',iroi>1)
        % Allocate variables for this ROI
        mean_activi_sess = zeros(1, n_sessions);
        basic_stats_mean = NaN(n_sessions, 1);
        basic_stats_rate = NaN(n_sessions, 1);
        basic_stats_quantiles = NaN(n_sessions, 5); 
        basic_stats_median = NaN(n_sessions, 1);
        event_attributes = struct();
        for iab = 1:n_attributes
            event_attributes.(event_attributes_name{iab}) = cell(n_sessions,1);
        end
        
        % Loop through days
        for i_sess = 1:n_sessions
            % Get and concatenate trials
            switch type_analysis
                case 'spikes'
                    data = SP_spikes{1, i_sess}{i_roi, 1};
                case 'traces'
                    data = SP_traces{1, i_sess}{i_roi, 1};
                    data = data - median (data);
            end
            
            % Patameters for re-convolution
            
            trial_duration = length(data); %max(reshape(cell2mat(cellfun(@length, spikes, 'UniformOutput',false)),[],1));
            frame_rate = GC.epifluorescence_downsample_to_frame_rate;
            t = 0 : 1/frame_rate : trial_duration/frame_rate - 1/frame_rate;
           
            switch type_analysis
                case 'spikes'
                    tau_growth = 0.1;
                    tau_decay = GC.decay_time_constant_calcium_reporter.GCaMP6f;
                    tau_growth_frames = round(frame_rate * tau_growth);
                    U = (1-exp(-t./tau_growth)) .* exp(-t./(tau_decay/1000));
                    U = (U-min(U)) ./ (max(U)-min(U)); % normalize in range [0 1]
                    
                    trace_conv = data .* 0 ;% cell (size(SP_spikes{1, i_sess},1),1);
                    % Get response and convolve it with unitary response
                    this_r = data;
                    this_r_convolved = conv(this_r, U, 'full');
                    % convolve spikes for this session
                    trace_conv = this_r_convolved(1+tau_growth_frames:length(this_r)+tau_growth_frames);
                    [ampl , Ca_events] = findpeaks(trace_conv, 'MinPeakDistance',5);
                case 'traces'
                    threshold_traces = median(data)  + 2.5 .* mad(data);
                    [ampl , Ca_events] = findpeaks(data, 'MinPeakHeight',threshold_traces, 'MinPeakWidth', 4 , 'MinPeakDistance',5);
            end
            
            % Find location of gaps between events
%             gap_location = Ca_events;
            n_ev = length(Ca_events);
            data_duration = length(data) ./ frame_rate;
            freq_this_cell = n_ev ./ data_duration;
            if analyze_interval
                event_attributes.interval{i_sess} = freq_this_cell;
            end
            
            % Get ampli_sessude of events
            if analyze_amplitude
                event_attributes.amplitude{i_sess} = mean(ampl);
            end
            
            % Compute basic statistics of firing: mean ampli_sessude event and
            % frequency
            basic_stats_mean(i_sess) = mean(ampl);
            basic_stats_rate(i_sess) = freq_this_cell;
            basic_stats_quantiles(i_sess, :) = quantile(ampl, [0, .25, .5, .75, 1]);
            basic_stats_median(i_sess, :) = median(ampl);
        end
        % When the cell was not active, replace NaN wi_sessh 0
        basic_stats_mean(isnan(basic_stats_mean)) = 0;
        basic_stats_median(isnan(basic_stats_median)) = 0;
        
%         % Compute probabili_sessy distributions functions (PDF) across all days
%         distributions = struct();
%         for iab = 1:n_attributes
%             attr_name = event_attributes_name{iab};
%             distributions.(attr_name) = zeros(length(PARAMETERS.(attr_name).bins), n_sessions);
%         end
%         for i_sess = 1:n_sessions
%             for iab = 1:n_attributes
%                 % Get parameters
%                 attr_name = event_attributes_name{iab};
%                 bins = PARAMETERS.(attr_name).bins;
%                 values = event_attributes.(attr_name){i_sess};
%                 % Compute PDF
%                 [~, ~, PDF]  = make_CDF(values, bins, true);
%                 distributions.(attr_name)(:,i_sess) = cumsum(PDF(:,1));
%             end
%         end
        
        % Perform statistics between 'before' and 'after' surgery
%         P_table = struct();
%         for iab = 1:n_attributes
%             P_table.(event_attributes_name{iab}) = [];
%         end
%         if n_sessions > 1
%             pairs = nchoosek(1:n_sessions, 2);
%             for iab = 1:n_attributes
%                 attr_name = event_attributes_name{iab};
%                 % Perform a Kolmogorov-Smirnov tests on all pairs of sessions
%                 P = ones(size(pairs,1), 4);
%                 for ipair = 1:size(pairs,1)
%                     [~, p] = kstest2(event_attributes.(attr_name){pairs(ipair, 1)}, event_attributes.(attr_name){pairs(ipair, 2)});
%                     P(ipair, 1:3) = [pairs(ipair,:), p];
%                 end
%                 % Correct p-values for multiple comparisons
%                 [~, P(:,4)] = adjust_Pvalues(P(:,3));
% 
%                 % Count number of significant tests in baseline, pre-post, and post comparisons
%                 timepoints_before = find(timepoints<0);
%                 timepoints_after = find(timepoints>0);
%                 % Before
%                 p = P(all(ismember(P(:,1:2), timepoints_before)'), 4);
%                 if isempty(p)
%                     P_table.(attr_name)(1) = 0;
%                 else
%                     P_table.(attr_name)(1) = sum(p<=0.05) / length(p) * 100;
%                 end
%                 % After
%                 p = P(all(ismember(P(:,1:2), timepoints_after)'), 4);
%                 if isempty(p)
%                     P_table.(attr_name)(end+1) = 0;
%                 else
%                     P_table.(attr_name)(end+1) = sum(p<=0.05) / length(p) * 100;
%                 end
%                 % Pre-post
%                 for ii = 1:length(timepoints_after)
%                     p = P((ismember(P(:,1),timepoints_before) & ismember(P(:,2), timepoints_after(ii))) | (ismember(P(:,2),timepoints_before) & ismember(P(:,1),timepoints_after(ii))),4);
%                     if isempty(p)
%                         P_table.(attr_name)(end+1) = 0;
%                     else
%                         P_table.(attr_name)(end+1) = sum(p<=0.05) / length(p) * 100;
%                     end
%                 end
%             end
%             
%         else
%             for iab = 1:n_attributes
%                 P_table.(event_attributes_name{iab}) = 0;
%             end
%         end
        
        % Store data
        Sbs = table();
        Sbs.ROI_id = i_roi;
        Sbs.mean_event_ampli_sessude = {basic_stats_mean};
        Sbs.mean_event_rate = {basic_stats_rate};
        Sbs.median_event_ampli_sessude = {basic_stats_median};
                
%         S = table();
%         S.ROI_id = i_roi;
%         S.mean_activi_sessy = {mean_activi_sess(:)'};
%         if analyze_amplitude
%             S.ampli_sessude_tests = {P_table.ampli_sessude(:)'};
%         end
%         if analyze_interval
%             S.interval_tests = {P_table.interval(:)'};
%         end

%         % Store distributions separately
%         timepoints_before = find(timepoints < 0);
%         timepoints_after = find(timepoints > 0);
%         for iab = 1:n_attributes
%             attr_name = event_attributes_name{iab};
%             % Compute mean distributions
%             S.([attr_name, '_dist_before']) = {distributions.(attr_name)(:, timepoints_before).'};
%             for ii = 1:length(timepoints_after)
%                 S.([attr_name, '_dist_after_', num2str(ii)]) = {distributions.(attr_name)(:, timepoints_after(ii)).'};
%             end
%         end
        % Concatenate wi_sessh output variable
        %SPONTANEOUS_ACTIVi_sessY = [SPONTANEOUS_ACTIVi_sessY; Sbs];
        SPONTANEOUS_ACTIVi_sessY_basic_stats = [SPONTANEOUS_ACTIVi_sessY_basic_stats; Sbs];
    end
    
%     % Assess whether attributes change between timepoints
%     LOGGER.trace('Computing change between time-points')
%     % Gather info on which cells have a stable baseline or do change
%     % distribution following surgery
%     attribute_changing = table();
%     attribute_changing.stable_baseline = true(n_ROIs, 1);
%     attribute_changing.change_after = false(n_ROIs, length(timepoints_after));
%     for i_roi = 1:n_ROIs
%         % Collect results from all tests
%         p = [];
%         for iab = 1:n_attributes
%             p = [p; SPONTANEOUS_ACTIVi_sessY.([event_attributes_name{iab}, '_tests']){i_roi};];
%         end
%         p = p.';
%         % Change between baseline sessions
%         if any(p(1,:) > 0)
%             attribute_changing.stable_baseline(i_roi) = false;
%         end
%         % Change between sessions after
%         p_others = p(3:end,:);
%         if ~isempty(p_others)
%             attribute_changing.change_after(i_roi, :) = sum(p_others > 0, 2);
%         end
%     end
% 
%     attribute_change = {};
%     stabili_sessy_outcomes = [true, false];
%     for i_sess = 1:length(timepoints_after)
%         timepoint_after_surgery = i_sess;
%         day_after_surgery = timepoints_after(timepoint_after_surgery);
% 
%         % Find out the direction of change in each cell and compute some sort of
%         % discrimination performance
%         for ii = 1:length(stabili_sessy_outcomes)
%             is_stable = stabili_sessy_outcomes(ii);
%             data = []; groups = [];
%             % TP (stable baseline & change) or FP (unstable baseline & change)
%             baseline_stabili_sessy = attribute_changing.stable_baseline == is_stable;
%             rows = baseline_stabili_sessy & attribute_changing.change_after(timepoint_after_surgery);
%             if any(rows)
%                 mean_activi_sess = cell2mat(SPONTANEOUS_ACTIVi_sessY.mean_activi_sessy(rows));
%                 activi_sessy = [mean(mean_activi_sess(:, timepoints_before), 2), mean(mean_activi_sess(:, day_after_surgery), 2)];
%                 activi_sessy = activi_sessy(:,2) - activi_sessy(:,1);
%                 type_of_change = attribute_changing.change_after(rows, timepoint_after_surgery);
%                 y = activi_sessy>0;  data=[data;activi_sessy(y)]; groups=[groups;[ones(sum(y),1)*1, type_of_change(y)]];
%                 y = activi_sessy<=0; data=[data;activi_sessy(y)]; groups=[groups;[ones(sum(y),1)*2, type_of_change(y)]];
%             end
%             % TN (stable baseline & no change) or FN (unstable baseline & no change)
%             rows = baseline_stabili_sessy & ~attribute_changing.change_after(timepoint_after_surgery);
%             if any(rows)
%                 mean_activi_sess = cell2mat(SPONTANEOUS_ACTIVi_sessY.mean_activi_sessy(rows));
%                 activi_sessy = [mean(mean_activi_sess(:, timepoints_before), 2), mean(mean_activi_sess(:, day_after_surgery), 2)];
%                 activi_sessy = activi_sessy(:,2) - activi_sessy(:,1);
%                 type_of_change = attribute_changing.change_after(rows, timepoint_after_surgery);
%                 y = activi_sessy>0;  data=[data;activi_sessy(y)]; groups=[groups;[ones(sum(y),1)*3, type_of_change(y)]];
%                 y = activi_sessy<=0; data=[data;activi_sessy(y)]; groups=[groups;[ones(sum(y),1)*4, type_of_change(y)]];
%             end
%             % Store data
%             if ~isempty(data)
%                 did_change_after = ismember(groups(:,1),[1,2]);
%                 attribute_change(end+1,:) = {i_sess, is_stable, hist(groups(data>=0, 2), 0:3), hist(groups(data<0, 2), 0:3), [data(did_change_after), groups(did_change_after, 2)]};
%             end
%         end
%     end
%     
%     % For each di_sessribution mark whether the baseline is stable and whether there
%     % was a change compared to after surgery
%     SPONTANEOUS_ACTIVi_sessY.stable_baseline = attribute_changing.stable_baseline;
%     SPONTANEOUS_ACTIVi_sessY.change_after = attribute_changing.change_after;
    
    % Store data
%     if ~isempty(attribute_change)
%         d = [repmat({experimental_group}, size(attribute_change,1),1), attribute_change];
%         d = cell2table(d, 'VariableNames',{'experimental_group', 'timepoint', 'stable_baseline', 'n_above_0', 'n_below_0', 'activi_sessy_change'});
%         STATS_spontaneous_activi_sessy.(animal_ID) = d;
%     else
%         STATS_spontaneous_activi_sessy.(animal_ID) = NaN;
%     end
%     DISTRIBUTIONS.(animal_ID) = SPONTANEOUS_ACTIVi_sessY;
    EVENT_STATS_spontaneous_activi_sessy.(animal_ID) = SPONTANEOUS_ACTIVi_sessY_basic_stats;
end

% Wri_sesse file to disk
LOGGER.info('Storing data to disk')
save(stats_filename, 'EVENT_STATS_spontaneous_activi_sessy', '-v7.3')
LOGGER.info('(ok)', 'append',true)

%% MLint exceptions
%#ok<*AGROW>
