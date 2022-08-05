function single_cell_response_modulation_AUROC(all_animal_IDs, do_overwrite)

% Get general_configs and LOGGER
global GC LOGGER

% Enable toolboxes used by this function
toolboxes_to_use = '_plotting';
toggle_toolbox(toolboxes_to_use, 'on')
previous_warning_state = warning('off', 'export_fig:transparency');

% List allowed stimuli and compounds for the analysis
allowed_stimuli = GC.analysis_allowed_stimuli.(GC.experiment_name);
allowed_compounds = GC.analysis_allowed_compounds.(GC.experiment_name);

% Load file containing single-cell response modulation statistics fot this animal
stats_filename = get_filename_of('response_modulation_stats', GC.experiment_name);
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');
all_animal_names = fieldnames(STATS_response_modulation);
if isempty(all_animal_IDs)
    all_animal_IDs = all_animal_names;
end
all_animal_IDs = natsort(all_animal_IDs);

% Loop through animals
for iid = 1:length(all_animal_IDs)
    % Skip if not completely analyzed
    animal_ID = all_animal_IDs{iid};
    if ~ismember(animal_ID, all_animal_names)
        LOGGER.warn([animal_ID, ' not yet analyzed. Skipping'])
        continue
    end
    
    % Skip if output file exists and user requested not to overwrite previous
    % runs of this script
    % Get name of folder where figures will be saved
    folder = [GC.data_root_path, GC.plots_path, GC.experiment_name, filesep, 'response_modulation', filesep];
    if ~exist(folder, 'dir'), mkdir(folder), end
    figure_filename = [folder, animal_ID, '_AUROC.pdf'];
    if exist(figure_filename, 'file') && ~do_overwrite
        LOGGER.info(['Plot for ' animal_ID, ' already exists. Skipping'])
        continue
    end
    
    LOGGER.info(['Processing ', animal_ID])
    % Get statistics
    STATS = STATS_response_modulation.(animal_ID);

    METADATA = SQL_database.read_table_where('trials', {'+','stimulus','experimental_condition','date'}, {animal_ID, GC.experiment_name}, {'animal_ID', 'experiment_name'});
    % Add column to mark unwanted trials
    METADATA(:, 'keep') = num2cell(true(height(METADATA),1));
    % Remove trials where a compound had been administered
    METADATA.keep(~ismember(METADATA.compound, allowed_compounds)) = false;
    if sum(METADATA.keep) < 1, continue, end
    % Get frame rate
    frame_rate = max(METADATA.frame_rate);
    % Load stimuli
    STIMULI = struct();
    for istim = 1:length(allowed_stimuli)
        n_frames = METADATA{ismember(table2cell(METADATA(:, 'stimulus')), allowed_stimuli{istim}), 'n_frames'}; n_frames = n_frames(1);
        STIMULI.(allowed_stimuli{istim}) = Metadata_workflow.load_stimuli(allowed_stimuli{istim}, 'frame_rate',frame_rate, 'pre_stimulus_window', GC.ROC_window, 'post_stimulus_window',GC.ROC_window, 'n_frames',n_frames);
    end
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

    LOGGER.trace('Loading spikes')
    spikes_filename = get_filename_of('spikes', animal_ID);
    spikes = load_variable(spikes_filename, 'spikes');
    % Get the number of ROIs
    n_ROIs = size(spikes, 1);
    all_ROIS = unique(STATS.ROI_id);
    % Get time constant for decay of calcium indicator
    calcium_reporter = SQL_database.read_table_where('experiments', 'Calcium_indicator', animal_ID,'animal_ID', 'return_as_table',false);
    decay_time_constant = GC.decay_time_constant_calcium_reporter.(calcium_reporter);

    % Make empty temporary folder to store each image
    temp_folder = [GC.temp_dir, animal_ID, '_response_modulation', filesep];
    if exist(temp_folder, 'dir'), rmdir(temp_folder, 's'), end
    mkdir(temp_folder)
    PDF_filenames = {};  % This variable will contain the list of printed PDFs
 
    % Loop through ROIS
    for iroi = 1:n_ROIs
        % Log current ROI
        LOGGER.trace(['ROI ', num2str(iroi), '/', num2str(n_ROIs)])
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
            N_SHOWN_TRIALS = NaN(1, 2);
            OFFSET = NaN(1, 2);

            % Make plot
            fig = figure('color','w', 'units','normalized','outerposition',[0 0 1 1], 'Visible','off');
            clf, ax = [];
            mb=.05; sh=.05; mt=.03;
            ax(1) = subaxis(1, 2, 1, 'mb',mb, 'ml',.05, 'sh',sh, 'mt',mt);
            ax(2) = subaxis(1, 2, 2, 'mb',mb, 'mr',.05, 'sh',sh, 'mt',mt);
            set(ax, 'Box','off', 'FontSize',12, 'YColor','w', 'TickDir','out', 'TickLength',[.005, 0], 'Clipping','off')
            
            % Get spontaneous data
            switch GC.experiment_name
                case 'ACC_CCI_anesth'
                    training_data = all_data(ismember(table2cell(conds_trial(:,'stimulus')),'SP') & ismember(conds_trial(:,'day_from_surgery'),conds(icond,'day_from_surgery')));
                case 'ACC_pain_LP-211'
                    training_data = all_data(ismember(table2cell(conds_trial(:,'stimulus')),'SP') & ismember(conds_trial(:,'compound_phase'),conds(icond,'compound_phase')));
            end
            SP_data2plot = cell2mat(convolve_spikes(training_data(:), frame_rate, decay_time_constant));
            N_SHOWN_TRIALS(1) = size(SP_data2plot, 1);
            % Get evoked data
            cond_data = all_data(trials_idx);
            EV_data2plot = cell2mat(convolve_spikes(cond_data(:), frame_rate, decay_time_constant));
            N_SHOWN_TRIALS(2) = size(EV_data2plot, 1);
            % Get optimal offset between traces of panel with more trials
            if N_SHOWN_TRIALS(1) > N_SHOWN_TRIALS(2)
                data = SP_data2plot;
                idx_offset = 1;
                trial_proportion = N_SHOWN_TRIALS(1) / N_SHOWN_TRIALS(2);
            else
                data = EV_data2plot;
                idx_offset = 2;
                trial_proportion = N_SHOWN_TRIALS(2) / N_SHOWN_TRIALS(1);
            end
            % Get maximum of data to plot
            m = max(data(:));
            if m == 0, m = 1; end  % Assign mock value to create offset between trials in case of no activity
            OFFSET(idx_offset) = m / 4;
            % Set a proportional offset on other panel
            OFFSET(setdiff(1:2, idx_offset)) = OFFSET(idx_offset) * trial_proportion;
            
            % Plot data
            % Get trial duration and make time axis
            t = linspace(0, size(SP_data2plot,2)./frame_rate, size(SP_data2plot,2));
            set(fig, 'CurrentAxes', ax(1)); cla, hold on
            for itrial = 1:size(SP_data2plot,1)
                offset = (size(SP_data2plot,1)-itrial) * OFFSET(1);
                % Plot outline
                y = SP_data2plot(itrial,:);
                if all(y==0), ls='--'; else, ls='-'; end
                plot(t, y+offset, 'color','k', 'linestyle',ls);
            end
            % Get ylim range
            limits = ylim();
            limits(1) = -.1;
            % Adjust axes
            xlim(t([1,end])), ylim([-.2, limits(2)])
            ylabel('Trial {\itno.}', 'Color','k', 'FontSize',16)
            xlabel('Time (s)', 'Color','k', 'FontSize',16)
            title('Spontaneous activity', 'Color','k', 'FontWeight','normal', 'FontSize',18)

            % Get trial duration and make time axis
            trial_duration = unique(cellfun(@length, cond_data));
            t = linspace(0, trial_duration./frame_rate, trial_duration) - (STIMULI.(stimulus_name).response_window(1) ./ frame_rate);
            set(fig, 'CurrentAxes', ax(2)); cla, hold on
            for itrial = 1:n_trials
                offset = (n_trials-itrial) * OFFSET(2);
                % Plot outline
                y = EV_data2plot(itrial,:);
                if all(y==0), ls='--'; else, ls='-'; end
                h = plot(t, y+offset, 'color','k', 'linestyle',ls);
                % Fill baseline
                bl_idx = [STIMULI.(stimulus_name).baseline_window, STIMULI.(stimulus_name).baseline_window(end)+1];
                x = t(bl_idx);
                y = EV_data2plot(itrial, bl_idx);
                fill([x, fliplr(x)], [y, zeros(size(y))]+offset, [0,.75,1], 'FaceAlpha',0.5, 'EdgeColor','none')
                % Fill response
                x = t(STIMULI.(stimulus_name).response_window);
                y = EV_data2plot(itrial, STIMULI.(stimulus_name).response_window);
                fill([x, fliplr(x)], [y, zeros(size(y))]+offset, 'r', 'FaceAlpha',0.75, 'EdgeColor','none')
                % Move outline above
                uistack(h,'top')
            end
            % Get ylim range
            limits = ylim();
            limits(1) = -.1;
            % Show baseline window
            bl = STIMULI.(stimulus_name).baseline_window - STIMULI.(stimulus_name).response_window(1);
            bl = bl(1) / frame_rate;
            h_bl = plot([bl, bl], limits, 'color',[0,.75,1], 'LineWidth',2);
            uistack(h_bl,'bottom')
            % Show response window
            ev = STIMULI.(stimulus_name).response_window - STIMULI.(stimulus_name).response_window(1);
            ev = ev(end) / frame_rate;
            h_ev = plot([ev, ev], limits, 'color','r', 'LineWidth',2);
            uistack(h_ev,'bottom')
            % Highlight stimulus onset time
            h_t0 = plot([0, 0], limits, 'color',[.7, .7, .7], 'LineWidth',2);
            uistack(h_t0,'bottom')
            % Adjust axes
            xlim(t([1,end])), ylim(limits)
            ylabel('Trial {\itno.}', 'Color','k', 'FontSize',16)
            xlabel('Time from stimulus onset (s)', 'Color','k', 'FontSize',16)

            % Get statistics for this condition
            row = ismember(STATS(:, ['stimulus', columns_experimental_condition]), conds(icond, :), 'rows') & ismember(STATS.ROI_id, all_ROIS(iroi));
            p = STATS.p(row);
            AUROC = STATS.AUROC(row);
            if isempty(AUROC) || isempty(p)
                title([stimulus_name, ': not enough spontaneous data'])
            else
                % Show statistics in title
                if p <= 0.05
                    activity = STATS.activity{row};
                    if sign(mean(activity(:,2) - activity(:,1))) == 1
                        str = '> than chance, ';
                    else
                         str = '< than chance, ';
                    end
                else
                    str = '';
                end
                p = [' = ', num2str(p)];
                title([stimulus_name, ': AUROC ', num2str(AUROC), ' (', str, '{\itp}', p, ')'], 'FontWeight','normal', 'FontSize',18)
            end
            % Add number of ROI to top left corner
            axes('Position',[0, 0, 1, 1], 'XLim',[0, 1], 'YLim',[0, 1], 'Box','off', 'Visible','off', 'Units','normalized', 'Clipping','off', 'HitTest','off');
            text(0, 1, ['ROI ', num2str(iroi), ' - ', cond_name], 'Hor','left', 'Ver','top', 'FontSize',20, 'FontWeight','bold', 'Interpreter','none')
            
            % Add scale bar to top rigth corner of the righthand graph
            xl = get(ax(2), 'XLim'); yl = get(ax(2), 'YLim');
            x_offset = diff(xl) * .05;
            text_offset = diff(xl) * .01;
            x_position = xl(2) + x_offset;
            m = max([SP_data2plot(:); EV_data2plot(:)]);
            if m > 0
                plot(ax(2), [x_position, x_position], [yl(2)-m/2, yl(2)], 'color','k', 'LineWidth',2);
                text(x_position - text_offset, yl(2), {'50%'; '\DeltaF/F_0'}, 'Parent',ax(2), 'Hor','right', 'Ver','top', 'FontSize',12);
            end
            
            % Export figure to .pdf
            filename = [temp_folder, num2str(length(PDF_filenames)+1, '%05i'), '.pdf'];
            PDF_filenames{end+1} = filename;
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters',fig)
            close(fig)
        end
    end

    % Merge all temporary PDFs in one on the remote server
    LOGGER.info('Concatenating PDF files on disk')
	if exist(figure_filename, 'file')
		delete(figure_filename)
	end
    append_pdfs(figure_filename, PDF_filenames{:})
    printFilePath(figure_filename, figure_filename, false)
    % Remove folder containing temporary PDF files from disk
    rmdir(temp_folder, 's')
end

% Disable toolboxes used by this function
toggle_toolbox(toolboxes_to_use, 'off')
warning(previous_warning_state)

%% MLint exceptions
%#ok<*AGROW>
