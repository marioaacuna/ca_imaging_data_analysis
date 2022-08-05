clear, clc, close all hidden

% Load general configs
global GC
if ~isfield(GC, 'repository_root_path'), load_repository(1); end
GC.experiment_name = 'ACC_CCI_anesth';
time_bin_activity_n_frames = GC.assemblies_time_bin_activity_n_frames;

toggle_toolbox('_plotting', 'on')
FP = figure_properties();

% Make folders to store figures
root_figure_path = os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name, 'ensembles');
if ~exist(root_figure_path, 'dir'), mkdir(root_figure_path), end
% Make subfolders
path_inside_vs_outside_ensemble       = os.path.join(root_figure_path, 'inside_vs_outside_ensemble');
path_cell_activity_ensembles          = os.path.join(root_figure_path, 'cell_activity_ensembles');
path_mean_activity_ensembles          = os.path.join(root_figure_path, 'mean_activity_ensembles');
path_mean_evoked_activity_to_SP       = os.path.join(root_figure_path, 'mean_evoked_activity_compared_to_SP');
path_mean_evoked_activity_to_baseline = os.path.join(root_figure_path, 'mean_evoked_activity_compared_to_baseline');
path_ensembles_activation             = os.path.join(root_figure_path, 'ensembles_activation');
path_mean_ensembles_activation        = os.path.join(root_figure_path, 'mean_ensemble_activation');
if ~exist(path_inside_vs_outside_ensemble, 'dir'),       mkdir(path_inside_vs_outside_ensemble), end
if ~exist(path_cell_activity_ensembles, 'dir'),          mkdir(path_cell_activity_ensembles), end
if ~exist(path_mean_activity_ensembles, 'dir'),          mkdir(path_mean_activity_ensembles), end
if ~exist(path_mean_evoked_activity_to_SP, 'dir'),       mkdir(path_mean_evoked_activity_to_SP), end
if ~exist(path_mean_evoked_activity_to_baseline, 'dir'), mkdir(path_mean_evoked_activity_to_baseline), end
if ~exist(path_ensembles_activation, 'dir'),             mkdir(path_ensembles_activation), end
if ~exist(path_mean_ensembles_activation, 'dir'),        mkdir(path_mean_ensembles_activation), end

% Load assemblies info
assemblies_info_filename = os.path.join(GC.repository_root_path, 'Figures_paper_FMP', '_data', 'assemblies_info.mat');
disp('Loading assemblies info')
CELL_ASSEMBLIES_INFO = load_variable(assemblies_info_filename, 'CELL_ASSEMBLIES_INFO');
ASSEMBLIES_INFO = load_variable(assemblies_info_filename, 'ASSEMBLIES_INFO');


% Set which figures to plot
overwrite_figures = false;
debug_mode = 0
plot_figure_inside_vs_outside_ensemble                = 1;
plot_figure_cell_activity_ensembles                   = 1;
plot_figure_mean_activity_ensembles                   = 1;
plot_figure_mean_evoked_activity_compared_to_SP       = 1;
plot_figure_mean_evoked_activity_compared_to_baseline = 1;
plot_figure_ensembles_activation                      = 1;
plot_figure_mean_ensemble_activation                  = 1;

% Set parameters for plotting
scatter_marker_size = 40;
plot_marker_size    = 8;

% Get list of animals
all_mice_list = animal_list();
n_mice = length(all_mice_list);


if debug_mode
    figure_visibility = 'on';
    print_figure = false;
else
    figure_visibility = 'off';
    print_figure = true;
end

%%
for i_mouse = 1:n_mice
    animal_ID = all_mice_list{i_mouse};
    disp(['Processing ', animal_ID])

    % Get data for this mouse
    cell_info = CELL_ASSEMBLIES_INFO(ismember(CELL_ASSEMBLIES_INFO.animal_ID, animal_ID), :);
    stimuli_analyzed = unique(cell_info.stimulus);
    n_stimuli = length(stimuli_analyzed);
    ensembles = unique(cell_info.ensemble);
    n_ensembles = length(ensembles);
    sessions = unique(cell_info.day_from_surgery);
    n_sessions = length(sessions);


    %% FIGURE: Activity of cells inside vs outside of ensemble
    filename = os.path.join(path_inside_vs_outside_ensemble, [animal_ID, '_inside_vs_outside_ensemble.pdf']);
    if plot_figure_inside_vs_outside_ensemble && (~exist(filename, 'file') || overwrite_figures)
        % Open figure
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 751, 1772, 603]);
        clf
        ax = NaN(1, n_stimuli);

        for i_stim = 1:n_stimuli
            rows = ismember(cell_info{:, 'stimulus'}, stimuli_analyzed{i_stim});
            data = cell_info{rows, {'day_from_surgery', 'ensemble', 'activity_in', 'activity_out'}};
            
            % Open axes
            ax(i_stim) = subaxis(1, length(stimuli_analyzed), i_stim, 'mb',.15, 'mt',.1, 'ml',.05, 'mr',.05);
            plot([-1, n_sessions + 1], [0, 0], 'color','k', 'XLimInclude','off')
            hold on
            for i_sess = 1:n_sessions
                % Get data for this session
                x = data(ismember(data(:, 1), sessions(i_sess)), 2:4);
                % Get number of ensembles and reset their index
                [~, ~, x(:, 1)] = unique(x(:, 1));
                n_ensembles = length(unique(x(:, 1)));
                % Set colors for each emsemble
                ensemble_colors = distinguishable_colors(n_ensembles, {'w', 'k'});

                % Plot activity of each cell
                for i_ens = 1:n_ensembles
                    ensemble_idx = ismember(x(:, 1), i_ens);
                    if any(ensemble_idx)
                        h = scatter(ones(sum(ensemble_idx), 1) * i_sess, x(ensemble_idx, 2) - x(ensemble_idx, 3), 'jitter','on', 'jitterAmount',0.3);
                        h.MarkerFaceColor = ensemble_colors(i_ens, :);
                        h.MarkerEdgeColor = 'w';
                        h.SizeData = scatter_marker_size;
                    end
                end
            end
        end
        % Set same limits
        lims = cell2mat(get(ax, 'YLim'));
        lims = [min(lims(:)), max(lims(:))];
        set(ax, 'XLim',[0, n_sessions + 1], 'YLim',lims, 'TickDir','out', 'FontSize',FP.font1.size.tick_labels, 'Box','off', ...
            'XTick',1:n_sessions, 'XTickLabel',sessions)

        % Set title
        for i_stim = 1:n_stimuli
            title(ax(i_stim), stimuli_analyzed{i_stim}, 'Interpreter','none', 'FontSize',FP.font1.size.title)
            xlabel(ax(i_stim), 'Days from surgery', 'FontSize',FP.font1.size.axes_labels)
        end
        ylabel(ax(1), 'In - not in an ensemble', 'FontSize',FP.font1.size.axes_labels)

        % Write animal ID in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end

    
    %% FIGURE: Activity of cells for different stimuli
    filename = os.path.join(path_cell_activity_ensembles, [animal_ID, '_cell_activity_ensembles.pdf']);
    if plot_figure_cell_activity_ensembles && (~exist(filename, 'file') || overwrite_figures)
        % Open figure
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 500, 1300, 800]);
        clf
        % Set position of each plot in the subplot grid
        n_subplot_rows = 2;  % pre and post surgery
        n_subplot_columns = max([sum(sessions < 0), sum(sessions >= 0)]);
        subplot_position = NaN(n_subplot_rows, n_subplot_columns);
        subplot_position(1, 1:sum(sessions < 0)) = sessions(sessions < 0);
        subplot_position(2, 1:sum(sessions >= 0)) = sessions(sessions >= 0);
        subplot_position_idx = reshape(1:n_subplot_rows * n_subplot_columns, n_subplot_columns, n_subplot_rows)';
        
        % Allocate variable for axes handles
        ax = NaN(1, n_sessions);

        for i_sess = 1:n_sessions
            data = cell_info(ismember(cell_info.day_from_surgery, sessions(i_sess)), :);
            
            % Get list of cells in ensembles
            ensembles = unique(data.ensemble);
            n_ensembles = length(ensembles);
            cells_in_ensembles = unique(data{:, 'cell'});
            n_cells_in_ensembles = length(cells_in_ensembles);
            % Set colors
            ensemble_colors = distinguishable_colors(n_ensembles, {'w', 'k'});

            % Make jitter array
            jitter = jitter_points(n_cells_in_ensembles, .2);

            % Plot each cell's activity
            subplot_location = subplot_position_idx(ismember(subplot_position, sessions(i_sess)));
            ax(i_sess) = subaxis(n_subplot_rows, n_subplot_columns, subplot_location, 'mb',.07, 'mt',.08, 'ml',.05, 'mr',.01, 'sv',.1);
            cla, hold on
            
            H = [];
            y_ens = cell(n_ensembles, n_stimuli);
            for i_cell = 1:n_cells_in_ensembles
                % Get data of this cell
                x_rows = ismember(data{:, 'cell'}, cells_in_ensembles(i_cell)) & ismember(data{:, 'epoch'}, 'evoked');
                x = data(x_rows, {'ensemble', 'stimulus', 'activity_in'});
                ensembles_this_cell = unique(x.ensemble);
                n_ensembles_this_cell = length(ensembles_this_cell);
                y = NaN(n_stimuli, n_ensembles_this_cell);
                
                for i_ens = 1:n_ensembles_this_cell
                    for i_stim = 1:n_stimuli
                        y_rows = ismember(x{:, 'stimulus'}, stimuli_analyzed{i_stim}) & ismember(x{:, 'ensemble'}, ensembles_this_cell(i_ens));
                        this_y = x{y_rows, {'activity_in'}};
                        if this_y > 0
                            y(i_stim, i_ens) = this_y;
                        end
                    end

                    % Plot data
                    plot((1:n_stimuli) + jitter(i_cell), y(:, i_ens), 'Color',[.7,.7,.7], 'LineWidth',.5);
                    h = plot((1:n_stimuli) + jitter(i_cell), y(:, i_ens), 'o', 'MarkerFaceColor',ensemble_colors(ensembles_this_cell(i_ens), :), 'MarkerEdgeColor','w', 'MarkerSize',plot_marker_size);
                    H(end + 1) = h;

                    % Keep data for later
                    for i_stim = 1:n_stimuli
                        y_ens{ensembles_this_cell(i_ens), i_stim} = [y_ens{ensembles_this_cell(i_ens), i_stim}; y(i_stim)];
                    end
                end
            end
            
            % Move all points to the top
            uistack(H, 'top')
        end
        
        % Set same limits
        % lims = cell2mat(get(ax, 'YLim')); lims = [0, max(lims(:))]; set(ax, 'YLim',lims)
        % Adjust axes appearance
        set(ax, 'TickDir','out', 'FontSize',FP.font1.size.tick_labels, 'Box','off', 'XColor','w')

        % Add title and labels
        for i_sess = 1:n_sessions
            title(ax(i_sess), ['session: ', num2str(sessions(i_sess), '%+i')], 'Interpreter','none', 'FontSize',FP.font1.size.title)
            for i_stim = 1:n_stimuli
                lims = get(ax(i_sess), 'YLim');
                text(i_stim, -lims(2) * .05, stimuli_analyzed{i_stim}, 'Hor','center', 'Ver','top', 'FontSize',FP.font1.size.axes_labels, 'Parent',ax(i_sess), 'Interpreter','none')
            end
        end
        for i_row = 1:n_subplot_rows
            if size(ax, 1) > i_row
                ylabel(ax(subplot_position_idx(i_row, 1)), 'Activity', 'FontSize',FP.font1.size.axes_labels)
            end
        end

            
        % Write animal ID and session in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end
        
    
    %% FIGURE: Mean activity of cells for different stimuli
    filename = os.path.join(path_mean_activity_ensembles, [animal_ID, '_mean_activity_ensembles.pdf']);
    if plot_figure_mean_activity_ensembles && (~exist(filename, 'file') || overwrite_figures)
        % Open figure
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 500, 1300, 800]);
        clf
        % Set position of each plot in the subplot grid
        n_subplot_rows = 2;  % pre and post surgery
        n_subplot_columns = max([sum(sessions < 0), sum(sessions >= 0)]);
        subplot_position = NaN(n_subplot_rows, n_subplot_columns);
        subplot_position(1, 1:sum(sessions < 0)) = sessions(sessions < 0);
        subplot_position(2, 1:sum(sessions >= 0)) = sessions(sessions >= 0);
        subplot_position_idx = reshape(1:n_subplot_rows * n_subplot_columns, n_subplot_columns, n_subplot_rows)';
        % Allocate variable for axes handles
        ax = NaN(1, n_sessions);

        for i_sess = 1:n_sessions
            data = cell_info(ismember(cell_info.day_from_surgery, sessions(i_sess)), :);

            % Get list of cells in ensembles
            ensembles = unique(data.ensemble);
            n_ensembles = length(ensembles);
            cells_in_ensembles = unique(data{:, 'cell'});
            n_cells_in_ensembles = length(cells_in_ensembles);
            % Set colors
            ensemble_colors = distinguishable_colors(n_ensembles, {'w', 'k'});

            y_ens = cell(n_ensembles, n_stimuli);
            for i_cell = 1:n_cells_in_ensembles
                % Get data of this cell
                x_rows = ismember(data{:, 'cell'}, cells_in_ensembles(i_cell)) & ismember(data{:, 'epoch'}, 'evoked');
                x = data(x_rows, {'ensemble', 'stimulus', 'activity_in'});
                ensembles_this_cell = unique(x.ensemble);
                n_ensembles_this_cell = length(ensembles_this_cell);
                y = NaN(n_stimuli, n_ensembles_this_cell);
                
                for i_ens = 1:n_ensembles_this_cell
                    for i_stim = 1:n_stimuli
                        y_rows = ismember(x{:, 'stimulus'}, stimuli_analyzed{i_stim}) & ismember(x{:, 'ensemble'}, ensembles_this_cell(i_ens));
                        this_y = x{y_rows, {'activity_in'}};
                        if this_y > 0
                            y(i_stim, i_ens) = this_y;
                        end
                    end

                    % Keep data for later
                    for i_stim = 1:n_stimuli
                        y_ens{ensembles_this_cell(i_ens), i_stim} = [y_ens{ensembles_this_cell(i_ens), i_stim}; y(i_stim)];
                    end
                end
            end

            % Plot mean cells activity
            subplot_location = subplot_position_idx(ismember(subplot_position, sessions(i_sess)));
            ax(i_sess) = subaxis(n_subplot_rows, n_subplot_columns, subplot_location, 'mb',.07, 'mt',.08, 'ml',.05, 'mr',.01, 'sv',.1);
            cla, hold on

            jitter = jitter_points(n_ensembles, .3);
            for i_ens = 1:n_ensembles
                y = cell2mat(y_ens(i_ens, :));
                avg = nanmean(y); avg(~isfinite(avg)) = 0;
                err = sem(y, 1); err(~isfinite(err)) = 0;
                errorbar((1:n_stimuli) + jitter(i_ens), avg, err, 'o', 'Color','k', 'MarkerFaceColor',ensemble_colors(i_ens, :), 'MarkerEdgeColor','w', 'MarkerSize',plot_marker_size, 'CapSize',0);
            end
        end
        
        % Set same y-limits
        % lims = cell2mat(get(ax, 'YLim')); lims = [0, max(lims(:))]; set(ax, 'YLim',lims)
        % Set same x-limits
        lims = get(ax, 'XLim'); if iscell(get(ax, 'XLim')); lims = cell2mat(lims); end; lims = [0, max(lims(:))]; set(ax, 'XLim',lims)
        % Adjust axes appearance
        set(ax, 'TickDir','out', 'FontSize',FP.font1.size.tick_labels, 'Box','off', 'XColor','w')

        % Add title and labels
        for i_sess = 1:n_sessions
            title(ax(i_sess), ['session: ', num2str(sessions(i_sess), '%+i')], 'Interpreter','none', 'FontSize',FP.font1.size.title)
            for i_stim = 1:n_stimuli
                lims = get(ax(i_sess), 'YLim');
                lims_range = lims(2) - lims(1);
                text(i_stim, lims(1) - lims_range * .05, stimuli_analyzed{i_stim}, 'Hor','center', 'Ver','top', 'FontSize',FP.font1.size.axes_labels, 'Parent',ax(i_sess), 'Interpreter','none')
            end
        end
        for i_row = 1:2
            if size(ax, 1) > i_row
                ylabel(ax(subplot_position_idx(i_row, 1)), 'Activity', 'FontSize',FP.font1.size.axes_labels)
            end
        end

        % Write animal ID and session in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end

    
    %% FIGURE: Mean activity of cells for different stimuli compared to SP
    filename = os.path.join(path_mean_evoked_activity_to_SP, [animal_ID, '_mean_evoked_activity_compared_to_SP.pdf']);
    if plot_figure_mean_evoked_activity_compared_to_SP && (~exist(filename, 'file') || overwrite_figures)
        % Open figure
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 500, 1300, 800]);
        clf
        H = [];
        % Set position of each plot in the subplot grid
        n_subplot_rows = 2;  % pre and post surgery
        n_subplot_columns = max([sum(sessions < 0), sum(sessions >= 0)]);
        subplot_position = NaN(n_subplot_rows, n_subplot_columns);
        subplot_position(1, 1:sum(sessions < 0)) = sessions(sessions < 0);
        subplot_position(2, 1:sum(sessions >= 0)) = sessions(sessions >= 0);
        subplot_position_idx = reshape(1:n_subplot_rows * n_subplot_columns, n_subplot_columns, n_subplot_rows)';
        % Allocate variable for axes handles
        ax = NaN(1, n_sessions);
        % Get column corresponding to SP
        SP_column = find(ismember(stimuli_analyzed, 'SP'));
        nonSP_columns = find(~ismember(stimuli_analyzed, 'SP'));
        n_nonSP_stimuli = length(nonSP_columns);
        % Set colors
        ensemble_colors = distinguishable_colors(n_nonSP_stimuli, {'w', 'k'});

        for i_sess = 1:n_sessions
            data = cell_info(ismember(cell_info.day_from_surgery, sessions(i_sess)), :);

            % Get list of cells in ensembles
            ensembles = unique(data.ensemble);
            n_ensembles = length(ensembles);
            cells_in_ensembles = unique(data{:, 'cell'});
            n_cells_in_ensembles = length(cells_in_ensembles);

            y_ens = cell(n_ensembles, n_stimuli);
            for i_cell = 1:n_cells_in_ensembles
                % Get data of this cell
                x_rows = ismember(data{:, 'cell'}, cells_in_ensembles(i_cell)) & ismember(data{:, 'epoch'}, 'evoked');
                x = data(x_rows, {'ensemble', 'stimulus', 'activity_in'});
                ensembles_this_cell = unique(x.ensemble);
                n_ensembles_this_cell = length(ensembles_this_cell);
                y = NaN(n_stimuli, n_ensembles_this_cell);
                
                for i_ens = 1:n_ensembles_this_cell
                    for i_stim = 1:n_stimuli
                        y_rows = ismember(x{:, 'stimulus'}, stimuli_analyzed{i_stim}) & ismember(x{:, 'ensemble'}, ensembles_this_cell(i_ens));
                        this_y = x{y_rows, {'activity_in'}};
                        if this_y > 0
                            y(i_stim, i_ens) = this_y;
                        end
                    end

                    % Keep data for later
                    for i_stim = 1:n_stimuli
                        y_ens{ensembles_this_cell(i_ens), i_stim} = [y_ens{ensembles_this_cell(i_ens), i_stim}; y(i_stim)];
                    end
                end
            end

            % Plot mean cells activity
            subplot_location = subplot_position_idx(ismember(subplot_position, sessions(i_sess)));
            ax(i_sess) = subaxis(n_subplot_rows, n_subplot_columns, subplot_location, 'mb',.07, 'mt',.08, 'ml',.05, 'mr',.01, 'sv',.1);
            cla, hold on

            jitter = jitter_points(n_nonSP_stimuli, .3);
            for i_ens = 1:n_ensembles
                y = cell2mat(y_ens(i_ens, :));
                y(~isfinite(y)) = 0;
                y = y(:, nonSP_columns) - y(:, SP_column);
                avg = mean(y); err = sem(y, 1);
                
                for i_stim = 1:n_nonSP_stimuli
                    H(i_stim) = bar(i_ens + jitter(i_stim), avg(i_stim), 'FaceColor',ensemble_colors(i_stim, :), 'EdgeColor','none', 'BarWidth',.3);
                    errorbar(i_ens + jitter(i_stim), avg(i_stim), err(i_stim), 'Marker','none', 'Color','k', 'CapSize',0);
                end
            end
            
            % Set x-limits
            y_limits = get(ax(i_sess), 'YLim'); if y_limits(1) > 0, y_limits(1) = 0; end; if y_limits(2) < 0, y_limits(2) = 0; end
            set(ax(i_sess), 'XLim',[.5, n_ensembles+.5], 'YLim',y_limits, 'XTick',1:n_ensembles)
        end
        % Adjust axes appearance
        set(ax, 'TickDir','out', 'FontSize',FP.font1.size.tick_labels, 'Box','off')

        % Add title and labels
        for i_sess = 1:n_sessions
            title(ax(i_sess), ['session: ', num2str(sessions(i_sess), '%+i')], 'Interpreter','none', 'FontSize',FP.font1.size.title)
        end
        for i_row = 1:2
            if size(ax, 1) > i_row
                ylabel(ax(subplot_position_idx(i_row, 1)), 'Activity compared to SP', 'FontSize',FP.font1.size.axes_labels)
            end
        end
        
        % Add legend
        legendflex(H, stimuli_analyzed(nonSP_columns), 'FontSize',FP.font1.size.axes_labels, 'ref',fig, ...
            'anchor',[6, 6], 'buffer',[0, 0], 'nrow',1, 'box','off', ...
            'fontsize',FP.font1.size.axes_labels, 'xscale',.5, 'nolisten',true, ...
            'Interpreter','none', 'padding',[0, 8, 20]);

        % Write animal ID and session in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end
    
    
    %% FIGURE: Mean activity of cells for different stimuli compared to baseline
    filename = os.path.join(path_mean_evoked_activity_to_baseline, [animal_ID, '_mean_evoked_activity_compared_to_baseline.pdf']);
    if plot_figure_mean_evoked_activity_compared_to_baseline && (~exist(filename, 'file') || overwrite_figures)
        % Open figure
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 500, 1300, 800]);
        clf
        H = [];
        % Set position of each plot in the subplot grid
        n_subplot_rows = 2;  % pre and post surgery
        n_subplot_columns = max([sum(sessions < 0), sum(sessions >= 0)]);
        subplot_position = NaN(n_subplot_rows, n_subplot_columns);
        subplot_position(1, 1:sum(sessions < 0)) = sessions(sessions < 0);
        subplot_position(2, 1:sum(sessions >= 0)) = sessions(sessions >= 0);
        subplot_position_idx = reshape(1:n_subplot_rows * n_subplot_columns, n_subplot_columns, n_subplot_rows)';
        % Allocate variable for axes handles
        ax = NaN(1, n_sessions);
        % Get column corresponding to SP
        nonSP_columns = find(~ismember(stimuli_analyzed, 'SP'));
        n_nonSP_stimuli = length(nonSP_columns);
        % Set colors
        ensemble_colors = distinguishable_colors(n_nonSP_stimuli, {'w', 'k'});

        for i_sess = 1:n_sessions
            data = cell_info(ismember(cell_info.day_from_surgery, sessions(i_sess)), :);

            % Get list of cells in ensembles
            ensembles = unique(data.ensemble);
            n_ensembles = length(ensembles);
            cells_in_ensembles = unique(data{:, 'cell'});
            n_cells_in_ensembles = length(cells_in_ensembles);

            y_ens = cell(n_ensembles, n_stimuli);
            for i_cell = 1:n_cells_in_ensembles
                % Get difference between evoked and baseline activity of this cell
                x_ev_rows = ismember(data{:, 'cell'}, cells_in_ensembles(i_cell)) & ismember(data{:, 'epoch'}, 'evoked');
                x = data(x_ev_rows, {'ensemble', 'stimulus', 'activity_in'});
                x_bl_rows = ismember(data{:, 'cell'}, cells_in_ensembles(i_cell)) & ismember(data{:, 'epoch'}, 'baseline');
                x_bl = data(x_bl_rows, {'ensemble', 'stimulus', 'activity_in'});
                x.activity_in = x.activity_in - x_bl.activity_in;
                ensembles_this_cell = unique(x.ensemble);
                n_ensembles_this_cell = length(ensembles_this_cell);
                y = NaN(n_stimuli, n_ensembles_this_cell);
                
                for i_ens = 1:n_ensembles_this_cell
                    for i_stim = 1:n_stimuli
                        y_rows = ismember(x{:, 'stimulus'}, stimuli_analyzed{i_stim}) & ismember(x{:, 'ensemble'}, ensembles_this_cell(i_ens));
                        this_y = x{y_rows, {'activity_in'}};
                        if this_y > 0
                            y(i_stim, i_ens) = this_y;
                        end
                    end

                    % Keep data for later
                    for i_stim = 1:n_stimuli
                        y_ens{ensembles_this_cell(i_ens), i_stim} = [y_ens{ensembles_this_cell(i_ens), i_stim}; y(i_stim)];
                    end
                end
            end

            % Plot mean cells activity
            subplot_location = subplot_position_idx(ismember(subplot_position, sessions(i_sess)));
            ax(i_sess) = subaxis(n_subplot_rows, n_subplot_columns, subplot_location, 'mb',.07, 'mt',.08, 'ml',.05, 'mr',.01, 'sv',.1);
            cla, hold on

            jitter = jitter_points(n_nonSP_stimuli, .3);
            for i_ens = 1:n_ensembles
                y = cell2mat(y_ens(i_ens, :));
                y(~isfinite(y)) = 0;
                y = y(:, nonSP_columns);
                avg = mean(y); err = sem(y, 1);
                
                for i_stim = 1:n_nonSP_stimuli
                    H(i_stim) = bar(i_ens + jitter(i_stim), avg(i_stim), 'FaceColor',ensemble_colors(i_stim, :), 'EdgeColor','none', 'BarWidth',.3);
                    errorbar(i_ens + jitter(i_stim), avg(i_stim), err(i_stim), 'Marker','none', 'Color','k', 'CapSize',0);
                end
            end
            
            % Set x-limits
            y_limits = get(ax(i_sess), 'YLim'); if y_limits(1) > 0, y_limits(1) = 0; end; if y_limits(2) < 0, y_limits(2) = 0; end
            set(ax(i_sess), 'XLim',[.5, n_ensembles+.5], 'YLim',y_limits, 'XTick',1:n_ensembles)
        end
        % Adjust axes appearance
        set(ax, 'TickDir','out', 'FontSize',FP.font1.size.tick_labels, 'Box','off', 'Layer','top')

        % Add title and labels
        for i_sess = 1:n_sessions
            title(ax(i_sess), ['session: ', num2str(sessions(i_sess), '%+i')], 'Interpreter','none', 'FontSize',FP.font1.size.title)
        end
        for i_row = 1:2
            if size(ax, 1) > i_row
                ylabel(ax(subplot_position_idx(i_row, 1)), 'Activity compared to baseline', 'FontSize',FP.font1.size.axes_labels)
            end
        end
        
        % Add legend
        legendflex(H, stimuli_analyzed(nonSP_columns), 'FontSize',FP.font1.size.axes_labels, 'ref',fig, ...
            'anchor',[6, 6], 'buffer',[0, 0], 'nrow',1, 'box','off', ...
            'fontsize',FP.font1.size.axes_labels, 'xscale',.5, 'nolisten',true, ...
            'Interpreter','none', 'padding',[0, 8, 20]);

        % Write animal ID and session in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end    
    
    
    %% FIGURE: Ensemble activation per stimulus - heatmap
    filename = os.path.join(path_ensembles_activation, [animal_ID, '_ensembles_activation.pdf']);
    if plot_figure_ensembles_activation && (~exist(filename, 'file') || overwrite_figures)
        assemblies_filename = get_filename_of('assemblies', animal_ID);
        ASSEMBLIES = load_variable(assemblies_filename, 'ASSEMBLIES');
        
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 6, 2500, 200 * n_sessions]);
        clf
        ax = NaN(n_sessions, n_stimuli);
        ax_subplots_idx = reshape(1:n_stimuli * n_sessions, n_stimuli, n_sessions)';
        colormap(1 - gray)

        for i_sess = 1:n_sessions
            session_fieldname = ['cond_', num2str(i_sess)];
            stimuli_fieldname = 'stim_1';
            
            % Extract info about ensemble activations
            try
                info = ASSEMBLIES.(session_fieldname).(stimuli_fieldname).info;
            catch ME
                if strcmp(ME.identifier, 'MATLAB:structRefFromNonStruct')
                    continue  % No assemblies on this day
                else
                    rethrow(ME)
                end
            end
            state_pks = ASSEMBLIES.(session_fieldname).(stimuli_fieldname).state_pks;
            core_svd = ASSEMBLIES.(session_fieldname).(stimuli_fieldname).core_svd;
            this_n_ensembles = size(core_svd, 1);
            
            % Initialize output variable
            for i_stim = 1:n_stimuli
                % Get info about this stimulus
                this_stimulus = stimuli_analyzed{i_stim};
                stim_idx = find(ismember(info.stimuli_to_analyze, this_stimulus));
                bin_idx = find(ismember(table2cell(info.bin_edges(:, 3)), this_stimulus));
                bin_edges = info.bin_edges(bin_idx, :);
                this_state_pks = state_pks(bin_idx);
                this_timestamps_bin = info.timestamps_bin(stim_idx);
                this_n_bins_per_trial = info.n_bins_per_trial(stim_idx);

                % Average network activation per trial
                raster = zeros(this_n_ensembles, length(this_state_pks));
                for state = 1:this_n_ensembles
                    raster(state, this_state_pks == state) = 1;
                end
                if ~isnan(this_timestamps_bin)
                    avg_raster = zeros(this_n_ensembles, this_n_bins_per_trial);
                    for state = 1:this_n_ensembles
                        avg_raster(state, :) = mean(reshape(raster(state, :), this_n_bins_per_trial, []), 2);
                    end
                else
                    avg_raster = raster;
                end

                % Make time axis
                time = size(avg_raster, 2) * time_bin_activity_n_frames / 8.49;
                time = linspace(0, time, size(avg_raster, 2));
                if ~isempty(info.all_timestamps{stim_idx})
                    time = time - info.all_timestamps{stim_idx}(1);
                end

                % Draw traces
                ax(i_sess, i_stim) = subaxis(n_sessions, n_stimuli, ax_subplots_idx(i_sess, stim_idx), 'mb',.05, 'mt',.07, 'mr',.02, 'ml',.05);
                hold on
                imagesc(time, 1:this_n_ensembles, avg_raster)
                % Plot timestamps
                if ~isempty(info.all_timestamps{stim_idx})
                    for i_ts = 1:length(info.all_timestamps{stim_idx})
                        plot(info.all_timestamps{stim_idx}(i_ts) * [1, 1] - info.all_timestamps{stim_idx}(1), [0, this_n_ensembles + 1], 'color','r', 'linewidth',2, 'YLimInclude','off')
                    end
                end
                
                % Adjust axes appearance
                set(ax(i_sess, i_stim), 'XLim',[time(1), time(end)], 'TickDir','out', 'YTick',1:n_ensembles, 'YLim',[.5, this_n_ensembles + .5], 'Layer','top', 'YDir','reverse', 'Box','off')
            end
        end
        
        % Adjust axes appearance
        ax_linear = ax(:);
        ax_linear(isnan(ax_linear)) = [];
        set(ax_linear, 'FontSize',FP.font1.size.tick_labels)
        for i_stim = 1:n_stimuli
            last_row = find(~isnan(ax(:, i_stim)), 1, 'last');
            xlabel(ax(last_row, i_stim), 'Time (s)', 'FontSize',FP.font1.size.axes_labels)
            first_row = find(~isnan(ax(:, i_stim)), 1, 'first');
            title(ax(first_row, i_stim), stimuli_analyzed{i_stim}, 'Interpreter','none', 'FontSize',FP.font1.size.axes_labels)
        end
        for i_sess = 1:n_sessions
            if ~isnan(ax(i_sess, 1))
                ylabel(ax(i_sess, 1), 'Ensembles', 'Color','k', 'FontSize',FP.font1.size.axes_labels)
            end
        end
            
        % Write animal ID and session in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end
    
    
    %% FIGURE: Mean ensemble activation per stimulus
    filename = os.path.join(path_mean_ensembles_activation, [animal_ID, '_mean_ensemble_activation.pdf']);
    if plot_figure_mean_ensemble_activation && (~exist(filename, 'file') || overwrite_figures)
        % Open figure
        fig = figure('Visible',figure_visibility, 'Color','w', 'Position',[61, 500, 1300, 800]);
        clf
        H = [];
        % Set position of each plot in the subplot grid
        n_subplot_rows = 2;  % pre and post surgery
        n_subplot_columns = max([sum(sessions < 0), sum(sessions >= 0)]);
        subplot_position = NaN(n_subplot_rows, n_subplot_columns);
        subplot_position(1, 1:sum(sessions < 0)) = sessions(sessions < 0);
        subplot_position(2, 1:sum(sessions >= 0)) = sessions(sessions >= 0);
        subplot_position_idx = reshape(1:n_subplot_rows * n_subplot_columns, n_subplot_columns, n_subplot_rows)';
        % Allocate variable for axes handles
        ax = NaN(1, n_sessions);

        nonSP_stimuli = stimuli_analyzed(~ismember(stimuli_analyzed, 'SP'));
        n_nonSP_stimuli = length(nonSP_stimuli);

        for i_sess = 1:n_sessions
            data = cell_info(ismember(cell_info.day_from_surgery, sessions(i_sess)), :);

            % Get list of cells in ensembles
            ensembles = unique(data.ensemble);
            n_ensembles = length(ensembles);

            y_ens = cell(n_ensembles, n_nonSP_stimuli);
            for i_ens = 1:n_ensembles
                x = data(data.ensemble == i_ens & ~ismember(data.stimulus, 'SP'), :);
                cells_in_ensembles = unique(x{:, 'cell'});
                n_cells_in_ensembles = length(cells_in_ensembles);
                for i_cell = 1:n_cells_in_ensembles
                    for i_stim = 1:n_nonSP_stimuli
                        this_stimulus = nonSP_stimuli{i_stim};
                        bl = x{x.cell == cells_in_ensembles(i_cell) & ismember(x.stimulus, this_stimulus) & ismember(x.epoch, 'baseline'), 'n_bins_in'};
                        ev = x{x.cell == cells_in_ensembles(i_cell) & ismember(x.stimulus, this_stimulus) & ismember(x.epoch, 'evoked'), 'n_bins_in'};
                        
                        y_ens{i_ens, i_stim} = [y_ens{i_ens, i_stim}; [bl, ev]];
                    end
                end
            end
            
            subplot_location = subplot_position_idx(ismember(subplot_position, sessions(i_sess)));
            ax(i_sess) = subaxis(n_subplot_rows, n_subplot_columns, subplot_location, 'mb',.07, 'mt',.08, 'ml',.05, 'mr',.01, 'sv',.1);
            cla, hold on
            % Set colors
            ensemble_colors = distinguishable_colors(n_ensembles, {'w', 'k'});
            
            % Plot mean ensemble activation
            jitter = jitter_points(n_ensembles, .3);
            for i_stim = 1:n_nonSP_stimuli
                for i_ens = 1:n_ensembles
                    y = nanmean(y_ens{i_ens, i_stim}, 1);
                    y(~isfinite(y)) = 0;
                    % Line
                    plot([i_stim, i_stim] + jitter(i_ens), y, '-', 'Color','k')
                    % Empty marker for baseline
                    plot(i_stim + jitter(i_ens), y(1), 'o', 'MarkerFaceColor','w', 'MarkerSize',plot_marker_size, 'MarkerEdgeColor',ensemble_colors(i_ens, :), 'LineWidth',2);
                    % Full marker for evoked
                    plot(i_stim + jitter(i_ens), y(2), 'o', 'MarkerFaceColor',ensemble_colors(i_ens, :), 'MarkerSize',plot_marker_size + 2, 'MarkerEdgeColor','w');
                end
            end
            
            % Set x-limits
            y_limits = get(ax(i_sess), 'YLim'); if y_limits(1) > 0, y_limits(1) = 0; end; if y_limits(2) < 0, y_limits(2) = 0; end
            set(ax(i_sess), 'XLim',[.5, n_nonSP_stimuli + .5], 'YLim',y_limits, 'XColor','k', 'Layer','top', 'XTick',[])
        end
        % Adjust axes appearance
        set(ax, 'TickDir','out', 'FontSize',FP.font1.size.tick_labels, 'Box','off')

        % Add title and labels
        for i_sess = 1:n_sessions
            title(ax(i_sess), ['session: ', num2str(sessions(i_sess), '%+i')], 'Interpreter','none', 'FontSize',FP.font1.size.title)
            for i_stim = 1:n_nonSP_stimuli
                lims = get(ax(i_sess), 'YLim');
                lims_range = lims(2) - lims(1);
                text(i_stim, lims(1) - lims_range * .05, nonSP_stimuli{i_stim}, 'Hor','center', 'Ver','top', 'FontSize',FP.font1.size.axes_labels, 'Parent',ax(i_sess), 'Interpreter','none')
            end
        end
        for i_row = 1:2
            if size(ax, 1) > i_row
                ylabel(ax(subplot_position_idx(i_row, 1)), 'Mean ensemble activation', 'FontSize',FP.font1.size.axes_labels)
            end
        end
        
        % Write animal ID and session in figure
        ax_bg = bgAxis();
        text(0.001, 0.999, animal_ID, 'HorizontalAlignment','left', 'VerticalAlignment','top', 'Interpreter','none', 'FontSize',FP.font1.size.suptitle, 'Parent',ax_bg)

        % Print figure
        if print_figure
            export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters', fig);
            close(fig)
        end
    end
    
    
end


%% MLint exceptions
%#ok<*UNRCH,*SAGROW,*NOPTS>
