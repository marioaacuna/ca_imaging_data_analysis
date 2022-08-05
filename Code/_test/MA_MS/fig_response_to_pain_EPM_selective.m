clear, clc

do_export_figure = 0;
global GC
GC.experiment_name = 'ACC_SNI_anxiety';
% load parametrers
animal_ID = 'MA_30epi';
session_to_take = 2;
n_cells = SQL_database.read_table_where('experiments',{'n_ROIs'}, animal_ID, 'animal_ID', 'return_as_table', false);
METADATA = SQL_database.read_epi_trials(animal_ID);
frame_rate = GC.epifluorescence_downsample_to_frame_rate;
dates_EPM = METADATA.date(ismemberCellRows(METADATA.experiment, {'EPM'}), :);
date = dates_EPM{session_to_take};

order_stimuli = {'heat', 'cold', 'pinprick', 'touch'};
% load traces
filename_traces = get_filename_of('dFF', animal_ID);
all_traces_this_animal = load_variable (filename_traces, 'dFF');
trials_idx = find(ismemberCellRows(METADATA.date, {date}) & ismemberCellRows(METADATA.type(:), {'evoked'})) ;
data = all_traces_this_animal(:, trials_idx);

% load session
day_from_surgery = cell2mat(unique (METADATA.day_from_surgery(trials_idx,: )));
selected_session =  make_variable_name(['session_',day_from_surgery ]);

% Load  pain data
pain_result_filename = get_filename_of('response_detection_epi', animal_ID);
RESULTS = load_variable(pain_result_filename, 'RESULTS');
session_fieldnames = fieldnames(RESULTS);
session_idx = find(ismember(session_fieldnames, selected_session));
stimuli = unique(METADATA.stimulus(trials_idx,:));
n_stim = length(stimuli);
excited_pain = cell(1, n_stim);

for i_stim = 1 : n_stim
    excited_pain{:, i_stim} = RESULTS.(selected_session).(stimuli{i_stim}).selectivity.excited;
end

% load EPM data
filename_EPM_data = os.path.join(GC.data_root_path, GC.aggregated_data_path,GC.experiment_name, 'EPM_results',animal_ID,['epm_idx_cells_open_',date,'_',animal_ID,'.json']);
txt = jsondecode(fileread(filename_EPM_data));
cells = txt.idxCellsRespondingInOpen;
cells_idx = fieldnames(cells);
n_EPM_cells = length (cells_idx);
cells_id_EPM = cell (n_EPM_cells,1);
for i_cell = 1:n_EPM_cells
    cells_id_EPM {i_cell, :}= str2num(cell2mat( strsplit((cells.(cell2mat(cells_idx(i_cell)))), 'cell_') ));  
end
%% get the intescet cells 
idx_cells_intersect = cell (n_stim,1);
for i_stim = 1 : n_stim
    idx_cells_intersect{i_stim} = intersect(find(excited_pain{i_stim}), cell2mat(cells_id_EPM(:,1)));
end



%% plot
addpath  C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP
toggle_toolbox ('_plotting', 'on')
FP = figure_properties();
modulations_colors = [FP.colors.groups.CCI; FP.colors.groups.sham; .7, .7, .7; .2, .2, .2];

YLIMITS = repmat([-0.5, 6], 4, 1);
subplot_n_rows = 1;
subplot_n_cols = n_stim;
mb = .1;
mt = .01;
fig = figure('color','k', 'pos',[100, 160, 900, 300]);
clf, ax = [];
hold on

for i_stim = 1: n_stim
    cells_idx_this_stim = excited_pain{i_stim};
    cells_intersect_this_stim = cell2mat(idx_cells_intersect(i_stim));
    name_stim = stimuli{i_stim};
    n_int_cells = length(cells_intersect_this_stim);
    trials_idx = ismember(METADATA.stimulus, stimuli(i_stim)) & ismember(METADATA.date, date);
    METADATA_this_stim = METADATA(trials_idx,:);
    data_this_stim = all_traces_this_animal(:, trials_idx);
    n_frames = cellfun(@length, data_this_stim);
    n_frames = unique(n_frames(:));
    time_axis = linspace(0,15,n_frames);
    data_to_plot = NaN(n_cells, n_frames);
    for i_cell = 1:n_cells
        data_this_cell = cell2mat(data_this_stim(i_cell, :).');
        data_this_cell_r = zeros (size(data_this_cell));
        for i_trial = 1 : size(data_this_cell,1)
            data_this_cell_r (i_trial,:) =  data_this_cell(i_trial,:) - median(data_this_cell(i_trial,:));
        end
        
        trial_duration = GC.miniscope_trial_duration.(name_stim);
        baseline_timestamp = 1: (frame_rate * trial_duration(1));
        evoked = (1 + (frame_rate * trial_duration (1)))  : ((trial_duration(2) + trial_duration(1))  * frame_rate);
        
        mean_baseline = mean(data_this_cell_r(:, baseline_timestamp), 2);
        std_baseline = std(data_this_cell_r(:, baseline_timestamp, :), 0, 2);
        data_this_cell_z = (data_this_cell_r - mean_baseline) ./ std_baseline;
        
        data_this_cell_z = mean(data_this_cell_z, 1);
        data_to_plot(i_cell, :) = data_this_cell_z;
%         data_to_plot(i_cell, :) = mean(data_this_cell_r, 1);


    end
    data_to_show = data_to_plot;
 
    ax(i_stim) = subaxis(subplot_n_rows, subplot_n_cols, find(ismember(order_stimuli, name_stim)), 'ml',.07, 'mt',mt, 'mr',.01, 'mb',mb, 'sh',.05);
    cla, hold on
    plot([-10, 100], [0, 0], 'Color',[.7, .7, .7], 'XLimInclude','off', 'Linestyle','--')
    hold on
    plot([5, 5], [-0.1, 5], 'Color',[.7, .7, .7], 'LineStyle','--', 'YLimInclude','on')
%     legend ((stimuli{i_stim}))
%     legend boxoff
%   
    y = cell(1,4);
    data_z = data_to_show(cells_idx_this_stim == 1, :); mod_idx = 2; % blue, selective
    if ~isempty(data_z)
        data_avg = smooth(mean(data_z, 1));
        data_err = repmat(sem(data_z, 1), 2, 1);
        [h_line_1, h_patch] = boundedline(time_axis', data_avg', data_err');
        set(h_patch, 'FaceAlpha',.1)
        set(h_line_1, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
        set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
        set(get(get(h_patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    y{1} = max (mean(data_z,1)); % to select max ylim
    
    cell_int_logic = zeros(n_cells,1);
    cell_int_logic(cells_intersect_this_stim) = 1;
    
    data_z = data_to_show(cell_int_logic == 1, :); mod_idx = 1; % red interserct (this_detection_results == 1, :)
    if ~isempty(data_z)
        data_avg = smooth(mean(data_z, 1));
        data_err = repmat(sem(data_z, 1), 2, 1);
        [h_line_2, h_patch] = boundedline(time_axis', data_avg', data_err');
        set(h_patch, 'FaceAlpha',.1)
        set(h_line_2, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
        set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
        set(get(get(h_patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    
    y{2} = max (mean(data_z,1));
    data_z = data_to_show(cells_idx_this_stim == 0, :); mod_idx = 3; % grey, non sel(this_detection_results == 1, :)
    if ~isempty(data_z)
        data_avg = smooth(mean(data_z, 1));
        data_err = repmat(sem(data_z, 1), 2, 1);
        [h_line_3, h_patch] = boundedline(time_axis', data_avg', data_err');
        set(h_patch, 'FaceAlpha',.1)
        set(h_line_3, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
        set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
        set(get(get(h_patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    y{3} = max (mean(data_z,1));

    end
    
    cell_EPM_logic = zeros(n_cells,1);
    cell_EPM_logic(cell2mat(cells_id_EPM(:,1))) = 1;

    data_z = data_to_show(cell_EPM_logic == 1, :); mod_idx = 4; % black, EPM (this_detection_results == 1, :)
    if ~isempty(data_z)
        data_avg = smooth(mean(data_z, 1));
        data_err = repmat(sem(data_z, 1), 2, 1);
        [h_line_4, h_patch] = boundedline(time_axis', data_avg', data_err');
        set(h_patch, 'FaceAlpha',.1)
        set(h_line_4, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
        set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
        set(get(get(h_patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    end
    y{4} = max (mean(data_z,1));
    
    q = quantile(cell2mat(y),4);
    set(ax(i_stim), 'Fontsize',FP.font1.size.tick_labels, 'Layer','top', 'Box','off', 'TickDir','out')
    xlim(ax(i_stim), [0, 15])
    ylim(ax(i_stim), [-0.5, max(cell2mat(y)) + q(3)]);%YLIMITS(i_stim, :))
    axis off

%    Add scalebar
    xlimits = xlim(ax(i_stim)); ylimits = YLIMITS(i_stim, :);%ylim(ax(i_stim));
    xrange = diff(xlimits); yrange = diff(ylimits);
    ax3 = axes('pos', getAxisPos(ax(i_stim)), 'XLim',xlimits, 'YLim',ylimits, 'Visible','off', 'Clipping','off', 'HitTest','off');
    h_scalebar = scalebar(ax3, 'Position',[xlimits(1)-xrange*.1, ylimits(1)], 'XUnit','s', 'YUnit','{\itz-score}', 'XLen',5, 'YLen',1);
    set(h_scalebar.hTextX, 'FontSize',FP.font1.size.axes_labels, 'FontName',FP.font1.family)
    set(h_scalebar.hTextY, 'FontSize',FP.font1.size.axes_labels, 'FontName',FP.font1.family)

end

if do_export_figure
    filename = [os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name), '\response_averages\response_average_pain_EPM_(', animal_ID, '_session', num2str(selected_session, '%+i') ').pdf'];
    export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters',fig); %close(fig)
end
    

toggle_toolbox ('_plotting', 'off')