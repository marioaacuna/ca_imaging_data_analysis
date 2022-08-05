clc
clear,

global GC
experiment_name = 'ACC_SNI_anxiety';
addpath  C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP
toggle_toolbox('_plotting', 'on')
FP = figure_properties();

do_export_figures = 0;

session_to_take = 1;
group_trials_by = GC.analysis.(experiment_name).group_trials_by_timepoint;
baseline_window = GC.detection_baseline_window;

activity_type = 'dFF';
% type_threshold = 'p_val';

do_print_figure = false;
do_separate_no_activated = 1; 


animal_ID = 'MA_30epi'; cells = [];
stimuli_to_show = {'heat', 'cold', 'pinprick', 'touch'};
stimuli_labels = {'heat', 'cold', 'pinprick', 'touch'};
activity_filename = get_filename_of(activity_type, animal_ID);
activity = load_variable(activity_filename, activity_type);
METADATA = SQL_database.read_epi_trials(animal_ID);

sessions = unique(METADATA.day_from_surgery, 'stable');
selected_session = sessions{session_to_take,:};
selected_session =  make_variable_name(['session_',selected_session ]);
result_filename = get_filename_of('response_detection_epi', animal_ID);
RESULTS = load_variable(result_filename, 'RESULTS');
session_fieldnames = fieldnames( RESULTS);
session_idx = find(ismember(session_fieldnames, selected_session));
n_cells = SQL_database.read_table_where('experiments', 'n_ROIs', animal_ID,'animal_ID', 'return_as_table',false);
cell_colors = [13, 122, 191 ; 255, 192, 0; 255, 0, 0; 255, 0, 255] ./ 255;
modulations_colors = [FP.colors.groups.CCI; FP.colors.groups.sham; .7, .7, .7];
do_downsample = 1;
baseline_window = GC.detection_baseline_window;
stimuls_profile = zeros(76,1);
stimuls_profile(25:end,:) =  1;


frame_rate = GC.epifluorescence_downsample_to_frame_rate;
fig = figure('color','w', 'pos',[100, 160, 900, 300]);

%%
clc
YLIMITS = repmat([0, 0.6], 4, 1);
subplot_n_rows = 1;
subplot_n_cols = 4;
mb = .1;
mt = .01;

clf, ax = [];
hold on
for i_stim = 1:length(stimuli_to_show)
    selected_stimulus = stimuli_to_show{i_stim};
    this_detection_results = RESULTS.(session_fieldnames{session_idx}).(selected_stimulus).selectivity;
    
    trials_idx = ismember(METADATA.stimulus, selected_stimulus) & ismember(METADATA{:, group_trials_by}, num2str( sessions{session_to_take,:}, '%+i'));
    METADATA_this_stim = METADATA(trials_idx,:);
    if ~any(trials_idx), continue, end
    data = activity(:, trials_idx);
    n_frames = cellfun(@length, data);
    n_frames = unique(n_frames(:));
    time_axis = linspace(0,15,n_frames);
    
    data_to_plot = NaN(n_cells, n_frames);

    % calculate z_score for each cell
    for i_cell = 1:n_cells
        data_this_cell = cell2mat(data(i_cell, :).');
        data_this_cell_r = zeros (size(data_this_cell));
        for i_trial = 1 : size(data_this_cell,1)
            data_this_cell_r (i_trial,:) =  data_this_cell(i_trial,:) - median(data_this_cell(i_trial,:));
        end
        trial_duration = GC.miniscope_trial_duration.(selected_stimulus);
        baseline_timestamp = 1: (frame_rate * trial_duration(1));
        evoked = (1 + (frame_rate * trial_duration (1)))  : ((trial_duration(2) + trial_duration(1))  * frame_rate);
        
        mean_baseline = mean(data_this_cell_r(:, baseline_timestamp), 2);
        std_baseline = std(data_this_cell_r(:, baseline_timestamp, :), 0, 2);
        data_this_cell_z = (data_this_cell_r - mean_baseline) ./ std_baseline;
        
        data_this_cell_z = mean(data_this_cell_z, 1);
        data_to_plot(i_cell, :) = data_this_cell_z;
    end
    
    data_to_plot(~isfinite(data_to_plot)) = 0;
    if do_downsample > 0
        data_to_show = data_to_plot;
        for i_cell = 1:n_cells
            data_to_show(i_cell, :) = smooth(data_to_plot(i_cell, :), do_downsample);
        end
    else
        data_to_show = data_to_plot;
    end
    ax(i_stim) = subaxis(subplot_n_rows, subplot_n_cols, i_stim, 'ml',.07, 'mt',mt, 'mr',.01, 'mb',mb, 'sh',.05);
    cla, hold on

    plot([-10, 100], [0, 0], 'Color',[.7, .7, .7], 'XLimInclude','off', 'Linestyle','--')
    plot([5, 5], [-0.1, 1], 'Color','k', 'LineStyle','--', 'YLimInclude','on')
%   
    if do_separate_no_activated
        data_z = data_to_show(this_detection_results.excited == 0, :); mod_idx = 3;
        if ~isempty(data_z)
            data_avg = smooth(mean(data_z, 1));
            data_err = repmat(sem(data_z, 1), 2, 1);
            [h_line, h_patch] = boundedline(time_axis', data_avg', data_err');
            set(h_patch, 'FaceAlpha',.1)
            set(h_line, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
            set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
        end
    end
    
    if do_separate_no_activated
        data_z = data_to_show(this_detection_results.inhibited == 1, :); mod_idx = 1;
        if ~isempty(data_z)
            data_avg = smooth(mean(data_z, 1));
            data_err = repmat(sem(data_z, 1), 2, 1);
            [h_line, h_patch] = boundedline(time_axis', data_avg', data_err');
            set(h_patch, 'FaceAlpha',.1)
            set(h_line, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
            set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
        end
    end

    
    if do_separate_no_activated, mod_idx = 2; else mod_idx = 1; end
    data_z = data_to_show(this_detection_results.excited == 1, :); %(this_detection_results == 1, :)
    if ~isempty(data_z)
        data_avg = smooth(mean(data_z, 1));
        data_err = repmat(sem(data_z, 1), 2, 1);
        [h_line, h_patch] = boundedline(time_axis', data_avg', data_err');
        set(h_patch, 'FaceAlpha',.1)
        set(h_line, 'color',tint(modulations_colors(mod_idx, :), 0), 'linewidth',2)
        set(h_patch, 'FaceColor',shade(modulations_colors(mod_idx, :), .75))
    end
    set(ax(i_stim), 'Fontsize',FP.font1.size.tick_labels, 'Layer','top', 'Box','off', 'TickDir','out')
    axis off

    % Add scalebar
    xlimits = xlim(ax(i_stim)); ylimits = ylim(ax(i_stim));
    xrange = diff(xlimits); yrange = diff(ylimits);
    ax3 = axes('pos', getAxisPos(ax(i_stim)), 'XLim',xlimits, 'YLim',ylimits, 'Visible','off', 'Clipping','off', 'HitTest','off');
    h_scalebar = scalebar(ax3, 'Position',[xlimits(1)-xrange*.1, ylimits(1)], 'XUnit','s', 'YUnit','{z-score}', 'XLen',5, 'YLen',1); %{\ita.u.}
    set(h_scalebar.hTextX, 'FontSize',FP.font1.size.axes_labels, 'FontName',FP.font1.family)
    set(h_scalebar.hTextY, 'FontSize',FP.font1.size.axes_labels, 'FontName',FP.font1.family)
    disp (sum(this_detection_results.inhibited))

end

if do_export_figures
    filename = [os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name), '\response_averages','\response_average_(', animal_ID, '_session', num2str(selected_session, '%+i') ').pdf'];
    export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters',fig);% close(fig)
end

    
toggle_toolbox('_plotting', 'off')


