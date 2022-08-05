clc
clear
close all

cd ('C:\Users\acuna\Documents\Two_photon_imaging_data_analysis')
addpath C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP

global GC
GC.experiment_name = 'ACC_SNI_anxiety';

toggle_toolbox('_plotting', 'on')
FP = figure_properties();
do_export_figure = 0;
session_to_take = 1;
activity_type = 'dFF';
do_print_figure = false;
animal_ID = 'MA_30epi'; cells = [];
stimuli_to_show = {'heat', 'cold', 'pinprick', 'touch'};
stimuli_labels = {'heat', 'cold', 'pinprick', 'touch'};
activity_filename = get_filename_of(activity_type, animal_ID);
activity = load_variable(activity_filename, activity_type);
METADATA = SQL_database.read_epi_trials(animal_ID);
sessions = (unique(METADATA.day_from_surgery, 'stable'));
selected_session = sessions{session_to_take,:};
baseline_window = GC.detection_baseline_window;
result_filename = get_filename_of('response_detection_epi', animal_ID);
RESULTS = load_variable(result_filename, 'RESULTS');
cells_to_show = [1,20, 49];
n_cells = size(cells_to_show,2);
cell_colors = [13, 122, 191 ; 255, 192, 0; 255, 0, 0; 255, 0, 255] ./ 255;

frame_rate = GC.epifluorescence_downsample_to_frame_rate;

%%
fig =figure('color','w', 'pos',[-1679, 160, 900, 973]);
do_downsample  = 4;
subplot_n_rows = 1;
subplot_n_cols = length (stimuli_to_show);
stimuls_profile = zeros(76,1);
stimuls_profile(25:end,:) =  1;
clf, ax = [];
did_sort_cells = false;
MIN_COLOR = []; MAX_COLOR = []; ZERO_COLOR = [];

%%
for selected_stimulus_idx = 1:length(stimuli_to_show)
    selected_stimulus = stimuli_to_show{selected_stimulus_idx};

    trials_idx =  find(ismember(METADATA.stimulus, {strcat(selected_stimulus)}) & ismember (METADATA.day_from_surgery, (selected_session)));
    METADATA_this_session = METADATA(trials_idx,:);
    if ~any(trials_idx), continue, end
    data = activity(:, trials_idx);
    
    n_cells = size(data, 1);
    n_frames = cellfun(@length, data);
    n_frames = unique(n_frames(:));
    data_to_plot = NaN(n_cells, n_frames);

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
    if ~did_sort_cells
        sum_activity = data_to_plot * stimuls_profile;
        [~, order] = sort(sum_activity, 'descend');
        cells_sorted_position = NaN(length(n_cells), 1);
        for i_cell = 1:length(cells)
            cells_sorted_position(i_cell) = find(order == n_cells(i_cell));
        end
        disp(cells_sorted_position)
        did_sort_cells = true;
    end
    data_to_plot = data_to_plot(order, :);
    
    if do_downsample > 0
        data_to_show = data_to_plot;
        for i_cell = 1:n_cells
            data_to_show(i_cell, :) = smooth(data_to_plot(i_cell, :), do_downsample);
        end        
    else
        data_to_show = data_to_plot;
    end

    %     data_to_show = data_to_plot;
    
    ax(end+1) = subaxis(subplot_n_rows, subplot_n_cols, selected_stimulus_idx, 'ml',.05, 'mt',.05, 'mb',.05, 'sh',.01, 'mr',.13);
    time_axis = linspace (0,15,n_frames);
    imagesc(data_to_show, 'Xdata',time_axis([1, end]))

    min_color = min(data_to_show(:));
    max_color = quantile(data_to_show(:), .99);
    zero_color = 0;
    MIN_COLOR = [MIN_COLOR, min_color];
    MAX_COLOR = [MAX_COLOR, max_color];
    ZERO_COLOR = [zero_color, zero_color];
    hold on
    plot([5, 5], [-1, n_cells+1], 'Color','k', 'LineStyle','--', 'YLimInclude','on')    
    if ax(selected_stimulus_idx) >1
        set(gca, 'ytick',[]) 
    end
end

MIN_COLOR = min(MIN_COLOR);
MAX_COLOR = max(MAX_COLOR);
ZERO_COLOR = mean(ZERO_COLOR);
colormap(divergent_colormap(MIN_COLOR, ZERO_COLOR, MAX_COLOR, 0, 1))
for iax = 1:length(ax)
    caxis(ax(iax), [MIN_COLOR, MAX_COLOR])
end
set(ax, 'Fontsize',FP.font1.size.tick_labels, 'Layer','top')

pos = getAxisPos(ax(1));
bgAxis([0, pos(2), .97, pos(4)])
hc = colorbar;
caxis([MIN_COLOR, MAX_COLOR])
set(hc.Label, 'String', 'Activity amplitude ({\itz}-score)', 'FontSize',FP.font1.size.axes_labels) % {\itz}-score, a.u., dF/F


xlabel(ax(1), 'Time (s)', 'Fontsize',FP.font1.size.axes_labels)
ylabel(ax(1), 'Cell ID', 'Fontsize',FP.font1.size.axes_labels)

for iax = 1:length(ax)
    title(ax(iax), stimuli_labels{iax}, 'Interpreter','none', 'FontSize',FP.font1.size.title)
end

if do_export_figure
    filename = [os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name), '\response_heatmaps\response_heatmap_(', animal_ID, '_session', num2str(selected_session, '%+i') ').pdf'];
    export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters',fig); close(fig)
end
    
toggle_toolbox('_plotting', 'off')

    

