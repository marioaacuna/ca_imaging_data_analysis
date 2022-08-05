%% Plot Proportion of neurons OPEN ARM selective
clear, clc
global GC

% parameters
do_export = 0;
PROPORTION_filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'proportions_EPM_responders.mat');
PROPORTION = load_variable(PROPORTION_filename, 'PERCENT_OPEN');
animals = fieldnames (PROPORTION); 
n_animals = length (animals);


% load data
proportion_all_animals = [];
for iid = 1:n_animals
    animal_ID = animals {iid};
    sessions = PROPORTION.(animal_ID).Properties.VariableNames;
    n_sessions = length(sessions);
    proportion_this_animal = cell(1,n_sessions);
    for i_sess = 1 : length (sessions)
        session_name = sessions{i_sess};
        proportion_this_animal{1,i_sess} = PROPORTION.(animal_ID).(session_name);
    end
    proportion_all_animals = [proportion_all_animals; proportion_this_animal ];
end
%%

mean_sess = zeros (1, n_sessions);
err_sess = mean_sess; 

for i_sess = 1 : n_sessions
    sessions = PROPORTION.(animal_ID).Properties.VariableNames;
    n_sessions = length(sessions);

    try
    mean_sess(1,i_sess) = mean(cell2mat (proportion_all_animals(:,i_sess)));
    err_sess(1,i_sess) = sem(cell2mat(proportion_all_animals(:, i_sess)));
    catch
        disp(['session ', {i_sess}, ' not analysed yet'])
    end
end


%%
close all
addpath ('C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP')
FP = figure_properties;
fig = figure('color','w', 'pos',[-0, 0, 600, 600]);
bar(mean_sess)
hold on
errorbar(mean_sess,err_sess, 'Color','k','LineStyle', 'none')
ylabel ('Proportion of Open Arm cells','FontName', FP.font1.family, 'FontSize', FP.font1.size.axes_labels)
xlabel ('Sessions', 'FontName', FP.font1.family, 'FontSize', FP.font1.size.axes_labels )
caxis ([0, 0.5])
ylim([0 , .5])

[p, h] = wilcoxon_signrank(cell2mat(proportion_all_animals(:,1)), cell2mat(proportion_all_animals(:,2)))









%% PLOT response AUC for EPM cells and responders 
clear,clc

analysis_filename = get_filename_of ('analysis_AUC');
ANALYSIS = load_variable(analysis_filename, 'AUC_ANALYSIS');
animals = fieldnames (ANALYSIS);

%% plot representative traces raw and z-scored

do_raw = 0;
data_this_cell = cell2mat(data_this_stim(80, :).');
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


figure('color','w', 'pos',[100, 160, 900, 300]);
time_axis = linspace(0,15,n_frames);
if do_raw
    plot (time_axis,data_this_cell_r.')
    hold on
    plot([5, 5], [min(data_this_cell_r(:)), max(data_this_cell_r(:))], 'Color',[.7, .7, .7], 'LineStyle','--', 'YLimInclude','on')
    
else
    plot (time_axis,data_this_cell_z.')
    hold on
    plot([5, 5], [min(data_this_cell_z(:)), max(data_this_cell_z(:))], 'Color',[.7, .7, .7], 'LineStyle','--', 'YLimInclude','on')
    
end

ylabel ('dF/F (%)')
xlabel ('time (s)')
% xticklabels (time_axis)

%% Plot representative traces for SP
clear, clc
do_export_figure = 0;
cells_to_take = 3:5;
type_traces = 'dFF';
% type_traces = 'spikes';
global GC
% global GC
GC.experiment_name = 'ACC_SNI_anxiety';
% load parametrers
animal_ID = 'MA_30epi';
session_to_take = 1; % of EPM sessions
n_cells = SQL_database.read_table_where('experiments',{'n_ROIs'}, animal_ID, 'animal_ID', 'return_as_table', false);
METADATA = SQL_database.read_epi_trials(animal_ID);
frame_rate = GC.epifluorescence_downsample_to_frame_rate;
dates_EPM = METADATA.date(ismemberCellRows(METADATA.experiment, {'EPM'}), :);
date = dates_EPM{session_to_take};
selected_session = date;


% load SP traces
disp (['loading spikes and traces session ', date, ' ',animal_ID])


filename_traces = get_filename_of ('dFF', animal_ID);
filename_spikes = get_filename_of ('spikes', animal_ID);
traces = load_variable (filename_traces, 'dFF');
spikes = load_variable (filename_spikes, 'spikes');
decay_time = GC.decay_time_constant_calcium_reporter;
spikes = convolve_spikes(spikes, frame_rate, decay_time.GCaMP6f);

trial_idx = ismember (METADATA.date, date) & ismember(METADATA.stimulus, 'SP');
traces_SP = cell2mat(traces (:, trial_idx));
spikes_SP = cell2mat(spikes (:, trial_idx));
n_frames = size (traces_SP,2);clear
% plot figure

time_axis = linspace(0,600,n_frames);
%%
toggle_toolbox ('_plotting', 'on')
addpath('C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP')
FP = figure_properties;

mb = .1;
mt = .01;
ax = [];

fig = figure('color','w', 'pos',[0, 0, 900, 500]);
plot (time_axis(1:2001), traces_SP(4, 1000:3000)'), hold on
plot (time_axis(1:2001), spikes_SP(4, 1000:3000)')
axis off
% Add scalebar
xlimits = xlim; ylimits = ylim;%ylim(ax(i_stim));
xrange = diff(xlimits); yrange = diff(ylimits);
% ax3 = axes('pos', 'XLim',xlimits, 'YLim',ylimits, 'Visible','off', 'Clipping','off', 'HitTest','off');
h_scalebar = scalebar( 'Position',[xlimits(1)-xrange*.1, ylimits(1)], 'XUnit','s', 'YUnit','{% dF/F}', 'XLen',20, 'YLen',1);
set(h_scalebar.hTextX, 'FontSize',FP.font1.size.axes_labels, 'FontName',FP.font1.family)
set(h_scalebar.hTextY, 'FontSize',FP.font1.size.axes_labels, 'FontName',FP.font1.family)

if do_export_figure
    filename = [os.path.join(GC.data_root_path, GC.plots_path, GC.experiment_name), '\material_for_figures\example_SP_traces_zoom_(', animal_ID, '_session', num2str(selected_session, '%+i') ').pdf'];
    export_fig(filename, '-pdf', '-q101', '-nocrop', '-painters',fig); %close(fig)
end

toggle_toolbox ('_plotting', 'off')
