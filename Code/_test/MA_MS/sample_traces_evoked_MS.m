%%%% plot shaded traces
clc, clear
global GC

animal_ID = 'MA_29epi';
METADATA = SQL_database.read_epi_trials(animal_ID);
sessions = (unique(METADATA.day_from_surgery));
i_sess = 3;
name_session = cell2mat(sessions(i_sess));
idx_sess = find(ismemberCellRows(METADATA.day_from_surgery, {name_session}));
METADATA_this_session = METADATA(idx_sess,:);

session_name = make_variable_name(['session_', GC.session_column_prefix]);

filename_traces = get_filename_of('dFF', animal_ID);
all_traces_this_animal = load_variable (filename_traces, 'dFF');

missing_only = 1;
result_filename = get_filename_of('response_detection_epi', animal_ID, missing_only);
RESULTS = load_variable(result_filename, 'RESULTS');
cells_to_show = [1,20, 49];
n_cells = size(cells_to_show,2);
cell_colors = [13, 122, 191 ; 255, 192, 0; 255, 0, 0; 255, 0, 255] ./ 255;
stimuli_to_show = {'cold', 'heat', 'pinprick', 'touch'};
% err_trace = std_trace_example ./ (sqrt(n_cells));

%plot the figure
toggle_toolbox('_plotting', 'on')
figure ('color', 'w')
ylabel ('Amplitude (dF/F)')
xlabel ('Time (s)')
title('example cells')
plot([5, 5], [-100, 1000], 'Color',[.7, .7, .7], 'LineStyle','--', 'YLimInclude','off')
hold on
for selected_stim_idx = 1:length (stimuli_to_show)
    selected_stim = stimuli_to_show{selected_stim_idx};
    selected_stim = make_variable_name(selected_stim);
    mean_trace_examples = RESULTS.(session_name).(selected_stim).mean_traces;
    std_trace_examples  = RESULTS.(session_name).(selected_stim).std_traces;
    stim_idx = find(ismember(METADATA.stimulus, {strcat(selected_stim)}) & ismember (METADATA.day_from_surgery, {name_session}));



    for i_cell = 1%: n_cells
        this_cell_trace = mean_trace_examples{cells_to_show(i_cell),:};
        this_cell_std = std_trace_examples{cells_to_show(i_cell),:};
        n_trials = size(stim_idx,1);
        error_this_cell  = (this_cell_std ./ sqrt(n_trials));
        x = linspace (0,15,size(this_cell_trace,2));
%         err_trace = this_cell_std ./ (sqrt(size(z_data, 1)));
        [l, p] = boundedline (x', this_cell_trace', error_this_cell');
        ax= subplot(1,1,1);
        outlinebounds(l,p, {'DisplayName',selected_stim} );
       
    end



    set(p, 'FaceAlpha',.05)
    set(l, 'color',tint(cell_colors(selected_stim_idx, :), 0), 'linewidth',2)
    set(p, 'FaceColor',shade(cell_colors(selected_stim_idx, :), .75))
    legend (l, stimuli_to_show{selected_stim_idx}) 


      
end


