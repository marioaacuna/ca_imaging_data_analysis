% thi script uses the 'results_manual_segmentation.mat', and deconvolve
% full-length traces.

cd('M:\Mario\MATLAB_scripts\Others\MovieAnalysis_B_Grewe\MovieAnalysis\segmentation')
load (['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\Session1\preprocessed\p.mat'])


% animal_ID = subject;
name_session = 'pain_1';
frame_rate = p.frameRate;


filename = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\jointExtraction\sorted\results_manual_segmentation.mat'];
load(filename, 'traces')

output_folder = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\',animal_ID,'\jointExtraction\sorted\'];

%%
global GC
GC = general_configs();

toggle_toolbox('OASIS_matlab', 'on')

n_cells = size(traces, 1);
spikes = traces .* 0;
for i_cell = 1:n_cells
    disp(['Cell ', num2str(i_cell)])
    [spikes(i_cell, :), ~] = deconvolve_spikes_cell(animal_ID, double(traces(i_cell, :)), frame_rate);
end

toggle_toolbox('OASIS_matlab', 'off')

%% save spikes
sessions = {'SP_1', 'pain_1', 'SP2_1', 'SP_2', 'pain_2', 'SP2_2'};

p_file = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\Session1\preprocessed\p.mat'];
load (p_file, 'p')
root_dir_on_server = [p.rootDir animal_ID '\Session1'];
% Load helper
local_saveDir = 'D:\miniscope_data\';
if ~exist(local_saveDir, 'dir'), mkdir(local_saveDir), end
helper_file = [local_saveDir, 'helper_', animal_ID, '_Session1.mat'];
if ~exist(helper_file, 'file')  % Copy helper file from server
    [~, helper_file_filename] = fileparts(helper_file);
    helper_filename_on_server = [root_dir_on_server, '\preprocessed\', helper_file_filename, '.mat'];
    load(helper_filename_on_server, 'helper')
else  % Load the local one
    load(helper_file, 'helper')
end

% Session 1 :  SP_1
frames_ses_1 = helper.preprocessed_frame_edges(1,:);
spikes_SP_1 = spikes (:,frames_ses_1(1):frames_ses_1(2));
filename = [output_folder, 'spikes_' sessions{1} '.mat'];
save (filename, 'spikes_SP_1')

% Session 2 : pain_1
frames_ses_2 = helper.preprocessed_frame_edges(2,:);
spikes_pain_1 = spikes (:,frames_ses_2(1):frames_ses_2(2));
filename = [output_folder, 'spikes_' sessions{2} '.mat'];
save (filename, 'spikes_pain_1')

% Session 3 : SP2_1
frames_ses_3 = helper.preprocessed_frame_edges(3,:);
spikes_SP2_1 = spikes (:,frames_ses_3(1):frames_ses_3(2));
filename = [output_folder, 'spikes_' sessions{3} '.mat'];
save (filename, 'spikes_SP2_1')

% Session 4 : SP_2
frames_ses_4 = helper.preprocessed_frame_edges(4,:);
spikes_SP_2 = spikes (:,frames_ses_4(1):frames_ses_4(2));
filename = [output_folder, 'spikes_' sessions{4} '.mat'];
save (filename, 'spikes_SP_2')

% Session 5 : pain_2
frames_ses_5 = helper.preprocessed_frame_edges(5,:);
spikes_pain_2 = spikes (:,frames_ses_5(1):frames_ses_5(2));
filename = [output_folder, 'spikes_' sessions{5} '.mat'];
save (filename, 'spikes_pain_2')

% Session 6 : SP2_2
frames_ses_6 = helper.preprocessed_frame_edges(6,:);
spikes_SP2_2 = spikes (:,frames_ses_6(1):frames_ses_6(2));
filename = [output_folder, 'spikes_' sessions{6} '.mat'];
save (filename, 'spikes_SP2_2')



%% plot examples

% figure('color','w')
% clf, hold on,
% plot(traces(2,:),'k'),
% fill([1:length(traces(2,:)), fliplr(1:length(traces(2,:)))], [spikes(2,:), zeros(1,length(s))],'b','EdgeColor','none')


