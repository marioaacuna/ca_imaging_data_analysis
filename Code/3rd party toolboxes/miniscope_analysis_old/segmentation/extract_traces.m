% This script is for extracting the full length traces after finishing the 
% manual segmentation. The output of this sctript is [output_folder, 'results_manual_segmentation.mat']
% and also the traces for each session. (Pain sessions are used for doing the time_stamps - in python)

% Write the animal_ID
cd('C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Figures_paper_FMP')
subject = 'I38_4';
% Get 'p' file
p_file = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', subject, '\Session1\preprocessed\p.mat'];
load (p_file, 'p')

% Load helper
local_saveDir = 'D:\miniscope_data\';
if ~exist(local_saveDir, 'dir'), mkdir(local_saveDir), end

root_dir_on_server = [p.rootDir subject '\Session1'];
helper_file = [local_saveDir, 'helper_', subject, '_Session1.mat'];
if ~exist(helper_file, 'file')  % Copy helper file from server
    [~, helper_file_filename] = fileparts(helper_file);
    helper_filename_on_server = [root_dir_on_server, '\preprocessed\', helper_file_filename, '.mat'];
    load(helper_filename_on_server, 'helper')

else  % Load the local one
    load(helper_file, 'helper')
end

%
animal_ID = p.subjects{1};

temp_dir = 'D:\miniscope_data\';
temp_dir_animal = [temp_dir, animal_ID, filesep];
filename_params = [temp_dir_animal, animal_ID, '_params.mat'];

load(filename_params, 'PARAMETERS')
load(PARAMETERS.filename_output, 'ROI_fluorescence')

% Normalize baseline
global GC
GC = general_configs();
Ca_indicator = 'GCaMP6f';
tau_decay = GC.decay_time_constant_calcium_reporter.(Ca_indicator);

frame_rate = p.frameRate;

% Make running-window
% Window is the maximum of 15s or 40x time decay
win_tau = 40 * tau_decay / (1000/frame_rate);
win_time = 15 * frame_rate;
win = max([win_tau, win_time]);
win = floor(win/2)*2 + 1;  % Round to nearest odd number
% Signal padders
left_pad = abs(ceil(-win/2));
right_pad = floor(win/2);

% Location of 'runline' function
cd('C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Code\3rd party toolboxes\Toolbox_Romano')

% Allocate output variable
traces = ROI_fluorescence .* 0;
n_ROIs = size(ROI_fluorescence, 2);

disp(['Normalizing baseline of ', num2str(n_ROIs), ' ROIs of mouse ', animal_ID])
for i_roi = 1:n_ROIs
    F = ROI_fluorescence(:, i_roi);
    % Get indices of each window
    all_windows = repmat(1:length(F),win,1).';
    all_windows = bsxfun(@plus, all_windows, 0:win-1);
    % Remove lines with out-of-bounds indices
    length_padded_F = length(F) + left_pad + right_pad;
    all_windows(any(all_windows.' > length_padded_F), :) = [];
    
    % Pad signal
    temp_F = [NaN(left_pad,1); F; NaN(right_pad,1)];
    % Get values in windows
    temp_F = temp_F(all_windows);
    % Compute 8th percentile in each window
    baseline = prctile(temp_F, 8, 2);
    F_smooth = runline(baseline, win, 1);
    
    % Normalize baseline
    traces(:, i_roi) = ROI_fluorescence(:, i_roi) - F_smooth;
end

% Transpose traces, so they are [ROI x time]
traces = traces.';
% filters corresponds to the ROI map
load(PARAMETERS.filename_output, 'ROIs')
filters = double(ROIs);

cellMap = max(filters, [], 3);


% Save file to disk
output_folder = ['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\jointExtraction\sorted\'];
if ~exist(output_folder, 'dir'), mkdir(output_folder), end

% Save data
filename = [output_folder, 'results_manual_segmentation.mat'];
save(filename, 'p','filters','traces')
disp(['Traces stored in ''', filename, ''''])

filename = [output_folder, 'cellMap.mat'];
save(filename, 'cellMap')
disp(['Cell map stored in ''', filename, ''''])

% Go to output folder
cd(output_folder)

%% Split traces by sessions
sessions = {'SP_1', 'pain_1','SP_2', 'pain_2', 'SP2_2'};%'SP2_1'

% do_resampling = true;
% 
% if do_resampling
%   sessions_to_rs = (sessions{4}, session{5}, session{6});
%   traces_ses

% Session 1 :  SP_1
frames_ses_1 = helper.preprocessed_frame_edges(1,:);
traces_SP_1 = traces (:,frames_ses_1(1):frames_ses_1(2));
filename = [output_folder, 'traces_' sessions{1} '.mat'];
save (filename, 'traces_SP_1')

% Session 2 : pain_1
frames_ses_2 = helper.preprocessed_frame_edges(2,:);
traces_pain_1 = traces (:,frames_ses_2(1):frames_ses_2(2));
filename = [output_folder, 'traces_' sessions{2} '.mat'];
save (filename, 'traces_pain_1')

% Session 3 : SP_2 % SP2_1
frames_ses_3 = helper.preprocessed_frame_edges(3,:);
traces_SP_2 = traces (:,frames_ses_3(1):frames_ses_3(2));
traces_SP_2 = resample_traces(traces_SP_2, 4, 3);
filename = [output_folder, 'traces_' sessions{3} '.mat'];
save (filename, 'traces_SP_2')

% Session 4 : pain_2
frames_ses_4 = helper.preprocessed_frame_edges(4,:);
traces_pain_2 = traces (:,frames_ses_4(1):frames_ses_4(2));
traces_pain_2 = resample_traces(traces_pain_2, 4, 3);
filename = [output_folder, 'traces_' sessions{4} '.mat'];
save (filename, 'traces_pain_2')

% Session 5 :pain_2
frames_ses_5 = helper.preprocessed_frame_edges(5,:);
traces_SP2_2 = traces (:,frames_ses_5(1):frames_ses_5(2));
traces_SP2_2 = resample_traces(traces_SP2_2 , 4, 3);
filename = [output_folder, 'traces_' sessions{5} '.mat'];
save (filename, 'traces_SP2_2')


% % Session 3 : SP_2 % SP2_1
% frames_ses_3 = helper.preprocessed_frame_edges(3,:);
% traces_SP2_1 = traces (:,frames_ses_3(1):frames_ses_3(2));
% %traces_SP2_1 = resample_traces(traces_SP_2, 4, 3);
% filename = [output_folder, 'traces_' sessions{3} '.mat'];
% save (filename, 'traces_SP2_1')
% 
% % Session 4 : SP_2
% frames_ses_4 = helper.preprocessed_frame_edges(4,:);
% traces_SP_2 = traces (:,frames_ses_4(1):frames_ses_4(2));
% %traces_SP_2 = resample_traces(traces_SP_2, 4, 3);
% filename = [output_folder, 'traces_' sessions{4} '.mat'];
% save (filename, 'traces_SP_2')
% 
% % Session 5 :pain_2
% frames_ses_5 = helper.preprocessed_frame_edges(5,:);
% traces_pain_2 = traces (:,frames_ses_5(1):frames_ses_5(2));
% %traces_pain_2 = resample_traces(traces_pain_2 , 4, 3);
% filename = [output_folder, 'traces_' sessions{5} '.mat'];
% save (filename, 'traces_pain_2')
% 
% % Session 6 : SP2_2
% frames_ses_6 = helper.preprocessed_frame_edges(6,:);
% traces_SP2_2 = traces (:,frames_ses_6(1):frames_ses_6(2));
% %traces_SP2_2 = resample_traces(traces_SP2_2 , 4, 3);
% filename = [output_folder, 'traces_' sessions{6} '.mat'];
% save (filename, 'traces_SP2_2')
% 




