
clear, clc
global GC


animal_ID = 'MA_35epi';
date_recording = '210528';
experiment = 'pain';
csv_filepath = 'T:\Mario\Behavior\ACC_pain\behavior cams\_stimuli_timestaps_csv_to_database\csvs';
table_filepath = 'T:\Mario\Behavior\ACC_pain\behavior cams\_stimuli_timestaps_csv_to_database\originals';

video_folder = os.path.join('T:\Mario\Behavior\', 'ACC_pain', GC.behavior_video_subfolder, animal_ID, date_recording, experiment);
filename = os.path.join(video_folder, 'stimuli_timestamps.csv');



original_table = readtable([table_filepath, '\',[date_recording, '_', animal_ID, '.xlsx']]);

% set parametres
time_start = num2cell(original_table.min .* 60 + original_table.sec + original_table.msec .*0.01);
time_end =  num2cell(NaN(length(time_start),1));
event = original_table.Event;
withdraw = original_table.Withdrawal;
aff = original_table.AffectiveResponse;

% Re set values to match data on database
withdraw = num2cell(strcmp(withdraw, 'yes'));
new_aff = aff;
new_aff(strcmp(aff, 'no')) = {''};
new_aff(strcmp(aff, 'other')) = {''};
% Outputs
var_names = {'Time_start', 'Time_end', 'Event', 'Withdraw', 'Affective_response'};
output_array = [time_start, time_end, event, withdraw, new_aff];
output_table = array2table(output_array,'VariableNames', var_names);

% Check if last value is nan
has_nan = isnan(cell2mat(output_table.Time_start));
if sum(has_nan) > 0
    output_table(has_nan,:)=[];
end

% Write csv
name_csv = [animal_ID, '_', date_recording, '_pain.csv'];
csv_filename = [csv_filepath,'\', name_csv];
writetable(output_table,csv_filename)
disp('done csv')
%% Copy files to original forlder to be next updated using C:\Users\acuna\Documents\Two_photon_imaging_data_analysis\Code\Superuser\import_behavior_timestamps_to_database.m

copyfile(csv_filename, filename)
disp('done copying')

