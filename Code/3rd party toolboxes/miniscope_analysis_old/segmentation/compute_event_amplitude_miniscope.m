clc
global GC
GC = general_configs();
%%
animals = {'I32_1', 'I32_3','I37_2','I38_3', 'I38_4', 'I39_1'};
n_animals = length(animals);
experimental_conditions = {'sham', 'cci', 'sham', 'sham', 'cci', 'cci'};
condition = 'SP'; % or 'SP2' (SP is the first SP of the sessioin, SP2 is after the pain protocol)
sessions = {'1', '2'};
n_sessions = length(sessions);
experiemental_group =[1,2,1,1,2,2];
experimental_gropus = unique(experimental_conditions);
experimental_groups = fliplr (experimental_gropus);
AMPLITUDE = cell(n_animals, n_sessions);
FREQUENCY = cell(n_animals, n_sessions);
TIMEPOINTS = cell(n_animals, 3);
TIMEPOINTS(:, 1) = animals;
TIMEPOINTS(:, 2) = experimental_conditions;
frame_rate = 5;
% frames_idx = 1:1:length(traces_SP_1);
tau_decay = cell2mat(GC.decay_time_constant_calcium_reporter.GCaMP6f);

%% Get data fromm animals and sessions
for i_session = 1:length(sessions)
    for i_animal = 1 : n_animals
        disp (['getting spikes from: ',animals{i_animal}])

        animal_ID = animals{i_animal};
        switch condition
            case 'SP'
                filename = (['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\jointExtraction\sorted\traces_', condition,'_',sessions{i_session}, '.mat']);
                traces = ['traces_SP_',sessions{i_session}];
            case 'SP2'
                filename = (['T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\', animal_ID, '\jointExtraction\sorted\traces_', condition,'_',sessions{i_session}, '.mat']);
                traces = ['traces_SP2_',sessions{i_session}];
        end
            
        try
            traces = load_variable(filename, traces);
            n_cells = size(traces,1);
            
            spike_folder =  'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\spikes\';
            file_name = [spike_folder, ['Spikes_', animal_ID,'_', condition,'_','S',sessions{i_session},'.mat' ]];
            if ~exist(file_name)
                toggle_toolbox('OASIS_matlab', 'on')
                Ca_events = cell(n_cells,1);
                spikes = (traces .* 0);
                for i_cell = 1:n_cells
                    disp(['Cell ', num2str(i_cell)])
                    [spikes(i_cell, :),Ca_events{i_cell}] = deconvolve_spikes_cell(animal_ID, double(traces(i_cell, :)), frame_rate);
                end
                toggle_toolbox('OASIS_matlab', 'off')
                if ~exist(spike_folder, 'dir'), mkdir(spike_folder), end
                % Save Spikes
                save (file_name, 'spikes');
                disp(['Spikes stored in ''', file_name, ''''])
                % Save Ca_events idx
                Ca_events_folder =  'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\Ca_events\';
                if ~exist(Ca_events_folder, 'dir'), mkdir(Ca_events_folder), end
                file_name = [Ca_events_folder, ['Ca_events_', animal_ID,'_', condition,'_','S',sessions{i_session},'.mat' ]];
                save (file_name, 'Ca_events');
                disp(['Ca_events stored in ''', file_name, ''''])
                
            else
                spike_folder =  'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\spikes\';
                file_name = [spike_folder, ['Spikes_', animal_ID,'_', condition,'_','S',sessions{i_session},'.mat' ]];
                spikes = load_variable (file_name, 'spikes');
                Ca_events_folder =  'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\Ca_events\';
                file_name = [Ca_events_folder, ['Ca_events_', animal_ID,'_', condition,'_','S',sessions{i_session},'.mat' ]];
                Ca_events = load_variable (file_name, 'Ca_events');
                disp (['spikes already saved, animal = ',animal_ID])
            end
            
            % Convolve spikes
            % Make time axis
            trial_duration = length (spikes);%max(reshape(cell2mat(cellfun(@length, spikes, 'UniformOutput',false)),[],1));
            t = 0 : 1/frame_rate : trial_duration/frame_rate - 1/frame_rate;
            tau_growth = 0.1;
            tau_growth_frames = round(frame_rate * tau_growth);
            U = (1-exp(-t./tau_growth)) .* exp(-t./(tau_decay/1000));
            U = (U-min(U)) ./ (max(U)-min(U)); % normalize in range [0 1]
            trace_conv = spikes .* 0;
            
            denoised_folder =  'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\denoised_traces\';
            file_name = [denoised_folder, ['Denoised_', animal_ID, '_',condition,'_','S',sessions{i_session},'.mat' ]];
            if ~exist(file_name)
                % Loop through trials
                for i_cell = 1:n_cells
                    % Get response and convolve it with unitary response
                    this_r = spikes(i_cell, :);
                    this_r_convolved = conv(this_r, U, 'full');
                    % Remove extra datapoints and store result
                    trace_conv (i_cell, :)= this_r_convolved(1+tau_growth_frames:length(this_r)+tau_growth_frames);
                end
                
                % Save denoised traces
                if ~exist(denoised_folder, 'dir'), mkdir(denoised_folder), end
                
                save (file_name, 'trace_conv');
                disp(['Denoised traces stored in ''', file_name, ''''])
            else
                trace_conv = load_variable(file_name, 'trace_conv');
                disp (['denoised traces already saved, animal = ',animal_ID])
            end
            
            
            % Calculate median ampl and freq
            ampl_this_animal = NaN(n_cells,1);
            freq_this_animal = ampl_this_animal;
            for i_cell = 1:n_cells
                idx = Ca_events{i_cell};
                n_ev = length(idx);
                ampl_this_cell = NaN (n_ev, 1);
                freq_this_cell = n_ev / (length(trace_conv)/frame_rate);
                for i_ev = 1:n_ev
                    ampl_this_cell(i_ev) = trace_conv(i_cell,idx (i_ev,2));
                end
                ampl_this_animal (i_cell) = nanmedian(ampl_this_cell,1);
                freq_this_animal (i_cell) = nanmedian(freq_this_cell,1);
            end
            AMPLITUDE{i_animal, i_session} = ampl_this_animal;
            FREQUENCY{i_animal, i_session} = freq_this_animal;
        catch
            disp(['no ', condition,'  for this animal, ', animal_ID])
            AMPLITUDE{i_animal, i_session} = NaN(n_cells, 1);
            FREQUENCY{i_animal, i_session} = NaN(n_cells, 1);

        end

    end
end
disp ('done computing spikes')

%% arrange data (amplitude and frequency)
% Day 1 is defined as the first session before surgery
disp('Re-arranging data for plotting')

n_columns = length(sessions); % only one timepoint (only after surgery)

% Mark groups
DATA_amplitude = NaN(0, n_columns + 1);
DATA_frequency = NaN(0, n_columns + 1);
sessions_after_surgery = str2double(sessions);


for i_mouse = 1:n_animals
    %         timepoints = 1; % ones(TIMEPOINTS{i_mouse, 3});
    %         if isempty(timepoints), continue, end
    
    amplitudes = NaN(size(AMPLITUDE{i_mouse}, 1), n_sessions);
    frequencies = amplitudes;
    
    try
        amplitudes(:, 1:2) = [AMPLITUDE{i_mouse,:}];
        frequencies(:,1:2) = [FREQUENCY{i_mouse,:}];
        %     group_idx = experimental_groups{ismember(experimental_groups(:, 1), TIMEPOINTS(i_mouse, 2)), 2};
        group_idx = experiemental_group(i_mouse);
    catch
        disp('animal with no SP2')
    end
    DATA_amplitude = [DATA_amplitude; [ones(size(amplitudes, 1), 1)  * group_idx, amplitudes]];
    DATA_frequency = [DATA_frequency; [ones(size(frequencies, 1), 1) * group_idx, frequencies]];
end


% Replace group index with name
column_names = {'Session_1', 'Session_2'};
DATA_amplitude_table = array2table(DATA_amplitude, 'VariableNames',['group', column_names]);
DATA_frequency_table = array2table(DATA_frequency, 'VariableNames',['group', column_names]);
DATA_amplitude_table.group = str2mat(experimental_groups(1,DATA_amplitude_table{:, 'group'})');
DATA_frequency_table.group = str2mat(experimental_groups(1,DATA_frequency_table{:, 'group'})');

% Replace group index with name
% DATA_amplitude_table.group = experimental_groups(DATA_amplitude_table{:, 'group'}, 1);
% DATA_frequency_table.group = experimental_groups(DATA_frequency_table{:, 'group'}, 1);


% Save data
output_folder = 'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results\';
file_name = [output_folder, [condition, '_analysis_all.mat']];
save(file_name,'DATA_amplitude_table', 'DATA_frequency_table')

%%
% ampl_sham = DATA_amplitude_table.Session_1(DATA_amplitude_table.group==1);
% ampl_cci = DATA_amplitude_table.Session_1(DATA_amplitude_table.group==2');
% freq_sham = DATA_frequency_table.Session_1(DATA_frequency_table.group==1);
% freq_cci = DATA_frequency_table.Session_1(DATA_frequency_table.group==1);


% disp ( ['ampl_sham: ',mat2str( nanmean(ampl_sham))])
% disp (['ampl_cci : ', mat2str(nanmean(ampl_cci))])
% [H, p] = ttest2 (ampl_sham,ampl_cci );
% disp (['p_value: ', mat2str(p)])
% disp ( ['freq_sham: ',mat2str( mean(freq_sham))])
% disp ( ['freq_cci: ',mat2str( mean(freq_cci))])



% Compute data range
amplitude_range = table2array(DATA_amplitude_table(:, 2:end));
amplitude_range = amplitude_range(:);
amplitude_range = [nanmin(amplitude_range), nanmax(amplitude_range)];
amplitude_range = [floor(amplitude_range(1) / 10) * 10, ceil(amplitude_range(2) / 10) * 5];
amplitude_range_str = ['"', strjoin(value2str(amplitude_range, '%.2f'), ', '), '"'];
frequency_range = table2array(DATA_frequency_table(:, 2:end));
frequency_range = frequency_range(:);
frequency_range = [nanmin(frequency_range), nanmax(frequency_range)];
frequency_range_str = ['"', strjoin(value2str(frequency_range, '%.5f'), ', '), '"'];


% Write data to disk so that it can be passed to the R script plotting

filename_amplitude = os.path.join(GC.temp_dir, 'amplitude.csv');
filename_frequency = os.path.join(GC.temp_dir, 'frequency.csv');
writetable(DATA_amplitude_table, filename_amplitude);
writetable(DATA_frequency_table, filename_frequency);
% Run analysis in R
output_folder = 'D:\';
% output_folder = [GC.data_root_path, GC.plots_path, 'miniscope_experiments'];
output_filename = os.path.join(output_folder, ['amplitude_',condition,'_MS.pdf']);
run_in_R(os.path.join('Figures_paper_FMP', 'amplitude_frequency_distributions_2days_MS'), ['--file_input ', filename_amplitude], ['--file_output ', output_filename], '--data_type amplitude', ['--data_limits ', amplitude_range_str], '--log_scale FALSE'); printFilePath(output_filename, [], 0)
output_filename = os.path.join(output_folder, ['frequency',condition,'_MS.pdf']);
run_in_R(os.path.join('Figures_paper_FMP', 'amplitude_frequency_distributions_2days_MS'), ['--file_input ', filename_frequency], ['--file_output ', output_filename], '--data_type frequency', ['--data_limits ', frequency_range_str]); printFilePath(output_filename, [], 0)

