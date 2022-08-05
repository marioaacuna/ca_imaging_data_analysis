% This file contains a set of general configurations. All variables should
% be loaded into a structure called GC.

function GC = general_configs()
    % Initialize structure
    GC = struct();
    
    % Set version of repository
    GC.version = '0.9.6';
    
    %% DATABASE
    GC.database_name = 'nevian2';  % Name of SQL database
    GC.database_table_prefix = 'ca_imaging';  % Name of SQL database
    GC.filename_match_metadata = {'stimulus' ,'experimental_condition'};  % Columns in INFO.experiments used to match the metadata to find unique sessions

    
    %% PATHS
    GroupNevian4_pc140_letter = 'V:\';  % This drive should register itself always with this letter
    
    % Folders
    root_folder                 = 'Ca_imaging_pain';
    GC.data_root_path           = fullfile(GroupNevian4_pc140_letter, root_folder);
    GC.log_path                 = '0_logs';
    GC.tiff_metadata_path       = '0_metadata';
    GC.movie_metadata_path      = '0_movie_metadata';
    GC.stimuli_path             = '0_stimuli';
    GC.tiff_raw_path            = '1_raw_tiff';
    GC.movie_raw_path           = '1_raw_movies';
    GC.miniscope_avi_path       = '\\pc200\GroupNevian2\Mario\miniscope data\Pain_behavior_miniscope';%'\\pc8\GroupNevian5\Pain_EPM_miniscope';
    GC.registered_images        = 'K:\Ca_imaging_pain\2_motion_corrected_movies';%'2_motion_corrected_movies';
    GC.segmentation_path        = '3_segmentation';
    GC.fluorescence_data_path   = '4_fluorescence_traces';
    GC.Ca_events_path           = '4_Ca_events';
    GC.Ca_events_path_amplitude = '5_Ca_events';
    GC.deconvolved_traces_path  = '5_deconvolved_traces';
    GC.plots_path               = '6_plots';
    GC.aggregated_data_path     = '6_data';
    
    % Get path of this file
    current_path = mfilename('fullpath');
    % Remove filename to get root path of the repository
    repository_root_path = regexp(current_path, filesep(), 'split');
    GC.repository_root_path = fullfile(repository_root_path{1:end-1});
    % Add path to 3rd-party toolboxes
    GC.toolboxes_root_path = fullfile(GC.repository_root_path, 'Code', '3rd party toolboxes');
    % List of folders to not add to MATLAB path
    GC.forbidden_folders = {'Superuser', '\.', '3rd party toolboxes', 'Documentation', 'python', '_test'};  % This list will be passed to regexp (make sure special characters are properly escaped)
	GC.figures_path_Mario = ['C:\Users\acuna\Documents\Figures_for_papers_Mario'];
    % The following is a local path where temporary files will be stored
    GC.temp_root_dir = 'D:\';
    if ~exist(GC.temp_root_dir, 'dir')
        GC.temp_root_dir = tempdir();
    end
    % Make sure there is the folder of interest at that location
    GC.temp_dir = fullfile(GC.temp_root_dir, '_MATLAB_CaImaging');
    try
        if ~exist(GC.temp_dir, 'dir')
            mkdir(GC.temp_dir)
        end
    catch
        GC.temp_root_dir = tempdir();
        GC.temp_dir = fullfile(GC.temp_root_dir, '_MATLAB_CaImaging');
        if ~exist(GC.temp_dir, 'dir')
            mkdir(GC.temp_dir)
        end
    end
    
    % Files suffix (files that are common to any experiment)
    GC.motion_corrected_file_suffix          = '.mat';
    GC.template_motion_corrected_file_suffix = '_motion_correction_template.mat';
    GC.crosscorrelation_file_suffix          = '_crosscorrelation.mat';
    GC.ROI_info_file_suffix                  = '_ROI_info.mat';
    GC.raw_deltaF_over_F_file_suffix         = '_raw_deltaF_over_F.mat';
    GC.Ca_events_file_suffix                 = '_Ca_events.mat';
    GC.Ca_events_amplitude_file_suffix       = '_inferred_spikes.mat';
    GC.inferred_spikes_file_suffix           = '_inferred_spikes.mat';
    GC.assemblies_file_suffix                = '_assemblies.mat';
    GC.visualization_traces_suffix           = '_traces.pdf';
    
    % Python
    GC.python = struct();
    GC.python.environment_name = 'env_ca_imaging';
    [~, msg] = system(sprintf('activate %s && python -c "import sys; print(sys.executable)"', GC.python.environment_name));
    GC.python.interpreter_path = msg(1:end-1);
    GC.python.scripts_path = fullfile(GC.repository_root_path, 'Code', 'python');
    
    % R
    GC.R = struct();
    GC.R.scripts_path = fullfile(GC.repository_root_path, 'Code', 'R');

    
    %% PREPROCESSING   
    % Motion correction
    GC.motion_correction_init_batch = 120;  % seconds of randomly selected frames to be taken for computing initial template
	GC.motion_correction_bin_width  = 10 ;  % length of bin (in seconds) over which the registered frames are averaged to update the template
    GC.crosscorrelation_time_smoothing_window = 5;  % number of seconds to smooth when computing crosscorrelation. Windows are all-but-1-frame overlapping with each other.
    % Remove temporary data file at the end of the preprocessing
    GC.preprocessing_delete_temp_file = true;
    
    % Set parameters to analyze miniscope data
    GC.epifluorescence_downsample_to_frame_rate = 5;  % frames / s
    % Set nPCs and nICs according to the numbers of expected neurons 
    GC.epifluorescence_PCAICA_nPCs = 100;  % ~ 2 x number of expected neurons
    GC.epifluorescence_PCAICA_nICs = 150;  % ~ 1.5 x number of expected neurons
    GC.epifluorescence_segmentation_n_frames_per_chunk = ceil(60 / (1 / GC.epifluorescence_downsample_to_frame_rate));  % number of frames to group together in each "trial". Combine frames to reach 60s
    GC.epifluorescence_skip_frames_with_zeros = .9;  % remove first frame if fraction of 0-pixels crosses this threshold    
    
    
    %% CACLIUM TRACE EXTRACTION
    GC.percentile_mean_luminance_subtraction = 1;  % E.g., 1% of mean luminance in each frame is subtracted from each frame.
    GC.neuropil_donut_size                   = [5, 22];  % Distances (in pixels) from the border of the cell to donut inner limit, and from the border of the cell to donut outer limit.
    GC.crosscorrelation_threshold            = 0.25;  % Threshold to determine that a given pixel in the maximal projection of the CC belongs to an active cell (which can be included in a ROI or not).
    GC.decay_time_constant_calcium_reporter = struct();  % ms
    GC.decay_time_constant_calcium_reporter.GCaMP6f = 380;
    GC.decay_time_constant_calcium_reporter.GCaMP6s = 380; % fix this value
    GC.decay_time_constant_calcium_reporter.GCaMP7f = 4000; % originally, 520 suggested by paper
    GC.rise_time_constant_calcium_reporter.GCaMP7f  = 1.2;
    GC.rise_time_constant_calcium_reporter.GCaMP6f  = 0.1;
    GC.calcium_transient_detection_SNR       = 2;  % This value is used to threshold the CWT tree of a fluorescence trace

    % Remove temporary data file at the end of the preprocessing
    GC.trace_extraction_delete_temp_file = false;

    % Set duration of baseline and evoked window for miniscope data
    GC.miniscope_trial_duration = struct();
    GC.miniscope_trial_duration.heat = [5, 10];
    GC.miniscope_trial_duration.cold = [5, 10];
    GC.miniscope_trial_duration.pinprick = [5, 10];
    GC.miniscope_trial_duration.touch = [5, 10];
    GC.miniscope_trial_duration.puff = [5, 10];
    GC.miniscope_trial_duration.TMT = [10, 15];
    GC.miniscope_trial_duration.water = [10, 15];
    
    
    %% ANALYSIS (set parameters per experiment)
     % --------------------------------------------------------------------------
        
    % --------------------------------------------------------------------------
    % Name: ACC_CCI_anesth
    % Aim: Effect of noxious and non-noxious stimuli onto ACC neuron activity of anesthetized mice.
    % --------------------------------------------------------------------------
    experiment_name = 'ACC_CCI_anesth';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
    GC.columns_experimental_condition.(experiment_name) = {'data_type=2p',  {'day_from_surgery', 'compound'};
                                                           'data_type=epi', {'day_from_surgery', 'experiment'}};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'SP', 'HPS', 'puff', 'temp_48', 'sound'};
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'day_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
    GC.analysis.(experiment_name).session_column = {'date'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'sham',[.05,.48,.75]; 'CCI',[1,0,0]};
    % Set which stimuli are shown in the heatmap that shows response modulation
    GC.response_modulation_heatmap_stimuli = {'HPS'};

    
    % --------------------------------------------------------------------------
    % Name: ACC_SNI_anxiety
    % Aim: Effect of noxious and non-noxious stimuli onto ACC neuron activity and anxiety of SNI mice.
    % --------------------------------------------------------------------------
    experiment_name = 'ACC_SNI_anxiety';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
    GC.columns_experimental_condition.(experiment_name) = {'data_type=epi',{'day_from_surgery', 'experiment'};
                                                           'data_type=epi', {'day_from_surgery', 'experiment'}};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'SP', 'pinprick', 'touch', 'heat', 'cold'};
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'day_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
    GC.analysis.(experiment_name).session_column = {'date'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};
    % Set number of frames to be analyzed
    GC.analysis.(experiment_name).n_frames = 76;
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'sham',[.05,.48,.75]; 'CCI',[1,0,.5]; 'SNI',[1,0,0]};
    

    % --------------------------------------------------------------------------
    % Name: CLA_pain
    % Aim: Effect of noxious and non-noxious stimuli onto CLA neuron activity of anesthetized mice.
    % --------------------------------------------------------------------------
    experiment_name = 'CLA_pain';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
     GC.columns_experimental_condition.(experiment_name) = {'data_type=2p',  {'day_from_surgery', 'compound'};
                                                           'data_type=epi', {'day_from_surgery', 'experiment'}};
%     GC.columns_experimental_condition.(experiment_name) = {'condition';
%                                                            'condition'};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli =  {'SPn', 'HPSn', 'puffn', 'soundn'};
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'cond_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
    GC.analysis.(experiment_name).session_column = {'date'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'SAL',[.05,.48,.75]; 'CFA',[1,0,0]};
    
    % --------------------------------------------------------------------------
    % Name: CLA-ACC_dreadd
    % Aim: Effect of noxious and non-noxious stimuli onto ACC neuronal
    % activity of anesthetized mice, while inhibiting the CLA by DREADDs
    % --------------------------------------------------------------------------
    experiment_name = 'CLA_ACC_dreadd';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
     GC.columns_experimental_condition.(experiment_name) = {'data_type=2p',  {'day_from_surgery', 'compound', 'compound_phase'};
                                                           'data_type=epi', {'day_from_surgery', 'experiment'}};
%     GC.columns_experimental_condition.(experiment_name) = {'condition';
%                                                            'condition'};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli =  {'SPn', 'HPSn', 'puffn', 'soundn'};
    % Set which column to use to determine how to group trials by timepoint
%     GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'experimental_condition'; % 'compound_phase';

    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'cond_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
%     GC.analysis.(experiment_name).session_column = {'date'};
%     GC.analysis.(experiment_name).session_column_prefix = 'session_';
%     
%     GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus','compound_phase'};
    GC.analysis.(experiment_name).session_column = {'date'};
%     GC.analysis.(experiment_name).session_column = {'compound_phase'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';

    
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {'SAL', 'CNO'}; % pre CNO, post CNO
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'SAL_dreadd',[.05,.48,.75]; 'CFA_dreadd',[1,0,0]; 'SAL_reporter',[.05,.48,.75]; 'CFA_reporter',[1,0,0]};

    
    % --------------------------------------------------------------------------
    % Name: Anaesth_check
    % Aim: Effect of anesthesia level on stimulus selectivity in ACC.
    % --------------------------------------------------------------------------
    experiment_name = 'Anaesth_check';
    GC.analysis.(experiment_name) = struct();
    GC.columns_experimental_condition.(experiment_name) = {'day_from_surgery', 'compound'};
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'SP', 'HPS', 'puff'};
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'experimental_condition';% anesthesia_level
    GC.analysis.(experiment_name).session_name_prefix = 'ISO_level_';
    GC.analysis.(experiment_name).group_trials_by_session = {'stimulus', 'experimental_condition'}; % 'date', 'stimulus'
    GC.analysis.(experiment_name).session_column = {'experimental_condition'}; % 'date'
    GC.analysis.(experiment_name).session_column_prefix = 'level_';
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};%{'SAL', ''};
    GC.experimental_groups_colors.(experiment_name) = {'Naive',[1,.65,0]; 'sham',[.05,.48,.75]; 'CCI',[1,0,0]};

    
    % --------------------------------------------------------------------------
    % Name: GBP_pain
    % Aim: Effect of Gabapentin (GBP) on the neuronal responses of of SNI mice.
    % --------------------------------------------------------------------------
    experiment_name = 'GBP_pain';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
    GC.columns_experimental_condition.(experiment_name) = {'data_type=epi',{'day_from_surgery', 'experiment'};
                                                           'data_type=epi', {'day_from_surgery', 'experiment'}};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'SP', 'pinprick', 'touch', 'heat', 'cold'};
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'day_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
    GC.analysis.(experiment_name).session_column = {'date'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};
    % Set number of frames to be analyzed
    GC.analysis.(experiment_name).n_frames = 76;
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'sham',[.05,.48,.75]; 'CCI',[1,0,.5]; 'SNI',[1,0,0]};
    

    % --------------------------------------------------------------------------
    % Name: ACC_pain_LP211
    % Aim: Effect of LP-211 on stimulus-evoked activity in ACC.
    % --------------------------------------------------------------------------
    experiment_name = 'ACC_pain_LP211';
    GC.analysis.(experiment_name) = struct();
    GC.columns_experimental_condition.(experiment_name) = {'day_from_surgery', 'compound', 'compound_phase'};
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'SP', 'HPS'};
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'compound_phase';
    GC.analysis.(experiment_name).session_name_prefix = '';
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus', 'compound_phase'};
    GC.analysis.(experiment_name).session_column = {'compound_phase'};
    GC.analysis.(experiment_name).session_column_prefix = '';
    GC.analysis.(experiment_name).analysis_allowed_compounds = {'SAL', 'LP-211'};
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'sham',[.05,.48,.75]; 'CCI',[1,0,0]};

    % --------------------------------------------------------------------------
    % Name: ACC_neuron_aversion
    % Aim: Effect of noxious and non-noxious stimuli onto ACC neuron activity and anxiety of SNI mice.
    % --------------------------------------------------------------------------
    experiment_name = 'ACC_neuron_aversion';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
    GC.columns_experimental_condition.(experiment_name) = {'data_type=epi',{'experiment'};
                                                           'data_type=epi', {'experiment'}};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'water', 'TMT'};
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'day_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
    GC.analysis.(experiment_name).session_column = {'date'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};
    % Set number of frames to be analyzed
    GC.analysis.(experiment_name).n_frames = 76;
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]};%; 'sham',[.05,.48,.75]; 'CCI',[1,0,.5]; 'SNI',[1,0,0]};

    % --------------------------------------------------------------------------
    % Name: ACC_neuron_aversion
    % Aim: Effect of noxious and non-noxious stimuli onto ACC neuron activity and anxiety of SNI mice.
    % --------------------------------------------------------------------------
    experiment_name = 'Aversion_puff';
    GC.analysis.(experiment_name) = struct();
    % Split the column 'experimental_condition' in the database into these
    GC.columns_experimental_condition.(experiment_name) = {'data_type=epi',{'day_from_surgery', 'experiment'};
                                                           'data_type=epi', {'day_from_surgery', 'experiment'}};
    % Only these stimuli will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_stimuli = {'TMT', 'water'};
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).group_trials_by_timepoint = 'day_from_surgery';
    % Set which column to use to determine how to group trials by timepoint
    GC.analysis.(experiment_name).session_name_prefix = 'day_';
    % Set which column to use to determine how to group trials by session
    GC.analysis.(experiment_name).group_trials_by_session = {'date', 'stimulus'};
    GC.analysis.(experiment_name).session_column = {'date'};
    GC.analysis.(experiment_name).session_column_prefix = 'session_';
    % Only these compounds will be analyzed
    GC.analysis.(experiment_name).analysis_allowed_compounds = {''};
    % Set number of frames to be analyzed
    GC.analysis.(experiment_name).n_frames = 126;
    % Visualization parameters
    GC.experimental_groups_colors.(experiment_name) = {'naive',[1,.65,0]; 'sham',[.05,.48,.75]; 'CCI',[1,0,.5]; 'SNI',[1,0,0]};





    % --------------------------------------------------------------------------
    % Response modulation (detection)
    GC.detection_baseline_window             = 5;  % s, before the first timestamp, of activity to consider as baseline
    GC.evoked_activity_max_latency           = 5;
    GC.response_detection_file_suffix        = ['_response_detection_', num2str(GC.evoked_activity_max_latency), 's_max_latency.mat'];
    GC.detection_n_permutations_significance = 1000;
    GC.response_detection_window_width       = 0.5;  % s
    GC.response_decoding_window_width        = 3;  % s
    GC.detection_k_fold_crossvalidation      = 5;
    % These stimuli don't have timestamps and are rather sustained. The algorithm 
    % to detect the response won't work on them
    GC.detection_no_timestamp_stimuli        = {'acet-', 'acet+', 'reward-', 'reward+', 'EPM'};
    
    % Response modulation (decoding)
    GC.response_decoding_stats_filename_all                         = 'response_decoding_stats.mat';
    GC.response_decoding_stats_filename_selective                   = 'response_decoding_stats_selective.mat';
    GC.response_decoding_stats_filename_stable                      = 'response_decoding_stats_stable.mat';
    GC.response_decoding_stats_filename_ensembles                   = 'response_decoding_stats_ensembles.mat';
    GC.response_decoding_stats_filename_specific                    = 'response_decoding_stats_specific.mat';
    GC.response_decoding_stats_filename_H_noci                      = 'response_decoding_stats_H_noci.mat';
    GC.response_decoding_stats_filename_H_noci_random               = 'response_decoding_stats_H_noci_random.mat';
    GC.response_decoding_stats_filename_H_sal                       = 'response_decoding_stats_H_sal.mat';
    GC.response_decoding_stats_filename_H_sal_random                = 'response_decoding_stats_H_sal_random.mat';
    GC.response_decoding_stats_filename_H_sal_categories            = 'response_decoding_stats_H_sal_categories.mat';
    GC.response_decoding_stats_filename_H_sal_random_categories     = 'response_decoding_stats_H_sal_random_categories.mat';
    GC.response_decoding_stats_filename_all_categories              = 'response_decoding_stats_all_cells_categories.mat';



    % Which stimuli do cells need to be selective to in the "selective"-only decoding?
    % Each cell of the array contains the list of stimuli to which cells should
    % be selective. If multiple stimuli are present, only cells that are
    % selective to all will be considered. Each cell of the array will be run
    % separately as a new analysis.
%     GC.response_decoding_selective_to         = {
%                                                  {'HPS'}, ...
%                                                  {'temp_48'}, ...
%                                                  {'HPS', 'temp_48'}, ...                                                 
%                                                  };
    GC.response_decoding_selective_to         = {
                                                     {'HPS'}, ...
                                                     {'temp_48'}, ...
                                                     {'puff'}, ...
                                                     {'temp_43'},...
                                                     {'sound'},...
                                                     {'pinprick'},...
                                                     {'touch'},...
                                                     {'FPS'}
                                                     };
    GC.response_decoding_baseline_window      = 0;
    GC.response_decoding_evoked_window        = 5;
    GC.decoding_n_permutations_significance   = 100;
    GC.decoding_CV_k_folds                    = 5;
    
    % Cell assemblies
    GC.assemblies_time_bin_activity_n_frames = 2;
    GC.assemblies_stimuli_to_combine         = {{'SP'},{'HPS'},{'FPS'},{'temp_48'}, {'pinprick'},{'puff'},{'sound'},{'temp_43'},{'touch'}};%{{'HPS', 'temp_48', 'SP'}}; % 'HPS', 'temp_48', 'SP'
    
    % Spotaneous activity (CDF)
    GC.spontaneous_activity_stats_filename_epi = 'spontaneous_activity_epi.mat';
	GC.spontaneous_activity_stats_filename     = 'spontaneous_activity_CDF.mat';
    GC.spontaneous_activity_attributes         = {'amplitude', 'interval'};
    GC.spontaneous_amplitude_bins              = linspace(0, 100, 10000);
    GC.spontaneous_interval_bins               = 1:5100;
    GC.spontaneous_groups_to_plot              = {'sham', 'CCI'};
    
    
    %% BEHAVIOR
    GC.behavior_raw_video_root_path = '\\pc200\GroupNevian2\Mario\Behavior\ACC_pain';
    GC.behavior_video_subfolder = 'behavior cams';
    GC.FFmpeg_exe = 'M:\Software\FFmpeg\bin\ffmpeg.exe';
    GC.freely_moving_stimuli = {'cold', 'heat', 'pinprick', 'touch', 'puff', 'TMT', 'water'};
    GC.freely_moving_affective_responses = {'other';
                                            'licking / biting / extended lifting / guarding / flinching'; 
                                            'escape / hyperlocomotion / rearing / jumping'};
    
    
    %% PLOTS
    % Set how cells are sorted in this heatmap. Options are:
    % 'none': no sorting.
    % 'session_1': sort cells in the session 1 from the most to the least
    % responsive.
    % 'all': sort each session from the most the least responsive cells
    % separately from each other.
    GC.response_modulation_heatmap_sort = 'session_1';
    % Set plotting options for graphs made in python
    GC.python.font_size_labels  = 16;
    GC.python.scatterplot_small = 7;
    GC.python.scatterplot_large = 10;
    
    
    
    %% TO BE FILLED IN GUI
    GC.experiment_name = '';
    

%% MLint exceptions
%#ok<*CTCH>
