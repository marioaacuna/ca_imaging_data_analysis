function filename = get_filename_of(filetype, varargin)

% Read general_configs
global GC
if length(varargin) > 1
    GC.experiment_name = varargin{2};
else
    GC.experiment_name = 'ACC_CCI_anesth';
end

switch filetype
    case 'stimulus'
        filename = os.path.join(GC.data_root_path, GC.stimuli_path, [varargin{1}, '.mat']);

    case 'motion_corrected'
        if length(varargin) == 1
            get_the_structural_channel = false;
        else
            get_the_structural_channel = varargin{2};
        end
        
        if ~get_the_structural_channel  % the functional channel
            filename = os.path.join(GC.registered_images, [varargin{1}, GC.motion_corrected_file_suffix]); % moved motion corrected to new server
        else  % the structural channel
            filename = os.path.join(GC.data_root_path, GC.registered_images, [varargin{1}, '_structural_channel', GC.motion_corrected_file_suffix]);
        end
                        
    case 'ROI_info'
        filename = os.path.join(GC.data_root_path, GC.segmentation_path, [varargin{1}, GC.ROI_info_file_suffix]);
    case 'source2D'
        filename = os.path.join(GC.data_root_path, GC.segmentation_path, [varargin{1}, '_Neuron_Source2D.mat']);
    case 'ROI_map_cnmfe'
        filename = os.path.join(GC.data_root_path, GC.plots_path, 'CNMFe_plots', [varargin{1}, '.png']);
        
    case 'dFF'
        filename = os.path.join(GC.data_root_path, GC.fluorescence_data_path, [varargin{1}, GC.raw_deltaF_over_F_file_suffix]);
        
    case 'Ca_events'
        filename = os.path.join(GC.data_root_path, GC.Ca_events_path, [varargin{1}, GC.Ca_events_file_suffix]);

    case 'Ca_events_amp'
        filename = os.path.join(GC.data_root_path, GC.Ca_events_path_amplitude, [varargin{1}, GC.Ca_events_amplitude_file_suffix]);

    case 'spikes'
        filename = os.path.join(GC.data_root_path, GC.deconvolved_traces_path, [varargin{1}, GC.inferred_spikes_file_suffix]);
    
    case 'assemblies'
        filename = os.path.join(GC.data_root_path, GC.deconvolved_traces_path, [varargin{1}, GC.assemblies_file_suffix]);
        
    case 'response_detection'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'response_detection', [varargin{1}, GC.response_detection_file_suffix]);
    
    case 'response_detection_z_score'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'response_detection',[varargin{1}, '_z_score.mat']);
   
    case 'response_detection_p_0_01'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'response_detection',[varargin{1}, '_p_0_01.mat']);

    case 'spontaneous_activity_stats'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.spontaneous_activity_stats_filename);
    
    case {'response_decoding_stats', 'response_decoding_stats_all'}
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_all);
    case {'response_decoding_stats_ACC', 'response_decoding_stats_ACC_all'}
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_all);

    case 'response_decoding_stats_selective'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_selective);
    case 'response_decoding_stats_stable'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_stable);
    case 'response_decoding_stats_ensembles'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_ensembles);
    case 'response_decoding_stats_specific'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{2}, GC.response_decoding_stats_filename_specific);
    case 'response_decoding_stats_H_noci'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_H_noci);
    case 'response_decoding_stats_H_noci_random'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_H_noci_random);
    case 'response_decoding_stats_H_sal'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_H_sal);
    case 'response_decoding_stats_H_sal_random'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_H_sal_random);
    case 'response_decoding_stats_H_sal_categories'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_H_sal_categories);
    case 'response_decoding_stats_H_sal_random_categories'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_H_sal_random_categories);

    case {'response_decoding_stats_all_categories'}
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.response_decoding_stats_filename_all_categories);

        
    %% Miniscope
    case 'miniscope_movie_avi'
        filename = os.path.join(GC.miniscope_avi_path, varargin{1}, varargin{2}, varargin{3}, ['msCam', varargin{4}, '.avi']);
        
    case 'miniscope_movie_parameters'
        filename = os.path.join(GC.data_root_path, GC.movie_raw_path, varargin{1}, 'p.mat');
        
    case 'miniscope_preprocessed_movie'
        filename = os.path.join(GC.registered_images, varargin{1}, [varargin{1}, '_', varargin{2}, '_', varargin{3}, '.h5']);
        
    case 'miniscope_cell_map'
        filename = os.path.join(GC.registered_images, varargin{1}, [varargin{1}, '_', varargin{2}, '_', varargin{3}, '_cellMap.mat']);
    
    case 'spontaneous_activity_stats_epi'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, GC.spontaneous_activity_stats_filename_epi);
   
    case 'response_detection_epi'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'response_detection', [varargin{1},'_',  'p_val', '.mat']);
        
    case 'analysis_AUC' 
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, GC.experiment_name, 'analysis_AUC_pain_EPM.mat');
    case 'analysis_p_EPM'
        filename = os.path.join(GC.data_root_path, GC.aggregated_data_path, varargin{1}, 'analysis_p_pain_EPM.mat');


        
        
end
