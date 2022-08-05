% This file contains a set of general configurations for the MVPA
% classification analysis.
% All variables should be loaded into a structure called GC.

function CFG_MVPA = general_configs_MVPA()
    % Initialize structure
    CFG_MVPA = struct();
    
    % Set version of repository
    CFG_MVPA.version = '0.0.1'; %#ok<STRNU>
    
    %% Set classifier configuration
    CFG_MVPA = [] ;
    CFG_MVPA.classifier                  =  'multiclass_lda'; % 'logreg', 'kernel_fda', multiclass_lda, naive_bayes
    CFG_MVPA.metric                      = {'accuracy' 'confusion'};
    CFG_MVPA.cv                          = 'kfold';
    CFG_MVPA.k                           = 5;
    CFG_MVPA.hyperparameter              = [];
    CFG_MVPA.feedback                    = 1;
    CFG_MVPA.preprocess                  = {};
    CFG_MVPA.preprocess_param            = {};
    CFG_MVPA.append                      = 0;
    CFG_MVPA.sample_dimension            = 1;
    CFG_MVPA.feature_dimension           = 2;
    CFG_MVPA.generalization_dimension    = [];
    CFG_MVPA.append                      = false;
    CFG_MVPA.dimension_names             = {'samples', 'cells'};
    
    % Other important parameters for organizing the data before decoding
    CFG_MVPA.time_bin_activity           = 0.4;
    CFG_MVPA.last_bin                    = 12;
    

    
end