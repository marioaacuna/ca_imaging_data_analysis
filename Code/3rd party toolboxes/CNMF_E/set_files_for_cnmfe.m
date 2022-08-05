% Pipeline for automatic segmentation using CNMF-e
% For 2-p imaging, so far
%% Concatenate all the TIFF files corresponding to an individual session.
% Transform the data into h5 and save it locally per session.

%% MAIN
function done = set_files_for_cnmfe(INFO)

    %% Run CNMF or CNMF-e cell extraction
    % It seems that files need to be concatenated in a single h5 file.
    % set paths
    animal_ID = 'MA_13b';
    prep_filename = ['D:\_MATLAB_CaImaging\test_CNMFE\', animal_ID, '.mat'];
    cnmfe_folder = 'C:\Users\acuna\Documents\CNMF_E';
      
    cnmfe_loaded = 0;
    if ~exist('cnmfe_loaded', 'var') || ~cnmfe_loaded
        addpath(fullfile(cnmfe_folder, 'ca_source_extraction'));
        addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'utilities')));
        addpath(genpath(fullfile(cnmfe_folder, 'ca_source_extraction', 'endoscope')));
        addpath(fullfile(cnmfe_folder, 'GUI'));
        addpath(fullfile(cnmfe_folder, 'GUI', 'gui_callbacks'));
        addpath(fullfile(cnmfe_folder, 'GUI', 'modules'));
        addpath(fullfile(cnmfe_folder, 'OASIS_matlab'));
        addpath(fullfile(cnmfe_folder, 'scripts'));
        cnmfe_loaded = true;
    end
    oasis_setup()
    
    % Set parameters
    cellWidth = 10;
    cnmfeOptions.gSiz = cellWidth;
    cnmfeOptions.gSig = ceil(cellWidth/4);
%     cnmfeOptions.batch_frames = 6000;
    cnmfeOptions.batch_frames = 500;

    cnmfeOptions.use_parallel = 1;
    cnmfeOptions.patch_dims = [128, 128]; % [64, 64]
%     cnmfeOptions.gSiz = 11;
	cnmfeOptions.Fs = 10;
    % Int: maximum number of neurons per patch. when K=[], take as many as possible.
    cnmfeOptions.K = 6;
    cnmfeOptions.min_corr = 0.6;
    cnmfeOptions.center_psf = false;  % set the value as true when the background fluctuation is large (usually 1p data)
	cnmfeOptions.min_corr_res = 0.3;
	% Float: stands for minimum peak-to-noise ratio to look for a cell
	cnmfeOptions.min_pnr_res = 4;
	cnmfeOptions.bg_model = 'ring'; %'ring';
    cnmfeOptions.ring_radius = 1; %18
    options.center_psf = 0;  % set the value as true when the background fluctuation is large (usually 1p data)
	% set the value as false when the background fluctuation is small (2p)
 	% ===COMPUTATION
	% Float: GB, memory space you allow to use in MATLAB
	cnmfeOptions.memory_size_to_use = 8; %
	% Float: GB, memory space you allow to use in MATLAB
	cnmfeOptions.memory_size_per_patch = 2; % 0.6
	% Int vector: patch size in pixels

    [cnmfeAnalysisOutput] = computeCnmfeSignalExtraction_batch(prep_filename,'options',cnmfeOptions);

    % Save outputs to NWB format
    if saveAnalysis==1
        % Save CNMF
        saveNeurodataWithoutBorders(cnmfAnalysisOutput.extractedImages,{cnmfAnalysisOutput.extractedSignals,cnmfAnalysisOutput.extractedSignalsEst},'cnmf',[nwbFilePath '_cnmf.nwb']);
        
        % Save CNMF-E
        saveNeurodataWithoutBorders(cnmfeAnalysisOutput.extractedImages,{cnmfeAnalysisOutput.extractedSignals,cnmfeAnalysisOutput.extractedSignalsEst},'cnmfe',[nwbFilePath '_cnmfe.nwb']);
    end
    [success] = cnmfVersionDirLoad('none');
    
    % =================================================
    %% Run EXTRACT cell extraction. Check each function with "edit" for options.
    % Load default configuration
    ciapkg.loadBatchFxns('loadEverything');
    extractConfig = get_defaults([]);
    
    % See https://github.com/schnitzer-lab/EXTRACT-public#configurations.
    extractConfig.avg_cell_radius = cellWidth;
    extractConfig.num_partitions_x = 2;
    extractConfig.num_partitions_y = 2;
    extractConfig.use_sparse_arrays = 0;
    
    outStruct = extractor(inputMovie3,extractConfig);
    
    % Grab outputs and put into standard format
    extractAnalysisOutput.filters = outStruct.spatial_weights;
    % permute so it is [nCells frames]
    extractAnalysisOutput.traces = permute(outStruct.temporal_weights, [2 1]);
    
    % Other run information if saving as a MAT-file.
    extractAnalysisOutput.info = outStruct.info;
    extractAnalysisOutput.config = outStruct.config;
    extractAnalysisOutput.info = outStruct.info;
    extractAnalysisOutput.userInputConfig = extractConfig;
    extractAnalysisOutput.opts = outStruct.config;
    
    % Save outputs to NWB format
    if saveAnalysis==1
        saveNeurodataWithoutBorders(extractAnalysisOutput.filters,{extractAnalysisOutput.traces},'extract',[nwbFilePath '_extract.nwb']);
    end
    
    % Remove EXTRACT from the path.
    ciapkg.loadBatchFxns();
    
    % =================================================
    %% USER INTERFACE Run cell sorting using matrix outputs from cell extraction.
    if guiEnabled==1
        [outImages, outSignals, choices] = signalSorter(pcaicaStruct.IcaFilters,pcaicaStruct.IcaTraces,'inputMovie',inputMovie3);
        
        % Plot results of sorting
        figure;
        subplot(1,2,1);imagesc(max(IcaFilters,[],3));axis equal tight; title('Raw filters')
        subplot(1,2,2);imagesc(max(outImages,[],3));axis equal tight; title('Sorted filters')
    end
    
    
    
    
    %% Run cross-session alignment of cells
    % Create input images, cell array of [x y nCells] matrices
    inputImages = cellfun(@(x) x.IcaFilters,pcaicaStructCell,'UniformOutput',false);
    
    % options to change
    opts.maxDistance = 5; % distance in pixels between centroids for them to be grouped
    opts.trialToAlign = 1; % which session to start alignment on
    opts.nCorrections = 1; %number of rounds to register session cell maps.
    opts.RegisTypeFinal = 2; % 3 = rotation/translation and iso scaling; 2 = rotation/translation, no iso scaling
    
    % Run alignment code
    [alignmentStruct] = matchObjBtwnTrials(inputImages,'options',opts);
    
    % Global IDs is a matrix of [globalID sessionID]
    % Each (globalID, sessionID) pair gives the within session ID for that particular global ID
    globalIDs = alignmentStruct.globalIDs;
    
    % View the cross-session matched cells, saved to `private\_tmpFiles` sub-folder.
    [success] = createMatchObjBtwnTrialsMaps(inputImages,alignmentStruct);


    %% Display cross-session matching movies
    
    
    
    
    
    
    
    
    
    
    % Disable toolboxes used in this section
    toggle_toolbox(toolboxes_to_use, 'off')
    
    % Increment counter
    current_action = current_action + 1;
end


%% MLint exceptions
%#ok<*AGROW,*NASGU,*STRNU>







