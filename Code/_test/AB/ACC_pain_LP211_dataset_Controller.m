%% Start Script
clear; clc
global GC
% Read general_configs
GC = general_configs();

% Disable warning
warning('off','MATLAB:rmpath:DirNotFound')
% Get list of folders in repository, and remove them all from MATLAB path
all_folders = genpath(GC.repository_root_path);
rmpath(all_folders)
% Split folder names
if ispc(), sep=';'; else, sep=':'; end
all_folders = regexp(all_folders, sep, 'split');
% Make sure that forbidden folders are not added to MATLAB's path
bad_folders_idx = cell2mat(cellfun(@(x) ~isempty(regexp(x,strjoin(GC.forbidden_folders,'|'),'once')), all_folders, 'UniformOutput',false));
good_folders = strjoin(all_folders(~bad_folders_idx), ';');
% Add good folders and remove bad ones
addpath(good_folders)
warning('on','MATLAB:rmpath:DirNotFound')

toggle_toolbox('_plotting', 'on');
% export_fig(filename, '-pdf', '-q101', '-nocrop', fig_handle)
%%
%-----------------%
% DATA EXTRACTION %
%-----------------%
%% Load data for SPONTANEOUS and EVOKED activity
   %-------------------------------------------%

% Load configs, names
clear, clc

global GC
GC = general_configs();
GC.experiment_name = 'ACC_pain_LP-211';
stats_filename = get_filename_of('response_modulation_stats', 'ACC_pain_LP-211');
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');

% Load dataset animals from the server
dataset_animal_names = SQL_database.read_table_where('sessions', {'animal_ID','experimental_group', 'experimental_condition'}, 'ACC_pain_LP-211','experiment_name');
dataset_animal_names = dataset_animal_names(:, {'animal_ID','experimental_group','compound'});
dataset_animal_names = unique(dataset_animal_names, 'rows', 'stable');
dataset_animal_names = natsortrows(dataset_animal_names,1);
all_animal_names = fieldnames(STATS_response_modulation);
n_datasets = size(dataset_animal_names,1);

% Find number of animals divided by 'experimental_group' and 'compound'
numbers = struct();
numbers.cells.HPS = struct();
numbers.cells.SP = struct();
numbers.animals = cell(2,2);
for a = 1:2
    for b = 1:2
        if a == 1 && b ==1
            idx = ismember(table2cell(dataset_animal_names(:,2:3)), {'CCI' 'LP-211'});
            idx(:,3) = idx(:,1).*idx(:,2);
            numbers.animals{a,b} = sum(idx(:,3));
        elseif a == 1 && b == 2
            idx = ismember(table2cell(dataset_animal_names(:,2:3)), {'CCI' 'SAL'});
            idx(:,3) = idx(:,1).*idx(:,2);
            numbers.animals{a,b} = sum(idx(:,3));
        elseif a == 2 && b == 1
            idx = ismember(table2cell(dataset_animal_names(:,2:3)), {'sham' 'LP-211'});
            idx(:,3) = idx(:,1).*idx(:,2);
            numbers.animals{a,b} = sum(idx(:,3));
        elseif a == 2 && b == 2
            idx = ismember(table2cell(dataset_animal_names(:,2:3)), {'sham' 'SAL'});
            idx(:,3) = idx(:,1).*idx(:,2);
            numbers.animals{a,b} = sum(idx(:,3));
        end
    end
end

% Load variables for spontaneous activity
stats_filename = get_filename_of('spontaneous_activity_stats', 'ACC_pain_LP-211');

% STATS_spontaneous_activity = load_variable(stats_filename, 'STATS_spontaneous_activity');
DISTRIBUTIONS = load_variable(stats_filename, 'DISTRIBUTIONS');

% Add parameters
phases_order = {'baseline', 'active_1'};
LEG_labels = {'non-selective, non-modulated', 'selective, non-modulated', 'only selective before', 'only selective after'};

toggle_toolbox('_plotting', 'on');

clearvars -except dataset_animal_names DISTRIBUTIONS GC LEG_labels numbers ...
                    phases_order STATS_response_modulation
%% Find cells AUROC/p and assign to results structure
compound_group_names = {'LP211','SAL'};
experiment_group_names = {'CCI','sham'};
results = struct();
results.HPS.DATA = cell(length(compound_group_names), length(experiment_group_names));
results.HPS.DATA_p = results.HPS.DATA;
results.SP.DATA = cell(length(compound_group_names), length(experiment_group_names));

% Results DATA structure and groups ID
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  1   | 2 | %
    % sham|  3   | 4 | %
    %------------------%
n_datasets = size(dataset_animal_names,1);

for iid = 1:n_datasets
    
    animal_ID = cell2mat(dataset_animal_names{iid,1});
    stats_HPS = STATS_response_modulation.(animal_ID);
    stats_SP = DISTRIBUTIONS.(animal_ID);
    
    [cells_ID, ~, rows_idx] = unique(stats_HPS.ROI_id);
    
    n_cells = length(cells_ID);
    
    HPS_DATA = NaN(n_cells, length(phases_order)+3);
    HPS_DATA_p = NaN(n_cells, length(phases_order)+1);
    SP_DATA = NaN(n_cells, length(phases_order)+3);
    
    % Get info on AUROC and p-value of each cell
    for icell = 1:n_cells
        data_HPS = stats_HPS(rows_idx == icell, {'ROI_id','compound_phase', 'AUROC', 'p'});
        data_SP = stats_SP(icell, {'ROI_id', 'mean_activity'});
        % Skip cells that were not active during all phases
        if height(data_HPS) ~= length(phases_order)
            continue
        end
        
        % Collect activity values 'baseline'/'active1' for HPS/SP
        for phi = 1:size(phases_order, 2)
            rows = ismember(table2cell(data_HPS(:, 'compound_phase')), phases_order{phi});
            if any(rows)
                AUROC = data_HPS{rows, 'AUROC'};
                HPS_DATA(icell, phi) = AUROC;
                
                MEAN_ACTIVITY = cell2mat(data_SP{1, 'mean_activity'});
                SP_DATA(icell, phi) = MEAN_ACTIVITY(phi);
                
                p = data_HPS{rows, 'p'};
                HPS_DATA_p(icell, phi) = p;
            end
        end
        
        % Assign group attribute to AUROC values
        %-----------------------------------------------------------%
        % GROUPS LEGEND |                                           %
        %---------------|                                           %
        % -- -- -- -- -- -- --  DATA  -- -- -- -- -- -- -- -- -- -- %
        % Group 1 = Decreased positive modulation by treatment      %
        % Group 2 = Increased positive modulation by treatment      %
        % Group 3 = Decreased negative modulation by treatment      %
        % Group 4 = Increased negative modulation by treatment      %
        % Group 0 = Any other                                       %
        % -- -- -- -- -- -- -- DATA_p -- -- -- -- -- -- -- -- -- -- %
        % Group 5 = Always significantly responsive to the stimulus %
        % Group 6 = Never significant responsiveness                %
        % Group 7 = Lost significant responsiveness                 %
        % Group 8 = Gain significant responsiveness                 %
        % Group 0 = Any other                                       %
        %-----------------------------------------------------------%
        
        % Assign group attribute to HPS_DATA
        if HPS_DATA(icell,1) > HPS_DATA(icell,2) && HPS_DATA(icell,2) > 0.5 % Decreased positive modulation by treatment
            HPS_DATA(icell, 3) = 1;
        elseif 0.5 <= HPS_DATA(icell,1)&& HPS_DATA(icell,1)<= HPS_DATA(icell,2) % Increased positive modulation by treatment
            HPS_DATA(icell, 3) = 2;
        elseif HPS_DATA(icell,1) <= HPS_DATA(icell,2)&& HPS_DATA(icell,2) <= 0.5 % Decreased negative modulation by treatment
            HPS_DATA(icell, 3) = 3;
        elseif 0.5 > HPS_DATA(icell,1) && HPS_DATA(icell,1) > HPS_DATA(icell,2) % Increased negative modulation by treatment
            HPS_DATA(icell, 3) = 4;
        else
            HPS_DATA(icell, 3) = 0;
        end
        
        % Assign group attribute to HPS_DATA_p values
        if HPS_DATA_p(icell,1)<=0.05 && HPS_DATA_p(icell,2)<=0.05 % Always significantly responsive to the stimulus
            HPS_DATA_p(icell, phi+1) = 5;
        elseif HPS_DATA_p(icell,phi-1)>0.05 && HPS_DATA_p(icell,phi)>0.05 % Never significant responsiveness
            HPS_DATA_p(icell, phi+1) = 6;
        elseif HPS_DATA_p(icell,phi-1)<=0.05 && HPS_DATA_p(icell,phi)>0.05 % Lost significant responsiveness
            HPS_DATA_p(icell, phi+1) = 7;
        elseif HPS_DATA_p(icell,1)>0.05 && HPS_DATA_p(icell,phi)<=0.05 % Gain significant responsiveness
            HPS_DATA_p(icell, phi+1) = 8;
        else
            HPS_DATA(icell, 3) = 0;
        end
        
        % Assign group attribute to SP_DATA
        if SP_DATA(icell,1) > SP_DATA(icell,phi) && SP_DATA(icell,phi) > 0.5 % Decreased positive modulation by treatment
            SP_DATA(icell, 3) = 1;
        elseif 0.5 <= SP_DATA(icell,1)&& SP_DATA(icell,1)<= SP_DATA(icell,phi) % Increased positive modulation by treatment
            SP_DATA(icell, 3) = 2;
        elseif SP_DATA(icell,1) <= SP_DATA(icell,phi)&& SP_DATA(icell,phi) <= 0.5 % Decreased negative modulation by treatment
            SP_DATA(icell, 3) = 3;
        elseif 0.5 > SP_DATA(icell,1) && SP_DATA(icell,1) > SP_DATA(icell,phi) % Increased negative modulation by treatment
            SP_DATA(icell, 3) = 4;
        else
            HPS_DATA(icell, 3) = 0;
        end

        % Compute the difference POST-PRE for each cell
        HPS_DATA(icell, 4) = HPS_DATA(icell, 2)-HPS_DATA(icell, 1);
%         HPS_DATA_p(icell, 4) = HPS_DATA_p(icell, 2)-HPS_DATA_p(icell, 1);
        SP_DATA(icell, 4) = SP_DATA(icell, 2)-SP_DATA(icell, 1);

        % Compute the normalized difference [(POST-PRE)/(POST+PRE)] for each cell
        HPS_DATA(icell, 5) = (HPS_DATA(icell, 2)-HPS_DATA(icell, 1))/(HPS_DATA(icell, 2)+HPS_DATA(icell, 1));
%         HPS_DATA_p(icell, 5) = (HPS_DATA_p(icell, 2)-HPS_DATA_p(icell, 1))/(HPS_DATA_p(icell, 2)+HPS_DATA_p(icell, 1));
        SP_DATA(icell, 5) = (SP_DATA(icell, 2)-SP_DATA(icell, 1))/(SP_DATA(icell, 2)+SP_DATA(icell, 1));
    end
  
    % Assign cell_id
    HPS_DATA = [cells_ID HPS_DATA];
    HPS_DATA_p = [cells_ID HPS_DATA_p];
    SP_DATA = [cells_ID SP_DATA];

    % Assign animal_id
    anim_ID = zeros(size(HPS_DATA,1),1);
    anim_ID(:,1) = iid; % Vector with the size of HPS_DATA containing the animal_ID
    HPS_DATA = [anim_ID HPS_DATA];
    HPS_DATA_p = [anim_ID HPS_DATA_p];
    SP_DATA = [anim_ID SP_DATA];
        
%     HPS_DATA(:,end+1) = iid;
%     HPS_DATA_p(:,end+1) = iid;
%     SP_DATA(:,end+1) = iid;
    
    % Remove NaNs HPS_DATA
    rows = any(isnan(HPS_DATA)');
    HPS_DATA(rows, :) = [];
    HPS_DATA_p(rows, :) = [];
    
    % Remove NaNs SP_DATA
    rows = any(isnan(SP_DATA)');
    SP_DATA(rows, :) = [];
    
    if strcmp(dataset_animal_names{iid, 'experimental_group'}, 'CCI')
        row = 1;
    elseif strcmp(dataset_animal_names{iid, 'experimental_group'}, 'sham')
        row = 2;
    else
        error('!!!')
    end
    if strcmp(dataset_animal_names{iid, 'compound'}, 'LP-211')
        col = 1;
    elseif strcmp(dataset_animal_names{iid, 'compound'}, 'SAL')
        col = 2;
    else
        error('!!!')
    end
    
    results.HPS.DATA{row, col} = cat(1, results.HPS.DATA{row, col}, HPS_DATA);
    results.HPS.DATA_p{row, col} = cat(1, results.HPS.DATA_p{row, col}, HPS_DATA_p);
    results.SP.DATA{row, col} = cat(1, results.SP.DATA{row, col}, SP_DATA);
end


% clearvars iid icell phi row col HPS_DATA HPS_DATA_p SP_DATA n_cells ...
%         p cells_ID AUROC animal_ID data_SP data_HPS rows rows_idx
    
clearvars -except compound_group_names dataset_animal_names DISTRIBUTIONS ...
                  experiment_group_names GC numbers results STATS_response_modulation

%____________________________________________________________________________
%|	  1    |   2   | 3 | 4  |    5     |       6       |	     7          |
%|Animal_ID|Cell_ID|PRE|POST|Cell_group|Diff_(POST-PRE)|Norm_diff_(POST-PRE)|
%% Isolate significant modulated cells during HPS and SP, Add cell ID to 'results.HPS.ID_' to keep track of cell behavior
results.HPS.DATA_sig = cell(2,2); % Significant selective cells in HPS

% Add 'zero-one' filter to discriminate significant modulated cells in DATA_p structure
for a = 1:2
    for b = 1:2
        for i = 1:size(results.HPS.DATA_p{a,b},1)
            if results.HPS.DATA_p{a,b}(i,3) <= 0.05;
                results.HPS.DATA_p{a,b}(i,6) = 1;
            else
                results.HPS.DATA_p{a,b}(i,6) = 0;
            end
            if results.HPS.DATA_p{a,b}(i,4) <= 0.05;
                results.HPS.DATA_p{a,b}(i,7) = 1;
            else
                results.HPS.DATA_p{a,b}(i,7) = 0;
            end
        end
    end
end

% Select only significant modulated cells in DATA structure with 'zero-one' vector filter
for a = 1:2
    for b = 1:2
        results.HPS.DATA{a,b}(:,8:9) = results.HPS.DATA_p{a,b}(:,6:7).*results.HPS.DATA{a,b}(:,3:4);
    end
end

% Keep significant cells during HPS
results.HPS.DATA_sig = results.HPS.DATA;
for a = 1:2
    for b = 1:2
        for i = 1:size(results.HPS.DATA_p{a,b},1)
            if results.HPS.DATA{a,b}(i,8) == 0 && results.HPS.DATA{a,b}(i,9) == 0 % Always non modulated
                results.HPS.DATA_sig{a,b}(i,10) = 0;
%                 results.HPS.ID_sig{a,b}(i,1) = 0;
%             elseif results.HPS.DATA{a,b}(i,8) > 0 && results.HPS.DATA{a,b}(i,9) > 0 % Always modulated
%                 results.HPS.DATA_sig{a,b}(i,10) = 1;
%                 results.HPS.ID_sig{a,b}(i,1) = 1;
            else
                results.HPS.DATA_sig{a,b}(i,10) = 1;
%                 results.HPS.ID_sig{a,b}(i,1) = 1;
            end
        end
        results.HPS.DATA_sig{a,b} = results.HPS.DATA_sig{a,b}(any(results.HPS.DATA_sig{a,b}(:,10),10),:); % Remove rows with zeros in the index column
        results.HPS.DATA_sig{a,b}(:,10) = []; % Remove index column
    end
end

% Keep significant cells during SP
results.SP.DATA_sig = cell(2,2);
for a = 1:2
    for b = 1:2
        idx_A = ismember(results.SP.DATA{a,b}(:,1:2), results.HPS.DATA_sig{a,b}(:,1:2), 'rows');
        results.SP.DATA_sig{a,b} = results.SP.DATA{a,b}(idx_A,:);
    end
end

% Find all and significant number of cells divided by 'experimental_group' and 'compound'
numbers.cells.HPS.all = cell(2,2);
numbers.cells.HPS.sig = cell(2,2);
numbers.cells.HPS.all = cellfun(@(x) size(x,1), results.HPS.DATA);
numbers.cells.HPS.sig = cellfun(@(x) size(x,1), results.HPS.DATA_sig);

clearvars -except compound_group_names dataset_animal_names DISTRIBUTIONS ...
                  experiment_group_names GC numbers results STATS_response_modulation
%% Values for baseline stability in SP and HPS

% Results DATA structure and groups ID
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  1   | 2 | %
    % sham|  3   | 4 | %
    %------------------%

% Form new structure
results.diff = struct();
results.diff.all = [];
results.diff.sig = [];

% Analysis for all cells
all_HPS = [];
all_HPS_cell = [];
all_SP = [];
all_SP_cell = [];

all_HPS = reshape(results.HPS.DATA',[4, 1]);
all_HPS_cell = cell2mat(all_HPS); % Concat all HPS data
all_SP = reshape(results.SP.DATA',[4, 1]);
all_SP_cell = cell2mat(all_SP); % Concat all SP data

all_groups = [ones(size(all_HPS{1,1},1),1); 2*ones(size(all_HPS{2,1},1),1); ...
            3*ones(size(all_HPS{3,1},1),1); 4*ones(size(all_HPS{4,1},1),1)]; % Keep track of experimental groups
all_HPS_diff = all_HPS_cell(:, 4) - all_HPS_cell(:, 3);
all_SP_diff = all_SP_cell(:, 4) - all_SP_cell(:, 3);
results.diff.all = [all_HPS_cell(:, 1:2) all_groups all_HPS_diff all_SP_diff];

% Analysis for significant cells
sig_HPS = [];
sig_HPS_cell = [];
sig_SP = [];
sig_SP_all = [];

sig_HPS = reshape(results.HPS.DATA_sig',[4, 1]);
sig_HPS_cell = cell2mat(sig_HPS); % Concat significant HPS data
sig_SP = reshape(results.SP.DATA_sig',[4, 1]);
sig_SP_all = cell2mat(sig_SP); % Concat significant SP data

sig_groups = [ones(size(sig_HPS{1,1},1),1); 2*ones(size(sig_HPS{2,1},1),1); ...
            3*ones(size(sig_HPS{3,1},1),1); 4*ones(size(sig_HPS{4,1},1),1)]; % Keep track of experimental groups
sig_HPS_diff = sig_HPS_cell(:, 4) - sig_HPS_cell(:, 3);
sig_SP_diff = sig_SP_all(:, 4) - sig_SP_all(:, 3);
results.diff.sig = [sig_HPS_cell(:, 1:2) sig_groups sig_HPS_diff sig_SP_diff];

%______________________________________________
%|	  1    |   2   |     3    |   4   |   5    |
%|Animal_ID|Cell_ID|exp_groups|Diff_SP|Diff_HPS|

clearvars -except compound_group_names dataset_animal_names DISTRIBUTIONS ...
                  experiment_group_names GC numbers results STATS_response_modulation   
%% Define data groups to plot

% Results DATA structure and groups ID
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  1   | 2 | %
    % sham|  3   | 4 | %
    %------------------%

% Specify the group from where data is used (results.HPS.DATA / results.HPS.DATA_sig_st)
state = 'HPS'; % 'SP' 'diff'
ds_group = 'DATA_sig'; % DATA 'all' 'sig'

clearvars -except compound_group_names datasets dataset_animal_names DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group numbers results state ...
            STATS_response_modulation
%% Define specific datasets

% Datasets structure
%__________________________________________________________________________
% datasets:            [1x1] struct
%   |--all:            [1x1] struct
%   |  |--CCIsham:     [1x2] cell
%   |  |--CCI:         [1x3] cell
%   |  |--sham:        [1x3] cell
%   |  '--All:         [1x6] cell
%   |--significant:    [1x1] struct
%   |  |--CCIsham:     [1x2] cell
%   |  |--CCI:         [1x3] cell
%   |  |--sham:        [1x3] cell
%   |  '--All:         [1x6] cell
%   '--baseline:       [1x1] struct
%      |--all:         [1x1] struct
%      |  |--all:      [1010x3] double
%      |  |--HPS:      [2x2] cell
%      |  |--SP:       [2x2] cell
%      |  |--HPS_CCI:  [641x1] double
%      |  |--HPS_sham: [369x1] double
%      |  |--SP_CCI:   [641x1] double
%      |  '--SP_sham:  [369x1] double
%      '--significant: [1x1] struct
%         |--all:      [330x3] double
%         |--HPS:      [2x2] cell
%         |--SP:       [2x2] cell
%         |--HPS_CCI:  [210x1] double
%         |--HPS_sham: [120x1] double
%         |--SP_CCI:   [210x1] double
%         '--SP_sham:  [120x1] double
%__________________________________________________________________________


% Define the dataset to plot
datasets.DATA_CCI_PRE = [results.(state).(ds_group){1,1}(:,3); results.(state).(ds_group){1,2}(:,3)]; % AUROC all CCI PRE treatment
datasets.DATA_sham_PRE = [results.(state).(ds_group){2,1}(:,3); results.(state).(ds_group){2,2}(:,3)]; % AUROC all sham PRE treatment

datasets.all = struct();
datasets.all.CCIsham = {datasets.DATA_CCI_PRE, datasets.DATA_sham_PRE};
datasets.all.CCI = {datasets.DATA_CCI_PRE, results.(state).(ds_group){1,1}(:,4), results.(state).(ds_group){1,2}(:,4)};
datasets.all.sham = {datasets.DATA_sham_PRE, results.(state).(ds_group){2,1}(:,4), results.(state).(ds_group){2,2}(:,4)};
datasets.all.all = {datasets.DATA_CCI_PRE, results.(state).(ds_group){1,1}(:,4), results.(state).(ds_group){1,2}(:,4), ...
                    datasets.DATA_sham_PRE, results.(state).(ds_group){2,1}(:,4), results.(state).(ds_group){2,2}(:,4)};

datasets.significant = struct();
datasets.significant.CCIsham = {datasets.DATA_CCI_PRE, datasets.DATA_sham_PRE};
datasets.significant.CCI = {datasets.DATA_CCI_PRE, results.(state).(ds_group){1,1}(:,4), results.(state).(ds_group){1,2}(:,4)};
datasets.significant.sham = {datasets.DATA_sham_PRE, results.(state).(ds_group){2,1}(:,4), results.(state).(ds_group){2,2}(:,4)};
datasets.significant.all = {datasets.DATA_CCI_PRE, results.(state).(ds_group){1,1}(:,4), results.(state).(ds_group){1,2}(:,4), ...
                    datasets.DATA_sham_PRE, results.(state).(ds_group){2,1}(:,4), results.(state).(ds_group){2,2}(:,4)};                

% Dataset for baseline CCI, sham, HPS, SP
datasets.baseline = struct();
datasets.baseline.all = struct();
datasets.baseline.significant = struct();

% All
a = reshape(results.HPS.DATA',[4,1]);
a1 = cell2mat(a);
HPS_baseline = a1(:,3);
b = reshape(results.SP.DATA',[4,1]);
b1 = cell2mat(b);
SP_baseline = b1(:,3);
groups = [ones(numbers.cells.HPS.all(1,1),1); 2*ones(numbers.cells.HPS.all(1,2),1); ...
        3*ones(numbers.cells.HPS.all(2,1),1); 4*ones(numbers.cells.HPS.all(2,2),1)];
datasets.baseline.all.all = [HPS_baseline SP_baseline groups];

datasets.baseline.all.HPS = cell(2,2);
datasets.baseline.all.SP = cell(2,2);

datasets.baseline.all.HPS{1,1} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==1, 1);
datasets.baseline.all.HPS{1,2} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==2, 1);
datasets.baseline.all.HPS{2,1} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==3, 1);
datasets.baseline.all.HPS{2,2} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==4, 1);

datasets.baseline.all.SP{1,1} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==1, 2);
datasets.baseline.all.SP{1,2} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==2, 2);
datasets.baseline.all.SP{2,1} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==3, 2);
datasets.baseline.all.SP{2,2} = datasets.baseline.all.all(datasets.baseline.all.all(:,3)==4, 2);

% Dataset for CCI and sham, HPS, SP
datasets.baseline.all.HPS_CCI = [datasets.baseline.all.all(datasets.baseline.all.all(:,3)==1, 1); ...
                            datasets.baseline.all.all(datasets.baseline.all.all(:,3)==2, 1)];
datasets.baseline.all.HPS_sham = [datasets.baseline.all.all(datasets.baseline.all.all(:,3)==3, 1); ...
                            datasets.baseline.all.all(datasets.baseline.all.all(:,3)==4, 1)];

datasets.baseline.all.SP_CCI = [datasets.baseline.all.all(datasets.baseline.all.all(:,3)==1, 2); ...
                            datasets.baseline.all.all(datasets.baseline.all.all(:,3)==2, 2)];
datasets.baseline.all.SP_sham = [datasets.baseline.all.all(datasets.baseline.all.all(:,3)==3, 2); ...
                            datasets.baseline.all.all(datasets.baseline.all.all(:,3)==4, 2)];

% Significant
a = [];
a1 = [];
b = [];
b1 = [];

a = reshape(results.HPS.DATA_sig',[4,1]);
a1 = cell2mat(a);
HPS_baseline_sig = a1(:,3);
b = reshape(results.SP.DATA_sig',[4,1]);
b1 = cell2mat(b);
SP_baseline_sig = b1(:,3);
groups_sig = [ones(size(a{1,1},1),1); 2*ones(size(a{2,1},1),1); ...
            3*ones(size(a{3,1},1),1); 4*ones(size(a{4,1},1),1)];
datasets.baseline.significant.all = [HPS_baseline_sig SP_baseline_sig groups_sig];

datasets.baseline.significant.HPS = cell(2,2);
datasets.baseline.significant.SP = cell(2,2);

datasets.baseline.significant.HPS{1,1} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==1, 1);
datasets.baseline.significant.HPS{1,2} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==2, 1);
datasets.baseline.significant.HPS{2,1} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==3, 1);
datasets.baseline.significant.HPS{2,2} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==4, 1);

datasets.baseline.significant.SP{1,1} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==1, 2);
datasets.baseline.significant.SP{1,2} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==2, 2);
datasets.baseline.significant.SP{2,1} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==3, 2);
datasets.baseline.significant.SP{2,2} = datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==4, 2);

% Dataset for CCI and sham, HPS, SP
datasets.baseline.significant.HPS_CCI = [datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==1, 1); ...
                            datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==2, 1)];
datasets.baseline.significant.HPS_sham = [datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==3, 1); ...
                            datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==4, 1)];

datasets.baseline.significant.SP_CCI = [datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==1, 2); ...
                            datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==2, 2)];
datasets.baseline.significant.SP_sham = [datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==3, 2); ...
                            datasets.baseline.significant.all(datasets.baseline.significant.all(:,3)==4, 2)];
                        
clearvars -except compound_group_names datasets dataset_animal_names DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group numbers results state ...
            STATS_response_modulation
%% Find SP, SP_norm, define animal_group, animal_ID, sessions

nel = sum(cell2mat(numbers.animals(:)));
ID = {fieldnames(DISTRIBUTIONS)};

% Find max number of imaging session and tot number of cells
sessions_cells_id = [];
for i = 1:nel
	animID = cell2mat(ID{1,1}(i));
	animal_session_cell_id = [i, size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),2), ...
		size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1)];
	sessions_cells_id = cat(1, sessions_cells_id, animal_session_cell_id);
end
%_________________________________________________________
%|	  1    |               2             |      3        |
% animal_ID|number of SP imaging sessions|number of cells|

max_sessions = max(sessions_cells_id(:,2),[], 1);
max_cells = sum(sessions_cells_id(:,3));

% Concatenate all SP/SP_normalized values for all animals
SP = [];
SP_tot = [];
SP_tot2 = [];
SP_tot3 = [];
normSP = [];
normSP_tot = [];

for i = 1:nel
	animID = cell2mat(ID{1,1}(i));

	% Spontaneous
	SP = NaN(size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1), max_sessions); % Preallocate matrix with NaN to load the SP values
	SP1 = cell2mat(DISTRIBUTIONS.(animID){:, 'mean_activity'}); % Select the SP values
	SP(:,[1:size(SP1,2)]) = SP1; % Substitute the SP values in the preallocated NaN matrix

	% Spontaneous normalized to the baseline (difference)
	normSP = NaN(size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1), max_sessions);
	normSP1 = cell2mat(DISTRIBUTIONS.(animID){:, 'mean_activity'});
	normSP(:,[1:size(normSP1,2)]) = normSP1 - normSP1(:,1); % Normalize activity on spontaneous baseline

	% Add animal ID
	animID = nan(size(SP, 1), 1);
	animID(:) = i;
	SP = [animID SP];
	normSP = [animID normSP];

	% Concatenate all the cells for all the animals
	SP_tot = cat(1,SP_tot, SP);
	normSP_tot = cat(1,normSP_tot, normSP);
end

% Order and find animal_ID values based on 'experimental_group'/'compound'
a = dataset_animal_names{ismember(dataset_animal_names.experimental_group, 'CCI') & ...
						ismember(dataset_animal_names.compound, 'LP-211'), 'animal_ID'}; % Find animal_ID which contains CCI AND LP-211
dataset_animal_names = natsortrows(dataset_animal_names, 1); % Sort rows based on ascending column 1

% Find same group / compound
%------------------%
%     |LP-211|SAL| %
%  CCI|  1   | 2 | %
% sham|  3   | 4 | %
%------------------%

group1 = [];
group = [];

for i = 1:nel % Assign group_ID
	if ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'LP-211')
		group1 = 1;
	elseif ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'SAL')
		group1 = 2;
	elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'LP-211')
		group1 = 3;
	elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'SAL')
		group1 = 4;
	end
	group = [group; group1];
end

TabGroup = [];
TabGroup = table(group,'VariableNames', {'group'});
dataset_animal_names = [dataset_animal_names TabGroup]; % Copy the group into the dataset

% Find group ID for each animal
animGroup = nan(size(SP_tot,1), 1);
for i = 1:length(SP_tot)
	animGroup(i,1) = dataset_animal_names.group(SP_tot(i,1));
end

animInfo = [SP_tot(:,1) animGroup]; % Cat animal_ID and group
SP_tot2 = [animInfo SP_tot(:,2:end)]; % Cat animal_ID and group in the SP file
normSP_tot = [animInfo normSP_tot(:,2:end)]; % Cat animal_ID and group the normalized SP file

% Find max number of imaging session per animal
SP_diff_tot = [];
evoked_diff_tot = [];
evoked_pre = [];
evoked_post = [];

for i = 1:nel
	animID = cell2mat(ID{1,1}(i));
	ROI_id = DISTRIBUTIONS.(animID){:,'ROI_id'};
	length(ROI_id);

	% Assign unique group value to surgery / treatment combination
	dataset_animal_names = natsortrows(dataset_animal_names, 1); % Sort rows based on ascending column 1
	if ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'LP-211')
		group1 = 1;
	elseif ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'SAL')
		group1 = 2;
	elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'LP-211')
		group1 = 3;
	elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'SAL')
		group1 = 4;
	end

	for ii = 1:length(ROI_id)
		cel_id = DISTRIBUTIONS.(animID){ii,'ROI_id'};
		
		% Calculate SP difference: POST - PRE and assign ID if positive or
		% negative
		SP = (cell2mat(DISTRIBUTIONS.(animID){cel_id,'mean_activity'}));
		SP_diff = SP(2) - SP(1);
		if SP_diff >= 0 % Increased activity after treatment
			SP_diff(1,2) = 1;
		else
			SP_diff(1,2) = -1; % Decreased activity after treatment
		end
		SP_diff = [SP_diff group1]; % Add animal group based on surgery / treatment
		
		% Calculate EVOKED difference: POST - PRE
		evoked_pre = STATS_response_modulation.(animID)((ismember(STATS_response_modulation.(animID){:,'ROI_id'}, cel_id) & ...
			ismember(STATS_response_modulation.(animID){:,'compound_phase'}, 'baseline')),:);
		evoked_pre = mean(cell2mat(evoked_pre{:,'activity'}));
		if isnan(evoked_pre)
			evoked_pre = 0;
		else
			evoked_pre = evoked_pre(1,2);
		end
		
		evoked_post = STATS_response_modulation.(animID)((ismember(STATS_response_modulation.(animID){:,'ROI_id'}, cel_id) & ...
			ismember(STATS_response_modulation.(animID){:,'compound_phase'}, 'active_1')),:);
		evoked_post = mean(cell2mat(evoked_post{:,'activity'}));
		
		if isnan(evoked_post)
			evoked_post = 0;
		else
			evoked_post = evoked_post(1,2);
		end
		evoked_diff = evoked_post - evoked_pre;
		evoked_diff = [evoked_diff group1]; % Add animal group based on surgery / treatment
		
		
		SP_diff_tot = cat(1, SP_diff_tot, SP_diff);
		evoked_diff_tot = cat(1, evoked_diff_tot, evoked_diff);
	end
end
diff_tot = cat(2, SP_diff_tot(:,1), evoked_diff_tot, SP_diff_tot(:,2)); % Create single variable containing SP_diff, HPS_diff, GROUP info for each cell

sessions_cells_id = array2table(sessions_cells_id, 'Variablenames', {'animal_ID', 'nsessions_SP', 'ncells'});
dataset_animal_names = [dataset_animal_names sessions_cells_id(:,2:3)];

clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_tot SP_tot2 STATS_response_modulation TabGroup
%% Collect data for ANOVA, save output as csv file
% Anova test script: 'N:\Alberto\Matlab Scripts Collection\ANOVArm.m'
% from Paolo De Luna
%______________________________________________________
%|       1        |    2    |    3    |   4   |   5   |
% Cell SP activity|animal_ID|treatment|surgery|session|

SP_ANOVA = [];
SP_ANOVA_animID = [];
SP_ANOVA_groupID = [];

% SP_ANOVA = reshape(SP_tot2(:,3:end), [], 1);
a = [];
b = [];
ID = [];
SP_ANOVA_tot = [];

for i = 1:size(dataset_animal_names, 1)
    a = dataset_animal_names.nsessions_SP(i);
    b = SP_tot2(SP_tot2(:, 1) == i, 3:2+a);
    ID = reshape(b, [], 1);
    SP_ANOVA_tot = [SP_ANOVA_tot; ID];
end
SP_ANOVA = SP_ANOVA_tot;

% Find animal_ID for ANOVA table
a = [];
b = [];
ID = [];
animal_ID_tot = [];
for i = 1:size(dataset_animal_names, 1)
    a = dataset_animal_names.nsessions_SP(i);
    b = repmat(a, [dataset_animal_names.ncells(i) 1]);
    ID = repmat(b, dataset_animal_names.nsessions_SP(i), 1);
    animal_ID_tot = [animal_ID_tot; ID];
end
SP_ANOVA_animID = animal_ID_tot;

% Find group_ID for ANOVA table
a = [];
b = [];
ID = [];
group_ID_tot = [];
for i = 1:size(dataset_animal_names, 1)
    a = dataset_animal_names.group(i);
    b = (repmat(a, [dataset_animal_names.ncells(i) 1]));
    ID = repmat(b, dataset_animal_names.nsessions_SP(i), 1);
    group_ID_tot = [group_ID_tot; ID];
end
SP_ANOVA_groupID = group_ID_tot;
SP_ANOVA_groupID = num2cell(SP_ANOVA_groupID);

for i = 1:length(SP_ANOVA_groupID)
    if isequal(SP_ANOVA_groupID(i,1), {1})
        SP_ANOVA_groupID(i, 2) = cellstr('CCI');
        SP_ANOVA_groupID(i, 3) = cellstr('LP-211');
    elseif isequal(SP_ANOVA_groupID(i,1), {2})
        SP_ANOVA_groupID(i, 2) = cellstr('CCI');
        SP_ANOVA_groupID(i, 3) = cellstr('SAL');
    elseif isequal(SP_ANOVA_groupID(i,1), {3})
        SP_ANOVA_groupID(i, 2) = cellstr('sham');
        SP_ANOVA_groupID(i, 3) = cellstr('LP-211');
    elseif isequal(SP_ANOVA_groupID(i,1), {4})
        SP_ANOVA_groupID(i, 2) = cellstr('sham');
        SP_ANOVA_groupID(i, 3) = cellstr('SAL');
    end
end

ID = [];
ID2 = [];
ID_tot = [];
ID_tot2 = [];

for i = 1:size(dataset_animal_names, 1)
    for ii = 1:dataset_animal_names.nsessions_SP(i);
        if ii == 1
            a = 'baseline';
            b = 'pre';
            ID = cellstr(repmat(a, [dataset_animal_names.ncells(i) 1]));
            ID2 = cellstr(repmat(b, [dataset_animal_names.ncells(i) 1]));
        elseif ii == 2
            a = 'active1';
            b = 'post';
            ID = cellstr(repmat(a, [dataset_animal_names.ncells(i) 1]));
            ID2 = cellstr(repmat(b, [dataset_animal_names.ncells(i) 1]));
        elseif ii == 3
            a = 'active2';
            b = 'post';
            ID = cellstr(repmat(a, [dataset_animal_names.ncells(i) 1]));
            ID2 = cellstr(repmat(b, [dataset_animal_names.ncells(i) 1]));
        elseif ii == 4
            a = 'active3';
            b = 'post';
            ID = cellstr(repmat(a, [dataset_animal_names.ncells(i) 1]));
            ID2 = cellstr(repmat(b, [dataset_animal_names.ncells(i) 1]));
        elseif ii == 5
            a = 'active4';
            b = 'post';
            ID = cellstr(repmat(a, [dataset_animal_names.ncells(i) 1]));
            ID2 = cellstr(repmat(b, [dataset_animal_names.ncells(i) 1]));
        end
        ID_tot = [ID_tot; ID];
        ID_tot2 = [ID_tot2; ID2];
    end
end

% Concatenate ID values to get vector of (5190, 1)
SP_ANOVA = num2cell(SP_ANOVA);
SP_ANOVA_animID = num2cell(SP_ANOVA_animID);

% Find unique cell_ID
a = [];
b = [];
c = [];
for i = 1:size(dataset_animal_names, 1)
    if i == 1
    a = (1:dataset_animal_names.ncells(i))';
    b = a;
    else
    b = (max(c) + 1:(max(c) + dataset_animal_names.ncells(i)))';
    end
    b = repmat(b, [dataset_animal_names.nsessions_SP(i) 1]);
    c = [c; b];
end

cell_ID = c;
cell_ID = num2cell(cell_ID);

% Concatenate
SP_ANOVA = [SP_ANOVA cell_ID SP_ANOVA_animID SP_ANOVA_groupID(:,2:3) ID_tot2 ID_tot];

SP_ANOVA_NaN = SP_ANOVA; % Cell with all values
rows = isnan(cell2mat(SP_ANOVA(:,1))); % Cell without NaNs
SP_ANOVA(rows, :) = [];

% Copy to table structure
SP_ANOVA2 = SP_ANOVA;
SP_ANOVA2 = cell2table(SP_ANOVA2, 'VariableNames',{'cell_SP', 'cell_ID', 'animal_ID', 'surgery', 'treatment', 'session', 'single_session'});
SP_ANOVA_NaN2 = SP_ANOVA_NaN;
SP_ANOVA_NaN2 = cell2table(SP_ANOVA_NaN2, 'VariableNames',{'cell_SP', 'cell_ID', 'animal_ID', 'surgery', 'treatment', 'session', 'single_session'});

% cd('N:\Alberto\Analysis\InVivo\2PLM\R ks density comparison groups');
% writetable(SP_ANOVA2, 'SP_ANOVA.csv');
% writetable(SP_ANOVA_NaN2, 'SP_ANOVA_NaN.csv');

clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_tot SP_tot2 ...
            STATS_response_modulation
%% Extract SP activity from significant modulated cells

% Adapt DISTRIBUTION structure erasing exceeding SP sessions from some animals
DISTRIBUTIONS2 = DISTRIBUTIONS;
DISTRIBUTIONS2.AB_10(:,'amplitude_dist_after_3') = [];
DISTRIBUTIONS2.AB_10(:,'interval_dist_after_3') = [];
DISTRIBUTIONS2.AB_10.change_after(:,end) = [];
DISTRIBUTIONS2.AB_13(:,'amplitude_dist_after_3') = [];
DISTRIBUTIONS2.AB_13(:,'amplitude_dist_after_4') = [];
DISTRIBUTIONS2.AB_13(:,'interval_dist_after_3') = [];
DISTRIBUTIONS2.AB_13(:,'interval_dist_after_4') = [];
DISTRIBUTIONS2.AB_13.change_after(:,3:4) = [];
DISTRIBUTIONS2.AB_14(:,'amplitude_dist_after_3') = [];
DISTRIBUTIONS2.AB_14(:,'interval_dist_after_3') = [];
DISTRIBUTIONS2.AB_14.change_after(:,3) = [];
DISTRIBUTIONS2.AB_15(:,'amplitude_dist_after_3') = [];
DISTRIBUTIONS2.AB_15(:,'interval_dist_after_3') = [];
DISTRIBUTIONS2.AB_15.change_after(:,3) = [];

% Create new structure to contain the SP_all activity
SP_all = cell(2, 2);
a = fieldnames(DISTRIBUTIONS2);

for i = 1:height(dataset_animal_names)
    gr = dataset_animal_names.group(i);
    an_ID = [];
    
    if gr == 1
        an_ID = i.*ones(height(DISTRIBUTIONS2.(a{i})), 1);
        SP_all{1, 1} = [SP_all{1, 1}; [array2table(an_ID,'VariableNames', {'anim_ID'})  DISTRIBUTIONS2.(a{i})]];
    elseif gr == 2
        an_ID = i.*ones(height(DISTRIBUTIONS2.(a{i})), 1);
        SP_all{1, 2} = [SP_all{1, 2}; [array2table(an_ID,'VariableNames', {'anim_ID'})  DISTRIBUTIONS2.(a{i})]];
    elseif gr == 3
        an_ID = i.*ones(height(DISTRIBUTIONS2.(a{i})), 1);
        SP_all{2, 1} = [SP_all{2, 1}; [array2table(an_ID,'VariableNames', {'anim_ID'})  DISTRIBUTIONS2.(a{i})]];
    elseif gr == 4
        an_ID = i.*ones(height(DISTRIBUTIONS2.(a{i})), 1);
        SP_all{2, 2} = [SP_all{2, 2}; [array2table(an_ID,'VariableNames', {'anim_ID'})  DISTRIBUTIONS2.(a{i})]];
    end
end

all_cellID = cell(2, 2);
for i = 1:2
    for ii = 1:2
        a = results.HPS.DATA_p{i, ii}(:, 6) + results.HPS.DATA_p{i, ii}(:, 7);
        b = a > 0;
        all_cellID{i, ii} = [results.HPS.DATA_p{i, ii}(:,1:2) b];
    end
end

% Create new structure to contain the SP_sig interval distribution activity

% sig_all = cell(2, 2);
% for i = 1:2
%     for ii = 1:2
%         a = results.HPS.DATA_p{i, ii}(:, 6) + results.HPS.DATA_p{i, ii}(:, 7);
%         b = a > 0;
%         sig_all{i, ii} = [results.HPS.DATA_p{i, ii}(:,1:2) b];
%     end
% end

SP_int = cell(2, 2);
for i = 1:2
    for ii = 1:2
        idx = [];
        idx = ismember(table2array(SP_all{i, ii}(:, 1:2)), all_cellID{i, ii}(:, 1:2), 'rows');
        SP_int{i, ii} = SP_all{i, ii}(idx, 9:11);
    end
end

SP_amp = cell(2, 2);
for i = 1:2
    for ii = 1:2
        idx = [];
        idx = ismember(table2array(SP_all{i, ii}(:, 1:2)), all_cellID{i, ii}(:, 1:2), 'rows');
        SP_amp{i, ii} = SP_all{i, ii}(idx, 6:8);
    end
end

% sig_cellID = cell(2, 2);
% for i = 1:2
%     for ii = 1:2
%         idx = [];
%         idx = ismember(table2array(SP_all{i, ii}(:, 1:2)), all_cellID{i, ii}(:, 1:2), 'rows');
%         idx = sig{i, ii}(any(sig{i, ii}(:,3), 3),:);
%         sig_cellID{i, ii} = SP_all{i, ii}(idx, 9:11)
%     end
% end

SP_int_sig = cell(2, 2);
for i = 1:2
    for ii = 1:2
        idx = [];
        idx = all_cellID{i, ii}(any(all_cellID{i, ii}(:,3), 3),:);
        idx = ismember(table2array(SP_all{i, ii}(:, 1:2)), idx(:, 1:2), 'rows');
        SP_int_sig{i, ii} = SP_all{i, ii}(idx, 9:11);
    end
end

SP_amp_sig = cell(2, 2);
for i = 1:2
    for ii = 1:2
        idx = [];
        idx = all_cellID{i, ii}(any(all_cellID{i, ii}(:,3), 3),:);
        idx = ismember(table2array(SP_all{i, ii}(:, 1:2)), idx(:, 1:2), 'rows');
        SP_amp_sig{i, ii} = SP_all{i, ii}(idx, 6:8);
    end
end

numbers.cells.SP.sig = cell(2, 2);
numbers.cells.SP.all = cell(2, 2);
numbers.cells.SP.sig = cellfun(@(x) size(x,1), SP_int_sig);
numbers.cells.SP.all = cellfun(@(x) size(x,1), SP_all);

clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            DISTRIBUTIONS2 experiment_group_names GC datasets ds_group normSP_tot numbers ...
            results state SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_all SP_sig SP_sig_int SP_tot ...
            SP_tot2 STATS_response_modulation
        
%%
%------------%
% STATISTICS %
%------------%
%% STATS 1: Wilcoxon statistics on the different conditions
% Wilcoxon rank sum test for unpaired distribution (CCI_PRE/sham_PRE)
% Wilcoxon signed rank test for paired distribution (CCI_PRE/LP | sham_PRE/LP)
numbers.stats = struct();

p_sign = [];
p_rank = []

% CCI_PRE vs sham_PRE
x = datasets.DATA_CCI_PRE; % CCI_LP PRE + CCI_SAL PRE
y = datasets.DATA_sham_PRE; % sham_LP PRE + sham_SAL PRE
p1 = ranksum(x,y);
p_sign = cat(1,p_sign, p1);
p_rank = cat(1,p_sign,p_rank);

% CCI_LP PRE vs CCI_LP POST
x = results.(state).(ds_group){1,1}(:,3); % CCI_LP PRE
y = results.(state).(ds_group){1,1}(:,4); % CCI_LP POST
p1 = signrank(x,y);
p_sign = cat(1,p_sign,p1);

% CCI_SAL PRE vs CCI_SAL POST
x = results.(state).(ds_group){1,2}(:,3); % CCI_SAL PRE
y = results.(state).(ds_group){1,2}(:,4); % CCI_SAL POST
p1 = signrank(x,y);
p_sign = cat(1,p_sign,p1);

% CCI_PRE vs CCI_LP POST
x = datasets.DATA_CCI_PRE;
y = results.(state).(ds_group){1,1}(:,4); % CCI_LP POST
p1 = ranksum(x,y);
p_rank = cat(1,p_rank,p1);

% sham_LP PRE vs sham_LP POST
x = results.(state).(ds_group){2,1}(:,3); % sham_LP PRE
y = results.(state).(ds_group){2,1}(:,4); % sham_LP POST
p1 = signrank(x,y);
p_sign = cat(1,p_sign,p1);

% sham_PRE vs sham_LP POST
x = datasets.DATA_sham_PRE;
y = results.(state).(ds_group){2,1}(:,4); % sham_LP POST
p1 = ranksum(x,y);
p_rank = cat(1,p_rank,p1);

% All cells baselines per treatment groups
%-----------------------------------------
numbers.stats.baseline.all = struct();

all_HPS = [];
all_HPS_cell = [];
all_SP = [];
all_SP_cell = [];

all_HPS = reshape(results.HPS.DATA',[4, 1]);
all_HPS_cell = cell2mat(all_HPS); % Concat all HPS data
all_SP = reshape(results.SP.DATA',[4, 1]);
all_SP_cell = cell2mat(all_SP); % Concat all SP data
all_groups = [ones(size(all_HPS{1,1},1),1); 2*ones(size(all_HPS{2,1},1),1); ...
            3*ones(size(all_HPS{3,1},1),1); 4*ones(size(all_HPS{4,1},1),1)]; % Keep track of experimental groups

numbers.stats.baseline.all.baseline_HPS_SP = [all_HPS_cell(:, 1:2) all_groups all_SP_cell(:,3) all_HPS_cell(:,3)] % Find baseline cell activity

%_______________________________________________________
%|	  1    |   2   |     3    |     4     |     5      |
%|Animal_ID|Cell_ID|exp_groups|Baseline_SP|Baseline_HPS|

% SP baseline

numbers.stats.baseline.all.p_baselineSP_CCI = ranksum((numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 1, 4)), ... % CCI LP-211
                (numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 2, 4))); % CCI SAL
numbers.stats.baseline.all.p_baselineSP_sham = ranksum((numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 3, 4)), ... % sham LP-211
                (numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 4, 4))); % sham SAL

numbers.stats.baseline.all.p_baselineSP_CCIsham = ranksum(numbers.stats.baseline.all.baseline_HPS_SP((numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 1 | numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 2), 4), ... % All CCI
                numbers.stats.baseline.all.baseline_HPS_SP((numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 3 | numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 4), 4)); % All sham
            

% HPS baseline            
numbers.stats.baseline.all.p_baselineHPS_CCI = ranksum((numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 1, 5)), ... % CCI LP-211
                (numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 2, 5))); % CCI SAL
numbers.stats.baseline.all.p_baselineHPS_sham = ranksum((numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 3, 5)), ... % sham LP-211
                (numbers.stats.baseline.all.baseline_HPS_SP(numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 4, 5))); % sham SAL
            
numbers.stats.baseline.all.p_baselineHPS_CCIsham = ranksum(numbers.stats.baseline.all.baseline_HPS_SP((numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 1 | numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 2), 5), ... % All CCI
                numbers.stats.baseline.all.baseline_HPS_SP((numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 3 | numbers.stats.baseline.all.baseline_HPS_SP(:,3) == 4), 5)); % All sham


% Significant cells baselines per treatment groups
%-------------------------------------------------
numbers.stats.baseline.sig = struct();

sig_HPS = [];
sig_HPS_cell = [];
sig_SP = [];
sig_SP_cell = [];

sig_HPS = reshape(results.HPS.DATA_sig',[4, 1]);
sig_HPS_cell = cell2mat(sig_HPS); % Concat sig HPS data
sig_SP = reshape(results.SP.DATA_sig',[4, 1]);
sig_SP_cell = cell2mat(sig_SP); % Concat sig SP data
sig_groups = [ones(size(sig_HPS{1,1},1),1); 2*ones(size(sig_HPS{2,1},1),1); ...
            3*ones(size(sig_HPS{3,1},1),1); 4*ones(size(sig_HPS{4,1},1),1)]; % Keep track of experimental groups

numbers.stats.baseline.sig.baseline_HPS_SP = [sig_HPS_cell(:, 1:2) sig_groups sig_SP_cell(:,3) sig_HPS_cell(:,3)] % Find baseline cell activity

%_______________________________________________________
%|	  1    |   2   |     3    |     4     |     5      |
%|Animal_ID|Cell_ID|exp_groups|Baseline_SP|Baseline_HPS|

% Sig SP baseline

numbers.stats.baseline.sig.p_baselineSP_CCI = ranksum((numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 1, 4)), ... % CCI LP-211
                (numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 2, 4))); % CCI SAL
numbers.stats.baseline.sig.p_baselineSP_sham = ranksum((numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 3, 4)), ... % sham LP-211
                (numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 4, 4))); % sham SAL

numbers.stats.baseline.sig.p_baselineSP_CCIsham = ranksum(numbers.stats.baseline.sig.baseline_HPS_SP((numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 1 | numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 2), 4), ... % All CCI
                numbers.stats.baseline.sig.baseline_HPS_SP((numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 3 | numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 4), 4)); % All sham
            

% Sig HPS baseline            
numbers.stats.baseline.sig.p_baselineHPS_CCI = ranksum((numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 1, 5)), ... % CCI LP-211
                (numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 2, 5))); % CCI SAL
numbers.stats.baseline.sig.p_baselineHPS_sham = ranksum((numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 3, 5)), ... % sham LP-211
                (numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 4, 5))); % sham SAL
            
numbers.stats.baseline.sig.p_baselineHPS_CCIsham = ranksum(numbers.stats.baseline.sig.baseline_HPS_SP((numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 1 | numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 2), 5), ... % All CCI
                numbers.stats.baseline.sig.baseline_HPS_SP((numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 3 | numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 4), 5)); % All sham
            
            
clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_tot SP_tot2 ...
            STATS_response_modulation 
%% STATS 2: Find mean/median AUROC values for each animal group based on significant / >=0 values
numbers.stats.general.all = struct();
numbers.stats.general.all.PrePost = struct();
numbers.stats.general.all.Post = struct();
numbers.stats.general.all.Pre = struct();
numbers.stats.general.all.Alw = struct();
numbers.stats.general.all.Nev = struct();

numbers.stats.general.all.PrePost.median = cell(2,2);
numbers.stats.general.all.PrePost.mean = cell(2,2);
numbers.stats.general.all.Post.median = cell(2,2);
numbers.stats.general.all.Post.mean = cell(2,2);
numbers.stats.general.all.Pre.median = cell(2,2);
numbers.stats.general.all.Pre.mean = cell(2,2);
numbers.stats.general.all.Alw.median = cell(2,2);
numbers.stats.general.all.Alw.mean = cell(2,2);
numbers.stats.general.all.Nev.median = cell(2,2);
numbers.stats.general.all.Nev.mean = cell(2,2);

numbers.stats.general.positiveAUROC = numbers.stats.general.all;

a = [];
b = [];
b1 = [];
c = [];
c1 = [];
d = [];

% Cells AUROC >= 0.5, significant PRE or POST
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 | results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
		c = b(any(b(:,3),3),:); % Keep only 1
		c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
		c1(:,3) = c1(:,1).*c1(:,2); % Find significant value PRE and POST
		%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
		d = c(any(c1(:,3),3),:); % Keep only 1
		d(:,3) = [];
		
		[~,p] = ttest(d(:, 1), d(:, 2));
		signr = signrank(d(:, 1), d(:, 2)); % Paired
		% ranksum(d(:,1), c(:,1))
		numbers.stats.general.positiveAUROC.PrePost.mean{i,ii} = mean(d);
		numbers.stats.general.positiveAUROC.PrePost.median{i,ii} = median(d);
		%         diff(median(d))
		%         diff(mean(d))
		med = median(d);
		me = mean(d);
		meddif = median(diff(d,1,2));
		medif = mean(diff(d,1,2));
		
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells AUROC >= 0.5, significant only during POST not PRE
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during POST
		c = b(any(b(:,3),3),:); % Keep only 1
		c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5 in PRE or POST
		c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
		%         c1(:,3) = c1(:,1)+c1(:,2); % Find any AUROC value >= 0.5 PRE or POST
		d = c(any(c1(:,3),3),:); % Keep only 1
		d(:,3) = [];
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.positiveAUROC.Post.mean{i,ii} = NaN;
			numbers.stats.general.positiveAUROC.Post.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.positiveAUROC.Post.mean{i,ii} = mean(d);
			numbers.stats.general.positiveAUROC.Post.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells AUROC >= 0.5, significant only during PRE not POST
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find significant values during PRE
		% b1 = results.(time).(ds_group){i,ii}(:,7)>0; % Find significant values during POST
		% b1(:,3) = b1(:,1).*b1(:,2);
		%         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
		c = b(any(b(:,3),3),:); % Keep only 1
		c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5
		c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
		%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
		d = c(any(c1(:,3),3),:); % Keep only 1
		d = c;
		d(:,3) = [];
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.positiveAUROC.Pre.mean{i,ii} = NaN;
			numbers.stats.general.positiveAUROC.Pre.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.positiveAUROC.Pre.mean{i,ii} = mean(d);
			numbers.stats.general.positiveAUROC.Pre.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells AUROC >= 0.5, always significant
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
		% b1 = results.(time).(ds_group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
		% b1(:,3) = b1(:,1).*b1(:,2);
		c = b(any(b(:,3),3),:); % Keep only 1
		c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
		c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
		d = c(any(c1(:,3),3),:); % Keep only 1
		d(:,3) = [];
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.positiveAUROC.Alw.mean{i,ii} = NaN;
			numbers.stats.general.positiveAUROC.Alw.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.positiveAUROC.Alw.mean{i,ii} = mean(d);
			numbers.stats.general.positiveAUROC.Alw.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells AUROC >= 0.5, never significant
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find never significant values
		c = b(any(b(:,3),3),:); % Keep only 1
		c1 = any(c(:,1:2)>=0.5,3);% Find only values >=0.5
		c1(:,3) = c1(:,1).*c1(:,2);
		d = c(any(c1(:,3),3),:); % Keep only 1
		d(:,3) = [];
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.positiveAUROC.Nev.mean{i,ii} = NaN;
			numbers.stats.general.positiveAUROC.Nev.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.positiveAUROC.Nev.mean{i,ii} = mean(d);
			numbers.stats.general.positiveAUROC.Nev.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end

%--------------------------------------------------------------------------
% Cells significant PRE or POST
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 | results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
		c = b(any(b(:,3),3),:); % Keep only 1
		c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
		%         c1(:,3) = c1(:,1).*c1(:,2); % Find significant value PRE and POST
		%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
		%         d = c(any(c1(:,3),3),:); % Keep only 1
		%         d(:,3) = [];
		d = c(:,1:2);
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.all.PrePost.mean{i,ii} = NaN;
			numbers.stats.general.all.PrePost.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.all.PrePost.mean{i,ii} = mean(d);
			numbers.stats.general.all.PrePost.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells significant only during POST not PRE
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during POST
		c = b(any(b(:,3),3),:); % Keep only 1
		%         c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5 in PRE or POST
		%         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
		%         c1(:,3) = c1(:,1)+c1(:,2); % Find any AUROC value >= 0.5 PRE or POST
		%         d = c(any(c1(:,3),3),:); % Keep only 1
		%         d(:,3) = [];
		d = c(:,1:2);
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.all.Post.mean{i,ii} = NaN;
			numbers.stats.general.all.Post.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.all.Post.mean{i,ii} = mean(d);
			numbers.stats.general.all.Post.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells significant only during PRE not POST
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find significant values during PRE
		% b1 = results.(time).(ds_group){i,ii}(:,7)>0; % Find significant values during POST
		% b1(:,3) = b1(:,1).*b1(:,2);
		%         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
		c = b(any(b(:,3),3),:); % Keep only 1
		%         c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5
		%         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
		%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
		%         d = c(any(c1(:,3),3),:); % Keep only 1
		%         d = c;
		%         d(:,3) = [];
		d = c(:,1:2);
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.all.Pre.mean{i,ii} = NaN;
			numbers.stats.general.all.Pre.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.all.Pre.mean{i,ii} = mean(d);
			numbers.stats.general.all.Pre.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells always significant
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
		% b1 = results.(time).(ds_group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
		% b1(:,3) = b1(:,1).*b1(:,2);
		c = b(any(b(:,3),3),:); % Keep only 1
		%         c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
		%         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
		%         d = c(any(c1(:,3),3),:); % Keep only 1
		%         d(:,3) = [];
		d = c(:,1:2);
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.all.Alw.mean{i,ii} = NaN;
			numbers.stats.general.all.Alw.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.all.Alw.mean{i,ii} = mean(d);
			numbers.stats.general.all.Alw.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end
% Cells never significant
for i = 1:2
	for ii = 1:2
		a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
		b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find never significant values
		c = b(any(b(:,3),3),:); % Keep only 1
		%         c1 = any(c(:,1:2)>=0.5,3);% Find only values >=0.5
		%         c1(:,3) = c1(:,1).*c1(:,2);
		%         d = c(any(c1(:,3),3),:); % Keep only 1
		%         d(:,3) = [];
		d = c(:,1:2);
		
		empty = isempty(d);
		if empty == 1
			numbers.stats.general.all.Nev.mean{i,ii} = NaN;
			numbers.stats.general.all.Nev.median{i,ii} = NaN;
		else
			[~,p] = ttest(d(:, 1), d(:, 2));
			signr = signrank(d(:, 1), d(:, 2)); % Paired
			% ranksum(d(:,1), c(:,1))
			numbers.stats.general.all.Nev.mean{i,ii} = mean(d);
			numbers.stats.general.all.Nev.median{i,ii} = median(d);
			%         diff(median(d))
			%         diff(mean(d))
			med = median(d);
			me = mean(d);
			meddif = median(diff(d,1,2));
			medif = mean(diff(d,1,2));
		end
		a = [];
		b = [];
		c = [];
		d = [];
	end
end


clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_tot SP_tot2 ...
            STATS_response_modulation 
%% STATS 3: R KSdensity Stats: selection of data pairs to compare

% Define the two groups to compare
a = numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 1, 5)
b = numbers.stats.baseline.sig.baseline_HPS_SP(numbers.stats.baseline.sig.baseline_HPS_SP(:,3) == 3, 5)
c = [a; b]; % Concatenate
c(:,2) = [ones(size(a)); 2*ones(size(b))]; % Add the grouping variable based on the number of items in each group
sig_baseline_HPS_CCIshamLP = array2table(c, 'VariableNames', {'data' 'group'}); % Save the variable

% csv save path
cd('N:\Alberto\Analysis\InVivo\2PLM\R ks density comparison groups');

% Save .csv file to use with R script
writetable(sig_baseline_SP_CCIsham, 'sig_baseline_SP_CCIsham.csv');
writetable(all_baseline_HPS_CCIshamLP, 'all_baseline_HPS_CCIshamLP.csv');
writetable(sig_baseline_HPS_CCIshamLP, 'sig_baseline_HPS_CCIshamLP.csv');
%% STATS 4: ANOVA test from script

data = cell2mat(SP_ANOVA(:,1));
groups = SP_ANOVA(:,2:end);

% Convert groups in boolean values for each column
[~,~,idx] = unique(groups(:,3), 'stable'); groups(:,3) = num2cell(idx);
[~,~,idx] = unique(groups(:,4), 'stable'); groups(:,4) = num2cell(idx);
[~,~,idx] = unique(groups(:,5), 'stable'); groups(:,5) = num2cell(idx);
[~,~,idx] = unique(groups(:,6), 'stable'); groups(:,6) = num2cell(idx);

% Define groups names and variables
groups = cell2mat(groups);
labels = {'cell','animal','surgery','LP_SAL','pre_post','phase'}
WITHINsubject = false(1,6);
WITHINsubject = logical([0, 0, 0, 0, 1, 1])
factorType = {'i', 'f', 'f', 'f', 'f', 'f'};

% Collect only pre-post data independently on phase order
idx = ismember(groups(:, 6), [1, 2]);
groups = groups(idx, :);
data = data(idx);

% % Remove phase column
% groups(:,end)=[];
% labels(end)=[]
% WITHINsubject(end)=[]
% factorType(end)=[]

% Run ANOVA test from script
% NB: in R open the script and remove the animal_ID variable in the aov 
% parameters to get rid of the error and allow the test to be performed on 
% the datasets. The output files (ANOVA_data.csv, ANOVA_results.csv, 
% ANOVArm_R.R) are located in:
% 
% C:\Users\bisco\AppData\Local\Temp

ANOVAtable = ANOVArm (data, groups, labels, WITHINsubject, factorType, true)
%% STATS 5: Stats of SP_tot based on grouping
SP_tot3 = array2table(SP_tot2,'VariableNames',{'animal_ID', 'animal_group', 'Baseline', 'Active1', 'Active2, 'Active3', 'Active4'});
grpstats(SP_tot3,'animal_group', 'median');

%%
%-------%
% PLOTS %
%-------%
%% Figure 1: Plot single cells based on treatment group
phases_order = {'baseline', 'active_1'};
LEG_labels = {'non-selective, non-modulated', 'selective, non-modulated', 'only selective before', 'only selective after'};

for irow = 1:2
    for icol = 1:2
        % Open figure and draw axes
        figure('color','w')
        clf; hold on;
        ax = subplot(1,1,1);
        plot([0, 1], [0, 1], 'k')
        plot([0.5, .5], [0, 1], 'k')
        plot([0, 1], [0.5, .5], 'k')
        
        % Plot each point with a different style to highlight the group they belong to
        LEG_H = NaN(1, 4);
        LEG_n = zeros(1, 4);
        % H = zeros(size(DATA,1), 1);
        H = [];
        
        for ix = 1:numbers.cells.HPS.all(irow,icol)
            x = results.HPS.DATA{irow, icol}(ix, 3);
            y = results.HPS.DATA{irow, icol}(ix, 4);
            p_pre = results.HPS.DATA_p{irow, icol}(ix, 3);
            p_post = results.HPS.DATA_p{irow, icol}(ix, 4);
            
            % Non-selective, non-modulated cells
            if p_pre >= 0.05 && p_post >= 0.05
                group = 1;
                mfc = [.7, .7, .7];
                mec = 'w';
                marker = 'o';
                markersize = 8;
                lw = 1;
                
            % Selective, non-modulated cells
            elseif p_pre < 0.05 && p_post < 0.05
                group = 2;
                mfc = 'k';
                mec = 'none';
                marker = 'p';
                markersize = 18;
                lw = 1;
                
            % Only selective before
            elseif p_pre < 0.05 && p_post >= 0.05
                group = 3;
                mfc = [.7, .7, .7];
                marker = 'o';
                markersize = 10;
                % Change color depending on sign of modulation
                if x < 0.5
                    mec = 'r';
                else
                    mec = [0, .67, 0];
                end
                lw = 2;
                
            % Only selective after
            elseif p_pre >= 0.05 && p_post < 0.05
                group = 4;
                mec = 'none';
                marker = 'p';
                markersize = 18;
                % Change color depending on sign of modulation
                if y < 0.5
                    mfc = 'r';
                else
                    mfc = [0, .67, 0];
                end
                lw = 1;
            end
            
            % Plot point
            H(ix) = plot(x, y, 'Marker',marker, 'MarkerFaceColor',mfc, 'MarkerEdgeColor',mec, 'MarkerSize',markersize, 'Linestyle','none', 'Linewidth',lw);
            
            % Store handle for legend
            if isnan(LEG_H(ds_group))
                LEG_H(ds_group) = H(ix);
            end
            LEG_n(ds_group) = LEG_n(ds_group) + 1;
        end
        
        axis square
        set(ax, 'FontSize',12, 'TickDir','out')
        xlabel(['Selectivity to  HPS before LP-211'], 'FontSize',13), ylabel(['Selectivity to HPS after LP-211'], 'FontSize',13)
        title([compound_group_names{irow},' + ',experiment_group_names{icol}, ' (n=', num2str(numbers.cells.HPS.all(irow, icol)), ')'], 'FontSize',16, 'Interpreter','none')
        
        % Add legend
        idx = isnan(LEG_H);
        LEG_H(idx) = [];
        leg_labels = LEG_labels(~idx);
        LEG_n(idx) = [];
        for igroup = 1:length(leg_labels)
            leg_labels{igroup} = [leg_labels{igroup}, ' (n=', num2str(LEG_n(igroup)),')'];
        end
        legend(LEG_H, leg_labels, 'box','off', 'Location','nw')
    end
end

clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_tot SP_tot2 ...
            STATS_response_modulation 
%% Figure 2: Plot specific dataset

% Datasets structure
%__________________________________________________________________________
% datasets:            [1x1] struct
%   |--all:            [1x1] struct
%   |  |--CCIsham:     [1x2] cell
%   |  |--CCI:         [1x3] cell
%   |  |--sham:        [1x3] cell
%   |  '--All:         [1x6] cell
%   |--significant:    [1x1] struct
%   |  |--CCIsham:     [1x2] cell
%   |  |--CCI:         [1x3] cell
%   |  |--sham:        [1x3] cell
%   |  '--All:         [1x6] cell
%   '--baseline:       [1x1] struct
%      |--all:         [1x1] struct
%      |  |--all:      [1010x3] double
%      |  |--HPS:      [2x2] cell
%      |  |--SP:       [2x2] cell
%      |  |--HPS_CCI:  [641x1] double
%      |  |--HPS_sham: [369x1] double
%      |  |--SP_CCI:   [641x1] double
%      |  '--SP_sham:  [369x1] double
%      '--significant: [1x1] struct
%         |--all:      [330x3] double
%         |--HPS:      [2x2] cell
%         |--SP:       [2x2] cell
%         |--HPS_CCI:  [210x1] double
%         |--HPS_sham: [120x1] double
%         |--SP_CCI:   [210x1] double
%         '--SP_sham:  [120x1] double
%__________________________________________________________________________

dataset_to_plot = {datasets.all.sham{1,1},  datasets.all.sham{1,2},  datasets.all.sham{1,3}};

dataset_to_plot = {datasets.baseline.all.SP{1,1}, datasets.baseline.all.SP{1,2}, ...
                    datasets.baseline.all.SP{2,1}, datasets.baseline.all.SP{2,2}}
data_type = 'SP'; % 'HPS'

% coord = struct(); % Keep peak x, y coordinates of ksdensity curves
% coord.dataset = [];


figure;
clf
hold on
for idata = 1:numel(dataset_to_plot)
    
    
    data = {dataset_to_plot{idata}};
    
    %     Settings for plotting normalized AUROC kdensity
    tf = strcmp(data_type, 'HPS');
    if tf == 1
        n_bins = 100;
        bin_edges = linspace(0, 1, n_bins);
        bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
        bin_centers = bin_centers';
        %         bin_edges_smooth = linspace(0, 1, 1000);
        %         bin_centers_smooth = bin_edges_smooth(1:end-1) + (bin_edges_smooth(2)-bin_edges_smooth(1))/2;
    else
        n_bins = 100;
        bin_edges = linspace(-0.0001, 0.03, n_bins);
        bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
        bin_centers = bin_centers';
        %         bin_edges_smooth = linspace(0, 1, 1000);
        %         bin_centers_smooth = bin_edges_smooth(1:end-1) + (bin_edges_smooth(2)-bin_edges_smooth(1))/2;
    end
    
    for i = 1:length(data)
        
        h1 = histcounts(data{1,i}, bin_edges);
        h1 = h1 / sum(h1);
        h1_smooth = ksdensity(data{1,i}, bin_centers, 'Support','positive');
        h1_smooth = h1_smooth / sum(h1_smooth);
        
        %         h1 = h1';
        %         join = [h1 bin_centers h1_smooth];
        %         join = [];
        %         h1 = [];
        %         h1_smooth = [];
        %
        %         bad_idx = h1_smooth < eps;
        %         h1_smooth(bad_idx) = [];
        %         x = bin_centers(~bad_idx);
        
        %         subplot(1,2,i),hold on,
        %         bar(bin_centers, h1),
        plot(bin_centers, h1_smooth, 'linewidth',3);
        
        % Peak h1_smooth coordinates
        %         [a,b] = max(h1_smooth);
        %         x = bin_centers(b,1);
        %         y = a;
        %         coord1 = [x y];
        %         coord.(dataset{idata}) = cat(1,coord.(dataset{idata}),coord1);
        %         plot([x x], [0 y]);
    end
end

% Legend and title
% if idata == 1
%     title('HPS')
%     legend('CCI', 'CCI-LP', 'CCI SAL')
% else
%     title('SP')
%     legend('CCI', 'sham')
% end

keyboard 

title('All HPS baseline'), legend('CCI LP','CCI SAL', 'sham LP', 'sham SAL')

title('All SP baseline'), legend('CCI LP','CCI SAL', 'sham LP', 'sham SAL')
title('PRE'),legend('CCI','sham')
title('CCI'),legend('PRE','LP','SAL')
title('sham'),legend('PRE','LP','SAL')
title('CCI/sham'),legend('CCI PRE','CCI LP','CCI SAL','sham PRE','sham LP','sham SAL')
title('PRE significant'),legend('CCI','sham')


clearvars -except compound_group_names datasets dataset_animal_names DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group numbers results state ...
            STATS_response_modulation
% name_of_script2
%% Figure 3: Boxplots
% Concat data
X = [datasets.DATA_CCI_PRE; datasets.DATA_sham_PRE; ...
    results.(state).(ds_group){1,1}(:,3); results.(state).(ds_group){1,1}(:,4); ...
    results.(state).(ds_group){2,1}(:,3); results.(state).(ds_group){2,1}(:,4); ...
    results.(state).(ds_group){1,2}(:,3); results.(state).(ds_group){1,2}(:,4); ...
    results.(state).(ds_group){2,2}(:,3); results.(state).(ds_group){2,2}(:,4)];
% Grouping variables
G = [ones(size(datasets.DATA_CCI_PRE)); 2*ones(size(datasets.DATA_sham_PRE)); ...
    3*ones(size(results.(state).(ds_group){1,1}(:,3))); 4*ones(size(results.(state).(ds_group){1,1}(:,4))); ...
    5*ones(size(results.(state).(ds_group){2,1}(:,3))); 6*ones(size(results.(state).(ds_group){2,1}(:,4))); ...
    7*ones(size(results.(state).(ds_group){1,2}(:,3))); 8*ones(size(results.(state).(ds_group){1,2}(:,4)));
    9*ones(size(results.(state).(ds_group){2,2}(:,3))); 10*ones(size(results.(state).(ds_group){2,2}(:,4)))];

% Print BoxPlot
figure;
clf
hold on
boxplot(X, G,'Notch','on','Labels',...
    {'CCI PRE','sham PRE',...
    'CCI_LP PRE', 'CCI_LP POST'...
    'sham_LP PRE', 'sham_LP POST', ...
    'CCI_SAL PRE', 'CCI_SAL POST', ...
    'sham_SAL PRE', 'sham_SAL POST'})
%% Figure 4: Violin Plots on grouped data
figure;
clf, hold on
x = 1; offset = 1;
violinDistribution(X(G==1), 'f',.02, 'Side',-1, 'Center',x) % DATA_CCI_PRE
violinDistribution(X(G==2), 'f',.02, 'Side',1, 'Center',x,'FaceColor','r') % DATA_sham_PRE
violinDistribution(X(G==3), 'f',.02, 'Side',-1, 'Center',x+offset) % CCI_LP PRE
violinDistribution(X(G==4), 'f',.02, 'Side',1, 'Center',x+offset,'FaceColor','r') % CCI_LP POST
violinDistribution(X(G==5), 'f',.02, 'Side',-1, 'Center',x+offset*2) % sham_LP PRE
violinDistribution(X(G==6), 'f',.02, 'Side',1, 'Center',x+offset*2,'FaceColor','r') % sham_LP POST
violinDistribution(X(G==7), 'f',.02, 'Side',-1, 'Center',x+offset*3) % CCI_SAL PRE
violinDistribution(X(G==8), 'f',.02, 'Side',1, 'Center',x+offset*3,'FaceColor','r') % CCI_SAL POST
violinDistribution(X(G==9), 'f',.02, 'Side',-1, 'Center',x+offset*4) % sham_SAL PRE
violinDistribution(X(G==10), 'f',.02, 'Side',1, 'Center',x+offset*4,'FaceColor','r') % sham_SAL POST
ylim([0 1])
%% Figure 5: Analyze all cells divided by subpopulation groups
state
group

a = [];
b = [];
c = [];

% Plot all cells of group before after
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)];
        
        figure
        hold on
        for icell = 1:size(results.(state).(ds_group){i,ii},1)
            if results.(state).(ds_group){i,ii}(icell,9) > 0 % results.(time).(ds_group){i,ii}(icell,8) == 0 && 
                plot(a(icell,:)','-or')
            else
                plot(a(icell,:)','-ob')
            end
        end
        xlim([.8, 2.2])
        ylim([0 1])
        
        if i == 1 && ii == 1
            title('CCI - LP')
        elseif i == 1 && ii == 2
            title('sham LP')
        elseif i == 2 && ii == 1
            title('CCI SAL')
        elseif i == 2 && ii == 2
            title('sham SAL')
        end
    end
end
%% Figure 6: ScatterPlot for all cells difference during SP and HPS
figure('color','w')
clf
subplot(221), hold on
    plot(diff_tot(diff_tot(:,3)==1, 2), diff_tot(diff_tot(:,3)==1, 1),'.k','MarkerSize',8), % lsline
    plot(sig_diff_tot(sig_diff_tot(:,3)==1, 2), sig_diff_tot(sig_diff_tot(:,3)==1, 1),'.r','MarkerSize',10),
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['CCI LP (n=', num2str(sum(diff_tot(:,3)==1)),'/',num2str(sum(sig_diff_tot(:,3)==1)), ')']),
%     xlim([-0.16 0.06]), ylim([-0.5 0.06]),
    axis equal square
    xlabel('SP'), ylabel('HPS')
subplot(222), hold on
    plot(diff_tot(diff_tot(:,3)==2, 2), diff_tot(diff_tot(:,3)==2, 1),'.k','MarkerSize',8),
    plot(sig_diff_tot(sig_diff_tot(:,3)==2, 2), sig_diff_tot(sig_diff_tot(:,3)==2, 1),'.r','MarkerSize',10),
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['CCI SAL (n=', num2str(sum(diff_tot(:,3)==2)),'/',num2str(sum(sig_diff_tot(:,3)==2)), ')']),
%     xlim([-0.16 0.06]), ylim([-0.5 0.06]),
    axis equal square
subplot(223), hold on
    plot(diff_tot(diff_tot(:,3)==3, 2), diff_tot(diff_tot(:,3)==3, 1),'.k','MarkerSize',8),
    plot(sig_diff_tot(sig_diff_tot(:,3)==3, 2), sig_diff_tot(sig_diff_tot(:,3)==3, 1),'.r','MarkerSize',10),
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['sham LP (n=', num2str(sum(diff_tot(:,3)==3)),'/',num2str(sum(sig_diff_tot(:,3)==3)), ')']),
%     xlim([-0.16 0.06]), ylim([-0.5 0.06]),
    axis equal square
subplot(224), hold on
    plot(diff_tot(diff_tot(:,3)==4, 2), diff_tot(diff_tot(:,3)==4, 1),'.k','MarkerSize',8),
    plot(sig_diff_tot(sig_diff_tot(:,3)==4, 2), sig_diff_tot(sig_diff_tot(:,3)==4, 1),'.r','MarkerSize',10),
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['sham SAL (n=', num2str(sum(diff_tot(:,3)==4)),'/',num2str(sum(sig_diff_tot(:,3)==4)), ')']),
%     xlim([-0.16 0.06]), ylim([-0.5 0.06]),
    axis equal square
%% Figure 7: Plot all cells for each animal

nel = sum(cell2mat(numbers.animals(:)));
ID = {fieldnames(DISTRIBUTIONS)};

% Find max number of imaging session and tot number of cells
sessions_cells_id = [];
for i = 1:nel
    animID = cell2mat(ID{1,1}(i));
    animal_session_cell_id = [i, size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),2), ...
        size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1)];
    sessions_cells_id = cat(1, sessions_cells_id, animal_session_cell_id);
end
%_________________________________________________________
%|	  1    |               2             |      3        |
% animal_ID|number of SP imaging sessions|number of cells|

max_sessions = max(sessions_cells_id(:,2),[], 1);
max_cells = sum(sessions_cells_id(:,3));

% Concatenate all SP/SP_normalized values for all animals
SP = [];
normSP = [];
SP_tot = [];
normSP_tot = [];

SP_tot2 = [];
SP_tot3 = [];

make_fig = 'off'; % 'off'
for i = 1:nel
    animID = cell2mat(ID{1,1}(i));
    
    % Spontaneous
    SP = NaN(size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1), max_sessions); % Preallocate matrix with NaN to load the SP values
    SP1 = cell2mat(DISTRIBUTIONS.(animID){:, 'mean_activity'}); % Select the SP values
    SP(:,[1:size(SP1,2)]) = SP1; % Substitute the SP values in the preallocated NaN matrix
    
    % Spontaneous normalized to the baseline (difference)
    normSP = NaN(size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1), max_sessions);
    normSP1 = cell2mat(DISTRIBUTIONS.(animID){:, 'mean_activity'});
    normSP(:,[1:size(normSP1,2)]) = normSP1 - normSP1(:,1); % Normalize activity on spontaneous baseline
    
    % Plot SP and normalized SP of all cells for each animal
    tf = strcmp(make_fig, 'on');
    if tf == 1;
        figure
        clf
        subplot(1,2,1),
        plot(SP', '-ok')
        xlim([0, (size(normSP,2)+1)])
        title(animID,'Interpreter','none')
        subplot(1,2,2),
        plot(normSP(:, 2:(size(normSP,2)))', '-ok')
        title({animID ' Normalized to baseline'},'Interpreter','none')
        xlim([0, (size(normSP,2))])
    else
    end
    
    % Add animal ID
    animID = nan(size(SP, 1), 1);
    animID(:) = i;
    SP = [animID SP];
    normSP = [animID normSP];
    
    % Concatenate all the cells for all the animals
    SP_tot = cat(1,SP_tot, SP);
    normSP_tot = cat(1,normSP_tot, normSP);
end

% Order and find animal_ID values based on 'experimental_group'/'compound'
a = dataset_animal_names{ismember(dataset_animal_names.experimental_group, 'CCI') & ...
                        ismember(dataset_animal_names.compound, 'LP-211'), 'animal_ID'}; % Find animal_ID which contains CCI AND LP-211
dataset_animal_names = natsortrows(dataset_animal_names, 1); % Sort rows based on ascending column 1

% Find same group / compound
%------------------%
%     |LP-211|SAL| %
%  CCI|  1   | 2 | %
% sham|  3   | 4 | %
%------------------%

group1 = [];
group = [];

for i = 1:nel % Assign group_ID
    if ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'LP-211')
        group1 = 1;
    elseif ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'SAL')
        group1 = 2;
    elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'LP-211')
        group1 = 3;
    elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'SAL')
        group1 = 4;
    end
    group = [group; group1];
end

TabGroup = table(ds_group);
dataset_animal_names = [dataset_animal_names TabGroup]; % Concatenate the group to the dataset

% Find group ID for each animal
animGroup = nan(size(SP_tot,1), 1);
for i = 1:length(SP_tot)
    animGroup(i,1) = dataset_animal_names.group(SP_tot(i,1));
end
animInfo = [SP_tot(:,1) animGroup]; % Cat animal_ID and group

SP_tot2 = [animInfo SP_tot(:,2:end)]; % Cat animal_ID and group in the SP file
normSP_tot = [animInfo normSP_tot(:,2:end)]; % Cat animal_ID and group the normalized SP file

% Make statistics based on grouping
SP_tot3 = array2table(SP_tot2,'VariableNames',{'animal_ID', 'animal_group', 'Baseline', 'Active1', 'Active2', 'Active3', 'Active4'});
grpstats(SP_tot3,'animal_group', 'median');

% Find max number of imaging session per animal
SP_diff_tot = [];
evoked_diff_tot = [];
evoked_pre = [];
evoked_post = [];

for i = 1:nel
    animID = cell2mat(ID{1,1}(i));
    ROI_id = DISTRIBUTIONS.(animID){:,'ROI_id'};
    length(ROI_id);
    
    % Assign unique group value to surgery / treatment combination
    dataset_animal_names = natsortrows(dataset_animal_names, 1); % Sort rows based on ascending column 1
    if ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'LP-211')
        group1 = 1;
    elseif ismember(dataset_animal_names.experimental_group(i), 'CCI') && ismember(dataset_animal_names.compound(i), 'SAL')
        group1 = 2;
    elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'LP-211')
        group1 = 3;
    elseif ismember(dataset_animal_names.experimental_group(i), 'sham') && ismember(dataset_animal_names.compound(i), 'SAL')
        group1 = 4;
    end
    
    for ii = 1:length(ROI_id)
        cel_id = DISTRIBUTIONS.(animID){ii,'ROI_id'};
        
        % Calculate SP difference: POST - PRE and assign ID if positive or
        % negative
        SP = (cell2mat(DISTRIBUTIONS.(animID){cel_id,'mean_activity'}));
        SP_diff = SP(2) - SP(1);
        if SP_diff >= 0 % Increased activity after treatment
            SP_diff(1,2) = 1;
        else
            SP_diff(1,2) = -1; % Decreased activity after treatment
        end
        SP_diff = [SP_diff group1]; % Add animal group based on surgery / treatment
        
        % Calculate EVOKED difference: POST - PRE
        evoked_pre = STATS_response_modulation.(animID)((ismember(STATS_response_modulation.(animID){:,'ROI_id'}, cel_id) & ...
            ismember(STATS_response_modulation.(animID){:,'compound_phase'}, 'baseline')),:);
        evoked_pre = mean(cell2mat(evoked_pre{:,'activity'}));
        if isnan(evoked_pre)
            evoked_pre = 0;
        else
            evoked_pre = evoked_pre(1,2);
        end
        
        evoked_post = STATS_response_modulation.(animID)((ismember(STATS_response_modulation.(animID){:,'ROI_id'}, cel_id) & ...
            ismember(STATS_response_modulation.(animID){:,'compound_phase'}, 'active_1')),:);
        evoked_post = mean(cell2mat(evoked_post{:,'activity'}));
        
        if isnan(evoked_post)
            evoked_post = 0;
        else
            evoked_post = evoked_post(1,2);
        end
        evoked_diff = evoked_post - evoked_pre;
        evoked_diff = [evoked_diff group1]; % Add animal group based on surgery / treatment
        
        
        SP_diff_tot = cat(1, SP_diff_tot, SP_diff);
        evoked_diff_tot = cat(1, evoked_diff_tot, evoked_diff);
    end
end
diff_tot = cat(2, SP_diff_tot(:,1), evoked_diff_tot, SP_diff_tot(:,2)); % Create single variable containing SP, EV, GROUP info for each cell

clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_tot SP_tot2 ...
            SP_tot3 STATS_response_modulation TabGroup
%% Figure 8: Plot SP and normalized SP of all cells for each animal
for i = 1: size(sessions_cells_id, 1)
	figure
	clf
	subplot(1,2,1),
	plot(((SP_tot((SP_tot(:, 1) == i), 2:sessions_cells_id(i, 2)+1)))', '-ok')
	xlim([0, (size(normSP,2)+1)])
	title(animID,'Interpreter','none')
	subplot(1,2,2),
	plot(((normSP_tot((normSP_tot(:, 1) == i), 3:sessions_cells_id(i, 2)+2)))', '-ok')
	title({animID ' Normalized to baseline'},'Interpreter','none')
	xlim([0, (size(normSP,2))])
end

clearvars i
%% Figure 9: Plot mean/median AUROC values for each animal group based on significant / >=0 values

a = [];
b = [];
b1 = [];
c = [];
c1 = [];
d = [];

make_fig = 'on'; % 'off' 'on'

% Cells AUROC >= 0.5, significant PRE or POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 | results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2); % Find significant value PRE and POST
        %         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP PRE/POST')
        elseif i == 1 && ii == 2
            title('sham LP PRE/POST')
        elseif i == 2 && ii == 1
            title('CCI SAL PRE/POST')
        elseif i == 2 && ii == 2
            title('sham SAL PRE/POST')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells AUROC >= 0.5, significant only during POST not PRE
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during POST
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5 in PRE or POST
        c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        %         c1(:,3) = c1(:,1)+c1(:,2); % Find any AUROC value >= 0.5 PRE or POST
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP POST not PRE')
        elseif i == 1 && ii == 2
            title('sham LP POST not PRE')
        elseif i == 2 && ii == 1
            title('CCI SAL POST not PRE')
        elseif i == 2 && ii == 2
            title('sham SAL POST not PRE')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells AUROC >= 0.5, significant only during PRE not POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find significant values during PRE
        % b1 = results.(time).(ds_group){i,ii}(:,7)>0; % Find significant values during POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        %         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        %         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
        d = c(any(c1(:,3),3),:); % Keep only 1
        d = c;
        d(:,3) = [];
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP PRE not POST')
        elseif i == 1 && ii == 2
            title('sham LP PRE not POST')
        elseif i == 2 && ii == 1
            title('CCI SAL PRE not POST')
        elseif i == 2 && ii == 2
            title('sham SAL PRE not POST')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells AUROC >= 0.5, always significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        % b1 = results.(time).(ds_group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP PRE/POST')
        elseif i == 1 && ii == 2
            title('sham LP PRE/POST')
        elseif i == 2 && ii == 1
            title('CCI SAL PRE/POST')
        elseif i == 2 && ii == 2
            title('sham SAL PRE/POST')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells AUROC >= 0.5, never significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find never significant values
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3);% Find only values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2);
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP never')
        elseif i == 1 && ii == 2
            title('sham LP never')
        elseif i == 2 && ii == 1
            title('CCI SAL never')
        elseif i == 2 && ii == 2
            title('sham SAL never')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end

%--------------------------------------------------------------------------

% Cells significant PRE or POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 | results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
        %         c1(:,3) = c1(:,1).*c1(:,2); % Find significant value PRE and POST
        %         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
        %         d = c(any(c1(:,3),3),:); % Keep only 1
        %         d(:,3) = [];
        d = c(:,1:2);

        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP PRE/POST')
        elseif i == 1 && ii == 2
            title('sham LP PRE/POST')
        elseif i == 2 && ii == 1
            title('CCI SAL PRE/POST')
        elseif i == 2 && ii == 2
            title('sham SAL PRE/POST')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells significant only during POST not PRE
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during POST
        c = b(any(b(:,3),3),:); % Keep only 1
        %         c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5 in PRE or POST
        %         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        %         c1(:,3) = c1(:,1)+c1(:,2); % Find any AUROC value >= 0.5 PRE or POST
        %         d = c(any(c1(:,3),3),:); % Keep only 1
        %         d(:,3) = [];
        d = c(:,1:2);
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP POST not PRE')
        elseif i == 1 && ii == 2
            title('sham LP POST not PRE')
        elseif i == 2 && ii == 1
            title('CCI SAL POST not PRE')
        elseif i == 2 && ii == 2
            title('sham SAL POST not PRE')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells significant only during PRE not POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find significant values during PRE
        % b1 = results.(time).(ds_group){i,ii}(:,7)>0; % Find significant values during POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        %         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
        c = b(any(b(:,3),3),:); % Keep only 1
        %         c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5
        %         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        %         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
        %         d = c(any(c1(:,3),3),:); % Keep only 1
        %         d = c;
        %         d(:,3) = [];
        d = c(:,1:2);
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP PRE not POST')
        elseif i == 1 && ii == 2
            title('sham LP PRE not POST')
        elseif i == 2 && ii == 1
            title('CCI SAL PRE not POST')
        elseif i == 2 && ii == 2
            title('sham SAL PRE not POST')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells always significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)>0 & results.(state).(ds_group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        % b1 = results.(time).(ds_group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        c = b(any(b(:,3),3),:); % Keep only 1
        %         c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
        %         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        %         d = c(any(c1(:,3),3),:); % Keep only 1
        %         d(:,3) = [];
        d = c(:,1:2);
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP PRE/POST')
        elseif i == 1 && ii == 2
            title('sham LP PRE/POST')
        elseif i == 2 && ii == 1
            title('CCI SAL PRE/POST')
        elseif i == 2 && ii == 2
            title('sham SAL PRE/POST')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end
% Cells never significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(ds_group){i,ii}(:,3) results.(state).(ds_group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(ds_group){i,ii}(:,8)==0 & results.(state).(ds_group){i,ii}(:,9)==0)); % Find never significant values
        c = b(any(b(:,3),3),:); % Keep only 1
        %         c1 = any(c(:,1:2)>=0.5,3);% Find only values >=0.5
        %         c1(:,3) = c1(:,1).*c1(:,2);
        %         d = c(any(c1(:,3),3),:); % Keep only 1
        %         d(:,3) = [];
        d = c(:,1:2);
        
        figure,
        clf,
        hold on,
        plot(d', '-or'),
        xlim([.8, 2.2]),
        ylim([0 1]),
        plot(mean(d),'-ok','MarkerFaceColor','k','MarkerSize',10,'LineWidth',2);
        plot(median(d),'-+g','MarkerFaceColor','k','LineWidth',2);
        txt = {'p',p,'signr',signr,'meddif',meddif,'medif',medif};
        text(1,0.2,txt)
        txt1 = {'med',med,'me',me};
        text(1.8,0.2,txt1)
        
        if i == 1 && ii == 1
            title('CCI - LP never')
        elseif i == 1 && ii == 2
            title('sham LP never')
        elseif i == 2 && ii == 1
            title('CCI SAL never')
        elseif i == 2 && ii == 2
            title('sham SAL never')
        end
        
        a = [];
        b = [];
        c = [];
        d = [];
    end
end


clearvars -except compound_group_names datasets dataset_animal_names diff_tot DISTRIBUTIONS ...
            experiment_group_names GC datasets ds_group normSP_tot numbers results state ...
            SP_ANOVA SP_ANOVA2 SP_ANOVA_NaN SP_tot SP_tot2 ...
            STATS_response_modulation TabGroup
%% Figure 10: Plot SP activity histogram over frames for each cell grouped by animal

% a = dataset_animal_names{ismember(dataset_animal_names.experimental_group, 'CCI') & ...
%                         ismember(dataset_animal_names.compound, 'LP-211'), 'animal_ID'}; % Find animal_ID which contains CCI AND LP-211


% plot(GC.spontaneous_interval_bins, conv([DISTRIBUTIONS.AB_1.interval_dist_before{1}(1), diff(DISTRIBUTIONS.AB_1.interval_dist_before{1})], ones(20,1)/20, 'same')), xlim([0,100])
% semilogx(GC.spontaneous_interval_bins, DISTRIBUTIONS.AB_1.interval_dist_before{2}), xlim([0,1000])

aname = dataset_animal_names.animal_ID;
    
for ii = 1:height(dataset_animal_names)
    a = [];
    b = [];
    a = cat(1, DISTRIBUTIONS.(aname{ii}).interval_dist_before{:});

    for i = 1:size(a, 1)
        b(i, :) = conv([a(i, 1), diff(a(i, :))], ones(20, 1)/20, 'same');
    end

    figure
    subplot(121)
    plot(GC.spontaneous_interval_bins, b);
    xlim([0 100]);
    title(['AB\_', num2str(ii)])
    subplot(122)
    semilogx(GC.spontaneous_interval_bins, a)
    xlim([0,1000]);
end

%-------------------------------------------------------------------------%
% Divide cell SP frequencies by animal group
freq = cell(2, 2);
aname = dataset_animal_names.animal_ID;
timepoint = 'interval_dist_before';

for i = 1:height(dataset_animal_names)
    if dataset_animal_names.group(i) == 1
        freq{1, 1} = cat(1, freq{1, 1}, DISTRIBUTIONS.(aname{i}).(timepoint){:});
    elseif dataset_animal_names.group(i) == 2
        freq{1, 2} = cat(1, freq{1, 2}, DISTRIBUTIONS.(aname{i}).(timepoint){:});
    elseif dataset_animal_names.group(i) == 3
        freq{2, 1} = cat(1, freq{2, 1}, DISTRIBUTIONS.(aname{i}).(timepoint){:});
    elseif dataset_animal_names.group(i) == 4
        freq{2, 2} = cat(1, freq{2, 2}, DISTRIBUTIONS.(aname{i}).(timepoint){:});
    end
end

results.SP.freq = freq;

% Plot cell SP frequencies by animal group
figure
hold on
for irow = 1:2
    for icol = 1:2
        a = mean(freq{irow, icol}, 1);
        b = [];
        b = conv([a(1), diff(a)], ones(20, 1)/20, 'same');
        
        plot(GC.spontaneous_interval_bins, b, 'LineWidth', 2)
    end
end
xlim([0 200]);
ylim([0 0.0055]);
legend('CCI + LP-211', 'CCI + SAL', 'sham + LP-211', 'sham + SAL');
%% Figure 11: Plot SP activity histogram of significantly modulated cells over frames for each cell grouped by group
% Spontaneous Interval of events
x = linspace(0, 600, 5100); % Convert frame into seconds
timepoint = 'interval_dist_before';

figure
hold on
for i = 1:2
    for ii = 1:2
        a = mean(cell2mat([SP_int{i, ii}.interval_dist_after_1; SP_int{i, ii}.interval_dist_after_2]), 1);
        b = [];
        b = conv([a(1), diff(a)], ones(20, 1)/20, 'same');
%         plot(GC.spontaneous_interval_bins, b, 'LineWidth', 2)
        plot(x, b, 'LineWidth', 2)
    end
end
set(gca, 'yscale', 'log')
xlim([0 50]);
ylim([10e-5 10e-2])
legend('CCI + LP-211', 'CCI + SAL', 'sham + LP-211', 'sham + SAL');
legend('CCI + before', 'CCI + LP-211');
legend('sham + before', 'sham + LP-211');
legend('CCI + before', 'sham + before');

title({'SP events intervall all'; timepoint});

% Collect different datasets to plot
c = [SP_int_sig{1, 1}.interval_dist_before; SP_int_sig{1, 2}.interval_dist_before]; a = mean(cell2mat(c), 1); % All CCI before
c = [SP_int_sig{1, 1}.interval_dist_after_1; SP_int_sig{1, 1}.interval_dist_after_2]; a = mean(cell2mat(c), 1); % All CCI+LP after
c = [SP_int_sig{2, 1}.interval_dist_before; SP_int_sig{2, 2}.interval_dist_before]; a = mean(cell2mat(c), 1); % All sham before
c = [SP_int_sig{2, 1}.interval_dist_after_1; SP_int_sig{2, 1}.interval_dist_after_2]; a = mean(cell2mat(c), 1); % All sham+LP after


%------------------
% Spontaneous Amplitude of events
x = linspace(0, 600, 10000); % Convert frame into seconds
timepoint = 'amplitude_dist_before'

figure
hold on
for i = 1:2
    for ii = 1:2
        a = mean(cell2mat(SP_amp{i, ii}.(timepoint)), 1);
        b = [];
        b = conv([a(1), diff(a)], ones(20, 1)/20, 'same');
%         plot(GC.spontaneous_amplitude_bins, b, 'LineWidth', 2)
        plot(x, b, 'LineWidth', 2)
    end
end
% set(gca, 'yscale', 'log')
xlim([0 15]);
legend('CCI + LP-211', 'CCI + SAL', 'sham + LP-211', 'sham + SAL');
title({'SP events amplitude all'; timepoint});

% Collect different datasets to plot
c = [SP_amp{1, 1}.amplitude_dist_before; SP_amp{1, 2}.amplitude_dist_before]; % All CCI before
c = [SP_amp{1, 1}.amplitude_dist_after_1; SP_amp{1, 1}.amplitude_dist_after_2]; % All CCI+LP after
c = [SP_amp{2, 1}.amplitude_dist_before; SP_amp{2, 2}.amplitude_dist_before]; % All sham before
c = [SP_amp{2, 1}.amplitude_dist_after_1; SP_amp{2, 1}.amplitude_dist_after_2]; % All sham+LP after

a = mean(cell2mat(c), 1);
