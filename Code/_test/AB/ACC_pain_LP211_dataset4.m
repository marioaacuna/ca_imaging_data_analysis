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

% Load variables for spontaneous activity
stats_filename = get_filename_of('spontaneous_activity_stats', 'ACC_pain_LP-211');

% STATS_spontaneous_activity = load_variable(stats_filename, 'STATS_spontaneous_activity');
DISTRIBUTIONS = load_variable(stats_filename, 'DISTRIBUTIONS');

% Add parameters
phases_order = {'baseline', 'active_1'};
LEG_labels = {'non-selective, non-modulated', 'selective, non-modulated', 'only selective before', 'only selective after'};

toggle_toolbox('_plotting', 'on');
%% Find cells AUROC/p and assign to results structure
compound_group_names = {'LP211','SAL'};
experiment_group_names = {'CCI','sham'};
results = struct();
results.HPS.DATA = cell(length(compound_group_names), length(experiment_group_names));
results.HPS.DATA_p = results.HPS.DATA;
results.SP.DATA = cell(length(compound_group_names), length(experiment_group_names));

for iid = 1:n_datasets
    
    animal_ID = cell2mat(dataset_animal_names{iid,1});
    stats_HPS = STATS_response_modulation.(animal_ID);
    stats_SP = DISTRIBUTIONS.(animal_ID);
    
    [cells_ID, ~, rows_idx] = unique(stats_HPS.ROI_id);
    
    n_cells = length(cells_ID);
    
    HPS_DATA = NaN(n_cells, length(phases_order)+1);
    HPS_DATA_p = NaN(n_cells, length(phases_order)+1);
    SP_DATA = NaN(n_cells, length(phases_order)+1);
    
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


clearvars iid icell phi row col HPS_DATA HPS_DATA_p SP_DATA n_cells ...
        p cells_ID AUROC animal_ID data_SP data_HPS rows rows_idx;
%____________________________________________________________________________
%|	  1    |   2   | 3 | 4  |    5     |       6       |	     7          |
%|Animal_ID|Cell_ID|PRE|POST|Cell_group|Diff_(POST-PRE)|Norm_diff_(POST-PRE)|
%% Isolate significant modulated cells | Add cell ID to 'results.HPS.ID_' to keep track of cell behavior during HPS
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

% Keep significant cells
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

clearvars a b anim_ID i n_datasets
%% Assign ID during SP to isolate 'increasing' / 'decreasing' activity cells

for a = 1:2
    for b = 1:2
        i = results.SP.DATA{a,b}(:,4) - results.SP.DATA{a,b}(:,3);
        ii = zeros(length(results.SP.DATA{a,b}),1);
        for id = 1:length(results.SP.DATA{a,b})
            if i(id) >= 0 % Increased activity after treatment
                ii(id) = 1;
            else
                ii(id) = -1; % Decreased activity after treatment
            end
        end
        results.SP.DATA{a,b} = cat(2, results.SP.DATA{a,b}, ii);
    end
end

clearvars a b i ii
%% Define data groups to plot

    % DATA/DATA_sig/DATA_sig_ groups
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  x   | x | %
    % sham|  x   | x | %
    %------------------%

% Specify the group from where data is used (results.HPS.DATA / results.HPS.DATA_sig_st)
state = 'HPS'; % 'SP'
group = 'DATA_sig'; % HPS.DATA; DATA_sig; DATA_sig_ch; DATA_sig_st, DATA_sig_post

DATA_CCI_PRE = [results.(state).(group){1,1}(:,3); results.(state).(group){1,2}(:,3)]; % AUROC CCI PRE treatment
DATA_sham_PRE = [results.(state).(group){2,1}(:,3); results.(state).(group){2,2}(:,3)]; % AUROC sham PRE treatment

%% Info for results dataset

    % Avg dataset
    % results. ... .DATA Table Structure:
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  x   | x | %
    % sham|  x   | x | %
    %------------------%
    
ncells = cellfun(@(x) size(x,1), results.HPS.DATA);
avg = cellfun(@(x) mean(x,1), results.HPS.DATA, 'UniformOutput',false);
avg_no_zeros = cellfun(@(x) (sum(x,1) ./ sum(x~=0,1)), results.HPS.DATA_sig, 'UniformOutput',false); % Average of non zeros values
cellfun(@(x) sum(x(:,6:7)~=0,1), results.HPS.DATA,'UniformOutput',false); % Number of positive values in column 6 and 7 (i.e number of cells responding PRE or POST)
%% Plot single cells based on treatment group
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
        
        for ix = 1:ncells(irow,icol)
            x = results.HPS.DATA{irow, icol}(ix, 3);
            y = results.HPS.DATA{irow, icol}(ix, 4);
            p_pre = results.HPS.DATA_p{irow, icol}(ix, 3);
            p_post = results.HPS.DATA_p{irow, icol}(ix, 4);
            
            
            % Non-selective, non-modulated cells
            if p_pre > 0.05 && p_post > 0.05
                group = 1;
                mfc = [.7, .7, .7];
                mec = 'w';
                marker = 'o';
                markersize = 8;
                lw = 1;
                
                % Selective, non-modulated cells
            elseif p_pre <= 0.05 && p_post <= 0.05
                group = 2;
                mfc = 'k';
                mec = 'none';
                marker = 'p';
                markersize = 18;
                lw = 1;
                
                % Only selective before
            elseif p_pre <= 0.05 && p_post > 0.05
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
            elseif p_pre > 0.05 && p_post <= 0.05
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
            if isnan(LEG_H(group))
                LEG_H(group) = H(ix);
            end
            LEG_n(group) = LEG_n(group) + 1;
        end
        
        axis square
        set(ax, 'FontSize',12, 'TickDir','out')
        xlabel('Selectivity to HPS before LP-211', 'FontSize',13), ylabel('Selectivity to HPS after LP-211', 'FontSize',13)
        title([compound_group_names{irow},' + ',experiment_group_names{icol}, ' (n=', num2str(ncells(irow, icol)), ')'], 'FontSize',16, 'Interpreter','none')
        
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
%% Plot specific datasets
% Define the dataset to plot
all_datasets = struct();
all_datasets.dataset1 = {DATA_CCI_PRE, DATA_sham_PRE};
all_datasets.dataset2 = {DATA_CCI_PRE, results.(state).(group){1,1}(:,4), results.(state).(group){1,2}(:,4)};
all_datasets.dataset3 = {DATA_sham_PRE, results.(state).(group){2,1}(:,4), results.(state).(group){2,2}(:,4)};
all_datasets.dataset4 = {DATA_CCI_PRE, results.(state).(group){1,1}(:,4), results.(state).(group){1,2}(:,4), ...
                        DATA_sham_PRE, results.(state).(group){2,1}(:,4), results.(state).(group){2,2}(:,4)};

HPS_basel = a1(:,3);
SP_basel = b1(:,3);
groups = [ones(size(a{1,1},1),1); 2*ones(size(a{2,1},1),1); 3*ones(size(a{3,1},1),1); 4*ones(size(a{4,1},1),1)];
baselines = [HPS_basel SP_basel groups];

% Dataset for CCI, sham, HPS, SP
HPS_basel_CCI_LP = baselines(baselines(:,3)==1, 1);
HPS_basel_CCI_SAL = baselines(baselines(:,3)==2, 1);
HPS_basel_sham_LP = baselines(baselines(:,3)==3, 1);
HPS_basel_sham_SAL = baselines(baselines(:,3)==4, 1);

SP_basel_CCI_LP = baselines(baselines(:,3)==1, 2);
SP_basel_CCI_SAL = baselines(baselines(:,3)==2, 2);
SP_basel_sham_LP = baselines(baselines(:,3)==3, 2);
SP_basel_sham_SAL = baselines(baselines(:,3)==4, 2);

all_datasets.dataset5 = {HPS_basel_CCI_LP, HPS_basel_CCI_SAL, HPS_basel_sham_LP, HPS_basel_sham_SAL};
all_datasets.dataset6 = {SP_basel_CCI_LP, SP_basel_CCI_SAL, SP_basel_sham_LP, SP_basel_sham_SAL};

% Dataset for CCI and sham, HPS, SP
HPS_basel_CCI = [baselines(baselines(:,3)==1, 1); baselines(baselines(:,3)==2, 1)];
HPS_basel_sham = [baselines(baselines(:,3)==3, 1); baselines(baselines(:,3)==4, 1)];

SP_basel_CCI = [baselines(baselines(:,3)==1, 2); baselines(baselines(:,3)==2, 2)];
SP_basel_sham = [baselines(baselines(:,3)==3, 2); baselines(baselines(:,3)==4, 2)];

all_datasets.dataset7a = {HPS_basel_CCI, HPS_basel_sham,}
all_datasets.dataset7b = {SP_basel_CCI, SP_basel_sham}
    
dataset = {'dataset1'; 'dataset2'; 'dataset3'; 'dataset4'};
dataset = {'dataset5', 'dataset6'};
dataset = {'dataset7a', 'dataset7b'}
coord = struct();
coord.dataset1 = [];
coord.dataset2 = [];
coord.dataset3 = [];
coord.dataset4 = [];
coord.dataset5 = [];
coord.dataset6 = [];
coord.dataset7a = [];
coord.dataset7b = [];

% Settings for plotting normalized AUROC kdensity
% n_bins = 100;
% bin_edges = linspace(-0.001, 0.03, n_bins);
% bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
% bin_centers = bin_centers';
% bin_edges_smooth = linspace(0, 1, 1000);
% bin_centers_smooth = bin_edges_smooth(1:end-1) + (bin_edges_smooth(2)-bin_edges_smooth(1))/2;

for idata = 1:numel(dataset)
    
    figure;
    clf
    hold on
    data = all_datasets.(dataset{idata});
    
    if idata == 1
        n_bins = 100;
        bin_edges = linspace(0, 1, n_bins);
        bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
        bin_centers = bin_centers';
    else
        n_bins = 100;
        bin_edges = linspace(-0.0001, 0.02, n_bins);
        bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
        bin_centers = bin_centers';
    end
    
    for i = 1:length(data)
        
        h1 = histcounts(data{1,i}, bin_edges);
        h1 = h1 / sum(h1);
        h1_smooth = ksdensity(data{1,i}, bin_centers, 'Support','positive');
        h1_smooth = h1_smooth / sum(h1_smooth);
        
        %     h1 = h1';
        %     join = [h1 bin_centers h1_smooth];
        
        %     join = [];
        %     h1 = [];
        %     h1_smooth = [];
        
        % bad_idx = h1_smooth < eps;
        % h1_smooth(bad_idx) = [];
        % x = bin_centers(~bad_idx);
        
        %subplot(1,2,i),hold on,
        %bar(bin_centers, h1),
        plot(bin_centers, h1_smooth, 'linewidth',3);
        
        % Peak h1_smooth coordinates
        [a,b] = max(h1_smooth);
        x = bin_centers(b,1);
        y = a;
        coord1 = [x y];
        coord.(dataset{idata}) = cat(1,coord.(dataset{idata}),coord1);
        % plot([x x], [0 y]);
    end
    if idata == 1
        title('HPS')
        legend('CCI', 'sham')
    else
        title('SP')
        legend('CCI', 'sham')
    end
%     if idata == 1
%         title('PRE'),legend('CCI','sham')
%     elseif idata == 2
%         title('CCI'),legend('PRE','LP','SAL')
%     elseif idata == 3
%         title('sham'),legend('PRE','LP','SAL')
%     elseif idata == 4
%         title('CCI/sham'),legend('CCI PRE','CCI LP','CCI SAL','sham PRE','sham LP','sham SAL')
%     end
end

% title('PRE'),legend('CCI','sham')
% title('CCI'),legend('PRE','LP','SAL')
% title('sham'),legend('PRE','LP','SAL')
% title('CCI/sham'),legend('CCI PRE','CCI LP','CCI SAL','sham PRE','sham LP','sham SAL')
% title('PRE significant'),legend('CCI','sham')
%% Summary values of grouped cells
% Cells numbers
results_values = struct();

results_values.ncells = struct();
results_values.ncells.tot = cellfun(@(x) size(x,1), results.(state).(group));
results_values.ncells.lost = cellfun(@(x) sum(x(:,6)~=0,1), results.(state).(group),'UniformOutput',false); % Number of positive values in column 6 (i.e number of cells responding during PRE, loosing selectivity during POST)
results_values.ncells.gain = cellfun(@(x) sum(x(:,7)~=0,1), results.(state).(group),'UniformOutput',false); % Number of positive values in column 7 (i.e number of cells responding during POST, gaining selectivity after PRE)

% Avg AUROC values
avg_AUROC = cellfun(@(x) (sum(x(:,6:7),1) ./ sum(x(:,6:7)~=0,1)), results.(state).(group), 'UniformOutput',false); % Average of non zeros values during PRE and POST (col 6 and 7)

results_values.ncells = ncells;
results_values.avg_AUROC = avg_AUROC;
%% Statistics on the distribution
% Wilcoxon rank sum test for unpaired distribution (CCI_PRE/sham_PRE)
% Wilcoxon signed rank test for paired distribution (CCI_PRE/LP | sham_PRE/LP)
p_sign = [];
p_rank = []

% CCI_PRE vs sham_PRE
x = DATA_CCI_PRE; % CCI_LP PRE + CCI_SAL PRE
y = DATA_sham_PRE; % sham_LP PRE + sham_SAL PRE
p1 = ranksum(x,y);
p_sign = cat(1,p_sign, p1);
p_rank = cat(1,p_sign,p_rank)

% CCI_LP PRE vs CCI_LP POST
x = results.(state).(group){1,1}(:,3); % CCI_LP PRE
y = results.(state).(group){1,1}(:,4); % CCI_LP POST
p1 = signrank(x,y);
p_sign = cat(1,p_sign,p1)

% CCI_SAL PRE vs CCI_SAL POST
x = results.(state).(group){1,2}(:,3); % CCI_SAL PRE
y = results.(state).(group){1,2}(:,4); % CCI_SAL POST
p1 = signrank(x,y);
p_sign = cat(1,p_sign,p1)

% CCI_PRE vs CCI_LP POST
x = DATA_CCI_PRE;
y = results.(state).(group){1,1}(:,4); % CCI_LP POST
p1 = ranksum(x,y);
p_rank = cat(1,p_rank,p1)

% sham_LP PRE vs sham_LP POST
x = results.(state).(group){2,1}(:,3); % sham_LP PRE
y = results.(state).(group){2,1}(:,4); % sham_LP POST
p1 = signrank(x,y);
p_sign = cat(1,p_sign,p1)

% sham_PRE vs sham_LP POST
x = DATA_sham_PRE;
y = results.(state).(group){2,1}(:,4); % sham_LP POST
p1 = ranksum(x,y);
p_rank = cat(1,p_rank,p1)
%% BOXPLOT
% Concat data
X = [DATA_CCI_PRE; DATA_sham_PRE; ...
    results.(state).(group){1,1}(:,3); results.(state).(group){1,1}(:,4); ...
    results.(state).(group){2,1}(:,3); results.(state).(group){2,1}(:,4); ...
    results.(state).(group){1,2}(:,3); results.(state).(group){1,2}(:,4); ...
    results.(state).(group){2,2}(:,3); results.(state).(group){2,2}(:,4)];
% Grouping variables
G = [ones(size(DATA_CCI_PRE)); 2*ones(size(DATA_sham_PRE)); ...
    3*ones(size(results.(state).(group){1,1}(:,3))); 4*ones(size(results.(state).(group){1,1}(:,4))); ...
    5*ones(size(results.(state).(group){2,1}(:,3))); 6*ones(size(results.(state).(group){2,1}(:,4))); ...
    7*ones(size(results.(state).(group){1,2}(:,3))); 8*ones(size(results.(state).(group){1,2}(:,4)));
    9*ones(size(results.(state).(group){2,2}(:,3))); 10*ones(size(results.(state).(group){2,2}(:,4)))];

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
%% Violin Plots on grouped data
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
%% Analyze subpopulation groups
state
group

a = [];
b = [];
c = [];

% Plot all cells of group before after
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)];
        
        figure
        hold on
        for icell = 1:size(results.(state).(group){i,ii},1)
            if results.(state).(group){i,ii}(icell,9) > 0 % results.(time).(group){i,ii}(icell,8) == 0 && 
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
%% Median statistics on defined group a (STABILITY)
% Prepare data to analyze
general_stats = struct();
general_stats.median = cell(2,2);
general_stats.mean = cell(2,2);

a = [];
b = [];
b1 = [];
c = [];
c1 = [];
d = [];

% Cells AUROC >= 0.5, significant PRE or POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)>0 | results.(state).(group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2); % Find significant value PRE and POST
%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
        
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, significant only during POST not PRE
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)==0 & results.(state).(group){i,ii}(:,9)>0)); % Find significant values during POST
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5 in PRE or POST
        c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
%         c1(:,3) = c1(:,1)+c1(:,2); % Find any AUROC value >= 0.5 PRE or POST
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, significant only during PRE not POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)>0 & results.(state).(group){i,ii}(:,9)==0)); % Find significant values during PRE
        % b1 = results.(time).(group){i,ii}(:,7)>0; % Find significant values during POST
        % b1(:,3) = b1(:,1).*b1(:,2);
%         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST         
        d = c(any(c1(:,3),3),:); % Keep only 1
        d = c;
        d(:,3) = [];
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, always significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)>0 & results.(state).(group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        % b1 = results.(time).(group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, never significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)==0 & results.(state).(group){i,ii}(:,9)==0)); % Find never significant values
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3);% Find only values >=0.5
        c1(:,3) = c1(:,1).*c1(:,2);
        d = c(any(c1(:,3),3),:); % Keep only 1
        d(:,3) = [];
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

%--------------------------------------------------------------------------
% Cells significant PRE or POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)>0 | results.(state).(group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
%         c1(:,3) = c1(:,1).*c1(:,2); % Find significant value PRE and POST
%         c1(:,3) = c1(:,1)+c1(:,2); % Find any significant value PRE or POST
%         d = c(any(c1(:,3),3),:); % Keep only 1
%         d(:,3) = [];
        d = c(:,1:2);
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, significant only during POST not PRE
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)==0 & results.(state).(group){i,ii}(:,9)>0)); % Find significant values during POST
        c = b(any(b(:,3),3),:); % Keep only 1
%         c1 = any(c(:,1:2)>=0.5,3);% Find only AUROC values >=0.5 in PRE or POST
%         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
%         c1(:,3) = c1(:,1)+c1(:,2); % Find any AUROC value >= 0.5 PRE or POST
%         d = c(any(c1(:,3),3),:); % Keep only 1
%         d(:,3) = [];
        d = c(:,1:2);
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, significant only during PRE not POST
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)>0 & results.(state).(group){i,ii}(:,9)==0)); % Find significant values during PRE
        % b1 = results.(time).(group){i,ii}(:,7)>0; % Find significant values during POST
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
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, always significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)>0 & results.(state).(group){i,ii}(:,9)>0)); % Find significant values during PRE/POST
        % b1 = results.(time).(group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        c = b(any(b(:,3),3),:); % Keep only 1
%         c1 = any(c(:,1:2)>=0.5,3); % Find only values >=0.5
%         c1(:,3) = c1(:,1).*c1(:,2); % Find AUROC values >= 0.5 in PRE and POST
%         d = c(any(c1(:,3),3),:); % Keep only 1
%         d(:,3) = [];
        d = c(:,1:2);
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end

% Cells AUROC >= 0.5, never significant
for i = 1:2
    for ii = 1:2
        a = [results.(state).(group){i,ii}(:,3) results.(state).(group){i,ii}(:,4)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(state).(group){i,ii}(:,8)==0 & results.(state).(group){i,ii}(:,9)==0)); % Find never significant values
        c = b(any(b(:,3),3),:); % Keep only 1
%         c1 = any(c(:,1:2)>=0.5,3);% Find only values >=0.5
%         c1(:,3) = c1(:,1).*c1(:,2);
%         d = c(any(c1(:,3),3),:); % Keep only 1
%         d(:,3) = [];
        d = c(:,1:2);
        
        [~,p] = ttest(d(:, 1), d(:, 2));
        signr = signrank(d(:, 1), d(:, 2)); % Paired
        % ranksum(d(:,1), c(:,1))
        general_stats.mean{i,ii} = mean(d);
        general_stats.median{i,ii} = median(d);
%         diff(median(d))
%         diff(mean(d))
        med = median(d);
        me = mean(d);
        meddif = median(diff(d,1,2));
        medif = mean(diff(d,1,2));
        
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
            
            a = [];
            b = [];
            c = [];
            d = [];
        end
    end
end


%% SPONTANEOUS ACTIVITY %
 %----------------------%
 
results.SP.DATA_HPS_sig = cell(2,2);

for a = 1:2
    for b = 1:2
        idx_A = ismember(results.SP.DATA{a,b}(:,1:2), results.HPS.DATA_sig{a,b}(:,1:2), 'rows');
        results.SP.DATA_HPS_sig{a,b} = results.SP.DATA{a,b}(idx_A,:);
    end
end
% Group of significant cells
a = [];
a1 = [];
b = [];
b1 = [];

a = reshape(results.HPS.DATA_sig',[4,1]);
a1 = cell2mat(a);
b = reshape(results.SP.DATA_HPS_sig',[4,1]);
b1 = cell2mat(b);
groups_sig = [ones(size(a{1,1},1),1); 2*ones(size(a{2,1},1),1); 3*ones(size(a{3,1},1),1); 4*ones(size(a{4,1},1),1)];
HPS_diff = a1(:,4)-a1(:,3);
SP_diff = b1(:,4)-b1(:,3);
sig_diff_tot = [HPS_diff SP_diff groups_sig];


% All cells
a = [];
a1 = [];
b = [];
b1 = [];
a = reshape(results.HPS.DATA',[4,1]);
a1 = cell2mat(a);
b = reshape(results.SP.DATA',[4,1]);
b1 = cell2mat(b);
groups = [ones(size(a{1,1},1),1); 2*ones(size(a{2,1},1),1); 3*ones(size(a{3,1},1),1); 4*ones(size(a{4,1},1),1)];
HPS_diff = a1(:,4)-a1(:,3);
SP_diff = b1(:,4)-b1(:,3);
diff_tot = [HPS_diff SP_diff groups];

% Plot for all cells
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

% Plot for only significant cells
figure('color','w')
clf
subplot(221), hold on
    plot(sig_diff_tot(sig_diff_tot(:,3)==1, 2), sig_diff_tot(sig_diff_tot(:,3)==1, 1),'.r','MarkerSize',8), lsline
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['CCI LP (n=', num2str(sum(sig_diff_tot(:,3)==1)), ')']),
    xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square
    xlabel('SP'), ylabel('HPS')
subplot(222), hold on
    plot(sig_diff_tot(sig_diff_tot(:,3)==2, 2), sig_diff_tot(sig_diff_tot(:,3)==2, 1),'.r','MarkerSize',8), lsline
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['CCI SAL (n=', num2str(sum(sig_diff_tot(:,3)==2)), ')']),
    xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square
subplot(223), hold on
    plot(sig_diff_tot(sig_diff_tot(:,3)==3, 2), sig_diff_tot(sig_diff_tot(:,3)==3, 1),'.r','MarkerSize',8), lsline
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['sham LP (n=', num2str(sum(sig_diff_tot(:,3)==3)), ')']),
    xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square
subplot(224), hold on
    plot(sig_diff_tot(sig_diff_tot(:,3)==4, 2), sig_diff_tot(sig_diff_tot(:,3)==4, 1),'.r','MarkerSize',8), lsline
    plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
    title(['sham SAL (n=', num2str(sum(sig_diff_tot(:,3)==4)), ')']),
    xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square