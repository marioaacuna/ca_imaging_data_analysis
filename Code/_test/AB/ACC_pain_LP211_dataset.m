%% EVOKED ACTIVITY %
 %-----------------%

% Load configs, names, plot toolbox
clear, clc

global GC
GC = general_configs();
GC.experiment_name = 'ACC_pain_LP-211';

stats_filename = get_filename_of('response_modulation_stats', 'ACC_pain_LP-211');
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');

% dataset = 'CCI' % 'sham', ''
% dataset_animal_names = SQL_database.read_table_where('2p_imaging_experiments','experimental_group',dataset);
% SQL_database.read_table_where('sessions', {'animal_ID','experimental_condition'}, {'CCI','ACC_pain_Lp-211'},{'experimental_group','experiment_name'})
% SQL_database.read_table_where('sessions', {'animal_ID','experimental_group','experimental_condition'}, 'ACC_pain_LP-211','experiment_name')
% SQL_database.read_table_where('sessions', {'animal_ID','experimental_condition'}, {'CCI','ACC_pain_Lp-211'},{'experimental_group','experiment_name'})
% dataset_animal_names = SQL_database.read_table_where('experiments', {'animal_ID','experimental_group'}, {'ACC_pain_LP-211',dataset},{'experiment_name','experimental_group'})

dataset_animal_names = SQL_database.read_table_where('sessions', {'animal_ID','experimental_group', 'experimental_condition'}, 'ACC_pain_LP-211','experiment_name');
dataset_animal_names = dataset_animal_names(:, {'animal_ID','experimental_group','compound'});
dataset_animal_names = unique(dataset_animal_names, 'rows', 'stable');


all_animal_names = fieldnames(STATS_response_modulation);
n_datasets = size(dataset_animal_names,1);

phases_order = {'baseline', 'active_1'};
LEG_labels = {'non-selective, non-modulated', 'selective, non-modulated', 'only selective before', 'only selective after'};

toggle_toolbox('_plotting', 'on')
%% Find cells AUROC/p and assign to results structure
close all

compound_group_names = {'LP211','SAL'};
experiment_group_names = {'CCI','sham'};
results = struct();
results.DATA = cell(length(compound_group_names), length(experiment_group_names));
results.DATA_p = results.DATA;

for iid = 1:n_datasets
    
    animal_ID = cell2mat(dataset_animal_names{iid,1});
    stats = STATS_response_modulation.(animal_ID);
    
    [cells_ID, ~, rows_idx] = unique(stats.ROI_id);
    n_cells = length(cells_ID);
    
    DATA = NaN(n_cells, length(phases_order)+1);
    DATA_p = NaN(n_cells, length(phases_order)+1);
    
    % Get info on AUROC and p-value of each cell
    for icell = 1:n_cells
        data = stats(rows_idx == icell, {'compound_phase', 'AUROC', 'p'});
        % Skip cells that were not active during all phases
        if height(data) ~= length(phases_order)
            continue
        end
        
        for phi = 1:length(phases_order)
            rows = ismember(table2cell(data(:, 'compound_phase')), phases_order{phi});
            if any(rows)
                AUROC = data{rows, 'AUROC'};
                DATA(icell, phi) = AUROC;
                
                p = data{rows, 'p'};
                DATA_p(icell, phi) = p;
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
        % -- -- -- -- -- -- -- DATA_p -- -- -- -- -- -- -- -- -- -- %
        % Group 5 = Always significantly responsive to the stimulus %
        % Group 6 = Never significant responsiveness                %
        % Group 7 = Lost significant responsiveness                 %
        % Group 8 = Gain significant responsiveness                 %
        %-----------------------------------------------------------%
        
        if DATA(icell,phi-1) > DATA(icell,phi) && DATA(icell,phi) > 0.5 % Decreased positive modulation by treatment
            DATA(icell, 3) = 1;
        elseif 0.5 <= DATA(icell,phi-1)&& DATA(icell,phi-1)<= DATA(icell,phi) % Increased positive modulation by treatment
            DATA(icell, 3) = 2;
        elseif DATA(icell,phi-1) <= DATA(icell,phi)&& DATA(icell,phi) <= 0.5 % Decreased negative modulation by treatment
            DATA(icell, 3) = 3;
        elseif 0.5 > DATA(icell,phi-1) && DATA(icell,phi-1) > DATA(icell,phi) % Increased negative modulation by treatment
            DATA(icell, 3) = 4;
        end
        % Assign group attribute to p values
        if DATA_p(icell,phi-1)<=0.05 && DATA_p(icell,phi)<=0.05 % Always significantly responsive to the stimulus
            DATA_p(icell, phi+1) = 5;
        elseif DATA_p(icell,phi-1)>0.05 && DATA_p(icell,phi)>0.05 % Never significant responsiveness
            DATA_p(icell, phi+1) = 6;
        elseif DATA_p(icell,phi-1)<=0.05 && DATA_p(icell,phi)>0.05 % Lost significant responsiveness
            DATA_p(icell, phi+1) = 7;
        elseif DATA_p(icell,phi-1)>0.05 && DATA_p(icell,phi)<=0.05 % Gain significant responsiveness
            DATA_p(icell, phi+1) = 8;
        end
        
        % Compute the difference POST-PRE for each cell
        DATA(icell, phi+2) = DATA(icell, phi)-DATA(icell, phi-1);
        DATA_p(icell, phi+2) = DATA_p(icell, phi)-DATA_p(icell, phi-1);
        
        % Compute the normalized difference [(POST-PRE)/(POST+PRE)] for each cell
        DATA(icell, phi+3) = (DATA(icell, phi)-DATA(icell, phi-1))/(DATA(icell, phi)+DATA(icell, phi-1));
        DATA_p(icell, phi+3) = (DATA_p(icell, phi)-DATA_p(icell, phi-1))/(DATA_p(icell, phi)+DATA_p(icell, phi-1));
    end
    
    % Remove NaNs
    rows = any(isnan(DATA)');
    DATA(rows, :) = [];
    DATA_p(rows, :) = [];
    
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
    results.DATA{row, col} = cat(1, results.DATA{row, col}, DATA);
    results.DATA_p{row, col} = cat(1, results.DATA_p{row, col}, DATA_p);
    
end
%% Info for results dataset

    % Avg dataset
    % results.DATA Table Structure:
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  x   | x | %
    % sham|  x   | x | %
    %------------------%
    
ncells = cellfun(@(x) size(x,1), results.DATA);
avg = cellfun(@(x) mean(x,1), results.DATA, 'UniformOutput',false);
avg_no_zeros = cellfun(@(x) (sum(x,1) ./ sum(x~=0,1)), results.DATA_sig, 'UniformOutput',false); % Average of non zeros values
cellfun(@(x) sum(x(:,6:7)~=0,1), results.DATA_sig,'UniformOutput',false); % Number of positive values in column 6 and 7 (i.e number of cells responding PRE or POST)
%% Isolate significant modulated cells
results.DATA_sig_ch = cell(2,2); % Significant changing selectivity before/after treatment
results.DATA_sig_st = cell(2,2); % Significant stable selectivity before/after treatment
results.DATA_sig = cell(2,2); % Significant selective cells
results.DATA_sig_post = cell(2,2); % Significant selective cells

% Add 'zeros-one' filter to discriminate significant modulated cells in DATA_p structure
for a = 1:2
    for b = 1:2
        for i = 1:length(results.DATA_p{a,b})
            if results.DATA_p{a,b}(i,1) <= 0.05;
                results.DATA_p{a,b}(i,6) = 1;
            else results.DATA_p{a,b}(i,1) > 0.05;
                results.DATA_p{a,b}(i,6) = 0;
            end
            if results.DATA_p{a,b}(i,2) <= 0.05;
                results.DATA_p{a,b}(i,7) = 1;
            else results.DATA_p{a,b}(i,2) > 0.05;
                results.DATA_p{a,b}(i,7) = 0;
            end
        end
    end
end

% Select only significant modulated cells in DATA structure with 'zeros-one' vector filter
for a = 1:2
    for b = 1:2
        results.DATA{a,b}(:,6:7) = results.DATA_p{a,b}(:,6:7).*results.DATA{a,b}(:,1:2);
    end
end

% Keep significant cells
results.DATA_sig = results.DATA;
sig_cell_ID = cell(2,2);
for a = 1:2
    for b = 1:2
        for i = 1:length(results.DATA_p{a,b})
            if results.DATA{a,b}(i,6) == 0 && results.DATA{a,b}(i,7) == 0 % Always non modulated
                results.DATA_sig{a,b}(i,8) = 0;
                sig_cell_ID{a,b}(i,1) = 0;
            elseif results.DATA{a,b}(i,6) > 0 && results.DATA{a,b}(i,7) > 0 % Always modulated
                results.DATA_sig{a,b}(i,8) = 1;
                sig_cell_ID{a,b}(i,1) = 1;
            else
                results.DATA_sig{a,b}(i,8) = 1;
                sig_cell_ID{a,b}(i,1) = 1;
            end
        end
        results.DATA_sig{a,b} = results.DATA_sig{a,b}(any(results.DATA_sig{a,b}(:,8),8),:); % Remove rows with zeros in the index column
        results.DATA_sig{a,b}(:,8) = []; % Remove index column
    end
end

% Keep significant unstable cells (gain/loss selectivity)
results.DATA_sig_ch = results.DATA;
a = [];
b = [];
for a = 1:2
    for b = 1:2
        for i = 1:length(results.DATA_p{a,b})
            if results.DATA{a,b}(i,6) == 0 && results.DATA{a,b}(i,7) == 0 % Always non modulated
                results.DATA_sig_ch{a,b}(i,8) = 0;
            elseif results.DATA{a,b}(i,6) > 0 && results.DATA{a,b}(i,7) > 0 % Always modulated
                results.DATA_sig_ch{a,b}(i,8) = 0;
            else
                results.DATA_sig_ch{a,b}(i,8) = 1;
            end
        end
        results.DATA_sig_ch{a,b} = results.DATA_sig_ch{a,b}(any(results.DATA_sig_ch{a,b}(:,8),8),:); % Remove rows with zeros in the index column
        results.DATA_sig_ch{a,b}(:,8) = []; % Remove index column
    end
end

% Keep significant stable cells (always significant/not significant)
results.DATA_sig_st = results.DATA;
a = [];
b = [];
for a = 1:2
    for b = 1:2
        for i = 1:length(results.DATA_p{a,b})
            if results.DATA{a,b}(i,6) == 0 && results.DATA{a,b}(i,7) == 0 % Always non modulated
                results.DATA_sig_st{a,b}(i,8) = 1;
            elseif results.DATA{a,b}(i,6) > 0 && results.DATA{a,b}(i,7) > 0 % Always modulated
                results.DATA_sig_st{a,b}(i,8) = 1;
            else
                results.DATA_sig_st{a,b}(i,8) = 0;
            end
        end
        results.DATA_sig_st{a,b} = results.DATA_sig_st{a,b}(any(results.DATA_sig_st{a,b}(:,8),8),:); % Remove rows with zeros in the index column
        results.DATA_sig_st{a,b}(:,8) = []; % Remove index column
    end
end

% Keep significant cells during POST
results.DATA_sig_post = results.DATA;
a = [];
b = [];
for a = 1:2
    for b = 1:2
        for i = 1:length(results.DATA_p{a,b})
            if results.DATA{a,b}(i,6) == 0 && results.DATA{a,b}(i,7) > 0 % Always non modulated
                results.DATA_sig_post{a,b}(i,8) = 1;
            else
                results.DATA_sig_post{a,b}(i,8) = 0;
            end
        end
        results.DATA_sig_post{a,b} = results.DATA_sig_post{a,b}(any(results.DATA_sig_post{a,b}(:,8),8),:); % Remove rows with zeros in the index column
        results.DATA_sig_post{a,b}(:,8) = []; % Remove index column
    end
end
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
        H = zeros(size(DATA,1), 1);
        
        for ix = 1:ncells(irow,icol)
            x = results.DATA{irow, icol}(ix, 1);
            y = results.DATA{irow, icol}(ix, 2);
            p_pre = results.DATA_p{irow, icol}(ix, 1);
            p_post = results.DATA_p{irow, icol}(ix, 2);
            
            
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

%             groups = table
%             if group == 1
%                 group1 = cat(1, group1, [x,y]);
%             elseif group == 2
%                 group2 = cat(1, group2, [x,y]);
%             elseif group == 3 && x <= 0.5
%                 group3 = cat(1, group3, [x,y]);
%             elseif group == 3 && x > 0.5
%                group3a = cat(1, group3a, [x,y]);
%             elseif group == 4 && y <= 0.5
%                 group4 = cat(1, group4, [x,y]);
%             elseif group == 4 && y > 0.5
%                 group4a = cat(1, group4a, [x,y]);
%             end
        
            % Store handle for legend
            if isnan(LEG_H(group))
                LEG_H(group) = H(ix);
            end
            LEG_n(group) = LEG_n(group) + 1;
        end
        
        % Fix axes appearance
        %     lims = DATA(:);
        %     lims = [min(lims), max(lims)];
        %     edge = max(abs(lims - 0.5));
        %     lims = [.5-edge, .5+edge];
        %     padding = 0.02;
        %     padding = range(lims) * padding;
        %     set(ax, 'Xlim',[lims(1)-padding, lims(2)+padding], 'ylim',[lims(1)-padding, lims(2)+padding])
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

keyboard
%% Define data groups to plot

    % DATA/DATA_sig/DATA_sig_ groups
    %------------------%
    %     |LP-211|SAL| %
    %  CCI|  x   | x | %
    % sham|  x   | x | %
    %------------------%

% Specify the group from where data is used (results.DATA / results.DATA_sig_st)
group = 'DATA_sig'; % DATA; DATA_sig; DATA_sig_ch; DATA_sig_st, DATA_sig_post

DATA_CCI_PrePost_norm = [results.DATA{1,1}(:,5); results.DATA{2,1}(:,5)]; % AUROC CCI (POST-PRE)/(POST+PRE)
DATA_sham_PrePost_norm = [results.DATA{1,2}(:,5); results.DATA{2,2}(:,5)]; % AUROC sham (POST-PRE)/(POST+PRE)
DATA_CCI_PRE = [results.(group){1,1}(:,1); results.(group){1,2}(:,1)]; % AUROC CCI PRE treatment
DATA_sham_PRE = [results.(group){2,1}(:,1); results.(group){2,2}(:,1)]; % AUROC sham PRE treatment


% Settings for plotting normalized AUROC kdensity
n_bins = 100;
bin_edges = linspace(0, 1, n_bins);
bin_centers = bin_edges(1:end-1) + (bin_edges(2)-bin_edges(1))/2;
bin_centers = bin_centers';
% bin_edges_smooth = linspace(0, 1, 1000);
% bin_centers_smooth = bin_edges_smooth(1:end-1) + (bin_edges_smooth(2)-bin_edges_smooth(1))/2;
%% Keep only significant responsive cells (OLD version)
% DATA_p_CCI_PRE = [results.DATA_p{1,1}(:,1) results.DATA_p{1,1}(:,2); results.DATA_p{1,2}(:,1) results.DATA_p{1,2}(:,2)];
% DATA_p_sham_PRE = [results.DATA_p{2,1}(:,1); results.DATA_p{2,2}(:,1)];

    % Keep only significant responsive cells based on DATA_p values
    DATA_CCI_PRE_sig = (group)_CCI_PRE;
    DATA_CCI_PRE_sig = [];

    for nel = 1:length(DATA_CCI_PRE)
       if DATA_p_CCI_PRE(nel,1)>0.05
           DATA_CCI_PRE_sig(nel,1) = DATA_CCI_PRE(nel);
       else
           DATA_CCI_PRE_sig(nel,1) = zeros;
       end
    end

    %rows = any(isnan(DATA_CCI_PRE_sig)');
    rows = ~any(DATA_CCI_PRE_sig)';
    %DATA_CCI_PRE_sig(rows, :) = [];
    DATA_CCI_PRE_sig = DATA_CCI_PRE_sig(DATA_CCI_PRE_sig>0);    

    DATA_sham_PRE_sig = DATA_sham_PRE;
    DATA_sham_PRE_sig = [];

    for nel = 1:length(DATA_sham_PRE)
       if DATA_p_sham_PRE(nel)>0.05
           DATA_sham_PRE_sig(nel,1) = DATA_sham_PRE(nel);
       else
           DATA_sham_PRE_sig(nel,1) = zeros;
       end
    end

    rows = any(isnan(DATA_sham_PRE_sig)');
    %DATA_sham_PRE_sig(rows, :) = [];
    DATA_sham_PRE_sig = DATA_sham_PRE_sig(DATA_sham_PRE_sig>0);
%% Plot specific datasets
% Define the dataset to plot
all_datasets = struct();
all_datasets.dataset1 = {DATA_CCI_PRE,DATA_sham_PRE};
all_datasets.dataset2 = {DATA_CCI_PRE,results.(group){1,1}(:,2),results.(group){1,2}(:,2)};
all_datasets.dataset3 = {DATA_sham_PRE,results.(group){2,1}(:,2),results.(group){2,2}(:,2)};
all_datasets.dataset4 = {DATA_CCI_PRE,results.(group){1,1}(:,2),results.(group){1,2}(:,2),DATA_sham_PRE,results.(group){2,1}(:,2),results.(group){2,2}(:,2)};

dataset = {'dataset1'; 'dataset2'; 'dataset3'; 'dataset4'};
coord = struct();
coord.dataset1 = [];
coord.dataset2 = [];
coord.dataset3 = [];
coord.dataset4 = [];

for idata = 1:numel(dataset)
    
    figure;
    clf
    hold on
    data = all_datasets.(dataset{idata});
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
        title('PRE'),legend('CCI','sham')
    elseif idata == 2
        title('CCI'),legend('PRE','LP','SAL')
    elseif idata == 3
        title('sham'),legend('PRE','LP','SAL')
    elseif idata == 4
        title('CCI/sham'),legend('CCI PRE','CCI LP','CCI SAL','sham PRE','sham LP','sham SAL')
    end
end

% title('PRE'),legend('CCI','sham')
% title('CCI'),legend('PRE','LP','SAL')
% title('sham'),legend('PRE','LP','SAL')
% title('CCI/sham'),legend('CCI PRE','CCI LP','CCI SAL','sham PRE','sham LP','sham SAL')
% title('PRE significant'),legend('CCI','sham')
%% Summary values of grouped cells
% Cells numbers
results_values = struct();

ncells = struct();
ncells.tot = cellfun(@(x) size(x,1), results.(group));
ncells.lost = cellfun(@(x) sum(x(:,6)~=0,1), results.(group),'UniformOutput',false); % Number of positive values in column 6 (i.e number of cells responding during PRE, loosing selectivity during POST)
ncells.gain = cellfun(@(x) sum(x(:,7)~=0,1), results.(group),'UniformOutput',false); % Number of positive values in column 7 (i.e number of cells responding during POST, gaining selectivity after PRE)

% Avg AUROC values
avg_AUROC = cellfun(@(x) (sum(x(:,6:7),1) ./ sum(x(:,6:7)~=0,1)), results.(group), 'UniformOutput',false); % Average of non zeros values during PRE and POST (col 6 and 7)

results_values.ncells = ncells;
results_values.avg_AUROC = avg_AUROC;
%% Statistics on the distribution
% Wilcoxon rank sum test for unpaired distribution (CCI_PRE/sham_PRE)
% Wilcoxon signed rank test for paired distribution (CCI_PRE/LP | sham_PRE/LP)
p = [];
p_all = []

% CCI_PRE vs sham_PRE
x = DATA_CCI_PRE; % CCI_LP PRE + CCI_SAL PRE
y = DATA_sham_PRE; % sham_LP PRE + sham_SAL PRE
p1 = ranksum(x,y);
p = cat(1,p, p1);
p_all = cat(1,p,p_all);

% CCI_LP PRE vs CCI_LP POST
x = results.(group){1,1}(:,1); % CCI_LP PRE
y = results.(group){1,1}(:,2); % CCI_LP POST
p1 = signrank(x,y);
p = cat(1,p,p1);

% CCI_SAL PRE vs CCI_SAL POST
x = results.(group){1,2}(:,1); % CCI_SAL PRE
y = results.(group){1,2}(:,2); % CCI_SAL POST
p1 = signrank(x,y);
p = cat(1,p,p1);

% CCI_PRE vs CCI_LP POST
x = DATA_CCI_PRE;
y = results.(group){1,1}(:,2); % CCI_LP POST
p1 = ranksum(x,y);
p_all = cat(1,p_all,p1);

% sham_LP PRE vs sham_LP POST
x = results.(group){2,1}(:,1); % sham_LP PRE
y = results.(group){2,1}(:,2); % sham_LP POST
p1 = signrank(x,y);
p = cat(1,p,p1);

% sham_PRE vs sham_LP POST
x = DATA_sham_PRE;
y = results.(group){2,1}(:,2); % sham_LP POST
p1 = ranksum(x,y);
p_all = cat(1,p_all,p1);

p = [p,p_all]; % Col 1 = p-values paired tests, col 2 = p-values unpaired tests

% Prepare BoxPlot data
% Concat data
X = [DATA_CCI_PRE; DATA_sham_PRE; ...
    results.(group){1,1}(:,1); results.(group){1,1}(:,2); ...
    results.(group){2,1}(:,1); results.(group){2,1}(:,2); ...
    results.(group){1,2}(:,1); results.(group){1,2}(:,2); ...
    results.(group){2,2}(:,1); results.(group){2,2}(:,2)];
% Grouping variables
G = [ones(size(DATA_CCI_PRE)); 2*ones(size(DATA_sham_PRE)); ...
    3*ones(size(results.(group){1,1}(:,1))); 4*ones(size(results.(group){1,1}(:,2))); ...
    5*ones(size(results.(group){2,1}(:,1))); 6*ones(size(results.(group){2,1}(:,2))); ...
    7*ones(size(results.(group){1,2}(:,1))); 8*ones(size(results.(group){1,2}(:,2)));
    9*ones(size(results.(group){2,2}(:,1))); 10*ones(size(results.(group){2,2}(:,2)))];

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
% Significant changing
a = [];
b = [];
c = [];

% Plot all cells of group before after
for i = 1:2
    for ii = 1:2
        a = [results.(group){i,ii}(:,1) results.(group){i,ii}(:,2)];
        
        figure
        hold on
        for icell = 1:size(results.(group){i,ii},1)
            if results.(group){i,ii}(icell,7) > 0 % results.(group){i,ii}(icell,6) == 0 && 
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

% Cells > 0.5, significant PRE/POST
for i = 1:2
    for ii = 1:2
        a = [results.(group){i,ii}(:,1) results.(group){i,ii}(:,2)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(group){i,ii}(:,6)>0 & results.(group){i,ii}(:,7)>0)); % Find significant values during PRE/POST
        % b1 = results.(group){i,ii}(:,6:7)>0; % Find significant values during PRE/POST
        % b1(:,3) = b1(:,1).*b1(:,2);
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>0.5,3);% Find only values >0.5
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

% Cells > 0.5, significant only during POST not PRE
for i = 1:2
    for ii = 1:2
        a = [results.(group){i,ii}(:,1) results.(group){i,ii}(:,2)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(group){i,ii}(:,6)==0 & results.(group){i,ii}(:,7)>0)); % Find significant values during POST
        % b1 = results.(group){i,ii}(:,7)>0; % Find significant values during POST
        % b1(:,3) = b1(:,1).*b1(:,2);
%         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>0.5,3);% Find only values >0.5
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

% Cells > 0.5, significant only during PRE not POST
for i = 1:2
    for ii = 1:2
        a = [results.(group){i,ii}(:,1) results.(group){i,ii}(:,2)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(group){i,ii}(:,6)>0 & results.(group){i,ii}(:,7)==0)); % Find significant values during PRE
        % b1 = results.(group){i,ii}(:,7)>0; % Find significant values during POST
        % b1(:,3) = b1(:,1).*b1(:,2);
%         b = cat(2,b1(:,1:2),b1(:,3).*b1(:,4));
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>0.5,3);% Find only values >0.5
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

% Cells > 0.5, never significant
for i = 1:2
    for ii = 1:2
        a = [results.(group){i,ii}(:,1) results.(group){i,ii}(:,2)]; % GROUP A = BEFORE/AFTER
        b = cat(2,a,(results.(group){i,ii}(:,6)==0 & results.(group){i,ii}(:,7)==0)); % Find never significant values
        c = b(any(b(:,3),3),:); % Keep only 1
        c1 = any(c(:,1:2)>0.5,3);% Find only values >0.5
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
%% SPONTANEOUS ACTIVITY %
 %----------------------%
 
% Load configs, names
global GC
GC = general_configs();
GC.experiment_name = 'ACC_pain_LP-211';
stats_filename = get_filename_of('response_modulation_stats', 'ACC_pain_LP-211');
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');

% Load dataset animals from the server
dataset_animal_names = SQL_database.read_table_where('sessions', {'animal_ID','experimental_group', 'experimental_condition'}, 'ACC_pain_LP-211','experiment_name');
dataset_animal_names = dataset_animal_names(:, {'animal_ID','experimental_group','compound'});
dataset_animal_names = unique(dataset_animal_names, 'rows', 'stable');
all_animal_names = fieldnames(STATS_response_modulation);
n_datasets = size(dataset_animal_names,1);

% Load variables for spontaneous activity
stats_filename = get_filename_of('spontaneous_activity_stats', 'ACC_pain_LP-211');
% STATS_spontaneous_activity = load_variable(stats_filename, 'STATS_spontaneous_activity');
DISTRIBUTIONS = load_variable(stats_filename, 'DISTRIBUTIONS');

    % a = DISTRIBUTIONS.AB_1(:, 'mean_activity')
    % a = DISTRIBUTIONS.AB_1{:, 'mean_activity'}
    % a = cell2mat(DISTRIBUTIONS.AB_1{:, 'mean_activity'})
    % SP = cell2mat(DISTRIBUTIONS.AB_1{:, 'mean_activity'})
    % SP = SP - SP(:,1);
    % 
    % plot(SP(:, 2:3)', '-ok')
    % xlim([0, 3])
    % 
    % SP = cell2mat(DISTRIBUTIONS.AB_1{:, 'mean_activity'});
    % plot(SP', '-ok')
    % xlim([0, 4])

% Prepare for collection spontaneous mean activity for each animal
nel = size(fieldnames(DISTRIBUTIONS),1);
ID = {fieldnames(DISTRIBUTIONS)};

% Find max number of imaging session and tot number of cells
sessions_cells_id = [];
for i = 1:nel
    animID = cell2mat(ID{1,1}(i));
    animal_session_cell_id = [size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),2), ...
        size(cell2mat(DISTRIBUTIONS.(animID){:,'mean_activity'}),1) ...
        , i];
    sessions_cells_id = cat(1, sessions_cells_id, animal_session_cell_id);
end
max_sessions = max(sessions_cells_id(:,1),[], 1);
max_cells = sum(sessions_cells_id(:,2));

% Concatenate all spontaneous/SP normalized values for all animals
SP = [];
normSP = [];
SP_tot = [];
normSP_tot = [];
% SP_tot = NaN(ncells_tot, max_sessions+1);
% normSP_tot = NaN(ncells_tot, max_sessions+1);

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
%     figure
%     clf
%     subplot(1,2,1),
%     plot(SP', '-ok')
%     xlim([0, (size(normSP,2)+1)])
%     title(animID,'Interpreter','none')
%     subplot(1,2,2),
%     plot(normSP(:, 2:(size(normSP,2)))', '-ok')
%     title({animID ' Normalized to baseline'},'Interpreter','none')
%     xlim([0, (size(normSP,2))])
    
    % Add animal ID
    SP(:,end+1) = i;
    normSP(:,end+1) = i;
    
    % Concatenate all the cells for all the animals
    SP_tot = cat(1,SP_tot,SP);
    normSP_tot = cat(1,normSP_tot,normSP);
end

% Order and find animal_ID values based on 'experimental_group' / 'compound'
a = dataset_animal_names{ismember(dataset_animal_names.experimental_group, 'CCI') & ...
    ismember(dataset_animal_names.compound, 'LP-211'), 'animal_ID'}; % Find animal_ID which are also CCI AND LP-211
dataset_animal_names = natsortrows(dataset_animal_names, 1); % Sort rows based on ascending column 1
dataset_animal_names = natsortrows(dataset_animal_names, [-2, -3]);

% Find same group / compound
group1 = [];
group = [];

for i = 1:nel
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

T = table(group);
dataset_animal_names = [dataset_animal_names T]; % Concatenate the group to the dataset

% Add group ID to SP_tot
for i = 1:length(SP_tot)
   SP_tot(i,7) = dataset_animal_names.group(SP_tot(i,6));
end

% Make statistics based on grouping
SP_tot = array2table(SP_tot,'VariableNames',{'Baseline','Active1','Active2','Active3','Active4','animal_ID','animal_group'});
grpstats(SP_tot,'animal_group', 'median')
%% Compare SP and evoked of all cells all animals

nel = size(fieldnames(DISTRIBUTIONS),1);
ID = {fieldnames(DISTRIBUTIONS)};

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
        
        % Calculate SP difference: POST - PRE
        SP = (cell2mat(DISTRIBUTIONS.(animID){cel_id,'mean_activity'}));
        SP_diff = SP(2) - SP(1);
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

diff_tot = cat(2, SP_diff_tot(:,1), evoked_diff_tot); % Create single variable containing SP, EV, GROUP info for each cell

    % Plot all cells SP vs EVOKED
    % plot(diff_tot(:,2),diff_tot(:,1),'.','MarkerSize',8)
    % plot(diff_tot(:,2),diff_tot(:,1),'o','MarkerSize',4)

    % Plot cells SP vs EVOKED by group
    % figure
    % clf
    % hold on
    % gscatter(diff_tot(:,2),diff_tot(:,1),diff_tot(:,3),'rycb','....',8) % Plot scatter by group in column 3
    % xlim([-0.1 0.06]), ylim([-0.1 0.06])
    % plot([0 0],[-0.1 0.06], 'k'), plot([-0.1 0.06],[0 0], 'k')

figure('color','w')
clf
subplot(221), hold on
plot(diff_tot(diff_tot(:,3)==1, 2), diff_tot(diff_tot(:,3)==1, 1),'.r','MarkerSize',8), lsline
plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
title(['CCI LP (n=', num2str(sum(diff_tot(:,3)==1)), ')']),
xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square
xlabel('SP'), ylabel('HPS')

subplot(222), hold on
plot(diff_tot(diff_tot(:,3)==2, 2), diff_tot(diff_tot(:,3)==2, 1),'.r','MarkerSize',8), lsline
plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
title(['CCI SAL (n=', num2str(sum(diff_tot(:,3)==2)), ')']),
xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square

subplot(223), hold on
plot(diff_tot(diff_tot(:,3)==3, 2), diff_tot(diff_tot(:,3)==3, 1),'.r','MarkerSize',8), lsline
plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
title(['sham LP (n=', num2str(sum(diff_tot(:,3)==3)), ')']),
xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square

subplot(224), hold on
plot(diff_tot(diff_tot(:,3)==4, 2), diff_tot(diff_tot(:,3)==4, 1),'.r','MarkerSize',8), lsline
plot([0 0],[-100 100], 'k', 'YLimInclude','off'), plot([-100 100],[0 0], 'k', 'XLimInclude','off')
title(['sham SAL (n=', num2str(sum(diff_tot(:,3)==4)), ')']),
xlim([-0.1 0.06]), ylim([-0.1 0.06]), axis equal square
