clear, clc

stats_filename = get_filename_of('response_modulation_stats', 'ACC_pain_LP-211');
STATS_response_modulation = load_variable(stats_filename, 'STATS_response_modulation');

all_animal_names = fieldnames(STATS_response_modulation);
n_datasets = length(all_animal_names);

phases_order = {'baseline', 'active_1'};
LEG_labels = {'non-selective, non-modulated', 'selective, non-modulated', 'only selective before', 'only selective after'};

toggle_toolbox('_plotting', 'on')

%%
close all

for iid = 1:n_datasets
    animal_ID = all_animal_names{iid};
    stats = STATS_response_modulation.(animal_ID);
    
    [cells_ID, ~, rows_idx] = unique(stats.ROI_id);
    n_cells = length(cells_ID);
    
    % Get info on AUROC and p-value of each cell
    DATA = NaN(n_cells, 2);
    DATA_p = NaN(n_cells, 2);
    for icell = 1:n_cells
        data = stats(rows_idx == icell, {'compound_phase', 'AUROC', 'p'});
        for phi = 1:length(phases_order)
            rows = ismember(table2cell(data(:, 'compound_phase')), phases_order{phi});
            if any(rows)
                AUROC = data{rows, 'AUROC'};
                DATA(icell, phi) = AUROC;
                
                p = data{rows, 'p'};
                DATA_p(icell, phi) = p;
            end
        end
    end
    % Remove NaNs
    rows = any(isnan(DATA)');
    DATA(rows, :) = [];
    DATA_p(rows, :) = [];

    % Open figure and draw axes
    figure('color','w')
    clf; hold on;
    ax = subplot(1,1,1);
    plot([0, 1], [0, 1], 'k')
    plot([0.5, .5], [0, 1], 'k')
    plot([0, 1], [0.5, .5], 'k')

    % Plot each point with a different style to highlight the group they belong to
    H = zeros(size(DATA,1), 1);
    LEG_H = NaN(1, 4);
    LEG_n = zeros(1, 4);
    for irow = 1:length(H)
        % Get coordinates and p-values
        x = DATA(irow, 1);
        y = DATA(irow, 2);
        p_pre = DATA_p(irow, 1);
        p_post = DATA_p(irow, 2);
        
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
        H(irow) = plot(x, y, 'Marker',marker, 'MarkerFaceColor',mfc, 'MarkerEdgeColor',mec, 'MarkerSize',markersize, 'Linestyle','none', 'Linewidth',lw);
        
        % Store handle for legend
        if isnan(LEG_H(group))
            LEG_H(group) = H(irow);
        end
        LEG_n(group) = LEG_n(group) + 1;
    end

    % Fix axes appearance
    lims = DATA(:);
    lims = [min(lims), max(lims)];
    edge = max(abs(lims - 0.5));
    lims = [.5-edge, .5+edge]; 
    padding = 0.02;
    padding = range(lims) * padding;
    set(ax, 'Xlim',[lims(1)-padding, lims(2)+padding], 'ylim',[lims(1)-padding, lims(2)+padding])
    axis square
    set(ax, 'FontSize',12, 'TickDir','out')
    xlabel('Selectivity to HPS before LP-211', 'FontSize',13), ylabel('Selectivity to HPS after LP-211', 'FontSize',13)
    title([animal_ID, ' (n=', num2str(size(DATA, 1)), ')'], 'FontSize',16, 'Interpreter','none')
    
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
