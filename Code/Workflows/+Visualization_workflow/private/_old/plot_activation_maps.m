function plot_activation_maps(METADATA, stimulus_filename, fig_filename)
    % This function plots calcium traces.

    %%
    % Read general configs
    GC = general_configs();
    
    % Load data
    deconvolved_traces_paths = METADATA.deconvolved_trace_path;
    data = load_deconvolved_traces(deconvolved_traces_paths, 'var_name','spikes', 'group_traces','by_ROI');
    n_ROIs = size(data,1);
    
    % Get frame rate
    frame_rate = unique(METADATA.frame_rate);
    
    % Get the different days for these data
    exp_sessions = METADATA.exp_session;
    [exp_sessions,~,exp_session_idx] = unique(exp_sessions, 'stable');
    n_sessions = length(exp_sessions);
    
    % Load stimulus
    STIMULUS = load_variable(stimulus_filename, 'STIMULUS');
    if strcmp(STIMULUS.code, 'SP')
        is_baseline = true;
    else
        is_baseline = false;
    end
    % Get time axis and timestamps
    T  = STIMULUS.time_axis;
    TS = STIMULUS.timestamps;
    if is_baseline
        EPOCHS = [1, length(T)];
    else
        % Get indices of timestamps in time axis
        TS_idx = zeros(size(TS));
        for its = 1:length(TS_idx)
            [~,TS_idx(its)] = min(abs(T-TS(its)));
        end
        EPOCHS = [1, TS_idx];
        EPOCHS = [EPOCHS(:), [TS_idx-1, TS_idx(1) + floor(GC.response_window * frame_rate)].'];
    end
    n_epochs = size(EPOCHS, 1);

    % Get firing rate before and after stimulus onset
    FR = zeros(n_ROIs, n_sessions);
    P = false(n_ROIs, n_sessions);
    for isess = 1:n_sessions
        for iroi = 1:n_ROIs
            % Get data
            d = data{iroi, 1}(:, exp_session_idx==isess);
            n_trials = size(d,2);
            % Get firing rate in epochs
            fr_this_ROI = zeros(n_trials, n_epochs);
            for iepoch = 1:n_epochs
                fr_this_ROI(:, iepoch) = sum(d(EPOCHS(iepoch,1):EPOCHS(iepoch,2),:), 1) ./ (EPOCHS(iepoch,2)-EPOCHS(iepoch,1)+1);
            end
            % Compute t-test
            if ~is_baseline
                p = ranksum(fr_this_ROI(:,1),fr_this_ROI(:,2));
                % Is p-value lower than 0.05?
                P(iroi, isess) = p < 0.05;
                % Keep FR change
                fr_mean = nanmean(fr_this_ROI,1);
                FR(iroi, isess) = log2(fr_mean(2) ./ fr_mean(1));
            else
                FR(iroi, isess) = nanmean(fr_this_ROI(:));
            end
        end
    end
    % Assign maximum value to Infs
    max_FR = nanmax(FR(isfinite(FR)));
    FR(isinf(FR)) = max_FR;
    FR(isnan(FR)) = 0;
    
    % Make map of ROIs
    ROI_info_path = unique(METADATA.ROI_info_path);
    ROI_info = load_variable(ROI_info_path{1}, 'ROIs');
    % Get all ROIs
    ROIs = cat(3,ROI_info.mask);
    frame_height = size(ROIs, 1);
    frame_width  = size(ROIs, 2);
    % Make a grid of pixels
    [columnsInImage, rowsInImage] = meshgrid(1:frame_height, 1:frame_width);
    % Get position of centroids
    CENTROIDS = zeros(n_ROIs, 2);
    for iroi = 1:n_ROIs
        % Get mask, its center and radius
        mask = ROIs(:, :, iroi);
        p = regionprops(mask, 'Centroid', 'Area');
        % Replace ROI with a circle?
        if strcmpi(GC.ROIs_activation_map_shape, 'approximate')
            if GC.ROIs_activation_map_radius < 1
                radius = sqrt(p.Area / pi);
            else
                radius = GC.ROIs_activation_map_radius;
            end
            ROIs(:, :, iroi) = (rowsInImage - p.Centroid(2)).^2 + (columnsInImage - p.Centroid(1)).^2 <= radius.^2;
        end
        % Keep position of centroids
        CENTROIDS(iroi, :) = p.Centroid;
    end
    
    % Get color limits
    CLIMITS = [nanmin(FR(:)), nanmax(FR(:))];
    if ~is_baseline
        if CLIMITS(1)>1, CLIMITS(1)=1; end
        if CLIMITS(2)<1, CLIMITS(2)=1; end
    end
    % Make colormap
    CMAP = divergent_colormap(CLIMITS(1), 0, CLIMITS(2), 0, 0, @(x) redblue(x));

    % Open figure
    FIG = figure('color','w', 'visible','off', 'units','normalized', 'outerposition',[0 0 1 1]);
    colormap(CMAP)
    
    ax = zeros(n_sessions,1);
    % Loop through sessions
    for isess = 1:n_sessions
        % Make image
        IM = NaN(frame_height, frame_width);
        for iroi = 1:n_ROIs
            IM(ROIs(:,:,iroi)==1) = FR(iroi, isess);
        end
        
        % Make axes
        ax(isess) = subaxis(1, n_sessions, isess, 1, 'm',.05, 's',0, 'sh',.02, 'Hold', 1);
        % Plot image (without shading)
        pcolor(IM);
        shading flat
        caxis(CLIMITS)
        
        % Make axis square and add title
        axis image
        title(['session ', num2str(exp_sessions(isess))], 'Fontsize',14)
        
        % Show ROI id
        for iroi = 1:n_ROIs
            if P(iroi, isess)
                str = ['\underline{\textbf{', num2str(iroi), '}}'];
                text_size = 13;
            else
                str = num2str(iroi);
                text_size = 11;
            end
            text(CENTROIDS(iroi,1), CENTROIDS(iroi,2), str, 'FontSize',text_size, 'Hor','center', 'Ver','middle', 'color','k', 'Interpreter','latex')
        end
    end
    % Change background color of axes
    set(ax, 'Color','k', 'Box','off', 'XTick',[], 'YTick',[], 'XColor','None', 'YColor','None', 'YDir','reverse')
        
    % Make colorbar
    pos = get(ax(end), 'position');
    CB = colorbar(ax(end));
    if is_baseline
        CB.Label.String = 'Spontaneous activity';
    else
        CB.Label.String = 'Before / After stimulus onset';
        ticks = CB.Ticks;
        ticks = unique([0:1:max(ticks), -(0:1:abs(min(ticks)))]);
        ticklabels = regexp(rats(2.^ticks), '\ +', 'split');
        ticklabels = ticklabels(2:end-1);
        ticklabels(ticks==0) = {'same'};
        CB.Ticks = ticks;
        CB.TickLabels = ticklabels;
    end
    drawnow; set(ax(end), 'Position', pos)

    % Add stimulus name at the top
    axes('Position',[0,0,1,1], 'XLim',[0 1], 'YLim',[0 1], 'Box','off', 'Visible','off', 'Units','normalized', 'Clipping','off', 'HitTest','off');
    text(.5, 1, STIMULUS.code, 'Interpreter','none', 'FontSize',18, 'FontWeight','bold', 'Hor','left', 'Ver','top');

    % Store figures in a multi-page pdf file
    export_fig(fig_filename, '-pdf', '-q101', FIG)
    
    % Close figures
    close(FIG)
