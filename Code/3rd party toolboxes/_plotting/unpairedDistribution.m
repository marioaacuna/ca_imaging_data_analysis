function [hData, hAvg, X, Y, M, E, I] = unpairedDistribution(data, groups, varargin)
%% TODO
% Fix 'rect' error when on y-axis
%

%% Parse inputs
% Make data a vector
data = data(:);

% Check inputs
if ~exist('groups','var') || isempty(groups)
    groups = ones(size(data));
elseif length(groups)==1
    groups = repmat(groups, size(data));
end

% Make groups a vector
groups = groups(:);

% Sort data by group
data = sortrows([data, groups], 2);
groups = data(:, 2);
data = data(:, 1);

% Remove NaNs
invalid_data = isnan(data) | isnan(groups);
data(invalid_data) = [];
groups(invalid_data) = [];
original_size = size(data);

% Find the number of groups
numPoints = length(data);
[groups,~,group_idx] = unique(groups,'stable'); % Use 1...n indices for rows labels
nGroups = length(groups);
if nGroups == 0
    groups = 1;
    group_idx = ones(size(data,1),1);
    nGroups = 1;
end

% Set defaul parameters
p = inputParser();
p.addParameter('Labels', '');
p.addParameter('Slices', (log(numPoints)/log(150)+3*(numPoints).^(1/3))');
p.addParameter('Width', 0.5);
p.addParameter('Compression', 4);
p.addParameter('Shape', 'v');  % '^' or 'v', i.e upward or downward
p.addParameter('Axis', 'x');
p.addParameter('Avg', 'median');  % 'mean' (mean+-SEM), or 'median' (+-quartiles), or 'none'
p.addParameter('AvgColor', {'k'});  % color of line for group average. 'same' is a special option that uses the color specified in 'MarkerFaceColor'
p.addParameter('Spread', 'iqr');  % 'iqr' (inter-quartile range), 'sem' (standard error of the mean), 'sd' (standard deviation), or '95ci' (95% confidence interval)
p.addParameter('ErrorBar', 'rect');  % 'line', 'rect', or 'none'
p.addParameter('ErrorBarColor', [.85, .85, .85]);  % color of error bars
p.addParameter('TickLength', 1);  % length ratio of mean and errorbar ticks, or width of 'rect'
p.addParameter('Marker', {'o'});
p.addParameter('MarkerSize', 6);
p.addParameter('MarkerFaceColor', {'k'});
p.addParameter('MarkerEdgeColor', '');
p.addParameter('MarkerFaceAlpha', 1);
p.addParameter('MarkerEdgeAlpha', 1);
p.addParameter('LineWidth', 1);
p.addParameter('Offset', 0);
p.addParameter('X_distance', 1);
p.addParameter('Visible', 'on');
p.addParameter('CollapseGroups', false);
p.addParameter('ShowPoints', true);
p.addParameter('ReshapeReturnedPoints', true);
parse(p, varargin{:})

% Unpack variables
Avg = lower(p.Results.Avg);
AvgColor = p.Results.AvgColor;
Axis = lower(p.Results.Axis);
CollapseGroups = p.Results.CollapseGroups;
Compression = p.Results.Compression;
ErrorBar = lower(p.Results.ErrorBar);
ErrorBarColor = p.Results.ErrorBarColor;
Labels = p.Results.Labels;
LineWidth = p.Results.LineWidth;
Marker = p.Results.Marker;
MarkerEdgeAlpha = p.Results.MarkerEdgeAlpha;
MarkerEdgeColor = p.Results.MarkerEdgeColor;
MarkerFaceAlpha = p.Results.MarkerFaceAlpha;
MarkerFaceColor = p.Results.MarkerFaceColor;
MarkerSize = p.Results.MarkerSize;
Offset = p.Results.Offset;
Shape = p.Results.Shape;
ShowPoints = p.Results.ShowPoints;
Slices = p.Results.Slices;
ReshapeReturnedPoints = p.Results.ReshapeReturnedPoints;
Spread = lower(p.Results.Spread);
TickLength = p.Results.TickLength;
Visible = p.Results.Visible;
Width = p.Results.Width;
X_distance = p.Results.X_distance;

% Check that some variables are of the correct data type
if ischar(Marker), Marker={Marker}; end
if ischar(MarkerFaceColor), MarkerFaceColor={MarkerFaceColor}; end
if isempty(MarkerEdgeColor)
    MarkerEdgeColor = MarkerFaceColor;
elseif ischar(MarkerEdgeColor)
    MarkerEdgeColor = {MarkerEdgeColor};
end
if ischar(Labels), Labels={Labels}; end

% Decide how to proceed for plotting / statistics
groups_to_plot = groups;
groups_to_plot_idx = group_idx;
nGroups_to_plot = nGroups;
if CollapseGroups
    groups = 1;
    group_idx = ones(size(data,1),1);  % Ignore groups for statistics
    nGroups = 1;
end

% Check that the type of data dispersion is of the right type according to the
% average type
if ~isempty(Spread)
    if strcmp(Avg,'mean')   &&  strcmp(Spread,'iqr'), warning('unpairedDistribution:mean','Your current choice of mean and inter-quartile range is not correct! Please use either SEM or SD. Errorbars will not be shown.'), end
    if strcmp(Avg,'median') && ~strcmp(Spread,'iqr'), warning('unpairedDistribution:median','Your current choice of median and SEM or SD is not correct! Please use the inter-quartile range (IQR). Errorbars will not be shown.'), end
end
    
% Expand style to other data groups
if length(Slices)<nGroups && length(Slices)==1
    Slices = repmat(Slices,1,nGroups);
end
if length(Width)<nGroups && length(Width)==1
    Width = repmat(Width,1,nGroups);
end
if length(Compression)<nGroups && length(Compression)==1
    Compression = repmat(Compression,1,nGroups);
end
if length(Marker)<nGroups_to_plot && length(Marker)==1
    Marker = repmat(Marker,1,nGroups_to_plot);
end
if length(MarkerSize)<nGroups_to_plot && length(MarkerSize)==1
    MarkerSize = repmat(MarkerSize,1,nGroups_to_plot);
end
if length(MarkerFaceColor)<nGroups_to_plot && length(MarkerFaceColor)==1
    MarkerFaceColor = repmat(MarkerFaceColor,1,nGroups_to_plot);
end
if length(MarkerEdgeColor)<nGroups_to_plot && length(MarkerEdgeColor)==1
    MarkerEdgeColor = repmat(MarkerEdgeColor,1,nGroups_to_plot);
end
if ~iscell(AvgColor)
    if strcmpi(AvgColor, 'same')
        AvgColor = MarkerFaceColor;
    else
        AvgColor = {AvgColor};
    end
end
if length(AvgColor)<nGroups_to_plot && length(AvgColor)==1
    AvgColor = repmat(AvgColor,1,nGroups_to_plot);
end
X_distance = 1 + linspace(0, X_distance * (nGroups_to_plot - 1), nGroups_to_plot);

% Assign the offset value to the correct dimension
switch Axis
    case 'x'
        Offset = [Offset, 0];
    case 'y'
        Offset = [0, Offset];
end

%% Estimate spread along dimensions
HOR = NaN(numPoints, 1);
VER = NaN(numPoints, 1);
IDX = NaN(numPoints, 2);
for g = 1:nGroups
    original_group_idx = find(group_idx == g);
    this_group_idx = original_group_idx;
    y = data(this_group_idx);
    if ~CollapseGroups
        switch Shape
            case '^'
                [y, order] = sort(y, 'ascend');
            case 'v'
                [y, order] = sort(y, 'descend');
        end
        this_group_idx = this_group_idx(order);
    end
    invalid_data_points = ~isfinite(y);
    n_invalid_data_points = sum(invalid_data_points);
    data_y = y(~invalid_data_points);
    x = scatterXYpositions(data_y, Slices(g), Width(g), Compression(g));
    x = [zeros(n_invalid_data_points,1); x];
    x = x + X_distance(g);
    this_group_idx = [this_group_idx(invalid_data_points); this_group_idx(~invalid_data_points)];

    HOR(original_group_idx) = x;
    VER(original_group_idx) = y;
    IDX(original_group_idx, :) = [this_group_idx, ones(size(this_group_idx)) * g];
end

% Get axes in which to plot
Ax1 = gca();
% Get and store current hold status
prevHold = ishold();
% Turn hold on
hold on

%% Show mean and errors
% Swap dimensions if necessary
switch Axis
    case 'x', x=HOR; y=VER;
    case 'y', x=VER; y=HOR;
end

% Get group dispersion
add_error_bars = ~strcmp(ErrorBar,'none');

M = zeros(nGroups,1);
E = NaN(nGroups,2);
switch Avg
    case 'mean'
        for g = 1:nGroups
            M(g) = nanmean(VER(group_idx==g));
            if add_error_bars
                err_data = y(group_idx==g);
                err_data(~isfinite(err_data)) = [];
                if isempty(err_data) || length(err_data)==1, continue, end
                switch Spread
                    case 'sd',   E(g,:) = nanstd(err_data);
                    case 'sem',  E(g,:) =    sem(err_data);
                    case '95ci', [~,~,E(g,:)] = normfit(err_data);
                end
            end
        end
        if add_error_bars && ~strcmpi(Spread,'95ci')
            E = [M-E(:,1) M+E(:,2)];
        end
    case 'median'
        for g = 1:nGroups
            M(g) = median(y(group_idx==g));
            if add_error_bars
                E(g,:) = quantile(y(group_idx==g), [.25 .75]);
            end
        end
end

% Get group average and plot it 
hAvg = [];
switch ErrorBar
    case 'none'
        if ~strcmpi(Avg,'none')
            for g = 1:nGroups
                if ~isfinite(M(g)), continue, end
                switch Axis
                    case 'x'
                        hAvg(g,1) = plot([X_distance(g)-Width(g)/2*TickLength, X_distance(g)+Width(g)/2*TickLength]+Offset(1), [M(g) M(g)]+Offset(2), 'Color',AvgColor{g}, 'LineWidth',2, 'Visible',Visible);
                    case 'y'
                        hAvg(g,1) = plot([M(g) M(g)]+Offset(1), [X_distance(g)-Width(g)/2*TickLength, X_distance(g)+Width(g)/2*TickLength]+Offset(2), 'Color',AvgColor{g}, 'LineWidth',2, 'Visible',Visible);
                end
            end
        end
        
    case 'line'
        for g = 1:nGroups
            switch Axis
                case 'x'
                    if ~strcmpi(Avg,'none') && isfinite(M(g))
                        hAvg(g,1) = plot([X_distance(g)-Width(g)/2*TickLength, X_distance(g)+Width(g)/2*TickLength]+Offset(1), [M(g) M(g)]+Offset(2), 'Color',AvgColor{g}, 'LineWidth',2, 'Visible',Visible);
                    end
                    if ~isempty(Spread) && all(isfinite(E(g,:)))
                        hAvg(g,2) = plot([X_distance(g)-Width(g)/Compression(g), X_distance(g)+Width(g)/Compression(g)]+Offset(1),[E(g,2) E(g,2)]+Offset(2), 'Color',ErrorBarColor, 'LineWidth',2, 'Visible',Visible);
                        hAvg(g,3) = plot([X_distance(g)-Width(g)/Compression(g), X_distance(g)+Width(g)/Compression(g)]+Offset(1),[E(g,1) E(g,1)]+Offset(2), 'Color',ErrorBarColor, 'LineWidth',2, 'Visible',Visible);
                        hAvg(g,4) = plot([X_distance(g)                          X_distance(g)]+Offset(1),                        [E(g,1) E(g,2)]+Offset(2), 'Color',ErrorBarColor, 'LineWidth',2, 'Visible',Visible);
                    end
                case 'y'
                    if ~strcmpi(Avg,'none') && isfinite(M(g))
                        hAvg(g,1) = plot([M(g) M(g)]+Offset(1), [X_distance(g)-Width(g)/2*TickLength, X_distance(g)+Width(g)/2*TickLength]+Offset(2), 'Color',AvgColor{g}, 'LineWidth',2, 'Visible',Visible);
                    end
                    if ~isempty(Spread) && all(isfinite(E(g,:)))
                        hAvg(g,2) = plot([E(g,2) E(g,2)]+Offset(1),[X_distance(g)-Width(g)/Compression(g)*TickLength, X_distance(g)+Width(g)/Compression(g)*TickLength]+Offset(2), 'Color',ErrorBarColor, 'LineWidth',2, 'Visible',Visible);
                        hAvg(g,3) = plot([E(g,1) E(g,1)]+Offset(1),[X_distance(g)-Width(g)/Compression(g)*TickLength, X_distance(g)+Width(g)/Compression(g)*TickLength]+Offset(2), 'Color',ErrorBarColor, 'LineWidth',2, 'Visible',Visible);
                        hAvg(g,4) = plot([E(g,1) E(g,2)]+Offset(1),[X_distance(g)                                     X_distance(g)]+Offset(2),                                    'Color',ErrorBarColor, 'LineWidth',2, 'Visible',Visible);
                    end
            end
        end
        
    case 'rect'
        for g = 1:nGroups
            if ~strcmpi(Avg,'none') && isfinite(M(g))
                if ~isempty(Spread) && all(isfinite(E(g,:)))
                    hAvg(g,2) = rectangle('Position', [X_distance(g)-Width(g)/Compression(g)*TickLength+Offset(1), E(g,1)+Offset(2), 2*Width(g)/Compression(g)*TickLength, E(g,2)-E(g,1)], 'FaceColor',ErrorBarColor, 'EdgeColor','none', 'Visible',Visible);
                end
                hAvg(g,1) = plot([X_distance(g)-Width(g)/Compression(g)*TickLength, X_distance(g)+Width(g)/Compression(g)*TickLength]+Offset(1),[M(g) M(g)]+Offset(2),'Color',AvgColor{g},'LineWidth',2, 'Visible',Visible);
            end
        end
end

%% Plot individual observations
hData = cell(nGroups_to_plot,1);
X = [];
Y = [];
I = [];
for g = 1:nGroups_to_plot
    idx = groups_to_plot_idx == g;
    xPoints = x(idx) + Offset(1);
    yPoints = y(idx) + Offset(2);
    iPoints = IDX(idx, :);
    % Remove points that cannot be shown
    bad_points = ~isfinite(yPoints);
    yPoints(bad_points) = [];
    xPoints(bad_points) = [];
    iPoints(bad_points) = [];
    % Skip to next group if no data left
    if isempty(yPoints), continue, end
    % Draw points
    if ShowPoints
        h = scatter(xPoints, yPoints, MarkerSize(groups_to_plot(g)));
        set(h, 'Marker',Marker{groups_to_plot(g)}, ...
               'MarkerFaceColor',MarkerFaceColor{groups_to_plot(g)}, ...
               'MarkerEdgeColor',MarkerEdgeColor{groups_to_plot(g)}, ...
               'MarkerFaceAlpha',MarkerFaceAlpha, ...
               'MarkerEdgeAlpha',MarkerEdgeAlpha, ...
               'LineWidth',LineWidth, ...
               'Visible',Visible);
        hData{g} = h;
    end
    % Store points
    X = [X; xPoints];
    Y = [Y; yPoints];
    I = [I; iPoints];
end

%% Fix axis appearance
if strcmpi(Visible,'on')
    % Make line errorbars more visible
    if strcmp(ErrorBar,'line'), uistack(hAvg(:),'top'), end
    set(Ax1, 'FontSize',10, 'TickDir','out', 'Layer','bottom', 'Box','off')

    % Set axis limits
    switch Axis
        case 'x'
            xL = xlim();
            newxL = [.5 nGroups+.5]+Offset(1);
            xlim(Ax1,[min([xL(1) newxL(1)]) max([xL(2) newxL(2)])])
            if ~all(cellfun(@isempty,Labels)) && length(Labels)==nGroups
                xTL = cellstr(get(gca,'XtickLabel'));
                ind = find(ismember(get(gca,'Xtick'), (1:nGroups)+Offset(1)));
                for ii = 1:length(ind)
                    if ~isempty(Labels{ii})
                        xTL{ind(ii)} = Labels{ii};
                    end
                end
                set(gca,'XTickLabel',xTL)
            end

        case 'y'
            yL = ylim();
            newyL = [.5 nGroups+.5]+Offset(2);
            ylim(Ax1,[min([yL(1) newyL(1)]) max([yL(2) newyL(2)])])
            if ~all(cellfun(@isempty,Labels)) && length(Labels)==nGroups
                yTL = cellstr(get(gca,'YtickLabel'));
                ind = find(ismember(get(gca,'Ytick'), (1:nGroups)+Offset(1)));
                for ii = 1:length(ind)
                    if ~isempty(Labels{ii})
                        yTL{ind(ii)} = Labels{ii};
                    end
                end
                set(gca,'YTickLabel',yTL)
            end
    end

    if ~prevHold, hold off, end

else
    % Set graph visibility
    try delete(hData); end
    try delete(hAvg); end
end

%% Return outputs
if nargout<2, clear hAvg, end
if nargout<1, clear hData, end

if ReshapeReturnedPoints
    % Sort point coordinates to original order
    if strcmp(Shape, 'v')
        sort_sign = +1;
    else
        sort_sign = -1;
    end
    [I, order] = sortrows(I, 1 * sort_sign);
    % Apply order
    X = X(order);
    Y = Y(order);
    
    % Reshape to original shape
    X = reshape(X, original_size);
    Y = reshape(Y, original_size);
end


function pos = scatterXYpositions (yValues, RangeCut, Width, Compression)
    range_val = norminv([0.05 0.95], mean(yValues), std(yValues));
    cuts = abs(range_val(1)-range_val(2)) / RangeCut;
    if isnan(cuts), cuts = mean(yValues)/2; end

    pos = zeros(length(yValues),1);
    for steps = [-1, +1]
        cutUp = mean(yValues) + cuts / 2;
        subsetind = yValues <= cutUp & yValues > cutUp - cuts;
        n = sum(subsetind);
        keep_going = true;
        tries = 1;
        max_tries = 10;
        while keep_going && tries<=max_tries
            if all(n ~= [1, 0])
                distmaker = Width ./ (exp((n^2 - 1) ./ (n*Compression)));
                xSubset = zeros(n,1);
                oddie = mod(n,2);
                xb = linspace(1 - Width + distmaker, 1 + Width - distmaker, n) - 1;
                xSubset(1:2:end-oddie) = xSubset(1:2:end-oddie) + xb(1:round(end/2)-oddie)';
                xSubset(2:2:end)       = xSubset(2:2:end)       - xb(1:round(end/2)-oddie)';
                pos(subsetind)= xSubset;    
            end
            keep_going = ~(cutUp > max(yValues) || cutUp < min(yValues));
            cutUp = cutUp + steps * cuts;
            subsetind = yValues < cutUp & yValues > cutUp - cuts;
            n = sum(subsetind);
            tries = tries + 1;
        end
    end

%% MLint exceptions
%#ok<*AGROW,*NASGU>
