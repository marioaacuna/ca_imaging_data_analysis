function offsetAxes(ax, do_add_listener)
% thanks to Pierre Morel, undocumented Matlab
% and https://stackoverflow.com/questions/38255048/separating-axes-from-plot-area-in-matlab
%
% by Anne Urai, 2016

if ~exist('ax', 'var'), ax = gca(); end
if ~exist('do_add_listener', 'var'), do_add_listener = true; end

% modify the x and y limits to below the data (by a small amount)
try
    ax.YLim(1) = ax.YLim(1)-(ax.YTick(2)-ax.YTick(1))/4;
    ax.XLim(1) = ax.XLim(1)-(ax.XTick(2)-ax.XTick(1))/4;
    ax.XLim(2) = ax.XLim(2)+(ax.XTick(2)-ax.XTick(1))/4;
catch ME
    if strcmp(ME.identifier, 'MATLAB:structRefFromNonStruct')  % old format axes
        YLims = get(ax, 'YLim');
        YTicks = get(ax, 'YTick');
        YLims(1) = YLims(1) - (YTicks(2)-YTicks(1)) / 4;
        YLims(2) = YLims(2) + (YTicks(2)-YTicks(1)) / 4;
        
        XLims = get(ax, 'XLim');
        XTicks = get(ax, 'XTick');
        XLims(1) = XLims(1) - (XTicks(2)-XTicks(1)) / 4;
        XLims(2) = XLims(2) + (XTicks(2)-XTicks(1)) / 4;

        set(ax, 'YLim',YLims, 'XLim',XLims)
        
    else
        rethrow(ME)
    end
end

% this will keep the changes constant even when resizing axes
if do_add_listener
    addlistener (ax, 'MarkedClean', @(obj,event)resetVertex(ax));
end

end

function resetVertex ( ax )
    XRuler = get(ax, 'XRuler');
    Axle = get(XRuler, 'Axle');
    VertexData = get(Axle, 'VertexData');
    VertexData(1,1) = min(get(ax, 'Xtick'));
    VertexData(1,2) = max(get(ax, 'Xtick'));
    set(Axle, 'VertexData',VertexData)

    YRuler = get(ax, 'YRuler');
    Axle = get(YRuler, 'Axle');
    VertexData = get(Axle, 'VertexData');
    VertexData(2,1) = min(get(ax, 'Ytick'));
    VertexData(2,2) = max(get(ax, 'Ytick'));
    set(Axle, 'VertexData',VertexData)

%     try
%         % extract the x axis vertext data
%         % X, Y and Z row of the start and end of the individual axle.
% %         ax.XRuler.Axle.VertexData(1,1) = min(get(ax, 'Xtick'));
% %         ax.XRuler.Axle.VertexData(1,2) = max(get(ax, 'Xtick'));
% %         % repeat for Y (set 2nd row)
% %         ax.YRuler.Axle.VertexData(2,1) = min(get(ax, 'Ytick'));
%     catch ME
%         if strcmp(ME.identifier, 'MATLAB:structRefFromNonStruct')  % old format axes
%             XRuler = get(ax, 'XRuler');
%             Axle = get(XRuler, 'Axle');
%             VertexData = get(Axle, 'VertexData');
%             VertexData(1,1) = min(get(ax, 'Xtick'));
%             VertexData(1,2) = max(get(ax, 'Xtick'));
%             set(Axle, 'VertexData', VertexData)
%             
%             YRuler = get(ax, 'YRuler');
%             Axle = get(YRuler, 'Axle');
%             VertexData = get(Axle, 'VertexData');
%             VertexData(2,1) = min(get(ax, 'Ytick'));
%             set(Axle, 'VertexData', VertexData)
%             
%         else
%             rethrow(ME)
%         end
%     end
end
