function pos = getAxisPos(h)
%PLOTBOXPOS Returns the position of the plotted axis region
%
% pos = plotboxpos(h)
%
% This function returns the position of the plotted region of an axis,
% which may differ from the actual axis position, depending on the axis
% limits, data aspect ratio, and plot box aspect ratio.  The position is
% returned in the same units as the those used to define the axis itself.
% This function can only be used for a 2D plot.  
%
% Input variables:
%
%   h:      axis handle of a 2D axis (if ommitted, current axis is used).
%
% Output variables:
%
%   pos:    four-element position vector, in same units as h

% Check input
h = h(ishandle(h(:)));
h = h(strcmp(get(h,'type'),'axes'));
if isempty(h), error('Input must be an axis handle'), end

% Get position of axis in pixels
pos = NaN(length(h),4);
for ax = 1:length(h)
    currunit = get(h(ax), 'units');
    set(h(ax), 'units', 'pixels')
    axisPos = get(h(ax), 'Position');
    set(h(ax), 'Units', currunit)

    % Calculate box position based axis limits and aspect ratios
    darismanual  = strcmpi(get(h(ax),'DataAspectRatioMode'),   'manual');
    pbarismanual = strcmpi(get(h(ax),'PlotBoxAspectRatioMode'),'manual');

    if ~darismanual && ~pbarismanual
        pos(ax,:) = axisPos;
    else
        dx = diff(get(h(ax), 'XLim'));
        dy = diff(get(h(ax), 'YLim'));
        dar  = get(h(ax), 'DataAspectRatio');
        pbar = get(h(ax), 'PlotBoxAspectRatio');

        limDarRatio = (dx/dar(1))/(dy/dar(2));
        pbarRatio = pbar(1)/pbar(2);
        axisRatio = axisPos(3)/axisPos(4);

        if darismanual
            if limDarRatio > axisRatio
                pos(ax,1) = axisPos(1);
                pos(ax,3) = axisPos(3);
                pos(ax,4) = axisPos(3)/limDarRatio;
                pos(ax,2) = (axisPos(4) - pos(ax,4))/2 + axisPos(2);
            else
                pos(ax,2) = axisPos(2);
                pos(ax,4) = axisPos(4);
                pos(ax,3) = axisPos(4) * limDarRatio;
                pos(ax,1) = (axisPos(3) - pos(ax,3))/2 + axisPos(1);
            end
        elseif pbarismanual
            if pbarRatio > axisRatio
                pos(ax,1) = axisPos(1);
                pos(ax,3) = axisPos(3);
                pos(ax,4) = axisPos(3)/pbarRatio;
                pos(ax,2) = (axisPos(4) - pos(ax,4))/2 + axisPos(2);
            else
                pos(ax,2) = axisPos(2);
                pos(ax,4) = axisPos(4);
                pos(ax,3) = axisPos(4) * pbarRatio;
                pos(ax,1) = (axisPos(3) - pos(ax,3))/2 + axisPos(1);
            end
        end
    end

    % Convert plot box position to the units used by the axis
    tmpAx = axes('Units','Pixels', 'Position',pos(ax,:), 'Visible','off', 'parent',get(h(ax),'parent'));
        set(tmpAx, 'Units',currunit)
        pos(ax,:) = get(tmpAx,'position');
        delete(tmpAx)
end
