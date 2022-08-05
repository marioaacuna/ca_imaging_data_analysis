function drawGrid(hAx, x, y, linewidth, linestyle, color, move_to)

if ~exist('linewidth', 'var') || isempty(linewidth), linewidth = 1; end
if ~exist('linestyle', 'var') || isempty(linestyle), linestyle = '--'; end
if ~exist('color', 'var') || isempty(color), color = [.7 .7 .7]; end
if ~exist('move_to', 'var') || isempty(move_to), move_to = 'bottom'; end


prevAx = get(gcf,'CurrentAxes');
set(gcf,'CurrentAxes',hAx)
prevHold = ishold; hold on

if exist('x', 'var') && ~isempty(x)
    % Vertical lines
    gr = hggroup();
    arrayfun(@(i) plot(gr, [i, i], ylim(hAx), 'linestyle',linestyle, 'LineWidth',linewidth, 'color',color, 'XLimInclude','off', 'YLimInclude','off'), x);
    uistack(gr,move_to)
end

if exist('y', 'var') && ~isempty(y)
    % Horizontal lines
    gr = hggroup();
    arrayfun(@(i) plot(gr, xlim(hAx), [i i], 'linestyle',linestyle, 'LineWidth',linewidth, 'color',color, 'XLimInclude','off', 'YLimInclude','off'), y);
    uistack(gr,move_to)
end

if ~prevHold, hold off, end
set(gcf,'CurrentAxes',prevAx)
