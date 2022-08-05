function [x, y] = correlogram_make_xy_outline(x, y)
    x = x(:);
    y = y(:);
    % Calculate the difference between x-values
    half_stepSize = diff(x) / 2;
    half_stepSize_before = [half_stepSize(1); half_stepSize];
    half_stepSize_after = [half_stepSize; half_stepSize(end)];
    % Calculate the values for x and y values of the outline
    x = reshape([x-half_stepSize_before x+half_stepSize_after]', [], 1);
    y = reshape(repmat(y, 1, 2)', [], 1);
