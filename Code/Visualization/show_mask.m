function show_mask(mask, invert, show_id)
%SHOWFRAME Shows either mean or max frame of data.

if ~exist('invert','var')
    invert = false;
end
if ~exist('show_id','var')
    show_id = true;
end

% Make colors
cmap = gray(256);
if invert
    cmap = flipud(cmap);
    text_color = [.9, .9, .9];
else
    text_color = [.1, .1, .1];
end
    
% Convert mask to labeled image
if islogical(mask)
    bwmask = mask;
    mask = bwlabel(mask);
else
    bwmask = mask > 0;
end

% If IDs have to be shown, calculate where to put them
if show_id
    % Get position of centroids
    CENTROIDS = regionprops(mask, 'Centroid');
    CENTROIDS = reshape([CENTROIDS.Centroid], 2, []).';
end

% Make figure and show data
figure('color','w');
clf()
if show_id
    imagesc(bwmask);
else
    imagesc(mask);
end
colormap(cmap)
axis image off

% Show ROI id
if show_id
    n_ROIs = size(CENTROIDS,1);
    for iroi = 1:n_ROIs
        text(CENTROIDS(iroi,1), CENTROIDS(iroi,2), num2str(iroi), 'FontSize',12, 'Hor','center', 'Ver','middle', 'color',text_color, 'Interpreter','latex')
    end
end
