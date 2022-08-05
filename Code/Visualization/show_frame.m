function show_frame(data, average_function, invert, show_cbar)
%SHOWFRAME Shows either mean or max frame of data.

if ~exist('average_function','var') || isempty(average_function)
    average_function = 'mean';
end
if ~exist('invert','var')
    invert = false;
end
if ~exist('show_cbar','var')
    show_cbar = false;
end

% Make data to show
switch average_function
    case 'mean'
        data_to_show = nanmean(data, 3);
    case 'max'
        data_to_show = nanmax(data, [], 3);
    case 'std'
        data_to_show = nanstd(data, 0, 3);
    case 'median'
        data_to_show = nanmedian(data, 3);
end

% Make figure and show data
figure('color','w');
imagesc(data_to_show);

% Adjust axes
cmap = gray(256);
if invert
    cmap = flipud(cmap);
end
colormap(cmap)
axis square off

% Show colorbar
if show_cbar
    colorbar();
end
