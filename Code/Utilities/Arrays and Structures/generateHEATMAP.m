function [heatmap, heatmap_HOR_values, heatmap_VER_values] = generateHEATMAP(data, varargin)

% Parse inputs
p = inputParser();
addParameter(p, 'range', [min(data, [], 1)', max(data, [], 1)']);
addParameter(p, 'bin_size', []);
addParameter(p, 'n_bins', [100, 100]);
addParameter(p, 'weights', ones(size(data,1), 1));
addParameter(p, 'occluder', []);
parse(p, varargin{:})
range = p.Results.range;
bin_size = p.Results.bin_size;
n_bins = p.Results.n_bins;
weights = p.Results.weights;
occluder = p.Results.occluder;
% Fix some inputs that depend on others
if length(n_bins) == 1, n_bins = [n_bins, n_bins]; end
if isempty(bin_size)
    bin_size = [(range(1, 2) - range(1, 1)) / n_bins(1), (range(2, 2) - range(2, 1)) / n_bins(2)];
end
if length(bin_size) == 1, bin_size = [bin_size, bin_size]; end

% Make grid --------------------------------------------------------------------
center_heatmap = mean(range,2);
range = range - repmat(center_heatmap,1,2);
heatmap_HOR_values = bin_size(1):bin_size(1):bin_size(1)*ceil((abs(diff(range(1,:))) /2) /bin_size(1));
heatmap_HOR_values = [-fliplr(heatmap_HOR_values) 0 heatmap_HOR_values];
heatmap_VER_values = bin_size(2):bin_size(2):bin_size(2)*ceil((abs(diff(range(2,:))) /2) /bin_size(2));
heatmap_VER_values = [-fliplr(heatmap_VER_values) 0 heatmap_VER_values];
heatmap_HOR_values = heatmap_HOR_values + center_heatmap(1);
heatmap_VER_values = heatmap_VER_values + center_heatmap(2);

% Calculate frequencies in each square -----------------------------------------
[~,idxX] = histc(data(:,1), heatmap_HOR_values);
[~,idxY] = histc(data(:,2), heatmap_VER_values);
OutOfRange = union(find(idxX==0|idxX==length(heatmap_HOR_values)),find(idxY==0|idxY==length(heatmap_VER_values)));
idxX(OutOfRange) = [];
idxY(OutOfRange) = [];
if length(weights)~=1; weights(OutOfRange) = []; end

% Sum up frequencies -----------------------------------------------------------
heatmap = accumarray([idxY idxX],weights,[length(heatmap_VER_values)-1 length(heatmap_HOR_values)-1],@sum,0);

% Remove occluded regions ------------------------------------------------------
if ~isempty(occluder)
    [heatmapGRIDx,heatmapGRIDy] = meshgrid(heatmap_HOR_values(1:end-1),heatmap_VER_values(1:end-1));
    for occ = 1:size(occluder,1)
        OCCvert = [occluder(occ,1) occluder(occ,2) occluder(occ,2) occluder(occ,1);...
                   occluder(occ,3) occluder(occ,3) occluder(occ,4) occluder(occ,4)];
               OCCvert = [OCCvert OCCvert(:,1)];
        heatmap(InPolygon(heatmapGRIDx,heatmapGRIDy,OCCvert(1,:)',OCCvert(2,:)')) = 0;
    end
end
