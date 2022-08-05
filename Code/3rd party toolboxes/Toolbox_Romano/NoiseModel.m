function [densityData, densityNoise, xev, yev] = NoiseModel(deltaFoF, sigma)

lambda_smoothing = 8;        
nevs = 1000;
sizeWinFilt = 8; % in sigmas


transfDataMatrix = deltaFoF ./ sigma;
points = [transfDataMatrix(1:end-1,:), transfDataMatrix(2:end,:)];
pointsNeg = points(all(points'<0), :);
Sigma = cov(pointsNeg(:,1), pointsNeg(:,2));
Mu = [0 0];
Sigma = [1, Sigma(1,2); Sigma(2,1), 1];
dataGaussianCov = mvnrnd(Mu,Sigma,4*length(pointsNeg));
clear pointsNeg

mn = floor(min([min(transfDataMatrix(:)), min(dataGaussianCov(:))]));
mx = floor(max([max(transfDataMatrix(:)), max(dataGaussianCov(:))]));
[xev, yev] = meshgrid(linspace(mn,mx,nevs),linspace(mn,mx,nevs));

binColIdx = xev(:);
binRowIdx = repmat(linspace(mn,mx,nevs), 1, nevs)';

[~, dist] = knnsearch(points, [binColIdx, binRowIdx], 'K',100);
binSpread = reshape(mean(dist,2), nevs, nevs);
[~, distNoise] = knnsearch(dataGaussianCov, [binColIdx, binRowIdx], 'K',100);
binSpreadNoise = reshape(mean(distNoise,2), nevs, nevs);
clear binColIdx binRowIdx dist distNoise

smoothParam      = 1/(sizeWinFilt * min(binSpread(:)));
smoothParamNoise = 1/(sizeWinFilt * min(binSpreadNoise(:)));


pointCol = interp1(xev(1,:), 1:nevs, points(:,1), 'nearest');
pointRow = interp1(xev(1,:), 1:nevs, points(:,2), 'nearest');
NaNs = isnan(pointCol) | isnan(pointRow);
pointCol(NaNs) = [];
pointRow(NaNs) = [];
hist2dim = accumarray([pointRow, pointCol], 1, [nevs, nevs]);

pointColNoise = interp1(xev(1,:), 1:nevs, dataGaussianCov(:,1), 'nearest');
pointRowNoise = interp1(xev(1,:), 1:nevs, dataGaussianCov(:,2), 'nearest');
NaNs = isnan(pointColNoise) | isnan(pointRowNoise);
pointColNoise(NaNs) = [];
pointRowNoise(NaNs) = [];
hist2dimNoise = accumarray([pointRowNoise, pointColNoise], 1, [nevs, nevs]);



densityData = zeros(size(xev));
n_rows = size(densityData, 1);
n_cols = size(densityData, 2);
[binRows, binCols] = find(hist2dim);
for irow = 1:length(binRows)
    sigmaFilt = smoothParam * binSpread(binRows(irow), binCols(irow));
    hsize = double([ceil(sizeWinFilt * sigmaFilt), ceil(sizeWinFilt*sigmaFilt)]);
    if mod(ceil(sizeWinFilt * sigmaFilt),2) == 0
        hsize = hsize + 1;
    end
    binFilter = fspecial('gaussian', hsize, double(sigmaFilt));
   
    widthRect = (size(binFilter,1)-1) / 2;
    centerRect = widthRect + 1;
         
    rectRows = binRows(irow) - min(binRows(irow)-1, widthRect):binRows(irow) + min(n_rows - binRows(irow), widthRect);
    rectCols = binCols(irow) - min(binCols(irow)-1, widthRect):binCols(irow) + min(n_cols - binCols(irow), widthRect);
    rectFiltRows = centerRect - min(binRows(irow)-1, widthRect):centerRect + min(n_rows - binRows(irow), widthRect); 
    rectFiltCols = centerRect - min(binCols(irow)-1, widthRect):centerRect + min(n_cols - binCols(irow), widthRect);
    
    densityData(rectRows, rectCols) = densityData(rectRows, rectCols) + hist2dim(binRows(irow), binCols(irow)) * binFilter(rectFiltRows, rectFiltCols);
end

densityNoise = zeros(size(xev));
n_rows = size(densityNoise, 1);
n_cols = size(densityNoise, 2);
[binRowsNoise, binColsNoise] = find(hist2dimNoise);
for irow = 1:length(binRowsNoise)
    sigmaFilt = smoothParamNoise * binSpreadNoise(binRowsNoise(irow), binColsNoise(irow));
    hsize = double([ceil(sizeWinFilt * sigmaFilt), ceil(sizeWinFilt * sigmaFilt)]);
    if mod(ceil(sizeWinFilt*sigmaFilt),2) == 0
        hsize = hsize + 1;
    end
    binFilter = fspecial('gaussian', hsize, double(sigmaFilt));
    
    widthRect = (size(binFilter,1)-1) / 2;
    centerRect = widthRect+1;
   
    rectRows = binRowsNoise(irow) - min(binRowsNoise(irow)-1, widthRect):binRowsNoise(irow) + min(n_rows - binRowsNoise(irow), widthRect);
    rectCols = binColsNoise(irow) - min(binColsNoise(irow)-1, widthRect):binColsNoise(irow) + min(n_cols - binColsNoise(irow), widthRect);
    rectFiltRows = centerRect - min(binRowsNoise(irow)-1, widthRect):centerRect + min(n_rows - binRowsNoise(irow), widthRect);
    rectFiltCols = centerRect - min(binColsNoise(irow)-1, widthRect):centerRect + min(n_cols - binColsNoise(irow), widthRect);
    
    densityNoise(rectRows, rectCols) = densityNoise(rectRows, rectCols) + hist2dimNoise(binRowsNoise(irow), binColsNoise(irow)) * binFilter(rectFiltRows, rectFiltCols);
end

densityData = Smooth1D(densityData,   lambda_smoothing);
densityData = Smooth1D(densityData.', lambda_smoothing)';
densityData = densityData ./ (sum(densityData(:)));
densityNoise = Smooth1D(densityNoise,   lambda_smoothing);
densityNoise = Smooth1D(densityNoise.', lambda_smoothing)';
densityNoise = densityNoise ./ (sum(densityNoise(:)));


function Z = Smooth1D(Y, lambda)
    m = size(Y, 1);
    E = eye(m);
    D1 = diff(E, 1);
    D2 = diff(D1, 1);
    P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
    Z = (E + P) \ Y;   
