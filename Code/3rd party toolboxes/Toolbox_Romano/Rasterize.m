function [raster, mapOfOddsJoint] = Rasterize(deltaFoF, sigma, mapOfOdds, xev, yev, fps, tauDecay)

transfDataMatrix = deltaFoF ./ sigma;

cc = bwconncomp(~mapOfOdds);
stats = regionprops(cc, 'PixelList');
indZero = find(xev(1,:)>0, 1 ,'first');

for iobj = 1:cc.NumObjects
    if ismember(indZero, stats(iobj).PixelList(:,1))
        break
    end
end
mapOfOddsCorrected = ones(size(mapOfOdds));
mapOfOddsCorrected(cc.PixelIdxList{iobj}) = 0;

%%
noiseBias = 1.5;
factorDecay = exp(-1 / (fps * tauDecay));
decayMap = ones(size(mapOfOdds));
rowsDecay = 1:size(decayMap,1);
for icol=1:size(decayMap,2)
    decayMap(rowsDecay(yev(:,1) < factorDecay * (xev(1,icol)-noiseBias)-noiseBias), icol) = 0;
end

riseMap = ones(size(mapOfOdds));
riseMap(end,:)=0; 

mapOfOddsJoint = mapOfOddsCorrected & riseMap & decayMap;

raster = zeros(size(transfDataMatrix,1), size(transfDataMatrix,2));
for numNeuron = 1:size(transfDataMatrix, 2)
    [~, bins] = histc(transfDataMatrix(:, numNeuron), xev(1,:));
    for numFrame = 3:size(transfDataMatrix,1)-2
        try
            optA = (mapOfOddsJoint(bins(numFrame+1), bins(numFrame))   & mapOfOddsJoint(bins(numFrame),   bins(numFrame-1)));
            optB = (mapOfOddsJoint(bins(numFrame),   bins(numFrame-1)) & mapOfOddsJoint(bins(numFrame-1), bins(numFrame-2)));
            optC = (mapOfOddsJoint(bins(numFrame+2), bins(numFrame+1)) & mapOfOddsJoint(bins(numFrame+1), bins(numFrame)));
            if optA || optB || optC
                raster(numFrame,numNeuron) = 1;
            end
        end
    end
end

%#ok<*TRYNC>
