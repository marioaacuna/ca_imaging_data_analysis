function [hubert_Gamma, hubert_GammaPC, Davies_Bouldin_value, Davies_Bouldin_value_PC] = compute_assembly_quality(cell_assemblies, data_zscore, lambda_max)

    % Compute distances of z-scored data
    dists = pdist(data_zscore', 'euclidean');
    distMat = squareform(dists);
    
    % Compute distances of PCs
    [PCs, ~, eigenvals] = pca(data_zscore, 'Algorithm','svd', 'Economy',false);
    QPCs = PCs(:,eigenvals >= lambda_max);
    distsPCs = pdist(QPCs, 'euclidean');
    distMatPC = squareform(distsPCs);

    hubert_Gamma   = hubert(distMat,   cell_assemblies);
    hubert_GammaPC = hubert(distMatPC, cell_assemblies);
    Davies_Bouldin_value    = dbIndex(data_zscore, cell_assemblies);
    Davies_Bouldin_value_PC = dbIndex(QPCs',       cell_assemblies);

    
function hubertGamma = hubert(distMat, cell_assemblies)
    Q = ones(size(distMat));
    for n = 1:length(cell_assemblies)
        for i = 1:length(cell_assemblies{n})
            cellsIn = cell_assemblies{n};
            for j = 1:length(cellsIn)
                Q(cellsIn(j),cellsIn) = 0;
            end
        end
    end
    s = size(distMat);
    idx = triu(true(s), 1);
    
    muP = mean(distMat(idx));
    muQ = mean(Q(idx));
    stdP = std(distMat(idx));
    stdQ = std(Q(idx));
    multiplied = ((distMat - muP) .* (Q - muQ));

    sumMultiplied = 0;
    for i = 1:size(distMat,1)-1
        sumMultiplied = sumMultiplied + sum(multiplied(i, i+1:size(distMat,1)));
    end

    pairs = s(1) * (s(1) - 1) / 2;
    hubertGamma = sumMultiplied / (stdP * stdQ * pairs);

    
function Davies_Bouldin_value = dbIndex(Q, cell_assemblies)
    sizeTemp = length([cell_assemblies{:}]);
    QTemp = nan(sizeTemp,size(Q,1));
    clust = nan(sizeTemp,1);

    cellsEns = unique(reshape([cell_assemblies{:}],[],1));
    count = 1;
    for n = 1:length(cell_assemblies)
        for j = 1:length(cell_assemblies{n})
            targetCell = cell_assemblies{n}(j);
            QTemp(count,:) = Q(:,targetCell)';
            clust(count) = n;
            count = count + 1;
        end
    end
    % Include cells not in assemblies
    cellsNotInEns = setdiff(1:size(Q, 2), cellsEns);
    QTemp = [QTemp; Q(:,cellsNotInEns)'];
    clust = [clust; (n+1) * ones(length(cellsNotInEns),1)];
    [Davies_Bouldin_value, ~] = db_index(QTemp, clust);

