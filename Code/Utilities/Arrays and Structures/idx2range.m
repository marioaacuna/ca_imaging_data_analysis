function dataIDX = idx2range(idx)

if isempty(idx)
    dataIDX = zeros(0,2);
    return
end

if numel(idx)~=1
    idx = idx(:);
    dataIDX = [1; find(diff(idx)>1)+1];
    if size(dataIDX,1)~=1
        dataIDX(:,2) = [dataIDX(2:end)-1; length(idx)];
        dataIDX = idx(dataIDX);
        dataIDX(:,3) = (dataIDX(:,2)-dataIDX(:,1))+1;
    else
        dataIDX(:,2) = length(idx);
        dataIDX = idx(dataIDX)';
        dataIDX(:,3) = ((dataIDX(:,2)-dataIDX(:,1))+1)';
    end
else
    dataIDX = [idx idx 1];
end

