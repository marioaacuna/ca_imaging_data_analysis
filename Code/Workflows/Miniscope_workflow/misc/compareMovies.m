function [ output_args ] = compareMovies(start,stop,FPS,varargin)

nMov = length(varargin);
x = ceil(nMov/3);
y = ceil(nMov/x);
f = figure('WindowStyle','normal','Position',[100,200,y*500,x*400]);

cvals = [];
for i = 1:nMov
    cvals(i,1) = min(varargin{i}(:));
    cvals(i,2) = max(varargin{i}(:));
end

for i = start:stop
    if ~ishandle(f)
        break
    end
    for j = 1:nMov
        subplot(x,y,j)
        imagesc(varargin{j}(:,:,i));
        caxis(cvals(j,:));
        colormap(gray)
    end
    pause(1/FPS-0.02);
end

end

