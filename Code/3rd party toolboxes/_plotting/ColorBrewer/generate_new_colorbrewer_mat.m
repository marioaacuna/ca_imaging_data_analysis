clear all;

warning('off');
% [num,txt] = xlsread('ColorBrewer_all_schemes_RGBonly3.XLS','','','basic');
[num,txt] = xlsread('ColorBrewer_all_schemes_RGBonly4_withPalette_and_Macro.xls','','','basic');
warning('on');

% remove header
header = txt(1,:);
num(1,:) = [];
txt(1,:) = [];
txt(size(num,1)+1:end,:) = [];

RGB = num(:,ismember(header,{'R','G','B'}));

ColorName = txt(:,ismember(header,{'ColorName','Type'}));
idx = find(cellfun(@(x)~isempty(x), ColorName(:,1)));
schemes = [ColorName(idx,[2 1]) num2cell(num(idx,ismember(header,'NumOfColors')))];


for s = 1:size(schemes,1)
    ctype = schemes{s,1};
    cname = schemes{s,2};
    ncol  = schemes{s,3};
    colorbrewer.(ctype).(cname){ncol} = RGB(1:ncol,:);
    RGB(1:ncol,:) = [];
end

% Sort fields alphabetically
fld = unique(schemes(:,1));
for f = 1:length(fld)
    colorbrewer.(fld{f}) = orderfields(colorbrewer.(fld{f}));
end
colorbrewer = orderfields(colorbrewer);

% Fix qualitative maps so they have only one field
QUALschemes = fieldnames(colorbrewer.qual);
for q = 1:length(QUALschemes)
    map = colorbrewer.qual.(QUALschemes{q}){end};
    if strcmpi(QUALschemes{q},'Paired')  % switch light/dark pairs (so dark ones come first)
        idx = flipdim(reshape(1:size(map,1),2,[]),1);
        idx = idx(:,[1:3 5 4 6:size(idx,2)]); % swap orange and purple
        idx = idx(:);
        map = map(idx,:);
    end
    colorbrewer.qual.(QUALschemes{q}) = map;
end   

% Store all colormaps on disk
save('colorbrewer.mat','colorbrewer');
