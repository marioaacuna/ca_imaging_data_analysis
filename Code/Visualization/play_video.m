function play_video(data, varargin)
%PLAYVIDEO Plays a data image and sets some defautls properties.

p = inputParser();
addParameter(p, 'scale', true)
addParameter(p, 'invert', false)
parse(p, varargin{:})
scale = p.Results.scale;
invert = p.Results.invert;

% Rescale data
if scale
    data = data -  nanmin(data(:));
    data = data ./ nanmax(data(:));
end 

% Open data to be played with a fast frame rate
h = implay(data, 100);

% Set colormap
cmap = gray(256);
if invert
    cmap = flipud(cmap);
end
h.Visual.ColorMap.Map = cmap;

% Enable fit to window
uitoggletool_Maintain = findall(0, 'tag','uimgr.uitoggletool_Maintain');
uitoggletool_Maintain(1).ClickedCallback();

% Make window larger 
h.Parent.Position = [100, 100, 700, 700];

