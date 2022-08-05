function [frame_rate, frame_width, frame_height, n_frames, bidirectional_scanning, fast_scanning, tiff_date, n_channels] = read_scanimage_file(filename, ~)

TIFF_file = ScanImageTiffReader(filename);
TIFF_header = strsplit(TIFF_file.metadata, '\n')';

% Get frame rate
frame_rate = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'SI.hRoiManager.scanFrameRate = ', 'once')), TIFF_header)), '=', 'split');
frame_rate = str2double(frame_rate{1}{2});

% Get frame size
frame_width = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'SI.hRoiManager.pixelsPerLine = ', 'once')), TIFF_header)), '=', 'split');
frame_width = str2double(frame_width{1}{2});
frame_height = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'SI.hRoiManager.linesPerFrame = ', 'once')), TIFF_header)), '=', 'split');
frame_height = str2double(frame_height{1}{2});

% Get whether image has been acquired by bidirectional scanning
bidirectional_scanning = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'SI.hScan2D.bidirectional = ', 'once')), TIFF_header)), '=', 'split');
switch strip(bidirectional_scanning{1}{2})
    case 'true', bidirectional_scanning = true;
    case 'false', bidirectional_scanning = false;
end

% Get whether image has been acquired in fast scanning mode
fast_scanning = false;

% Get exact date of experiment
descriptions = strsplit(TIFF_file.descriptions{1}, '\n')';
epoch = regexp(descriptions(cellfun(@(x) ~isempty(regexp(x, 'epoch = ', 'once')), descriptions)), '=', 'split');
epoch = strip(epoch{1}{2});
try
    tiff_date = datenum(epoch, '[yyyy,mm,dd,HH,MM,SS.FFF]');
catch
     tiff_date = datenum(epoch, '[yyyy,mm,dd,HH,MM,SS]');
end
% Get number of channels imaged
n_channels = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'SI.hChannels.channelsActive = ', 'once')), TIFF_header)), '=', 'split');
n_channels = str2double(n_channels{1}{2});

% Load whole file
data = TIFF_file.data;
n_frames = size(data, 3);

TIFF_file.close();

