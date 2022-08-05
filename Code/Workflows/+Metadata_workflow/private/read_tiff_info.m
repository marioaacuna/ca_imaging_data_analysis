function [frame_rate, frame_width, frame_height, n_frames, bidirectional_scanning, fast_scanning, tiff_date, n_channels] = read_tiff_info(filename, skip_check_file)
    if ~exist('skip_check_file','var')
        skip_check_file = false;
    end

    % Remove harmless warning
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning')
    
    % Open file and read header
    TIFF_file = Tiff(filename, 'r');
    TIFF_header = regexp(TIFF_file.getTag('ImageDescription'), '\r', 'split')';
    
    % Get frame rate
    try
        frame_rate = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.frameRate=', 'once')), TIFF_header)), '=', 'split');
        frame_rate = str2double(frame_rate{1}{2});
        
    catch ME
        if strcmp(ME.identifier, 'MATLAB:badsubscript')
            TIFF_file.close()
            [frame_rate, frame_width, frame_height, n_frames, bidirectional_scanning, fast_scanning, tiff_date, n_channels] = read_scanimage_file(filename, skip_check_file);
            return
        else
            rethrow(ME)
        end
    end
        
    % Get frame size
    frame_width = TIFF_file.getTag('ImageWidth');
    frame_height = TIFF_file.getTag('ImageLength');

    % Get number of frames
    n_frames = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.numberOfFrames=', 'once')), TIFF_header)), '=', 'split');
    n_frames = str2double(n_frames{1}{2});
    
    % Get whether image has been acquired by bidirectional scanning
    bidirectional_scanning = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.bidirectionalScan=', 'once')), TIFF_header)), '=', 'split');
    bidirectional_scanning = logical(str2double(bidirectional_scanning{1}{2}));

    % Get whether image has been acquired in fast scanning mode
    fast_scanning_X = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.fastScanningX=', 'once')), TIFF_header)), '=', 'split');
    fast_scanning_X = logical(str2double(fast_scanning_X{1}{2}));
    fast_scanning_Y = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.fastScanningY=', 'once')), TIFF_header)), '=', 'split');
    fast_scanning_Y = logical(str2double(fast_scanning_Y{1}{2}));
    fast_scanning = fast_scanning_X | fast_scanning_Y;
    
    % Get exact date of experiment
    tiff_date = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.internal.triggerTimeString=', 'once')), TIFF_header)), '=', 'split');
    tiff_date = datenum(tiff_date{1}{2}(2:end-1), 'mm/dd/yyyy HH:MM:SS.FFF');

    % Get number of channels imaged
    n_channels = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.numberOfChannelsSave=', 'once')), TIFF_header)), '=', 'split');
    n_channels = str2double(n_channels{1}{2});

    % Close handle
    TIFF_file.close()
    
    if ~skip_check_file
        % Load whole file
        data = load_tiff_stack(filename);
        real_n_frames = size(data, 3);
        if n_frames * n_channels ~= real_n_frames
            error('TIFFinfo:CorruptedFile', 'The file "%s" might be corrupted because the number of frames preesnt in the file (%i) is not equal to the number of frames in the metadata of the file (%i) * the number of imaged channels (%i).', filename, real_n_frames, n_frames, n_channels)
        end
    end
    
    % Re-enable warning
    warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning')
