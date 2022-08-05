function n_channels = read_tiff_n_channels(filename)
    % Remove harmless warning
%     warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning')
    
    % Open file and read header
    TIFF_file = Tiff(filename, 'r');
    TIFF_header = regexp(TIFF_file.getTag('ImageDescription'), '\r', 'split')';
    
    % Get number of channels imaged
    n_channels = regexp(TIFF_header(cellfun(@(x) ~isempty(regexp(x, 'state.acq.numberOfChannelsSave=', 'once')), TIFF_header)), '=', 'split');
    n_channels = str2double(n_channels{1}{2});

    % Close handle
    TIFF_file.close()
    
%     if ~skip_check_file
%         % Load whole file
%         data = load_tiff_stack(filename);
%         real_n_frames = size(data, 3);
%         if n_frames * n_channels ~= real_n_frames
%             error('TIFFinfo:CorruptedFile', 'The file "%s" might be corrupted because the number of frames preesnt in the file (%i) is not equal to the number of frames in the metadata of the file (%i) * the number of imaged channels (%i).', filename, real_n_frames, n_frames, n_channels)
%         end
%     end
    
    % Re-enable warning
%     warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning')
