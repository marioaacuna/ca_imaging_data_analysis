function METADATA = read_preprocessing_metadata(metadata_file, ignore_duplicates)

    if ~exist(metadata_file, 'file')
        error('csv_metadata:MetadataCsvDoesNotExist', [metadata_file, ' does not exist'])
    end
    
    if ~exist('ignore_duplicates','var')
        ignore_duplicates = false;
    end
    
    % Open file to read
    fileID = fopen(metadata_file);

    % Read lines
    lines = textscan(fileID, '%s', 'delimiter','\n');
    % Close file
    fclose(fileID);
    
    % The first element contains the text
    lines = lines{1};
    if isempty(lines)
        error('csv_metadata:MetadataCsvIsEmpty', 'Metadata file ''%s'' is empty',metadata_file)
    end
    
    % Remove blanks
    lines = strrep(lines, ' ','');
        
    % The first line contains the header
    header = regexp(lines{1}, ',','split');
    % Allocate output variable
    METADATA = repmat({''}, [length(lines)-1, length(header)]);
    % Fill in table
    for irow = 2:length(lines)
        METADATA(irow-1, :) = regexp(lines{irow}, ',', 'split');
    end
    
    if ~ignore_duplicates
        % Check for duplicates
        if size(METADATA,1) ~= length(unique(METADATA(:,1)))
            error('csv_metadata:MetadataCsvContainsDuplicates', 'Metadata file ''%s'' contains duplicates', metadata_file)
        end
    end
    
    % Convert to a MATLAB table
    METADATA = cell2table(METADATA, 'VariableNames',header);
