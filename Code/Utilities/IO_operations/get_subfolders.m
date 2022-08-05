% This function looks returns subfolder names, and it can also match
% against a pattern.

function SUBFOLDERS = get_subfolders(path, pattern, full_path)
    % Provide default values
    if ~exist('full_path', 'var'), full_path=true; end
    if ~exist('pattern', 'var'), pattern=''; end
    
    % Check inputs
    p = inputParser;
    addRequired(p, 'path', @ischar);
    addRequired(p, 'full_path', @islogical);
    addRequired(p, 'pattern', @ischar);
    parse(p, path,full_path,pattern);
    
    % Get list of files in this folder
    file_list = dir(p.Results.path);
    % Remove '.' and '..'
    file_list = file_list(~ismember({file_list.name}, {'.','..','settings'}));
    % Keep only folders
    SUBFOLDERS = file_list(cell2mat({file_list.isdir}));
    % Get names
    SUBFOLDERS = {SUBFOLDERS.name}';

    % If user has to check a pattern
    if ~strcmp(p.Results.pattern, '')
        SUBFOLDERS = SUBFOLDERS(cell2mat(cellfun(@(x) ~isempty(strfind(x,p.Results.pattern)), SUBFOLDERS, 'UniformOutput',false)));
    end
    
    % If user wants the full path or just the name
    if p.Results.full_path
        SUBFOLDERS = strcat(p.Results.path, filesep, SUBFOLDERS);
    end
