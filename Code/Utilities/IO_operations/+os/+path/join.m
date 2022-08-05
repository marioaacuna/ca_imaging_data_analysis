function joined_path = join(varargin)

% Keep only strings
if ~iscell(varargin), varargin = {varargin}; end
varargin(~cellfun(@ischar, varargin)) = [];

% Add path separator at the end of each element, except the last one
path_separator = filesep();
for i_path = 1:length(varargin) - 1
    if ~endsWith(varargin{i_path}, path_separator)
        varargin{i_path} = [varargin{i_path}, path_separator];
    end
end

% Join path
joined_path = strjoin(varargin, '');
