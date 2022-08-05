function varargout = run_in_python(python_script_name, varargin)

global GC LOGGER

% Parse additional arguments
if ~isempty(varargin)
    % Convert all arguments to strings
    non_char_cells = cellfun(@(x) ~ischar(x), varargin);
    if any(non_char_cells)
        idx = find(non_char_cells);
        for ii = 1:length(idx)
            value = varargin{idx(ii)};
            if isnumeric(value)
                value = num2str(value);
            elseif islogical(value)
                if value
                    value = 'True';
                else
                    value = 'False';
                end
            end
            varargin{idx(ii)} = value;
        end
    end
    
    % Add double quotes around each parameter
    args = cell(1, length(varargin));
    for ii = 1:length(args)
        args{ii} = ['"', varargin{ii}, '"'];
    end
    % Join all arguments in one string separated by a single space
    args = strjoin(args, ' ');
    % Escape slashes in paths
    args = strrep(args, '\','\\');
    
else
    % Make an empty string
    args = '';
end

% Make sure that filename ends with extension '.py'
if ~endsWith(python_script_name, '.py')
    python_script_name = [python_script_name, '.py'];
end
% Get full name of python script
python_script_path = os.path.join(GC.python.scripts_path, python_script_name);
% Check that script exists
if ~exist(python_script_path, 'file')
    msg = 'Python script does not exist';
    LOGGER.critical([python_script_path, ':', msg], 'contains_path',true)
    error([GC.python.environment_name, ':python'], msg)
end

% Make command string
if ispc
    command_string = ['conda activate ', GC.python.environment_name, ' && "' GC.python.interpreter_path, '" "', python_script_path, '" ', args];
else
    command_string = ['source ', GC.python.enviroment_path, ' activate' ,' && "' GC.python.interpreter_path, '" "', python_script_path, '" ', args];
end
% Run command and print shell's stdout in MATLAB's Command Window
[status, msg] = system(command_string, '-echo');
% Catch errors
if status ~= 0
    if ~isempty(LOGGER)
        LOGGER.critical([python_script_name, ':', msg])
    end
    error([GC.python.environment_name, ':python'], strrep(msg, '\', '\\'))
end

if nargout > 0
    varargout{1} = command_string;
end
