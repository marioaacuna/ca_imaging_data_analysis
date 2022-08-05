function varargout = run_in_R(R_script_name, varargin)

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
    
    % Add double quotes around paths
    path_args_idx = cellfun(@(x) contains(x, '\'), varargin);
    if any(path_args_idx)
        path_args_idx = find(path_args_idx);
        for ii = 1:length(path_args_idx)
            path_arg = varargin{path_args_idx(ii)};
            path_arg = regexp(path_arg, '\s', 'split');
            k = find(cellfun(@(x) contains(x, '\'), path_arg));
            % Add quotes and escape slashes
            for jj = 1:length(k)
                path_arg{k(jj)} =  ['"', strrep(path_arg{k(jj)}, '\','\\'), '"'];
            end
            varargin{ii} = strjoin(path_arg, ' ');
        end
    end
    % Join all arguments in one string
    args = strjoin(varargin, ' ');
    
else
    % Make an empty string
    args = '';
end

% Make sure that filename ends with extension '.R'
if ~endsWith(R_script_name, '.R')
    R_script_name = [R_script_name, '.R'];
end
% Get full name of R script
R_script_path = os.path.join(GC.R.scripts_path, R_script_name);
% Check that script exists
if ~exist(R_script_path, 'file')
    msg = 'R script does not exist';
    LOGGER.critical([R_script_path, ':', msg], 'contains_path',true)
    error('MATLAB:R', msg)
end

% Make command string
R_exe = find_path_of_R();
command_string = ['"', R_exe, '" --no-save "', R_script_path, '" ', args];
% Check whether to print output or not
if ~contains(args, {'-v', '--verbose'})
    R_log = [GC.temp_dir, 'R.log'];
    command_string = [command_string, ' > "', R_log, '" 2>&1'];
    has_log = true;
else
    has_log = false;
end
% Run command and print shell's stdout in MATLAB's Command Window
[status, msg] = system(command_string, '-echo');
flag_error = false;
% Catch errors
if status ~= 0
    if isempty(LOGGER)
        disp(msg)
        flag_error = 0;
    else
        %         LOGGER.critical([R_script_name, ':', msg], 'contains_path',true)
        flag_error = 1;
    end
    if has_log,
        %         edit(R_log),
    end
    %     error('MATLAB:R', ['ERROR encountered with following command:\n\n', strrep(command_string, '\','\\'), '\n\n'])
end
% Delete log, if it exists
if has_log
    delete(R_log)
end

varargout{2} = flag_error;
if nargout > 0
    varargout{1} = command_string;
    varargout{2} = flag_error;
end
