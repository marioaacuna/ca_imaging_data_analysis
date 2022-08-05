function varargout = load_variable(fname, varargin)

if ~exist(fname,'file')
    error('Unable to read file %s. No such file or directory.',fname)
end % Let MATLAB raise an error
% Get access to file as a matfile object
matObj = matfile(fname, 'Writable',false);
vars = fieldnames(matObj);
vars(1) = []; % The first one is "Properties"

if isempty(varargin)
    if length(vars)>2, allVars=natsort(vars); allVars=sprintf('%s, ',allVars{:});
        error('File "%s" contains multiple variables. Please indicate which one should be loaded.\nVariables are: %s',fname,allVars(1:end-2))
    end
    varargout{1} = matObj.(vars{1}); % If there is only one variable, load it
else
    for v = 1:length(varargin)
        varname = varargin{v}; % For each variable
        
        if ismember('.',varname) % User inputed a structure field
            varIsStructure = true;
        else
            varIsStructure = false;
        end
        if varIsStructure
            idx = find(ismember(varname,'.'),1,'first');
            fieldname = varname(idx+1:end);
            varname = varname(1:idx-1);
        else
            fieldname = '';
        end

        idx = ismember(vars, varname);
        if all(~idx)
            error(['Variable "', varname, '" is not present in the file "', fname, '".\nTry one among ', strjoin(vars, ', ')])
        end
        var = matObj.(vars{idx});

        % Find substructure
        if varIsStructure
            try
                eval(['var = var.' fieldname ';']);
            catch
                warning('It is not possible to access field "%s" in variable "%s" in the file "%s".\nThe output contains the full structure.',fieldname,varname,fname)
            end
        end

        if nargout == 0 % Copy variable to workspace
            assignin('caller',varname,var)
        else
            varargout{v} = var;
        end
    end
end

if ~exist('varargout','var'), varargout={}; end

%#ok<*ERTAG,*AGROW,*WNTAG,*CTCH>
