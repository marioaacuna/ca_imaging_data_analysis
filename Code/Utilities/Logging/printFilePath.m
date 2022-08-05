function str = printFilePath(filePath, string2show, printToVar)

if ~exist('string2show','var') || isempty(string2show), string2show = filePath; end
if ~exist('printToVar','var') || isempty(printToVar), printToVar = true; end

% Add the command to open the file
if ispc
    % Replace slashes
    command = ['winopen(''' filePath ''')'];
else
    command = ['unix(''open "' filePath '" &'')'];
end

% Make the string
str0 = sprintf('<a href="matlab:%s">%s</a>',command,string2show);

% Output string
if printToVar
    str = str0;
else
    fprintf('%s\n',str0)
end
