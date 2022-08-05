function str = make_variable_name(str, hyphen_to_minus, convert_back)
% MAKE_VARIABLE_NAME Convert a string to a valid variable name.

% Whether hyphen has to be translated to a minus or a space
if ~exist('hyphen_to_minus','var') || isempty(hyphen_to_minus)
    hyphen_to_minus = true;
end

if ~exist('convert_back','var') || isempty(convert_back)
    convert_back = false;
end

if ~convert_back
    % Replace "-" with "minus"
    if hyphen_to_minus
        str = strrep(str, '-', ' minus');
    else
        str = strrep(str, '-', ' ');
    end
    % Replace "+" with "plus"
    str = strrep(str, '+', ' plus');

    % Replace "%" with "percent"
    str = strrep(str, '%', ' percent');

    % Replace " " with "_"
    str = strrep(str, ' ', '_');

else
    if hyphen_to_minus
        str = strrep(str, 'minus', '-');
    end
    % Replace "+" with "plus"
    str = strrep(str, 'plus', '+');

    % Replace "-" with "minus"
    remove_underscore = contains(str, '-') | contains(str, '+');

    % Replace "%" with "percent"
    str = strrep(str, '_percent', '%');
    
    % Replace "_" with ""
    if remove_underscore
        str = strrep(str, '_', '');
    end
end
