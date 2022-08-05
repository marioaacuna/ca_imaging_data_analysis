function str = value2str (value, format, force_cell_output)

if ~exist('format','var') || isempty(format)
    format = '%i';
end

if ~exist('force_cell_output','var')
    force_cell_output = false;
end

% Convert values to matrix
if iscell(value), value=cell2mat(value); end

% If format contains 'g', convert to 'f'
if contains(format, 'g')
    format(ismember(format,'g')) = 'f';
    remove_zeros = true;
else
    remove_zeros = false;
end
 
% Convert values to strings
str = reshape(strrep(cellstr(num2str(value(:),format)),' ',''),size(value));

if remove_zeros
    for ival = 1:length(str)
        if contains(str{ival}, '.')
            % Get all the digits after the point
            after_the_point = str{ival}(find(ismember(str{ival}, '.'))+1:end);
            % Invert and remove all digits up to the last non-zero digit
            after_the_point = fliplr(after_the_point);
            idx_nonzero = find(~ismember(after_the_point, '0'), 1, 'first');
            if isempty(idx_nonzero)  % all digits are 0s
                str{ival}(find(ismember(str{ival}, '.')):end) = [];
            else
                after_the_point = fliplr(after_the_point(idx_nonzero:end));
                str{ival}(find(ismember(str{ival}, '.')):end) = [];
                str{ival} = [str{ival}, '.', after_the_point];
            end
        end
    end
end

if length(value) == 1 && ~force_cell_output
    str = str{1};
end
