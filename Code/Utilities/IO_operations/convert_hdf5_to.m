function convert_hdf5_to(filename_in, filename_out, target_dtype, varargin)

% Unpack inputs
p = inputParser();
addParameter(p, 'scale', true);  % Whether to scale the input data to cover whole range of target type
parse(p, varargin{:});
scale = p.Results.scale;

% Check whether conversion will up- or down-cast the input data
target_bytesize = cast(0, target_dtype);
target_bytesize = whos('target_bytesize');
target_bytesize = target_bytesize.bytes;

file_info = h5info(filename_in);
file_bytesize = file_info.Datasets.Datatype.Size;

is_downcasting = file_bytesize > target_bytesize;

% Create output file
array_size = file_info.Datasets.Dataspace.Size;
dataset_name = ['/', file_info.Datasets.Name];
h5create(filename_out, dataset_name, array_size, 'Datatype',target_dtype);

% Here, we need to do something different in each case
if is_downcasting
    % Load data
    array = h5read(filename_in, dataset_name);
    if any(isnan(array))
        disp('CHECK ME!!!')
        % Check whether there are NaNs. If so, assign the lowest value to them
        keyboard
    end   
    
    if scale
        % Get data type of array
        source_dtype = whos('array');
        source_dtype = source_dtype.class;

        % Get range of target type
        if contains(target_dtype, 'int')
            max_value = intmax(target_dtype);
            min_value = intmin(target_dtype);
        else  % single or double
            max_value = realmax(target_dtype);
            min_value = realmin(target_dtype);
        end
        % Cast values to source data type, so we can use them in next operations
        max_value = cast(max_value, source_dtype);
        min_value = cast(min_value, source_dtype);
    
        % Normalize in range [0, 1]
        array = array - nanmin(array(:));
        array = array ./ nanmax(array(:));
        
        % Normalize to range of target type
        array = array * (max_value - min_value - 1);
        array = array + min_value + 1;
    end
    
    % Write to disk
    h5write(filename_out, dataset_name, array);

else
    keyboard
    
end


%#ok<*NASGU>
