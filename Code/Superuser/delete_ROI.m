function delete_ROI(animal_ID, ROIs_to_delete)
%DELETE_ROI Deletes ROIs from dataset.

% Clear Command Window
clc
disp(['Deleting ROI ', num2str(ROIs_to_delete), ' from ', animal_ID])

% Delete ROIs from ROI_info
filename = get_filename_of('ROI_info', animal_ID);
if exist(filename, 'file')
    disp('Loading ROI info')
    load(filename);
    disp('Deleting data')
    ROIs(:, :, ROIs_to_delete) = [];
    ROI_fluorescence(:, ROIs_to_delete) = [];
    ROI_info(ROIs_to_delete, :) = [];
    disp('Saving new file')
    save(filename, 'MAPPING_TABLE', 'modified_MAPPING', 'ROI_fluorescence', 'ROI_info', 'ROIs', 'TRANSFORMATION_IDX', 'TRANSFORMATION_MATRICES', '-v6')
else
    disp(['File "', filename, '" does not exist'])
end

% Delete ROIs from deconvolved traces
filename = get_filename_of('dFF', animal_ID);
if exist(filename, 'file')
    disp('Loading dFF')
    load(filename);
    disp('Deleting data')
    dFF(ROIs_to_delete, :) = [];
    disp('Saving new file')
    save(filename, 'dFF', '-v7.3')
else
    disp(['File "', filename, '" does not exist'])
end

% Delete ROIs from deconvolved traces
filename = get_filename_of('spikes', animal_ID);
if exist(filename, 'file')
    disp('Loading spikes')
    load(filename);
    disp('Deleting data')
    Ca_events(ROIs_to_delete) = [];
    spikes(ROIs_to_delete, :) = [];
    disp('Saving new file')
    save(filename, 'Ca_events','spikes', '-v7.3')
else
    disp(['File "', filename, '" does not exist'])
end

disp('done')

%#ok<*LOAD,*NASGU>
