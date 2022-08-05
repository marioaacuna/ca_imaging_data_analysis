function f = open_motion_corrected_file(animal_ID, the_structural_channel)

% Get general_configs
global GC

if ~exist('the_structural_channel','var') || isempty(the_structural_channel)
    the_structural_channel = false;
end

% Open output of NoRMCorre
if ~the_structural_channel
    local_normcorre_filename = [GC.temp_dir, animal_ID, '.mat'];
    if exist(local_normcorre_filename, 'file')
        % Get access to file if not done before or if we have to read a different file
        f = matfile(local_normcorre_filename, 'Writable',false);

    else
        % Open remote file for reading
        remote_motioncorrected_file = get_filename_of('motion_corrected', animal_ID);
        f = matfile(remote_motioncorrected_file, 'Writable',false);
    end

else
    local_normcorre_filename = [GC.temp_dir, animal_ID, '_structural_channel.mat'];
    if exist(local_normcorre_filename, 'file')
        % Get access to file if not done before or if we have to read a different file
        f = matfile(local_normcorre_filename, 'Writable',false);

    else
        % Open remote file for reading
        remote_motioncorrected_file = get_filename_of('motion_corrected', animal_ID, the_structural_channel);
        f = matfile(remote_motioncorrected_file, 'Writable',false);
    end
end
