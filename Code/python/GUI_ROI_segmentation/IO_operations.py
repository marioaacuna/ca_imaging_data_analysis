# System packages
import sys
import os
# Numerical packages
import numpy as np
import tables as tb


def log(message):
    sys.stdout.write(message + '\n')
    sys.stdout.flush()


def prepare_data_for_GUI(PARAMETERS, n_frames_per_chunk=1000):
    # Calculate number of chunks
    n_chunks = int(np.ceil(np.float64(PARAMETERS['n_frames']) / n_frames_per_chunk))
    if n_chunks == 1:
        n_frames_per_chunk = PARAMETERS['n_frames']

    # Get path of python executable and script
    python_exe = sys.executable
    file_path = os.path.realpath(__file__)
    python_script = os.path.join(os.path.split(file_path)[0], 'write_data_to_file.py')
    base_cmd = 'conda activate env_ca_imaging && "%s" "%s" ' % (python_exe, python_script)

    # Check which file should be used as input
    if os.path.exists(PARAMETERS['filename_local']):
        filename = PARAMETERS['filename_local']
    else:
        filename = PARAMETERS['filename_remote']

    # Get name of hdf5 dataset to read
    variable_name = PARAMETERS['variable_name']

    # Get last modification date of data file
    data_file_modification_time = os.path.getmtime(filename)

    # Check whether files need conversion
    do_convert_data_frame_file = False  # Initialize flag
    if not os.path.exists(PARAMETERS['filename_data_frame']):
        do_convert_data_frame_file = True
    else:
        # Compare last modification time of this file to the original data
        # file: If data file is later, reconvert this file
        data_frame_modification_time = os.path.getmtime(PARAMETERS['filename_data_frame'])
        if data_file_modification_time > data_frame_modification_time:
            do_convert_data_frame_file = True

    do_convert_data_time_file = False  # Initialize switch
    if not os.path.exists(PARAMETERS['filename_data_time']):
        do_convert_data_time_file = True
    else:
        # Compare last modification time of this file to the original data
        # file: If data file is later, reconvert this file
        data_time_modification_time = os.path.getmtime(PARAMETERS['filename_data_time'])
        if data_file_modification_time > data_time_modification_time:
            do_convert_data_time_file = True

    if do_convert_data_frame_file:
        # Copy entire file in a memory-mapped numpy array
        log('Copying whole dataset in \'%s\'' % PARAMETERS['filename_data_frame'])
        destination_shape = (PARAMETERS['n_frames'], PARAMETERS['frame_height'], PARAMETERS['frame_width'])
        file_data_frame = np.memmap(PARAMETERS['filename_data_frame'], dtype=np.float32, mode="w+", shape=destination_shape)
        # Loop through chunks
        with tb.File(filename, 'r') as file_in:
            node = file_in.get_node('/%s' % variable_name)

            for ichunk in range(n_chunks):
                log('Copying chunk %i/%i' % (ichunk + 1, n_chunks))
                # Get edges of chunk
                start_sample = int(ichunk * n_frames_per_chunk)
                end_sample = int(np.clip(start_sample + n_frames_per_chunk, a_min=start_sample + 1, a_max=PARAMETERS['n_frames']))
                # Make slice objects
                source_slice = (slice(start_sample, end_sample), slice(0, PARAMETERS['frame_height']), slice(0, PARAMETERS['frame_width']))
                # Swap x and y axes
                source_slice = (source_slice[0], source_slice[2], source_slice[1])
                destination_slice = (slice(start_sample, end_sample), slice(None), slice(None))

                # Get frames
                frames = node.read(start=source_slice[0].start, stop=source_slice[0].stop, step=source_slice[0].step)
                # Replace NaNs with 0s
                frames[np.isnan(frames)] = 0
                # Transpose dimensions
                frames = frames.transpose((0, 2, 1))
                # Store frames
                file_data_frame[destination_slice] = frames

    if do_convert_data_time_file:
        # Open file for writing
        log('Copying time-data in %s' % PARAMETERS['filename_data_time'])
        destination_shape = (PARAMETERS['n_pixels'], PARAMETERS['n_frames'])
        # file_data_time = np.memmap(PARAMETERS['filename_data_time'], dtype=np.float32, mode="w+", shape=destination_shape)
        data_frame = np.memmap(PARAMETERS['filename_data_frame'], dtype=np.float32, mode='r', offset=0, shape=(PARAMETERS['n_frames'], PARAMETERS['frame_height'], PARAMETERS['frame_width']))
        # Create file
        array_name = 'time_first'
        with tb.open_file(PARAMETERS['filename_data_time'], 'w', title=PARAMETERS['animal_ID']) as hdf5_file:
            # Remove node if it exists
            if hdf5_file.__contains__('/%s' % array_name):
                hdf5_file.remove_node('/', name=array_name, recursive=False)
            # Create array
            hdf5_array = hdf5_file.create_array('/', array_name, title=array_name, atom=tb.Float32Col(), shape=destination_shape)

            # Loop through chunks
            for ichunk in range(n_chunks):
                log('Reshaping chunk %i/%i' % (ichunk + 1, n_chunks))
                # Get edges of chunk
                start_sample = int(ichunk * n_frames_per_chunk)
                end_sample = int(np.clip(start_sample + n_frames_per_chunk, a_min=start_sample + 1, a_max=PARAMETERS['n_frames']))
                # Make slice objects
                source_slice = (slice(start_sample, end_sample), slice(0, PARAMETERS['frame_height']), slice(0, PARAMETERS['frame_width']))
                # Swap x and y axes
                # source_slice = (source_slice[0], source_slice[2], source_slice[1])
                destination_slice = (slice(0, PARAMETERS['n_pixels']), slice(start_sample, end_sample))

                # Reading
                frames = data_frame[source_slice]
                n_frames_read = frames.shape[0]
                frames = frames.reshape((n_frames_read, -1)).transpose()

                # Copy data
                hdf5_array[destination_slice] = frames

    if do_convert_data_frame_file or do_convert_data_time_file:
        log('Finished file conversion')

    return PARAMETERS

