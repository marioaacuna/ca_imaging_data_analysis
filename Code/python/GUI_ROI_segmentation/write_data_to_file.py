import sys
import numpy as np
from time import time
import tables as tb



def main(source_filename, source_variable_name, destination_filename, destination_shape, source_slice, destination_slice, file_format):
    """

    :param source_filename:
    :param source_variable_name:
    :param destination_filename:
    :param destination_shape:
    :param source_slice:
    :param destination_slice:
    :return:
    """
    # Open .mat file for reading (this context will automatically close the file on errors)
    t0 = time()
    with tb.File(source_filename, 'r') as file_in:
        node = file_in.get_node('/%s' % source_variable_name)
        # Get frames
        frames = node.read(start=source_slice[0].start, stop=source_slice[0].stop, step=source_slice[0].step)

    # Replace NaNs with 0s
    frames[np.isnan(frames)] = 0
    # Transpose dimensions
    frames = frames.transpose((0, 2, 1))

    print('\t%.4fs to read frames' % (time() - t0))

    # Open destination file for reading and writing
    file_out = np.memmap(destination_filename, dtype=np.float32, mode='r+', shape=destination_shape)

    # Copy frames
    t0 = time()
    if file_format == 'time':
        # Copy data collapsing all pixels in a row
        n_frames_read = frames.shape[0]
        file_out[destination_slice] = frames.reshape((n_frames_read, -1)).transpose()

    elif file_format == 'frame':
        file_out[destination_slice] = frames

    print('\t%.4fs to copy data' % (time() - t0))


if __name__ == '__main__':
    if len(sys.argv) > 1:  # Command line
        # Unpack arguments
        args = dict(arg.split('=') for arg in sys.argv[1:])

    else:
        argv = 'source_filename=D:\\_MATLAB_CaImaging\\MA_21epi\\MA_21epi_joint.h5" "source_variable_name=1" "destination_filename=D:\\_MATLAB_CaImaging\\MA_21epi\\MA_21epi_data_frame.npy" "destination_shape=(32606, 450, 550)" "source_slice=(slice(0, 10000, None), slice(0, 550, None), slice(0, 450, None))" "destination_slice=(slice(0, 10000, None), slice(None, None, None), slice(None, None, None))" "file_format=frame'
        argv = argv.split('" "')
        args = dict(arg.split('=') for arg in argv)

    # Convert some arguments to objects
    args['destination_shape'] = eval(args['destination_shape'])
    args['source_slice'] = eval(args['source_slice'])
    args['destination_slice'] = eval(args['destination_slice'])

    main(**args)
