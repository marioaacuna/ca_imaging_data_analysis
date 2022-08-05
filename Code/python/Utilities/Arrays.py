# Imports
import numpy as np


def idx2range(idx):
    # Convert to numpy array
    if not type(idx) is np.ndarray:
        idx = np.array([idx], dtype=int).ravel()

    if idx.shape[0] > 1:
        # Find discontinuities in index
        dataIDX = np.atleast_2d(np.unique(np.hstack((0, np.where(np.diff(idx) > 1)[0]+1)))).transpose()
        dataIDX = np.hstack((dataIDX, np.atleast_2d(np.hstack((dataIDX[1:,0]-1, idx.shape[0]-1))).transpose()))
        # Get original values
        dataIDX = idx[dataIDX]

        # Add column for duration
        dataIDX = np.hstack((dataIDX, np.atleast_2d(dataIDX[:,1] - dataIDX[:,0] + 1).transpose()))

    else:
        dataIDX = np.empty((0, 3), dtype=int)

    return dataIDX
