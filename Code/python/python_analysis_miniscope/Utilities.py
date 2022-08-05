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
        dataIDX = np.hstack((dataIDX, np.atleast_2d(dataIDX[:, 1] - dataIDX[:, 0] + 1).transpose()))

    else:
        dataIDX = np.empty((0, 3), dtype=int)

    return dataIDX


def rotate(origin, points, angle):
    """
    Rotate a point clockwise by a given angle around a given origin.

    The angle should be given in degrees.
    """
    ox, oy = origin
    angle_rad = angle * np.pi / 180

    qx = ox + np.cos(angle_rad) * (points[:, 0] - ox) - np.sin(angle_rad) * (points[:, 1] - oy)
    qy = oy + np.sin(angle_rad) * (points[:, 0] - ox) + np.cos(angle_rad) * (points[:, 1] - oy)
    return qx, qy


def segment_epm_maze(x, y, open_arms):
    left_arm = np.where((x <= -10) & (x >= -40) & (y <= 10) & (y >= -5))[0]
    right_arm = np.where((x <= 40) & (x >= 5) & (y <= 10) & (y >= -5))[0]
    top_arm = np.where((x <= 5) & (x >= -5) & (y <= 40) & (y >= 7))[0]
    bottom_arm = np.where((x <= 5) & (x >= -5) & (y <= -3) & (y >= -40))[0]
    center_arm = np.where((x < 5) & (x > -10) & (y < 7) & (y > -3))[0]

    if open_arms == 'top_bottom':
        open_arms_pos = [top_arm, bottom_arm]
        closed_arms_pos = [left_arm, right_arm]
    elif open_arms == 'left_right':
        open_arms_pos = [left_arm, right_arm]
        closed_arms_pos = [top_arm, bottom_arm]

    return open_arms_pos, closed_arms_pos, center_arm
