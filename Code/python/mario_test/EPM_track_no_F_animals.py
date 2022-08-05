import os
import sys
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
from pyparsing import col
import seaborn as sns
import numpy as np
from scipy.signal import resample
import matplotlib.pyplot as plt
import h5py
from openpyxl import load_workbook
from matplotlib import cm
from mpl_toolkits import mplot3d
from Utilities_miniscope import idx2range, rotate, segment_epm_maze
from scipy.signal import find_peaks
from scipy.stats import stats
import matplotlib.backends.backend_pdf
from astropy.convolution import convolve
from astropy.convolution.kernels import Gaussian2DKernel
import scipy.io

################################################################################
np.set_printoptions(suppress=True)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
################################################################################
date_exp = '190916'
animal_id = 'MA_28epi'
# dFF_trial = 191


root_path_data = r"V:\Ca_imaging_pain\6_data\ACC_SNI_anxiety\EPM_results\%s" % (
    animal_id)

excel_file_path = r"T:\Mario\Behavior\EPM\raw_files"
filename = os.path.join(excel_file_path, date_exp + '_' + animal_id + '.xlsx')
wb = load_workbook(filename, read_only=True)
ws = wb['Track-Arena 1-Subject 1']
df = pd.DataFrame(ws.values)
skip_rows = int(df.loc[0, 1])
header = df.loc[1:skip_rows - 4, 0:1].copy()
header.reset_index(drop=True, inplace=True)

info = dict()
for row in range(header.shape[0]):
    key = header.loc[row, 0]
    value = header.loc[row, 1]
    info[key] = value

DATA = df.loc[skip_rows:, :].copy()
column_names = df.loc[skip_rows - 2].values
DATA.columns = column_names

DATA.reset_index(drop=True, inplace=True)

epm_position = DATA[
    ['Recording time', 'X center', 'Y center', 'Velocity']].copy()
epm_position.columns = ['time', 'x', 'y', 'vel']
epm_position.replace("-", np.nan,
                     inplace=True)  # replace the '-' values for NaN
epm_position['time'] = epm_position['time'].map(float)
epm_position['x'] = epm_position['x'].map(float)
epm_position['vel'] = epm_position['vel'].map(float)

if np.isnan(epm_position.loc[0, 'x']):
    idx = idx2range(np.where(np.isnan(epm_position['x']))[0])
    epm_position.loc[0, 'x'] = epm_position.loc[idx[0, 1] + 1, 'x']
epm_position['x'] = epm_position['x'].interpolate()
if np.isnan(epm_position.loc[0, 'y']):
    idx = idx2range(np.where(np.isnan(epm_position['y']))[0])
    epm_position.loc[0, 'y'] = epm_position.loc[idx[0, 1] + 1, 'y']
epm_position['y'] = epm_position['y'].interpolate()
if np.isnan(epm_position.loc[0, 'vel']):
    idx = idx2range(np.where(np.isnan(epm_position['vel']))[0])
    epm_position.loc[0, 'vel'] = epm_position.loc[idx[0, 1] + 1, 'vel']
epm_position['vel'] = epm_position['vel'].interpolate()

qx, qy = rotate(origin=[0, 0], points=epm_position[['x', 'y']].values, angle=45)
epm_position['x'] = qx
epm_position['y'] = qy

n_samples_fluorescence = 3000

time = np.linspace(epm_position.loc[0, 'time'],
                   epm_position.loc[epm_position.shape[0] - 1, 'time'],
                   n_samples_fluorescence)

epm_position_and_F_resampled = pd.DataFrame(columns=['time', 'x', 'y', 'vel'])
epm_position_and_F_resampled['time'] = time
epm_position_and_F_resampled['x'] = resample(epm_position["x"],
                                             num=n_samples_fluorescence)
epm_position_and_F_resampled['y'] = resample(epm_position["y"],
                                             num=n_samples_fluorescence)
epm_position_and_F_resampled['vel'] = resample(epm_position["vel"],
                                               num=n_samples_fluorescence)

X = epm_position_and_F_resampled['x'].values
Y = epm_position_and_F_resampled['y'].values
VEL = epm_position_and_F_resampled['vel'].values

open_idx, closed_idx, center_idx = segment_epm_maze(X, Y,
                                                    open_arms='top_bottom')
open_idx = np.hstack(open_idx)
closed_idx = np.hstack(closed_idx)
epm_position_and_F_resampled['open_arm'] = False
epm_position_and_F_resampled.loc[open_idx, 'open_arm'] = True
epm_position_and_F_resampled['closed_arm'] = False
epm_position_and_F_resampled.loc[closed_idx, 'closed_arm'] = True
epm_position_and_F_resampled['center_arm'] = False
epm_position_and_F_resampled.loc[center_idx, 'center_arm'] = True

# Velocity per arm
arms = ['open', 'closed', 'center']
vel_columns = ['vel']

Vel_open = epm_position_and_F_resampled.loc[
    epm_position_and_F_resampled['open_arm'], vel_columns]
Vel_closed = epm_position_and_F_resampled.loc[
    epm_position_and_F_resampled['closed_arm'], vel_columns]
Vel_center = epm_position_and_F_resampled.loc[
    epm_position_and_F_resampled['center_arm'], vel_columns]
mean_vel_arms = np.mean(Vel_open.values), np.mean(Vel_closed.values), np.mean(
    Vel_center.values)
mean_vel_open = mean_vel_arms[0]
mean_vel_closed = mean_vel_arms[1]
mean_vel_center = mean_vel_arms[2]

duration_frame = time[1] - time[0]
n_frames_spent_in_arms, _ = np.histogram(epm_position_and_F_resampled['open_arm'], np.arange(0, 4))

save_csv_metadata_path = os.path.join(root_path_data,
                                      r"epm_metadata_%s_%s.csv") % (
                         date_exp, animal_id)
data_from_analysis = pd.DataFrame(
    columns=['Rel time spent in arms', 'mean vel open', 'mean vel close',
             'mean vel center'])

data_from_analysis['Rel time spent in arms'] = relative_time_in_arms
data_from_analysis['mean vel open'] = mean_vel_arms[0]
data_from_analysis['mean vel close'] = mean_vel_arms[1]
data_from_analysis['mean vel center'] = mean_vel_arms[2]
# data_from_analysis['Cells resp in open'] = cell_responding_in_open
# data_from_analysis['Cells resp in center']= cell_responding_in_center
# data_from_analysis['prop of cells responding in open'] =  prop_of_cells_responding_in_open
# data_from_analysis['prop of cells responding in center'] = prop_of_cells_responding_in_center
# data_from_analysis['prop cells overlap open and center'] = prop_cell_overlap_open_center
# data_from_analysis['total cells']= n_cells
save_metadata = input("Save metadata too?")
data_from_analysis.to_csv(save_csv_metadata_path)




