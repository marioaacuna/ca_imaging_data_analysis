from time import time
import numpy as np
from Utilities import matlab_file


def get_mean_luminance(f, indices, time_first=True):
    t0 = time();
    if time_first:
        np.mean(f[indices, :], axis=0)
    else:
        np.mean(f[:, indices[0], indices[1]], axis=0);
    t1 = time()

    return t1 - t0


use_time_first_file = False
tries = 5


filename_data_frame = r"D:\_MATLAB_2PI\FK_11_data_frame.npy"
filename_data_time = r"D:\_MATLAB_2PI\FK_11_data_frame.npy"
params_filename = r'D:\_MATLAB_2PI\FK_11_params.mat'
PARAMETERS = matlab_file.load(params_filename, variable_names='PARAMETERS')

# Get indices of pixels in ROIs
ROI_pixels_small = np.array([32083, 32084, 32085, 32338, 32339, 32340, 32341, 32342, 32593, 32594, 32595, 32596, 32597, 32598, 32599, 32600, 32849, 32850, 32851, 32852, 32853, 32854, 32855, 33104, 33105, 33106, 33107, 33108, 33109, 33110, 33360, 33361, 33362, 33363, 33364, 33365, 33366, 33618, 33619, 33620], dtype=np.int64)
ROI_pixels_large = np.array([21366, 21367, 21619, 21620, 21621, 21622, 21623, 21624, 21625, 21872, 21873, 21874, 21875, 21876, 21877, 21878, 21879, 21880, 21881, 21882, 22125, 22126, 22127, 22128, 22129, 22130, 22131, 22132, 22133, 22134, 22135, 22136, 22137, 22138, 22139, 22378, 22379, 22380, 22381, 22382, 22383, 22384, 22385, 22386, 22387, 22388, 22389, 22390, 22391, 22392, 22393, 22394, 22395, 22396, 22631, 22632, 22633, 22634, 22635, 22636, 22637, 22638, 22639, 22640, 22641, 22642, 22643, 22644, 22645, 22646, 22647, 22648, 22649, 22650, 22651, 22652, 22653, 22654, 22885, 22886, 22887, 22888, 22889, 22890, 22891, 22892, 22893, 22894, 22895, 22896, 22897, 22898, 22899, 22900, 22901, 22902, 22903, 22904, 22905, 22906, 22907, 22908, 22909, 22910, 22911, 23141, 23142, 23143, 23144, 23145, 23146, 23147, 23148, 23149, 23150, 23151, 23152, 23153, 23154, 23155, 23156, 23157, 23158, 23159, 23160, 23161, 23162, 23163, 23164, 23165, 23166, 23167, 23398, 23399, 23400, 23401, 23402, 23403, 23404, 23405, 23406, 23407, 23408, 23409, 23410, 23411, 23412, 23413, 23414, 23415, 23416, 23417, 23418, 23419, 23420, 23421, 23422, 23654, 23655, 23656, 23657, 23658, 23659, 23660, 23661, 23662, 23663, 23664, 23665, 23666, 23667, 23668, 23669, 23670, 23671, 23672, 23673, 23674, 23675, 23676, 23677, 23678, 23910, 23911, 23912, 23913, 23914, 23915, 23916, 23917, 23918, 23919, 23920, 23921, 23922, 23923, 23924, 23925, 23926, 23927, 23928, 23929, 23930, 23931, 23932, 23933, 24167, 24168, 24169, 24170, 24171, 24172, 24173, 24174, 24175, 24176, 24177, 24178, 24179, 24180, 24181, 24182, 24183, 24184, 24185, 24186, 24187, 24188, 24189, 24423, 24424, 24425, 24426, 24427, 24428, 24429, 24430, 24431, 24432, 24433, 24434, 24435, 24436, 24437, 24438, 24439, 24440, 24441, 24442, 24443, 24444, 24679, 24680, 24681, 24682, 24683, 24684, 24685, 24686, 24687, 24688, 24689, 24690, 24691, 24692, 24693, 24694, 24695, 24696, 24697, 24698, 24699, 24700, 24937, 24938, 24939, 24940, 24941, 24942, 24943, 24944, 24945, 24946, 24947, 24948, 24949, 24950, 24951, 24952, 24953, 24954, 24955, 25195, 25196, 25197, 25198, 25199, 25200, 25201, 25202, 25203, 25204, 25205, 25206, 25207, 25208, 25454, 25455, 25456, 25457, 25458, 25459, 25460, 25461, 25712, 25713], dtype=np.int64)
ROI_pixels_random = np.arange(PARAMETERS['n_pixels'])
np.random.shuffle(ROI_pixels_random)
ROI_pixels_random_small = ROI_pixels_random[:ROI_pixels_small.shape[0]]
ROI_pixels_random_large = ROI_pixels_random[:ROI_pixels_large.shape[0]]
# Convert linear to 2d indices
if not use_time_first_file:
    ROI_pixels_small = np.unravel_index(ROI_pixels_small, (PARAMETERS['frame_height'], PARAMETERS['frame_width']), order='C')
    ROI_pixels_large = np.unravel_index(ROI_pixels_large, (PARAMETERS['frame_height'], PARAMETERS['frame_width']), order='C')
    ROI_pixels_random_small = np.unravel_index(ROI_pixels_random_small, (PARAMETERS['frame_height'], PARAMETERS['frame_width']), order='C')
    ROI_pixels_random_large = np.unravel_index(ROI_pixels_random_large, (PARAMETERS['frame_height'], PARAMETERS['frame_width']), order='C')


# Open files
if use_time_first_file:
    f = np.memmap(PARAMETERS['filename_data_time'], dtype=np.float32, mode='r', offset=0).reshape((PARAMETERS['n_pixels'], PARAMETERS['n_frames']))
    print('Testing speed of time-first file')
else:
    f = np.memmap(PARAMETERS['filename_data_frame'], dtype=np.float32, mode='r', offset=0).reshape((PARAMETERS['n_frames'], PARAMETERS['frame_height'], PARAMETERS['frame_width']))
    print('Testing speed of frame-first file')

# 1 pixel
elapsed_time = np.zeros((tries, ), dtype=np.float64)
for trial in range(tries):
    elapsed_time[trial] = get_mean_luminance(f, ROI_pixels_small[0], time_first=use_time_first_file)
print('Mean 1 pixel: %.5fs' % elapsed_time.mean())

# Small ROI
elapsed_time = np.zeros((tries, ), dtype=np.float64)
for trial in range(tries):
    elapsed_time[trial] = get_mean_luminance(f, ROI_pixels_small, time_first=use_time_first_file)
print('Mean small ROI (40 pixels): %.5fs' % elapsed_time.mean())

# Small random ROI
elapsed_time = np.zeros((tries, ), dtype=np.float64)
for trial in range(tries):
    elapsed_time[trial] = get_mean_luminance(f, ROI_pixels_random_small, time_first=use_time_first_file)
print('Mean 40 random pixels: %.5fs' % elapsed_time.mean())

# Large ROI
elapsed_time = np.zeros((tries, ), dtype=np.float64)
for trial in range(tries):
    elapsed_time[trial] = get_mean_luminance(f, ROI_pixels_large, time_first=use_time_first_file)
print('Mean large ROI (300 pixels): %.5fs' % elapsed_time.mean())

# Large random ROI
elapsed_time = np.zeros((tries, ), dtype=np.float64)
for trial in range(tries):
    elapsed_time[trial] = get_mean_luminance(f, ROI_pixels_random_large, time_first=use_time_first_file)
print('Mean 300 random pixels: %.5fs' % elapsed_time.mean())
