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

np.set_printoptions(suppress=True)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

#######################################################################
############# I . Load the EPM track and the F traces ############
########################################################################


## Load the excel file with the coordinates and timestamp of the behavior
#########################################################################



# date_exp = input ("What date do you want to take?")
# trial_id = input ("What trial do you want to take?")

def epm_track(date_exp , animal_id, dFF_trial):
    root_path_data = r"V:\Ca_imaging_pain\6_data\ACC_SNI_anxiety\EPM_results\%s" % (
        animal_id)
    root_path_plots = r"V:\Ca_imaging_pain\6_plots\ACC_SNI_anxiety\EPM_figures\%s" % (
        animal_id)
    if not os.path.exists(root_path_data):
        os.makedirs(root_path_data)
    if not os.path.exists(root_path_plots):
        os.makedirs(root_path_plots)

    dFF_filename = r"V:\Ca_imaging_pain\4_fluorescence_traces\%s_raw_deltaF_over_F.mat" % animal_id
    excel_file_path = r"T:\Mario\Behavior\EPM\raw_files"
    filename = os.path.join(excel_file_path, date_exp + '_' + animal_id + '.xlsx')
    wb = load_workbook(filename, read_only=True)
    ws = wb['Track-Arena 1-Subject 1']
    df = pd.DataFrame(ws.values)
    skip_rows = int(df.loc[0, 1])
    header = df.loc[1:skip_rows-4, 0:1].copy()
    header.reset_index(drop=True, inplace=True)


    ## Assign a dict with key and values
    info = dict()
    for row in range(header.shape[0]):
        key = header.loc[row, 0]
        value = header.loc[row, 1]
        info[key] = value

    # animal_id = info['animal ID']
    DATA = df.loc[skip_rows:, :].copy()
    column_names = df.loc[skip_rows - 2].values
    DATA.columns = column_names

    DATA.reset_index(drop=True, inplace=True)

    epm_position = DATA[['Recording time', 'X center', 'Y center', 'Velocity']].copy()
    epm_position.columns = ['time', 'x', 'y', 'vel']
    epm_position.replace("-", np.nan, inplace=True) # replace the '-' values for NaN
    epm_position['time'] = epm_position['time'].map(float)
    epm_position['x'] = epm_position['x'].map(float)
    epm_position['vel'] = epm_position['vel'].map(float)

    ## Interpolate the NaN values from the EPM data (x,y).

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

    ####################################################
    # plotting scatter plot for checking a given cell
    ####################################################

    # #plotting_a_cell = input("Do you also want to plot a cell?")
    # #plotting_a_cell == "yes":
    # #choose_a_cell = input("What cell do you want tpo plot?")
    # #x = int(choose_a_cell)-1
    # #cell_x = ("cell_"+str(x))
    # fig = plt.figure()
    # ax = fig.add_subplot(111)#, projection = '3d' )
    # #ax = Axes3D(fig)
    # #plt.title(cell_x)
    # ax.scatter(epm_position.loc[idx, 'x'], epm_position.loc[idx, 'y'])#, epm_position.loc[idx, cell_x])
    # plt.show()

    ## Rotate coordinates around (0, 0)

    qx, qy = rotate(origin=[0, 0], points=epm_position[['x', 'y']].values, angle=45)
    epm_position['x'] = qx
    epm_position['y'] = qy



    ## load the mat file and get the df/f values per ROI
    #####################################################

    # For miniPipe data use this line:
    #f = h5py.File('msCam_data_processed_' + str(date_exp) + '_'+ "I"+str(animal_id)+'.mat')
    #df_f_data = f.get("sigfn")[:]

    # For Grewe's pipeline use this line.
    # miniscope_file = ('N:/Mario/miniscope data/' +str (date_exp)+ '_' +'I'+str(animal_id)+ '/Session1/sorted/PCAICAsorted.mat')
    # f = scipy.io.loadmat(miniscope_file)
    # df_f_data = f.get("traces")

    print('reading dFF')
    with h5py.File(dFF_filename, mode='r') as f:
        h5_references = f.get('dFF')[:]
        n_cells = h5_references.shape[1]
        # Gather names of actual data where references point to
        h5names = [h5py.h5r.get_name(h5_references[dFF_trial - 1, i_cell], f.id) for i_cell in range(n_cells)]
        # Read actual data
        data = np.hstack([f[name].value for name in h5names]).T
    print('done')
    df_f_data_t = data.copy()

    #f.close()
    #f = h5py.File(miniscope_file)
    #f = h5py.File('PCAICAsorted.mat')


    # with h5py.File('msCam_data_processed_190130_I26_3.mat') as f: # another option to open a file
    #     data = f.get("sigfn")



    # df_f_data_t = df_f_data.transpose()
    # df_f_data_t = df_f_data#df_f_data.transpose()
    # !!!!! for Grewe's pipeline do not transpose!!!

    all_normalized_amp = df_f_data_t / np.max(df_f_data_t, axis=1).reshape((-1, 1))
    all_F_raw_amp = df_f_data_t
    n_samples_fluorescence = all_normalized_amp.shape[1]

    ###  Transform the time of the EPM relative to the nr frames of the df_f data.
    time = np.linspace(epm_position.loc[0, 'time'], epm_position.loc[epm_position.shape[0] - 1, 'time'], n_samples_fluorescence)



    ##  Convert the data from position in DataFrame
    ###############################################

    epm_position_and_F_resampled = pd.DataFrame(columns=['time', 'x', 'y', 'vel'])
    epm_position_and_F_resampled['time'] = time
    epm_position_and_F_resampled['x'] = resample(epm_position["x"], num=n_samples_fluorescence)
    epm_position_and_F_resampled['y'] = resample(epm_position["y"], num=n_samples_fluorescence)
    epm_position_and_F_resampled['vel'] = resample(epm_position["vel"], num=n_samples_fluorescence)

    #epm_position_resampled['cell'] = cell_x


    ## Add to the DataFrame epm_position_and_F_resampled the values of normalized df/f of each cell.


    #'%s_cell_%i_%s' % ('MA_1', 2121, '1000')
    for i in range(len(all_normalized_amp[:,0])):
        epm_position_and_F_resampled["cell_%i" % (i+1)] = all_normalized_amp[i,:]

    # print(epm_position_and_F_resampled)
    print("Animal id: "+ str(animal_id))
    print('Date : ' + str(date_exp))
    # plt.clf(); plt.plot(epm_position['time'], epm_position['x'], '-or', label='original'), plt.plot(time, epm_position_x_resampled, '-ok', label='resampled'); plt.legend()


    ## Save the DataFrame epm_position_and_F_resampled in a folder
    ################################################
    # save_csv_path = r"V:\Ca_imaging_pain\6_data\ACC_SNI_anxiety\EPM_results\%s\epm_data_%s_%s.csv" % (animal_id, date_exp, animal_id)
    # save_option = input("do you want to save the data as csv? ")
    # if save_option == "yes":
    save_csv_path = os.path.join(root_path_data, r"EPM_position_and_F_resampled_%s_%s.csv") % (date_exp, animal_id)
    epm_position_and_F_resampled.to_csv(save_csv_path)

    # save_figures_option = input("do you want to save figures in pdf?")
    # if save_figures_option == 'yes':
    save_pdf_path = os.path.join(root_path_plots,r"%s_%s.pdf") % (date_exp, animal_id)
    pdf = matplotlib.backends.backend_pdf.PdfPages (save_pdf_path)




    ########################################################################
    ################  II   Compute values of peak amplitude and F rate per arm
    ########################################################################
    # t0 = 0
    # t1 = 600
    # idx = np.where((epm_position_and_F_resampled['time'] >= t0) & (epm_position_and_F_resampled['time'] <= t1))[0]


    X = epm_position_and_F_resampled['x'].values
    Y = epm_position_and_F_resampled['y'].values
    VEL = epm_position_and_F_resampled['vel'].values

    open_idx, closed_idx, center_idx = segment_epm_maze(X, Y, open_arms='top_bottom')
    open_idx = np.hstack(open_idx)
    closed_idx = np.hstack(closed_idx)
    epm_position_and_F_resampled['open_arm'] = False
    epm_position_and_F_resampled.loc[open_idx, 'open_arm'] = True
    epm_position_and_F_resampled['closed_arm'] = False
    epm_position_and_F_resampled.loc[closed_idx, 'closed_arm'] = True
    epm_position_and_F_resampled['center_arm'] = False
    epm_position_and_F_resampled.loc[center_idx, 'center_arm'] = True


    #    Plot the EPM trajectory of this mouse
    # plotting_EPM_areas = input("do you want to print the epm areas?")
    # if plotting_EPM_areas == 'yes':
    #     print("Plotting trajectory of animal :" + str(animal_id))
    marker_size = 3
    fig = plt.figure()
    ax_1 = fig.add_subplot(1, 1, 1)
    ax_1.cla()
    ax_1.plot(X[center_idx], Y[center_idx], 'o', color=[.7, .7, .7], label='center', ms=marker_size)
    ax_1.plot(X[open_idx], Y[open_idx], 'or', label='open', ms=marker_size)
    ax_1.plot(X[closed_idx], Y[closed_idx], 'ok', label='closed', ms=marker_size)
    ax_1.legend()
    # if save_figures_option == 'yes':
    pdf.savefig(fig)


    # Compute amplitude of events in maze arms
    ###########################################################
    cell_columns = epm_position_and_F_resampled.columns
    cell_columns = [i for i in cell_columns if i.startswith('cell_')]
    n_cells = len(cell_columns)
    # Fluorescence per arm
    F_open = epm_position_and_F_resampled.loc[epm_position_and_F_resampled['open_arm'], cell_columns]
    F_closed = epm_position_and_F_resampled.loc[epm_position_and_F_resampled['closed_arm'], cell_columns]
    F_center = epm_position_and_F_resampled.loc[epm_position_and_F_resampled['center_arm'], cell_columns]

    # Velocity per arm
    arms = ['open', 'closed', 'center']
    vel_columns =['vel']

    Vel_open = epm_position_and_F_resampled.loc[epm_position_and_F_resampled['open_arm'], vel_columns]
    Vel_closed = epm_position_and_F_resampled.loc[epm_position_and_F_resampled['closed_arm'], vel_columns]
    Vel_center = epm_position_and_F_resampled.loc[epm_position_and_F_resampled['center_arm'], vel_columns]
    mean_vel_arms = np.mean(Vel_open.values), np.mean(Vel_closed.values), np.mean(Vel_center.values)
    mean_vel_open = mean_vel_arms[0]
    mean_vel_closed = mean_vel_arms[1]
    mean_vel_center = mean_vel_arms[2]



    # fig = plt.figure() # plot an eventual anxiety cell
    # plt.title('Speed / Cell in different arms')
    # plt.xlabel('velocity (cm/s)')
    # plt.ylabel('Normalized Amplitude')
    # plt.scatter(Vel_open, F_open['cell_30'])
    # plt.scatter(Vel_closed, F_closed['cell_30'])
    # plt.legend()
    # plt.xlim(0, 20)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)


    #plot the mean velocity per arm
    fig = plt.figure()
    plt.bar(np.arange(len(arms)), mean_vel_arms,)
    plt.title("Mean velocity per arm.\n " + "Animal: " + str(animal_id))
    plt.ylabel('Average Speed cm/s')
    plt.xticks(np.arange(len(arms)), arms)
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    # Look for peaks in signal
    F = pd.concat((F_open, F_closed, F_center), ignore_index=True)
    threshold = F.median() + 3 * F.mad()
    F_transposed = F.transpose()
    # Add column to mark arms
    F['arm'] = np.hstack([[0] * F_open.shape[0]] + [[1] * F_closed.shape[0]] + [[2] * F_center.shape[0]])
    n_frames_spent_in_arms, _ = np.histogram(F['arm'], np.arange(0, 4))
    duration_frame = time[1] - time[0]
    time_spent_in_arms = n_frames_spent_in_arms * duration_frame
    relative_time_in_arms =(time_spent_in_arms/max(time))

    # Plot pie chart, time spent in each arm
    fig = plt.figure()
    plt.pie(time_spent_in_arms, labels=['open', 'closed', 'center'], autopct='%1.1f%%',  colors= ['r','g', 'y'])
    plt.title("Relative Time Spent in Each Arm.\n Animal:" + str(animal_id))
    plt.axis('square')
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    #####################################################
    # Plot amplitude vs Speed:
    fig = plt.figure()
    plt.title("Mean Amplitude per velocity")
    plt.scatter (VEL[: F.shape[0]],np.mean(F_transposed))
    plt.xlabel('Velocity (cm/s)')
    plt.xlim(0,20)
    plt.ylabel('Normalized amplitude')
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    ## Find cross correlation amplitude/speed and Plot bar plot for each cell and speed
    v = VEL[: F.shape[0]]
    idx = np.where(np.abs(v) <= 20)[0]
    cc = [np.corrcoef(v[idx], F.loc[idx, 'cell_%i' % (i+1)])[0, 1] for i in range(n_cells)]
    fig = plt.figure()
    plt.clf(); plt.bar(np.arange(n_cells), cc)
    plt.title("Correlation coefficients Amplitude and Speed")
    plt.ylim(-0.5, 0.5)
    plt.xlabel("Cell id")
    plt.ylabel("Correlation coefficient")
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    time_here = time
    time_here = time[: F.shape[0]]
    fig = plt.figure()
    fig,axx1 = plt.subplots()
    y1= v
    y2 = F['cell_1']*20
    axx1.plot(time_here,y1, 'b-')
    axx2 = axx1.twinx ()
    axx2.plot(time_here,y2,  'r-')

    axx1.set_xlabel("Time (s)")
    axx1.set_ylabel("Velocity (cm/s)",  color='b')
    axx1.set_ylim(-7, 20)
    axx2.set_ylabel ("dF/F",  color='r')
    axx2.set_ylim(-20, 50)
    #########################################
    ##  III Find Peaks
    ###########################################

    # make an average of all peaks found in a given area
    n_cells = len(cell_columns)
    min_distance_s = 1
    peak_distance = np.where(time <= min_distance_s)[0].shape[0]
    avgF = np.zeros((n_cells, 3), dtype=float)
    avgRate = np.zeros((n_cells, 3), dtype=float)
    total_avgF = np.zeros((n_cells, 3), dtype=float)
    avg_peak_cells = np.zeros((n_cells, 3), dtype=float)
    for i_cell in range(n_cells):
        F_this_cell = F[cell_columns[i_cell]]
        peaks_this_cell, _ = find_peaks(F_this_cell, height=threshold[cell_columns[i_cell]], distance=peak_distance)
        arms_this_cell = np.array([F.loc[i, 'arm'] for i in peaks_this_cell])
        avg_peak_cells [i_cell] = np.mean(F_this_cell.loc[peaks_this_cell])
        for i_area in range(len(arms)):
            idx_peaks_this_area = peaks_this_cell[arms_this_cell == i_area]
            F_peaks_this_area = F_this_cell.loc[idx_peaks_this_area]
            mean_F_peaks_this_area = np.nanmean(F_peaks_this_area)
            avgF[i_cell, i_area] = mean_F_peaks_this_area
            avgRate[i_cell, i_area] = idx_peaks_this_area.shape[0] / time_spent_in_arms[i_area]
            total_avgF [i_cell] =np.mean(avgF[i_cell])



    #####
    ## Calculate the frequency rate and the difference between open and close arms
    ###############################################################################

    avgF[np.isnan(avgF)] = 0
    avgRate[np.isnan(avgRate)] = 0
    total_avgF[np.isnan(total_avgF)] = 0 #not used anymore, used instad : avg_peak_cells

    #cell_name = "cell_3"
    # look_at_signal = input("do you want to look at the signal?" )
    # if look_at_signal == 'yes':
    #cell_name = input("what cell do you wnt to plot?")
    #cell_name = 'cell_1'; plt.clf(); plt.plot(F[cell_name], 'k'); plt.plot(peaks_this_cell[cell_name], 'or')


    ## Get  the difference of Rates and amplitudes (Open - Close) and plot them
    ################################################
    diff_Rate = avgRate[:, 0] - avgRate[:, 1] # 0 is open; 1 is closed; 2 is center
    diff_amp = avgF[:, 0] - avgF[:, 1]
    diff_F_times_Rate = avgF [:,0] * avgRate[:,0] - avgF [:,1] * avgRate[:,1]
    diff_amp_open_arm_total = avgF[:, 0] - avg_peak_cells[:, 0]
    diff_amp_close_arm_total = avgF[:, 1] - avg_peak_cells[:, 0]
    diff_amp_center_arm_total = avgF[:, 2] - avg_peak_cells[:, 0]
    diff_open_vs_close_rel_total_amp = diff_amp_open_arm_total - diff_amp_close_arm_total
    diff_center_vs_close_rel_total_amp = diff_amp_center_arm_total - diff_amp_close_arm_total
    [cell_responding_in_open] = np.where(diff_open_vs_close_rel_total_amp>0) # This is the idx of cells with positive score values in open
    cell_responding_in_center = np.where(diff_center_vs_close_rel_total_amp>0) # This is the idx of cells with positive score values in center
    relative_avgF_time_open = diff_amp_open_arm_total / relative_time_in_arms [0]# This is the relative amplitude of peaks relative to time spent in a location
    relative_avgF_time_close = diff_amp_close_arm_total / relative_time_in_arms [1]# This is the relative amplitude of peaks relative to time spent in a location
    relative_avgF_time_center = diff_amp_center_arm_total / relative_time_in_arms [2]# This is the relative amplitude of peaks relative to time spent in a location
    relative_open_minus_close = relative_avgF_time_open - relative_avgF_time_close
    prop_of_cells_responding_in_open =np.count_nonzero(cell_responding_in_open)/(n_cells)
    prop_of_cells_non_responding_in_open = 1-prop_of_cells_responding_in_open
    prop_of_cells_responding_in_center =np.count_nonzero(cell_responding_in_center)/(n_cells)
    prop_of_cells_non_responding_in_center = 1-prop_of_cells_responding_in_center

    overlap_open_center = np.intersect1d(cell_responding_in_open, cell_responding_in_center)
    prop_cell_overlap_open_center = np.count_nonzero(overlap_open_center)/n_cells

    ## Save Csv with data analysed ---- that can be changed when more variables added
    save_csv_metadata_path = os.path.join (root_path_data, r"epm_metadata_%s_%s.csv") % (date_exp, animal_id)
    data_from_analysis = pd.DataFrame(columns= ['total cells','Rel time spent in arms', 'mean vel open','mean vel close',
                                                'mean vel center',
                                                'prop of cells responding in open','prop of cells responding in center',
                                                'prop cells overlap open and center'])


    data_from_analysis['Rel time spent in arms']=relative_time_in_arms
    data_from_analysis['mean vel open'] =mean_vel_arms [0]
    data_from_analysis['mean vel close'] =mean_vel_arms [1]
    data_from_analysis['mean vel center'] =mean_vel_arms [2]
    #data_from_analysis['Cells resp in open'] = cell_responding_in_open
    #data_from_analysis['Cells resp in center']= cell_responding_in_center
    data_from_analysis['prop of cells responding in open'] =  prop_of_cells_responding_in_open
    data_from_analysis['prop of cells responding in center'] = prop_of_cells_responding_in_center
    data_from_analysis['prop cells overlap open and center'] = prop_cell_overlap_open_center

    data_from_analysis['total cells']= n_cells

    # save_metadata = input("Save metadata too?")
    # if save_metadata == 'yes':
    data_from_analysis.to_csv(save_csv_metadata_path)

    save_json_responding_idx_open_path = os.path.join(root_path_data,r"epm_idx_cells_open_%s_%s.json") % (date_exp, animal_id)
    cells_idx_responding_in_open =  pd.DataFrame (columns=['idx cells responding in open'])
    cells_idx_responding_in_open['idx cells responding in open'] = F_open.columns[cell_responding_in_open]
    cells_idx_responding_in_open.reset_index(drop=True, inplace=True)
    cells_idx_responding_in_open.to_json(save_json_responding_idx_open_path)
    save_json_responding_idx_center_path = os.path.join(root_path_data, r"epm_idx_cells_center_%s_%s.json") % (
    date_exp, animal_id)
    cells_idx_responding_in_center = pd.DataFrame(columns=['idx cells responding in center'])
    cells_idx_responding_in_center['idx cells responding in center'] = F_center.columns[cell_responding_in_center]
    cells_idx_responding_in_center.to_json(save_json_responding_idx_center_path)
    # plt.figure()
    # plt.subplot(3,1, 2)
    # plt.bar(np.arange(n_cells), diff_amp_close_arm_total)
    # plt.title("Avg close - Avg Total ")
    # plt.subplot(3,1, 1)
    # plt.bar(np.arange(n_cells), diff_amp_open_arm_total)
    # plt.title("Avg open - Avg Total ")
    # plt.subplot(3, 1, 3)
    # plt.bar(np.arange(n_cells), diff_amp_center_arm_total)
    # plt.title("Avg center - Avg Total ")

    # plot figure for the difference of amplitude in the open - close relative to the total amplitude average

    ### plot the center and open cells
    fig =plt.figure()
    plt.suptitle("Open and center activation. " + str(animal_id))
    ax =fig.add_subplot(2,1,1)
    ax.bar(np.arange(n_cells), diff_open_vs_close_rel_total_amp)
    #ax.get_children()[4].set_color('r')
    lop = np.zeros((len(overlap_open_center), 1), dtype=int)
    lop.tolist()
    for i in range(len(overlap_open_center)):
        lop = overlap_open_center[i]
        ax.get_children()[lop].set_color('r')
    ax.legend(["Overlap"])
    leg = ax.get_legend()
    leg.legendHandles[0].set_color('red')
    plt.ylabel("Score (Open/Avg - Closed/Avg)")
    ay = fig.add_subplot(2,1,2)
    plt.ylabel("Score Center")
    ay.bar(np.arange(n_cells), diff_center_vs_close_rel_total_amp)
    for i in range(len(overlap_open_center)):
        lop = overlap_open_center[i]
        ay.get_children()[lop].set_color('r')
    #plt.xticks(np.arange(n_cells), (np.arange(0,n_cells)+1))
    plt.xlabel("Cells")
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    # Plot the proportion of cells responding in both open and center
    fig = plt.figure()
    plt.pie([prop_cell_overlap_open_center, (1-prop_cell_overlap_open_center)],labels= ["Overlap resp in open and center", "no overlap"], autopct='%1.1f%%')
    plt.title("Overlap Open and Center responding cells")
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    # plot the difference of amplitude open-close (over the average amplitude) relatice to the time spent in each arm

    # fig = plt.figure()
    # plt.bar(np.arange(n_cells), relative_open_minus_close)
    # plt.title("difference amplitude Open vs Close \n normalized to the average amplitude in the entire session\n and relative to the time spent in each arm " + str(animal_id))
    # plt.xlabel("Cells")
    # plt.ylabel("Preferential Score")
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)


    # # Plot the difference of amplitude sorted
    # #diff_amp.sort(); diff_amp = diff_amp[::-1]  # Sort ascending order
    # fig = plt.figure()
    # plt.clf(); plt.bar(np.arange(n_cells), diff_amp)
    # plt.title("Amplitude Difference (Open - Close)" + str(animal_id))
    # plt.xlabel("Cell id")
    # plt.ylabel("Difference score")
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)

    # # Plot the difference of rates sorted
    # #diff_Rate.sort(); diff_Rate = diff_Rate[::-1]  # Sort ascending order
    # fig = plt.figure()
    # plt.clf(); plt.bar(np.arange(n_cells), diff_Rate)
    # plt.title("Event Rate Difference (Open - Close)"+ str(animal_id))
    # plt.xlabel("Cell sorted")
    # plt.ylabel("Difference score")
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)


    # ## plot amplitude and rate Scores  together for each cell
    # fig = plt.figure()
    # #plt.set_tickslabel= cell_columns
    # plt.bar(np.arange(n_cells), diff_amp)
    # plt.bar(np.arange(n_cells), diff_Rate)
    # plt.legend(('Amp', 'Freq'))
    # plt.title("Amplitude & Frequency Difference (Open - Close)"+ str(animal_id))
    # plt.xlabel("Cell id")
    # plt.ylabel("Difference score")
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)

    ## plot amplitude*freq Scores
    fig = plt.figure()
    #plt.set_tickslabel= cell_columns
    plt.bar(np.arange(n_cells), diff_F_times_Rate)
    plt.title("Amplitude * Frequency Difference (Open - Close)"+ str(animal_id))
    plt.xlabel("Cell id")
    plt.ylabel("A*F Difference score")
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    ##  Calculate the proportion of cells being activated in open arms

    fig = plt.figure()
    plt.title("Proportion of cells responding in Open - closed \n "
              "relative to average amplitude \n"
              "Animal :" +str(animal_id))
    plt.pie([prop_of_cells_responding_in_open, prop_of_cells_non_responding_in_open], labels=["responding cells", "nonresponding"], autopct='%1.1f%%')
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    fig = plt.figure()
    plt.title("Proportion of cells responding in Center - closed \n "
              "relative to average amplitude \n"
              "Animal :" +str(animal_id))
    plt.pie([prop_of_cells_responding_in_center, prop_of_cells_non_responding_in_center], labels=["responding cells", "nonresponding"], autopct='%1.1f%%')
    # if save_figures_option == 'yes':
    pdf.savefig(fig)

    ###
    # calculate the difference in amp and freq score:
    # interesting_cells = abs(diff_amp - diff_Rate)
    # #interesting_cells.sort(), interesting_cells = interesting_cells[::-1]
    # plt.figure()
    # plt.bar(np.arange(n_cells), interesting_cells)


    #make a DataFrame with the average fluorescence per cell and the areas
    avgF = pd.DataFrame(avgF.ravel().reshape(-1, 1), columns=['F'])
    avgF['area'] = np.array(['1_open', '2_closed', '3_center'])[np.tile(np.arange(3), (n_cells, 1)).ravel()]
    avgF['cell'] = np.tile(np.arange(n_cells).reshape(-1, 1), (1, 3)).ravel()
    # make a DataFrame with the frequency rates per cell per area
    f_Rate = pd.DataFrame(avgRate.ravel().reshape(-1, 1), columns=['F'])
    f_Rate['area'] = np.array(['1_open', '2_closed', '3_center'])[np.tile(np.arange(3), (n_cells, 1)).ravel()]
    f_Rate['cell'] = np.tile(np.arange(n_cells).reshape(-1, 1), (1, 3)).ravel()
    # Make a Dataframe with the Ampl*Rate
    F_times_Rate = avgF ['F'] * f_Rate['F']
    F_times_Rate = pd.DataFrame(F_times_Rate.ravel().reshape(-1, 1), columns=['F'])
    F_times_Rate['area'] = np.array(['1_open', '2_closed', '3_center'])[np.tile(np.arange(3), (n_cells, 1)).ravel()]
    F_times_Rate['cell'] = np.tile(np.arange(n_cells).reshape(-1, 1), (1, 3)).ravel()


    ## Plot figures which depict the difference by cell and by area
    ##########################################################
    # # Mean frequency scatter plot
    # fig = plt.figure()
    # sns.swarmplot(x='area', y='F', data=f_Rate, dodge=False, color=None, size=5, edgecolor='None', linewidth=0); plt.ylabel("Mean Frequency"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # # Mean amplitude scatter plot
    # fig = plt.figure()
    # sns.swarmplot(x='area', y='F', data=avgF, dodge=False, color=None, size=5, edgecolor='None', linewidth=0); plt.ylabel("Mean amplitude"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # # Mean Frequency bar plot
    # fig = plt.figure()
    # sns.barplot(x='cell', y='F', hue='area', data=f_Rate, edgecolor='None', linewidth=0); plt.ylabel("Mean Frequency"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # # Mean amplitude bar plot
    # fig = plt.figure()
    # sns.barplot(x='cell', y='F', hue='area', data=avgF, edgecolor='None', linewidth=0); plt.ylabel("Mean amplitude"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # # Mean amplitude line plot
    # fig = plt.figure()
    # sns.lineplot(x='area', y='F', hue='cell', data=avgF, linewidth=1); plt.ylabel("Normalized average df/F peaks"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # # plotting the product of Amplitude and frequency
    # fig = plt.figure()
    # sns.swarmplot(x='area', y='F', data=F_times_Rate, dodge=False, color=None, size=5, edgecolor='None', linewidth=0); plt.ylabel("amplitude * frequency"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # fig = plt.figure()
    # sns.barplot(x='cell', y='F', hue='area', data=F_times_Rate, edgecolor='None', linewidth=0); plt.ylabel("amplitude * frequency"); plt.title("Animal :" + animal_id)
    # if save_figures_option == 'yes':
    #     pdf.savefig(fig)
    # #ax_2 = fig.add_subplot(1, 2, 2)



    # epm_position_resampled['arm_idx'] = np.nan
    # epm_position_resampled.loc[epm_position_resampled['open_arm'], 'arm_idx'] = 1
    # epm_position_resampled.loc[epm_position_resampled['closed_arm'], 'arm_idx'] = 2
    # epm_position_resampled.loc[epm_position_resampled['center_arm'], 'arm_idx'] = 3


    # ax_2.axhline(0, color='k')
    # ax_2.bar(1, avgF_open, facecolor='r')
    # ax_2.bar(2, avgF_closed, facecolor='k')
    # ax_2.bar(3, avgF_center, facecolor=[.7, .7, .7])


    #####################################################################
    ##   III Create heatmap
    ######################################################################
    n_bins =256#64
    ### Make the heatmap for the behavior:

    heatmap_behavior, xedges, yedges = np.histogram2d(epm_position_and_F_resampled['x'], epm_position_and_F_resampled['y'], bins=(n_bins, n_bins))
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    edges = [xedges, yedges]
    sample = epm_position_and_F_resampled[['x', 'y']].values
    bin_position = np.vstack([np.searchsorted(edges[i], sample[:, i], side='right') for i in range(2)]).transpose()
    bin_position[bin_position < 0] = 0
    bin_position[bin_position > n_bins - 1] = n_bins - 1

    print("printing all cells position / fluorescence heatmap")

    #cell_id = cell_columns[0]
    # For the normalized amplitude heatmap:
    for n in range(n_cells):
        cell_id = cell_columns [n]
        F = epm_position_and_F_resampled[cell_id].values
        heatmap_F = np.zeros_like(heatmap_behavior, dtype=float)

        for i_sample in range(F.shape[0]):
            heatmap_F[bin_position[i_sample, 0], bin_position[i_sample, 1]] = F[i_sample]

        fig = plt.figure()
        plt.suptitle("Animal: " + animal_id)
        plt.subplot(1,2,1); plt.imshow(convolve(heatmap_behavior, Gaussian2DKernel(stddev=1))); plt.title('EPM Position'); plt.clim(-1, 1)
        plt.xlabel("X coordinate")
        plt.ylabel("Y coordinate")
        #plt.text (x=o)
        plt.subplot(1,2,2); plt.imshow(convolve(heatmap_F, Gaussian2DKernel(stddev=1))); plt.title('Normalized dF/F')#; plt.clim(0, .4)
        plt.title(cell_id)
        plt.xlabel("X coordinate")
        plt.ylabel("Y coordinate")
        #plt.colorbar()
        # if save_figures_option == 'yes':
        pdf.savefig(fig)
    #
    # if save_figures_option == 'yes':
    pdf.close()

    # if save_figures_option != 'yes':
    # plt.show()


if __name__ == '__main__':
    # Get user inputs
    if len(sys.argv) > 1:
        date_exp = sys.argv[1]
        animal_id = sys.argv[2]
        dFF_trial = int(sys.argv[3])
    else:
        date_exp = '190917'
        animal_id = 'MA_30epi'
        dFF_trial = 191
    # run
    epm_track(date_exp, animal_id, dFF_trial)
