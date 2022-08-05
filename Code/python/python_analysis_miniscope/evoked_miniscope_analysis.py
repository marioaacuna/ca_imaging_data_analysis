import scipy
import numpy as np
import pandas as pd
from scipy.signal import resample
import scipy.io
import matplotlib.pyplot as plt
from scipy.signal import resample
from scipy.signal import find_peaks
from scipy import stats
import seaborn as sns
import pylab
import os

# plt.clf()
# plt.close()
plotting = True
## GENERAL STATS OF EVOKED ACTIVITY
session = '2'# '1', '2'
miniscope_path = r'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results'
miniscope_file = os.path.join(miniscope_path, str('all_results' + '_ses_'+ str(session)+'.mat'))
f = scipy.io.loadmat(miniscope_file)
data_cci = f.get('all_data_cci')
data_sham = f.get('all_data_sham')
n_stimuli = np.size(data_cci,1)
# p_values= np.zeros((1,n_stimuli))
# p_values[0,0] = f.get('p_cold')
# p_values[0,1] = f.get('p_heat')
# p_values [0,2]  =  f.get('p_pin')
# p_values [0,3] = f.get('p_touch')
p_cold = f.get('p_cold')
p_heat = f.get('p_heat')
p_pin = f.get('p_pin')
p_touch = f.get('p_touch')
prop_responders_cci = f.get('prop_resp_cci')
prop_responders_sham = f.get('prop_resp_sham')
df_data_cci = pd.DataFrame()
df_data_sham = pd.DataFrame()
df_data_cci['cold'] = (data_cci[:,0])
df_data_cci['heat'] = (data_cci[:,1])
df_data_cci['pinprick'] = (data_cci[:,2])
df_data_cci['touch'] = (data_cci[:,3])

df_data_sham['cold'] = (data_sham[:,0])
df_data_sham['heat'] = (data_sham[:,1])
df_data_sham['pinprick'] = (data_sham[:,2])
df_data_sham['touch'] = (data_sham[:,3])

n_responding_cci = df_data_cci[df_data_cci.notnull()].count()
n_responding_sham = df_data_sham[df_data_sham.notnull()].count()
print('p_cold ', p_cold)
print('p_heat ', p_heat)
print('p_prinprick ', p_pin)
print('p_touch ', p_touch)
if p_cold > 0.05:
    p_cold = 'ns'
else:
    p_cold = '*'

if p_heat > 0.05:
    p_heat = 'ns'
else:
    p_heat = '*'

if p_pin > 0.05:
    p_pin = 'ns'
else:
    p_pin = '*'
if p_touch > 0.05:
    p_touch = 'ns'
else:
    p_touch = '*'








if plotting:
    axes = ['CCI', 'Sham']
    means_cold = [np.nanmean(df_data_cci['cold']), np.nanmean(df_data_sham['cold'])]
    errors_cold = [np.nanstd(df_data_cci['cold'])/np.sqrt(n_responding_cci['cold']), np.nanstd(df_data_sham['cold'])/np.sqrt(n_responding_sham['cold'])]
    means_heat = [np.nanmean(df_data_cci['heat']), np.nanmean(df_data_sham['heat'])]
    errors_heat = [np.nanstd(df_data_cci['heat'])/np.sqrt(n_responding_cci['heat']), np.nanstd(df_data_sham['heat'])/np.sqrt(n_responding_sham['heat'])]

    means_pin = [np.nanmean(df_data_cci['pinprick']), np.nanmean(df_data_sham['pinprick'])]
    errors_pin = [np.nanstd(df_data_cci['pinprick'])/np.sqrt(n_responding_cci['pinprick']), np.nanstd(df_data_sham['pinprick'])/np.sqrt(n_responding_sham['pinprick'])]


    means_touch = [np.nanmean(df_data_cci['touch']), np.nanmean(df_data_sham['touch'])]
    errors_touch = [np.nanstd(df_data_cci['touch'])/np.sqrt(n_responding_cci['touch']), np.nanstd(df_data_sham['touch'])/np.sqrt(n_responding_sham['touch'])]
    y_label= 'Amplitude (dF/F)' # 'Amplitude (dF/F)', 'z_score'

   # plotting bars
    barlist = plt.bar (axes, means_cold, yerr=errors_cold);plt.title(str('cold ' + 'Session ' + session))
    barlist[0].set_color('r');plt.text(-0.2, means_cold[0], p_cold);plt.ylabel(y_label);plt.show()# ***
    barlist = plt.bar (axes, means_heat, yerr=errors_heat);plt.title(str('heat ' + 'Session ' + session));plt.text(-0.2, means_heat[0], p_heat)
    barlist[0].set_color('r')
    plt.ylabel(y_label);plt.show() #ns
    barlist = plt.bar (axes, means_pin, yerr=errors_pin)
    barlist[0].set_color('r')
    plt.title(str('pinprick ' + 'Session ' + session));plt.text(-0.2, means_pin[0],p_pin);plt.ylabel(y_label);plt.show() #***
    barlist = plt.bar (axes, means_touch, yerr=errors_touch)
    barlist[0].set_color('r')
    plt.title(str('touch ' + 'Session ' + session));plt.text(-0.2, means_touch[0], p_touch);plt.ylabel(y_label);plt.show() #***

    ## plotting boxplot
    data_cold_sham = df_data_sham['cold'][~np.isnan(df_data_sham['cold'])]
    data_heat_sham = df_data_sham['heat'][~np.isnan(df_data_sham['heat'])]
    data_pin_sham = df_data_sham['pinprick'][~np.isnan(df_data_sham['pinprick'])]
    data_touch_sham = df_data_sham['touch'][~np.isnan(df_data_sham['touch'])]

    data_cold_cci = df_data_cci['cold'][~np.isnan(df_data_cci['cold'])]
    data_heat_cci = df_data_cci['heat'][~np.isnan(df_data_cci['heat'])]
    data_pin_cci = df_data_cci['pinprick'][~np.isnan(df_data_cci['pinprick'])]
    data_touch_cci= df_data_cci['touch'][~np.isnan(df_data_cci['touch'])]

    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.boxplot([data_cold_cci, data_cold_sham])
    # ax.set_xticks(range(1,len(axes)+1))
    # ax.set_xticklabels(axes)
    # ax.set_ylim(0,5) ; ax.set_ylabel(y_label) ; ax.set_title((str('cold ' + 'Session ' + session)))

    ## Pie Charts
    non_resp_cci = 1- sum(prop_responders_cci[0],0)
    non_resp_sham = 1- sum(prop_responders_sham[0],0)
    plt.pie(np.concatenate((prop_responders_cci[0],np.array([non_resp_cci]))), autopct='%1.1f%%',labels=['cold', 'heat', 'pin', 'touch', 'non-resp'],  colors= ['b','r', 'g', 'c', 'y'], shadow=True)#autopct='%1.1f%%'
    plt.title(str('Sham ' + 'Session ' + session));plt.show()
    plt.pie(np.concatenate((prop_responders_sham[0],np.array([non_resp_sham]))), autopct='%1.1f%%', labels=['cold', 'heat', 'pin', 'touch', 'non-resp'],  colors= ['b','r', 'g', 'c', 'y'], shadow=True)#autopct='%1.1f%%'
    plt.title(str('Sham ' + 'Session ' + session));plt.show()

#################################################################################
## SPONTANEOUS ACTIVITY OF RESPONDING CELLS FOR SESSION 1
###############################################################################
animal_ID = [ 'I37_2', 'I38_3', 'I38_4', 'I39_1']#,'I32_1', 'I32_3',
is_CCI = [ False, False, True, True] #False, True,
n_animals_cci = 2
n_animals_sham = 2
for animal in range(len(animal_ID)):
    if is_CCI[animal]:
        traces_path = 'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\%s\jointExtraction\sorted' % animal_ID [animal]
        traces_file = os.path.join (traces_path, 'traces_SP_1.mat' )
        f_traces = scipy.io.loadmat(traces_file)
        df_f_data_cci = f_traces.get("traces_SP_1")
        dF_F = pd.DataFrame()
        for i in range(len(df_f_data_cci[:, 0])):
            dF_F["cell_%i" % (i + 1)] = df_f_data_cci[i, :]
        # Load Selectivity variable
        select_path= r'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results'
        select_file = os.path.join(select_path, str('results_analyis_'+ str(animal_ID [animal])+'_ses_'+session+ '.mat')) #%  animal_ID [animal],% session
        f_select = scipy.io.loadmat(select_file)
        selective_cells =f_select.get("selective_cells")
        selective_cells = np.logical_or(selective_cells,0)
        # parameters
        n_stimuli =len(selective_cells[0])
        n_cells = len(selective_cells)
        threshold = dF_F.median() + 3 * dF_F.mad()
        min_distance = 50
        mean_amplitude_per_animal = np.zeros((2, n_stimuli))
        stats_amplitude_per_animal = np.zeros((1, n_stimuli))
        # get the traces of selective cells per stim
        for i_stim in range (len(selective_cells[0])):
            idx_select = np.where(selective_cells[:, i_stim]== True)
            # print('animal: ', animal_ID[animal], ', stim: ', i_stim, ', selective cells: ', str(idx_select))
            idx_nonselect = np.where(selective_cells[:, i_stim] == False)
            selective_traces = np.zeros((len(df_f_data_cci[0]), len(idx_select[0])))
            n_selec_cells = len(idx_select[0])
            nonselect_traces = np.zeros((len(df_f_data_cci[0]), len(idx_nonselect[0])))
            n_nonselect_cells = len(idx_nonselect[0])
            for i_cell in range(n_selec_cells):
                selective_traces[:,i_cell] = df_f_data_cci[[idx_select[0][i_cell]]]
            for i_cell in range(n_nonselect_cells):
                nonselect_traces [:,i_cell] =  df_f_data_cci[[idx_nonselect[0][i_cell]]]
            # Amplitude  selective cells this stim
            amp_peak_sel_cells = np.zeros((n_selec_cells, 1), dtype=float)
            for i_cell in range(n_selec_cells):
                df_select_cells = pd.DataFrame()
                for i in range(n_selec_cells):
                    df_select_cells["cell_%i" % (i + 1)] = selective_traces[:,i]
                cell_columns = df_select_cells.columns
                F_this_cell = df_select_cells[cell_columns[i_cell]]
                peaks_this_cell, _ = find_peaks(F_this_cell,
                                                height=threshold[i_cell],
                                                distance=min_distance)
                amp_peak_sel_cells[i_cell] = np.mean(
                        F_this_cell.loc[peaks_this_cell])

            # Amplitude NON SELECTIVE cells this stim
            amp_peak_nonsel_cells = np.zeros((n_nonselect_cells, 1), dtype=float)
            for i_cell in range(n_nonselect_cells):
                df_nonselect_cells = pd.DataFrame()
                for i in range(n_nonselect_cells):
                    df_nonselect_cells["cell_%i" % (i + 1)] = nonselect_traces[:,i]
                cell_columns = df_nonselect_cells.columns
                F_this_cell = df_nonselect_cells[cell_columns[i_cell]]
                peaks_this_cell, _ = find_peaks(F_this_cell,
                                                height=threshold[i_cell],
                                                distance=min_distance)
                amp_peak_nonsel_cells[i_cell] = np.mean(
                        F_this_cell.loc[peaks_this_cell])


            mean_amplitude_per_animal[0, i_stim] = np.nanmean (amp_peak_sel_cells)
            mean_amplitude_per_animal[1, i_stim] = np.nanmean(amp_peak_nonsel_cells)
            [_ , stats_amplitude_per_animal[0,i_stim]] = stats.ttest_ind(amp_peak_sel_cells,amp_peak_nonsel_cells, nan_policy='omit')

        print('animal: ', str(animal_ID[animal]), mean_amplitude_per_animal)
        print('p_value: ', str(animal_ID[animal]), stats_amplitude_per_animal)

    else:
        miniscope_path = 'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\%s\jointExtraction\sorted' % animal_ID [animal]
        miniscope_file = os.path.join (miniscope_path, 'traces_SP_1.mat' )
        f = scipy.io.loadmat(miniscope_file)
        df_f_data_sham = f.get("traces_SP_1")
        dF_F = pd.DataFrame()
        for i in range(len(df_f_data_sham[:, 0])):
            dF_F["cell_%i" % (i + 1)] = df_f_data_sham[i, :]
        # Load Selectivity variable
        select_path= r'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results'
        select_file = os.path.join(select_path, str('results_analyis_'+ str(animal_ID [animal])+'_ses_'+session+ '.mat'))#%  animal_ID [animal]
        f_select = scipy.io.loadmat(select_file)
        selective_cells =f_select.get("selective_cells")
        selective_cells = np.logical_or(selective_cells,0)
        # parameters
        n_stimuli =len(selective_cells[0])
        n_cells = len(selective_cells)
        threshold = dF_F.median() + 3 * dF_F.mad()
        min_distance = 50
        mean_amplitude_per_animal = np.zeros((2, n_stimuli))
        stats_amplitude_per_animal = np.zeros((1, n_stimuli))
        # get the traces of selective cells per stim
        for i_stim in range (len(selective_cells[0])):
            idx_select = np.where(selective_cells[:, i_stim]== True)
            idx_nonselect = np.where(selective_cells[:, i_stim] == False)
            # print('animal: ', animal_ID[animal], ', stim: ', i_stim,
            #       ', selective cells: ', str(idx_select))
            selective_traces = np.zeros((len(df_f_data_sham[0]), len(idx_select[0])))
            n_selec_cells = len(idx_select[0])
            nonselect_traces = np.zeros((len(df_f_data_sham[0]), len(idx_nonselect[0])))
            n_nonselect_cells = len(idx_nonselect[0])
            for i_cell in range(n_selec_cells):
                selective_traces[:,i_cell] = df_f_data_sham[[idx_select[0][i_cell]]]
            for i_cell in range(n_nonselect_cells):
                nonselect_traces [:,i_cell] =  df_f_data_sham[[idx_nonselect[0][i_cell]]]
            # Amplitude  selective cells this stim
            amp_peak_sel_cells = np.zeros((n_selec_cells, 1), dtype=float)
            for i_cell in range(n_selec_cells):
                df_select_cells = pd.DataFrame()
                for i in range(n_selec_cells):
                    df_select_cells["cell_%i" % (i + 1)] = selective_traces[:,i]
                cell_columns = df_select_cells.columns
                F_this_cell = df_select_cells[cell_columns[i_cell]]
                peaks_this_cell, _ = find_peaks(F_this_cell,
                                                height=threshold[i_cell],
                                                distance=min_distance)
                amp_peak_sel_cells[i_cell] = np.mean(
                        F_this_cell.loc[peaks_this_cell])

            # Amplitude NON SELECTIVE cells this stim
            amp_peak_nonsel_cells = np.zeros((n_nonselect_cells, 1), dtype=float)
            for i_cell in range(n_nonselect_cells):
                df_nonselect_cells = pd.DataFrame()
                for i in range(n_nonselect_cells):
                    df_nonselect_cells["cell_%i" % (i + 1)] = nonselect_traces[:,i]
                cell_columns = df_nonselect_cells.columns
                F_this_cell = df_nonselect_cells[cell_columns[i_cell]]
                peaks_this_cell, _ = find_peaks(F_this_cell,
                                                height=threshold[i_cell],
                                                distance=min_distance)
                amp_peak_nonsel_cells[i_cell] = np.mean(
                        F_this_cell.loc[peaks_this_cell])


            mean_amplitude_per_animal[0, i_stim] = np.nanmean (amp_peak_sel_cells)
            mean_amplitude_per_animal[1, i_stim] = np.nanmean(amp_peak_nonsel_cells)
            [_ , stats_amplitude_per_animal[0,i_stim]] = stats.ttest_ind(amp_peak_sel_cells,amp_peak_nonsel_cells, nan_policy='omit')
        print('animal: ', str(animal_ID[animal]), mean_amplitude_per_animal)
        print('p_value: ', str(animal_ID[animal]), stats_amplitude_per_animal)



################################################################################
# ## checking selective cells across sessions
# for animal in range (n_animals):
#     sessions = ['1', '2']
#     n_sess = len(sessions)
#
#     for i_sess in range (n_sessions):
#         select_path= r'T:\Mario\miniscope data\Pain_behavior_miniscope\Analysis_full\all_results'
#         select_file = os.path.join(select_path, str('results_analyis_'+ str(animal_ID [animal])+'_ses_'+sessions[i_sess]+ '.mat')) #%  animal_ID [animal],% session)
#         f_select = scipy.io.loadmat(select_file)
#         selective_cells = f_select.get("selective_cells")
#         selective_cells = np.logical_or(selective_cells,0)

