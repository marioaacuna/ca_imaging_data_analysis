import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
# animal_id = '26_1'#, 26_2, 26_3, 27_2, 27_3
# date_exp = '190130'#, 190131
# path_one_animal = 'C:/Users/acuna/Desktop/fun_python/epm_data/epm_data_' + str(date_exp)+ '_'+str(animal_id)+'_metadata.csv'
path = 'C:/Users/acuna/Desktop/fun_python/epm_data/' # epm_data_' + str(date_exp)+ '_'+str(animal_id)+'_metadata.csv'
list_of_files = os.listdir(path)
all_metadata = [os.path.join(path, i) for i in list_of_files if i.endswith("_metadata.csv")]
n_animals = len(all_metadata)

save_figures_option = input("do you want to save figures in pdf?")
if save_figures_option == 'yes':
    save_pdf_path = (
                "C:/Users/acuna/Desktop/fun_python/epm_data/figures/_global_data"+".pdf")
    pdf = matplotlib.backends.backend_pdf.PdfPages (save_pdf_path)



#i_animal = 0
n_cells = np.zeros(n_animals, dtype=float)
rel_time_open = np.zeros(n_animals, dtype=float)
rel_time_closed = np.zeros(n_animals, dtype=float)
rel_time_center = np.zeros(n_animals, dtype=float)
vel_open = np.zeros(n_animals, dtype=float)
vel_closed = np.zeros(n_animals, dtype=float)
vel_center = np.zeros(n_animals, dtype=float)
prop_cells_in_open = np.zeros(n_animals, dtype=float)
prop_cells_in_center = np.zeros(n_animals, dtype=float)
prop_cells_overlap_op_cen = np.zeros(n_animals, dtype=float)
for i_animal in range(n_animals):
    a = pd.read_csv(all_metadata[i_animal], index_col=0)
    n_cells_this_animal = np.array(a['total cells'])
    n_cells [i_animal] = n_cells_this_animal [0]
    rel_time_in_arms_this_animal = np.array(a['Rel time spent in arms'])
    rel_time_open [i_animal] = rel_time_in_arms_this_animal[0]
    rel_time_closed [i_animal] = rel_time_in_arms_this_animal[1]
    rel_time_center [i_animal] = rel_time_in_arms_this_animal[2]

    vel_open_this_animal =  np.array(a['mean vel open'])
    vel_closed_this_animal = np.array(a['mean vel close'])
    vel_center_this_animal = np.array(a['mean vel center'])
    vel_open [i_animal] = vel_open_this_animal[0]
    vel_closed [i_animal] = vel_closed_this_animal[0]
    vel_center [i_animal] = vel_center_this_animal[0]

    prop_cells_in_open_this_animal = np.array(a["prop of cells responding in open"])
    prop_cells_in_center_this_animal = np.array(a["prop of cells responding in center"])
    prop_cells_in_open [i_animal] = prop_cells_in_open_this_animal[0]
    prop_cells_in_center [i_animal] = prop_cells_in_center_this_animal[0]

    prop_cell_overlap_this_animal = np.array(a["prop cells overlap open and center"])
    prop_cells_overlap_op_cen[i_animal] = prop_cell_overlap_this_animal[0]
total_cells = int(np.sum(n_cells))
## Calculate the means
# Time spent in each compartment:
mean_rel_time_open = np.mean(rel_time_open)
mean_rel_time_closed = np.mean(rel_time_closed)
mean_rel_time_center = np.mean(rel_time_center)

# Velocities :
mean_vel_open = np.mean(vel_open)
mean_vel_closed = np.mean(vel_closed)
mean_vel_center =np.mean (vel_center)

# Proportion of cells in each compartment and overlap
mean_prop_cells_open = np.mean(prop_cells_in_open)
mean_prop_cells_center = np.mean(prop_cells_in_center)
mean_overlap_cells = np.mean(prop_cells_overlap_op_cen)
## Plotting:
n_arms = 3
arms = ["Open", "Closed", "Center"]

# Time spent in each arm:
fig = plt.figure ()
plt.pie ([mean_rel_time_open, mean_rel_time_closed,mean_rel_time_center], labels=["Open", "closed", "center"], autopct='%1.1f%%')
plt.title("Proportion of time spent in each compartment \n  "
          + "n = " + str(n_animals) + " mice")
if save_figures_option == 'yes':
    pdf.savefig(fig)


# Mean velocity in each arm
means_v = (mean_vel_open, mean_vel_closed, mean_vel_center)
err_v = (np.std(vel_open)/np.sqrt(n_animals), np.std(vel_closed)/np.sqrt(n_animals), np.std(vel_center)/np.sqrt(n_animals))

fig = plt.figure()
plt.bar(np.arange(n_arms),means_v, yerr=err_v)
plt.title("Velocity in arms\n N = "+ str (n_animals) + " mice")
plt.ylabel("Velocity (cm/s)")
plt.xticks(np.arange(n_arms), arms)
if save_figures_option == 'yes':
    pdf.savefig(fig)

# Nr. of cells with high activity in open and center arms and overlap
fig_mean_overlap = (mean_overlap_cells,mean_overlap_cells)
prop_cells_in_arms = (mean_prop_cells_open-mean_overlap_cells, mean_prop_cells_center-mean_overlap_cells)
fig_err_overlap = ( np.std(prop_cells_overlap_op_cen)/np.sqrt(n_animals),  np.std(prop_cells_overlap_op_cen)/np.sqrt(n_animals))
err_prop = (np.std(prop_cells_in_open)/np.sqrt(n_animals), np.std(prop_cells_in_center)/np.sqrt(n_animals))

fig = plt.figure()
width = 0.8
p1 = plt.bar(np.arange(n_arms-1), fig_mean_overlap, width, yerr = fig_err_overlap)
p2 =plt.bar(np.arange(n_arms-1), prop_cells_in_arms,width, yerr= err_prop, bottom = fig_mean_overlap)
plt.xticks(np.arange(n_arms-1), (arms[0], arms[2]))
plt.legend((p1[0], p2 [0]), ("Overlap", "Specific"))
plt.title("Proportion of neurons activated in open arm and center\n"
          "N = " + str (n_animals) + " mice\n"
                                     "n = " + str(total_cells) + " cells")
plt.ylabel("Proportion of cells")
if save_figures_option == 'yes':
    pdf.savefig(fig)

# prop cells with overlap activity in open and center
plt.show()
