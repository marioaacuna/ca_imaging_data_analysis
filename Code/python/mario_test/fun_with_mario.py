## Tutorials downloaded from
## https://matplotlib.org/tutorials/introductory/pyplot.html#sphx-glr-tutorials-introductory-pyplot-py

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# y_1  = (1, 2, 4, 16, 32)
# x_1 = (0, 1, 2, 3 , 4)
#
# plt.figure(1, figsize=(9, 3))
# plt.subplot(131)
# plt.plot(x_1,y_1)
# plt.ylabel("Some numbers")
# plt.xlabel("samples")
# plt.subplot(132)
# plt.bar(x_1,y_1)
# plt.ylabel("Some numbers")
# plt.xlabel("samples")
# plt.subplot(133)
# plt.scatter(x_1,y_1)
# plt.ylabel("Some numbers")
# plt.xlabel("samples")
#
#
# plt.show()

# N = 5
# menMeans = (20, 35, 30, 35, 27)
# womenMeans = (25, 32, 34, 20, 25)
# menStd = (2, 3, 4, 1, 2)
# womenStd = (3, 5, 2, 3, 3)
# ind = np.arange(len(menMeans))    # the x locations for the groups
# width = 0.35       # the width of the bars: can also be len(x) sequence
#
# p1 = plt.bar(ind, menMeans, width, yerr=menStd)
# p2 = plt.bar(ind, womenMeans, width,
#              bottom=menMeans, yerr=womenStd)
#
# plt.ylabel('Scores')
# plt.title('Scores by group and gender')
# plt.xticks(ind, ('G1', 'G2', 'G3', 'G4', 'G5'))
# plt.yticks(np.arange(0, 81, 10))
# plt.legend((p1[0], p2[0]), ('Men', 'Women'))
#
# plt.show()
#
# import numpy as np
# import matplotlib.pyplot as plt
#
#
# # Fixing random state for reproducibility
# np.random.seed(19680801)
#
# # Compute pie slices
# N = 20
# theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
# radii = 10 * np.random.rand(N)
# width = np.pi / 4 * np.random.rand(N)
#
# ax = plt.subplot(111, projection='polar')
# bars = ax.bar(theta, radii, width=width, bottom=0.0)
#
# # Use custom colors and opacity
# for r, bar in zip(radii, bars):
#     bar.set_facecolor(plt.cm.viridis(r / 10.))
#     bar.set_alpha(0.5)
#
# plt.show()

#


## -- - importing from CSV file and plotting bar graphs----

males=pd.read_csv ("test.csv",usecols=[0])
females =pd.read_csv ("test.csv",usecols=[1])
#print(file)
#pd.read_excel("test_1.xls")  # reading file
#print(males)
#print(females)
mean_males =np.mean(males)
std_men = np.std(males)/np.sqrt(len(males))
mean_females =np.mean(females)
std_females = np.std(females)/np.sqrt(len(females))
from scipy import stats
p_value=stats.ttest_ind(males, females)
p_value= p_value[1]
print(p_value)
print(mean_males[0])
print( mean_females[0])
data = (mean_males[0], mean_females[0])
ind= np.arange(2)
print(data, std_men[0], std_females[0])
#print(ind)
#data=str(data)
plt.annotate("p_value = "+ str(p_value), xy=(0.75, 25))
p1 =plt.bar(0,mean_males, yerr=std_men[0])
p2 =plt.bar(1, mean_females, yerr=std_females[0])

plt.xticks(ind, ("Males", "Females"))
#plt.bar(ind,data)
plt.show()#file.head()


