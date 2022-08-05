#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("cpm", "MASS")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load libraries
require(cpm)
require(MASS)

# Unpack inputs
# file_input = 'C:\\Users\\deluna\\AppData\\Local\\Temp\\cpm.csv'
file_input = args[1]

cpmType = 'Kolmogorov-Smirnov'
ARL0 = 370
startup = 1

# Load trace
print(paste("Reading data from: ",file_input, sep=""))
data <- read.table(file_input, sep =",", header=FALSE);
y = as.matrix(data)

# Compute change-points
cpm_results1 = processStream(y, cpmType=cpmType, ARL0=ARL0, startup=startup)
change_points1 = as.vector(cpm_results1[['changePoints']])

change_points = change_points1 + 1  # Timestamps now mark beginning of new section

if (length(change_points1) > 0) {
  write.matrix(change_points, file_input)
} else {
  writeLines('None', file_input)
}



# The code below has been commented because it refers to the case when trials can be concatenated in time. However, in most cases, it cannot be done because there is a slight (and sometimes large) delay between them.

# if (length(change_points1) > 0) {
#   # Take first change point and make it time 0 of the sequence
#   central_change_point = cpm_results1[["changePoints"]][1]
#   z = as.vector(y)
#   z = z[(central_change_point):(n_elements-(n_rows-central_change_point)-1)]
#   z = matrix(z, nrow = n_rows)
#   z = rowMeans(z)
#   cpm_results2 = processStream(z, cpmType=cpmType, ARL0=ARL0, startup=startup)
#   change_points2 = as.vector(cpm_results2[['changePoints']])
#
#   # Concatenate results
#   change_points = change_points1 + 1  # Timestamps now mark beginning of new section
#
#   # Re-align timestamps of second
#   if (length(change_points2) > 0) {
#     # Add offset
#     change_points2 = change_points2 + central_change_point
#     # Timestamps now mark beginning of new section
#     change_points2 = change_points2 + 1
#     # Wrap around max number of rows
#     change_points2[change_points2 > n_rows] = change_points2[change_points2 > n_rows] - n_rows
#
#     # Add to final output
#     change_points = sort(unique(c(change_points, change_points2)))
#   }
#   write.matrix(change_points, file_input)
#
# } else {
#   writeLines('None', file_input)
# }

print("done!")
