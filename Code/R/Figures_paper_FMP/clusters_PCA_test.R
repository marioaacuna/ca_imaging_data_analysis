#!/usr/bin/env Rscript
rm(list = ls())
cat("\014")
set.seed(17)

## dependencies
# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("ggplot2", "ggbeeswarm", "tidyr", "optparse", "grid", "tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://stat.ethz.ch/CRAN/")
# Installs pacman ("package manager") if needed
if (!require("pacman")) install.packages("pacman")

# Use pacman to load add-on packages as desired
pacman::p_load(pacman, tidyverse) 
library(cluster)
library(factoextra)

file_input = "/Users/marioacuna/Desktop/trial_CCI.csv"
file_input_sham = "/Users/marioacuna/Desktop/trial_sham.csv"
# file_output = "D:\\_MATLAB_CaImaging\\test.pdf"
# data_type = "amplitude"
# data_limits = "0, 100.0"
# verbose = TRUE
data <- read.table(file_input, header = F, sep = ",")
#head(data) 
data_sham <- read.table(file_input_sham, header = F, sep = ",")
#head(data_sham)

# Make dataframes
col_names = c("Painful", "aversive", "innocuous")
df_s = as.data.frame(data_sham)
colnames(df_s) <- col_names
df_c = as.data.frame(data)
colnames(df_c) <- col_names

# Do hierarchical clustering for sham and CCI
hc <- t(data)   %>%  # Get data
  dist   %>%  # Compute distance/dissimilarity matrix
  hclust      # Computer hierarchical clusters
hc_s <- t(data_sham)   %>%  
  dist   %>%  # Compute distance/dissimilarity matrix
  hclust # Computer hierarchical clusters
par(mfrow=c(1,2)) # Select 1 row and 2 columns to plot
# Plot hierarchical clustering
plot(hc, main = "Dendrogram CCI")
plot(hc_s,main = "Dendrogram sham")

# Cluster data into k clusters (clustering of large applications)
k = length(col_names)
cl_s = clara(df_s, k, samples = 100)
cl_cci = clara(df_c, k, samples= 100)
par(mfrow=c(1,2))
plot(cl_s, main = "Sham")
plot(cl_cci,main = "CCI")

# cl_agnes_cci = cluster::agnes(data)
# plot(cl_agnes_cci)

# Do PCA on the data
pc_cci = prcomp(df_c, center = T, scale = T)
pc_sham = prcomp(df_s, center = T, scale = T)
#  Plot principal componets
biplot(pc_cci, main = "PCA CCI", xlim = c(-.3, .2),  # CCI
       ylim = c(-0.2, 0.2), 
       xlabs = rep("●", nrow(df_c)))

biplot(pc_sham,main = "PCA sham", xlim = c(-.3, .2), # SHAM
       ylim = c(-0.2, 0.2),
       xlabs = rep("●", nrow(df_s)))

# PAM clustering
# CCI
# n_k  =fviz_nbclust(data, pam, method = "silhouette", k.max = 10) # Function to determine size of clusters
pam_c =pam(data, 3, metric = "euclidean")
fviz_cluster(pam_c, ellipse.type = "norm",  main = "CCI clustering",
             xlim = c(-8, 3),
             ylim = c(-5, 3),
             show.clust.cent = T,
             labelsize = 0, 
             pointsize = 1,
             ggtheme = theme_classic()
             )
# Sham
# fviz_nbclust(data_sham, pam, method = "silhouette", k.max = 10) # Function to determine size of clusters
pam_s = pam(data_sham, 3, metric = "euclidean")
fviz_cluster(pam_s, ellipse.type = "norm", main = "Sham clustering", 
             xlim = c(-8, 3),
             ylim = c(-5, 3),
             show.clust.cent = T,
             labelsize = 0, 
             pointsize = 1,
             ggtheme = theme_classic()
            )

# CLEAN UP #################################################

# Clear environment
rm(list = ls())

# Clear packages
p_unload(all)  # Remove all add-ons


# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)

