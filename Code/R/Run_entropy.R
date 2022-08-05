#!/usr/bin/env Rscript
rm(list = ls())
cat("\014")
set.seed(17)

## dependencies
# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("optparse", "entropy", "MASS", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://stat.ethz.ch/CRAN/")
library("entropy")
suppressMessages(require(MASS))
suppressMessages(require(optparse))
############################################
# Unpack inputs
option_list = list(
  make_option("--file_input", action = "store", type = "character", default = NULL),
  make_option("--file_output", action = "store", type = "character", default = NULL),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE)
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

file_input = opt$file_input
file_output = opt$file_output
verbose = opt$verbose

##############
# file_input = "C:/Users/acuna/Desktop/Entropy_analysis/R_entropy/data_entropy.csv"
# file_output = "C:/Users/acuna/Desktop/Entropy_analysis/R_entropy/result_entropy.csv"
# verbose = FALSE
#############


DATA = read.table(file_input, sep = ",", header = FALSE)
data = as.matrix(DATA)
disc_data = discretize(data, numBins = 10)
# View(disc_data)
# use for the discretized version
data = disc_data 
## Get the entropy values
# E_empirical = entropy.empirical(data)
# E_Dirichlet = entropy.Dirichlet(data,a=sqrt(sum(data))/length(data))
# E_MillerMadow = entropy.MillerMadow(data)
# E_ChaoShen = entropy.ChaoShen(data)
# E_plugin= entropy.plugin(data)
# E_shrink=entropy.shrink(data)


# Get the mutual information values
# mi_empirical = mi.empirical(data)
# mi_Dirichlet = mi.Dirichlet(data, a=sqrt(sum(data))/length(data))

## entropy values after using the frequency of each value
# freqs_data = freqs.shrink(data) 
# E_Shrink = entropy(data, freqs_data, method = "shrink", unit=c("log2"),verbose = FALSE)
E_Shrink = entropy(data, method = "shrink", unit=c("log2"),verbose = verbose)
# E_MillerMadow = entropy.MillerMadow(data , unit = "log2")


###################
write.matrix(E_Shrink[1], file_output)
# write.matrix(E_MillerMadow, file_output)





# #############################################
# x1 = runif(10000)
# y = c(4,3,2,5,8,4,2,0,0,0,2,4,1)
# y_2 = c(2,1,5,4,6,3,8,0,0,0,1,2,7)
# y_3 = c(0,0,0,2,5,2,0,6,0,1,0,4,0)
# hist(x1, xlim=c(0,1), freq = FALSE)
# 
# y1 = discretize(x1, numBins = 100, r= c(0,1))
# y.mat = matrix(y2d, ncol = 3)
# entropy.empirical(y.mat) # it can be done for a 2D matrix
# entropy.Dirichlet(y.mat, a = 0)
# mi.Dirichlet(y.mat, a=1)
# 
# entropy(y, method= "MM", unit="log2") # entropy in bits
# entropy.ChaoShen(y)
# entropy.empirical(y)
# 
# ## Dirichet
# freqs.Dirichlet(y, a=1/2)
# 
# entropy.Dirichlet(y, a=1)
# entropy.Dirichlet(y, a=sqrt(sum(y))/length(y))
# 
# # Bayesian estimate of Kullback-Leibler divergence (a=1/6)
# KL.Dirichlet(y1, y_2, y_3,a1=1/6, a2=1/6, a3=1/6)
# 
# ###
# # two independent random variables
# x1 = runif(10000)
# x2 = runif(10000)
# y2d = discretize2d(x1, x2, numBins1 = 10, numBins2 = 10)
# View(y2d)
# sum(y2d)
# # joint entropy
# H12 = entropy(y2d )
# H12
# log(100)
# 
# y.mat = matrix(y2d, ncol = 3)
# mi = mi.empirical(y.mat)
# mi
# 
# # empirical chi-squared divergence of independence
# cs.indep = chi2indep.empirical(y.mat)
# cs.indep
