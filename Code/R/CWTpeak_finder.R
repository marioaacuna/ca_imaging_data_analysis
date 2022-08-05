#!/usr/bin/env Rscript
rm(list = ls())

# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("wmtsa", "MASS", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://stat.ethz.ch/CRAN/")

# Load libraries
suppressMessages(require(wmtsa))
suppressMessages(require(MASS))
suppressMessages(require(optparse))

# Unpack inputs
option_list = list(
  make_option("--file_input", action = "store", type = "character", default = NULL),
  make_option("--file_output", action = "store", type = "character", default = NULL),
  make_option("--frame_rate", action = "store", type = "double", default = NULL),
  make_option("--SNR", action = "store", type = "double", default = NULL),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE)
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

file_input = opt$file_input
file_output = opt$file_output
frame_rate = opt$frame_rate
SNR = opt$SNR
verbose = opt$verbose

# For test ####################################################################
# file_input = "D:\\_MATLAB_2PI\\peak_detection_data_input_I32_1.csv"
# file_output = "D:\\_MATLAB_2PI\\peak_detection_data_output_I32_1.csv"
# frame_rate = 5.0
# SNR = 2
# verbose = TRUE
###############################################################################


# Load trace
if (verbose) {
    print(paste("Reading data from: ",file_input, sep = ""))
}
data <- read.table(file_input, sep = ",", header = FALSE)
if (verbose) {
    print(paste("Making time series at ",frame_rate,"Hz", sep = ""))
}
series = ts(c(t(data)), start = 0, frequency = frame_rate)

# Compute good scale range to detect peaks
if (verbose) {
    print("Computing optimal scale range to detect calcium peaks")
}
gaus2_central_frequency = .3
lo_freq = 100   # number of peaks in 1 s
hi_freq = 1/30  # 1 peak in 30 s
scale_range = gaus2_central_frequency / (c(lo_freq, hi_freq) * (1/frame_rate))
scale_range[1] = max(floor(scale_range[1]) - 1, 1/frame_rate)
scale_range[2] = ceiling(scale_range[2])

# Compute CWT and its tree
if (verbose) {
    print("Computing CWT and its tree")
}
wav_trans = wavCWT(series, wavelet = "gaussian2", n.scale = 200, scale.range = scale_range)
wav_tree = wavCWTTree(wav_trans, n.octave.min = 0.8, tolerance = 0.0, type = "maxima")

# Compute peaks
if (verbose) {
    print("Computing peaks")
}
result <- wavCWTPeaks(wav_tree, snr.min = SNR, scale.range = scale_range, length.min = 5, noise.span = NULL, noise.fun = "quantile", noise.min = NULL)
# Get the list of good peaks
peaks = list(result[['x']])

# Store data to disk
if (verbose) {
    print("Writing to disk")
}
write.matrix(peaks, file_output)

print("done!")
