#!/usr/bin/env Rscript
rm(list = ls())
cat("\014")
set.seed(17)

## dependencies
# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("ggplot2", "ggbeeswarm", "tidyr", "optparse", "grid", "tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://stat.ethz.ch/CRAN/")

# Load libraries
suppressMessages(library(tidyr))
suppressMessages(library(tibble))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(optparse))

# Unpack inputs
option_list = list(
  make_option("--file_input", action = "store", type = "character", default = NULL),
  make_option("--file_output", action = "store", type = "character", default = NULL),
  make_option("--data_type", action = "store", type = "character", default = NULL),
  make_option("--data_limits", action = "store", type = "character", default = NULL),
  make_option("--log_scale", action = "store", default = FALSE),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE)
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

file_input = opt$file_input
file_output = opt$file_output
data_type = opt$data_type
print(data_type)
data_limits= opt$data_limits
print(data_limits)
log_scale = opt$log_scale
verbose = opt$verbose


# For test ####################################################################
 # file_input = "D:\\\amplitude.csv"
 # file_output = "D:\\\test.pdf"
 # data_type = "amplitude"
 # data_limits = "0, 5"
 # verbose = TRUE
###############################################################################

# Read data
if (verbose) {
  print("Reading data")
}
data <- read.table(file_input, header = TRUE, sep = ",")
data <- na.omit(data)
group_names = unique(data$group)

# Get limits data
data_limits = as.numeric(unlist(lapply(strsplit(data_limits, split = ","), trimws)))

# Convert the variable day_names to a list of strings
day_names = colnames(data)
day_names = day_names[day_names != "group"]
# Replace underscores with spaces
day_names_str = unlist(lapply(day_names, sub, pattern = "_", replacement = " "))

# Assign datapoints in each group to a quantile based on day 1
data.wide <- data
data.wide$qnt_idx <- 0
for (group in group_names) {
  group_idx = which(data.wide$group == group)
  x <- data.wide$Session_1[group_idx]
  data.wide$qnt_idx[group_idx] <- as.numeric(cut(x, breaks = quantile(x, probs = c(0, .25, .50, .75, 1.0), type = 8, na.rm = TRUE), include.lowest = TRUE))
}
# Convert data to long (i.e., tall) format
day_1 <- data$Session_1
day_2 <- data$Session_2
day_3 <- data$Session_3
day_4 <- data$Session_4
day_5 <- data$Session_5
day_6 <- data$Session_6
day_7 <- data$Session_7
day_8 <- data$Session_8



#day_1_post <- data$day_1_post
#day_2_post <- data$day_2_post
data.long <- tidyr::gather(data = data.wide, key = "day", value = "value", Session_1, Session_2, Session_3, Session_4, Session_5, Session_6, Session_7, Session_8)#, day_2, day_1_post, day_2_post)
data.long$group <- as.factor(data.long$group)
data.long$day <- factor(data.long$day, levels = day_names)

# 1D scatterplots + superimposed deciles
if (verbose) {
  print("Plotting distributions with deciles")
}
if (data_type == "amplitude") {
  ylabel = expression(paste("Event amplitude (% ", Delta, plain(F/F)[0], ")", sep = ""))
} else {
  ylabel = expression(paste("Event frequency (Hz)", sep = ""))
}
# Remove 0s from plot
data.toplot <- data.long[data.long$value != 0,]
zero.position = 10^-5
n_offset = .001

if (verbose) {
  print("Plotting distributions per group")
}
# Set position for dodging distributions
pjd <- position_jitterdodge(seed = 17)
pd <- position_dodge(width = pjd$dodge.width)
# Draw plot
graph3 <- ggplot(data.long, aes(x = .data$group, y = .data$value, group = .data$day)) +
  geom_quasirandom(aes(colour = .data$qnt_idx), method = "pseudorandom", size = 1.5, alpha = .75, varwidth = TRUE, dodge.width = .75) +
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, geom = "crossbar", width = .75, position = pd, size = .5) +
  stat_summary(fun.y = function(x) {quantile(x, 0.25)}, fun.ymin = function(x) {quantile(x, 0.25)}, fun.ymax = function(x) {quantile(x, 0.25)}, geom = "crossbar", width = 0.75, position = pd, size = .1) +
  stat_summary(fun.y = function(x) {quantile(x, 0.75)}, fun.ymin = function(x) {quantile(x, 0.75)}, fun.ymax = function(x) {quantile(x, 0.75)}, geom = "crossbar", width = 0.75, position = pd, size = .1) +
  scale_colour_gradient2(midpoint = 2.5, low = "blue", mid = "white", high = "red", breaks = c(4, 3, 2, 1), labels = c("Q4", "Q3", "Q2", "Q1")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(),
        panel.background = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 4))) + 
  ylab(ylabel) 
  # scale_y_continuous(limits = data_limits, expand = expand_scale(mult = c(0, 0)), breaks = round(seq(data_limits[1], data_limits[2], length.out = 11), 2))
# Apply logarithmic scale
if (log_scale) {
  if (data_limits[1] == 0) {
    data_limits[1] = 1e-2
  }
  if (data_type == "amplitude") {
    ylabel = expression(paste("Event amplitude (log % ", Delta, plain(F/F)[0], ")", sep = ""))
  } else {
    ylabel = expression(paste("Event frequency (log Hz)", sep = ""))
  }
  graph3 <- graph3 + 
    scale_y_log10(breaks = exp(log(10) * seq(-5, 1)), labels = function(x) sprintf("%.1f", x), limits = data_limits, expand = expand_scale(mult = c(0, 0))) + 
    annotation_logticks(base = 10, sides = "l", short = unit(-0.1, "cm"), mid = unit(-0.2, "cm"), long = unit(-0.3, "cm")) + 
    theme(axis.text.y = element_text(margin = margin(r = 10))) +
    ylab(ylabel)
}

# Print plot
if (verbose) {
  print("Printing figures")
}
pdf(file = file_output, width = 10, height = 10)
print(graph3)
dev.off()

# Log end of script
if (verbose) {
  print('done')
}
