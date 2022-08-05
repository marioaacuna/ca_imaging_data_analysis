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
suppressMessages(library(plyr)) 
suppressMessages(library(Rmisc)) 
suppressMessages(library(dplyr))

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
data_limits = opt$data_limits
print(data_limits)
log_scale = opt$log_scale
verbose = opt$verbose


# For test ####################################################################
# file_input = "D:\\_MATLAB_CaImaging\\amplitude.csv"
# file_output = "D:\\_MATLAB_CaImaging\\test.pdf"
# data_type = "amplitude"
# data_limits = "0, 3"
# verbose = TRUE
################################################################################

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
  x <- data.wide$day_1[group_idx]
  data.wide$qnt_idx[group_idx] <- as.numeric(cut(x, breaks = quantile(x, probs = c(0, .25, .50, .75, 1.0), type = 8, na.rm = TRUE), include.lowest = TRUE))
}
# Assing ID to each row
data.wide$ID <- 0
id = 0
for (row in 1: sum(lengths(data.wide$group))) {
  ID = id + row
  data.wide$ID[row] <- ID
}

# Convert data to long (i.e., tall) format
day_1 <- data$day_1
# day_2 <- data$day_2
day_1_post <- data$day_1_post
day_2_post <- data$day_2_post
data.long <- tidyr::gather(data = data.wide, key = "day", value = "value", day_1, day_1_post, day_2_post)
data.long$group <- as.factor(data.long$group)
data.long$day <- factor(data.long$day, levels = day_names)

# 1D scatterplots + superimposed deciles
if (verbose) {
  print("Plotting distributions with deciles")
}
if (data_type == "amplitude") {
  ylabel = expression(paste("Event amplitude (% ", Delta, plain(F/F)[0], ")", sep = ""))
} else if (data_type == "amplitude_times_frequency") {
  ylabel = expression(paste("Event A*F (% ", Delta, plain(F/F)[0], ")", sep = ""))
} else {
  
  ylabel = expression(paste("Event frequency (Hz)", sep = ""))
}
# Remove 0s from plot
# data.toplot <- data.long[data.long$value != 0,]
zero.position = 10^-5
n_offset = .001

if (verbose) {
  print("Plotting distributions per group")
}
# Set position for dodging distributions
pjd <- position_jitterdodge(seed = 17)
pd <- position_dodge(width = pjd$dodge.width)

###################################################
# Plot data all together

# data.long %>%
#   group_by(day, group) %>% 
#   summarise (
#     mean_day = mean(day),
#     sd_day = sd(day),
#     mean_group= mean(group),
#     sd_group = sd(group)
#     )

stats_data = summarySE(data.long, measurevar="value", groupvars=c("day","group"))
gg.base <- ggplot(data = data.long, aes(x = day, y = value,  group = day))
gg.idline = gg.base + geom_line(aes(color = group, group = ID), size = 0.5, alpha = 0.3)
gg.Gline = gg.idline + geom_point(aes(color = group), alpha = 0.3)
graph3 <- gg.Gline + 
  stat_summary(aes(group = group),
                         geom = "line",
                         fun.y = mean,
                         size = 3) +
  stat_summary(
    aes(group = group, color = group),
    geom = "point",
    fun.y = mean,
    shape = 16,
    size = 6
  ) +
  geom_errorbar(
    data = stats_data,
    aes(
      ymin = value - se,
      ymax = value + se,
      color = group
    ),
    width = 0.2,
    size = 2
  ) +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    panel.background = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 4))) +
  ylab(ylabel) +
  scale_y_continuous(
    limits = data_limits,
    expand = expand_scale(mult = c(0, 0)),
    breaks = round(seq(data_limits[1], data_limits[2], length.out = 11), 2)
  )

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

