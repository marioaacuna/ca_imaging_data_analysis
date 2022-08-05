#!/usr/bin/env Rscript
rm(list = ls())
cat("\014")
set.seed(17)

## dependencies
# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("ggplot2", "ggbeeswarm", "tidyr", "optparse", "grid", "tibble", "plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://stat.ethz.ch/CRAN/")
library("reshape2")
# library(ggplot2)
# theme_set(
#   theme_classic() +
#     theme(legend.position = "top")
# )
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
  # make_option("--data_type", action = "store", type = "character", default = NULL),
  # make_option("--data_limits", action = "store", type = "character", default = NULL),
  # make_option("--log_scale", action = "store", default = FALSE),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE)
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)


file_input = opt$file_input
file_output = opt$file_output
verbose = opt$verbose


# For test ####################################################################
# file_input = "D:\\_MATLAB_CaImaging\\Proportion_EPM.csv"
# file_output = "D:\\_MATLAB_CaImaging\\test.pdf"
# verbose = TRUE
################################################################################


DATA <- read.csv(file_input, header = FALSE)

sessions = c(1:length(DATA))
session_names = paste('S', sessions,sep = "")
# DATA_df = data.frame(DATA)
colnames(DATA) <- session_names
if (verbose) {
  print("Plotting Proportion of EPM cells per session")
}

# Prepare data for plotting
melted <- melt(DATA)
stats_data = summarySE(melted, measurevar="value", groupvars= "variable")
ylabel = "Proportion of EPM neurons"
data_limits = c(0, 1)



# Plot scatter and mean+-SEM
gg.base <- ggplot(melted, aes(variable, value))
graph <- gg.base +
           geom_jitter(position = position_jitter(0.1), color = "darkgray", size = 7) + 
           geom_pointrange(aes(ymin = value-se, ymax = value+se), data = stats_data, size = 2) +
                    guides(colour = guide_legend(override.aes = list(size = 7))) +
                    ylab(ylabel) +
                    scale_y_continuous(
                      limits = data_limits,
                      expand = expand_scale(mult = c(0, 0)),
                      breaks = round(seq(data_limits[1], data_limits[2], length.out = 11), 2)
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
            )
# Print plot
if (verbose) {
  print("Printing figures")
}
pdf(file = file_output, width = 10, height = 10)
print(graph)
dev.off()

# Log end of script
if (verbose) {
  print('done')
}




