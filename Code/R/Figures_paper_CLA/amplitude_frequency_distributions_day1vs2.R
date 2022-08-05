#!/usr/bin/env Rscript
rm(list = ls())
cat("\014")  
set.seed(17)

## dependencies
if (suppressMessages(!require(rogme))) {
  if (suppressMessages(!require(devtools))) {
    install.packages("devtools", repos = "https://stat.ethz.ch/CRAN/")}
  suppressMessages(library(devtools))
  devtools::install_github("GRousselet/rogme")}

# See: http://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
list.of.packages <- c("ggplot2", "ggbeeswarm", "tidyr", "optparse", "grid", "tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages, repos = "https://stat.ethz.ch/CRAN/")

# Load libraries
suppressMessages(library(rogme))
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
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE)
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

file_input = opt$file_input
file_output = opt$file_output
data_type = opt$data_type
data_limits = opt$data_limits
verbose = opt$verbose

###############################################################################
# Split violin
###############################################################################
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             # Original function by Jan Gleixner (@jan-glx)
                             # Adjustments by Wouter van der Bijl (@Axeman)
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
                               quantiles <- create_quantile_segment_frame(data, draw_quantiles, split = TRUE, grp = grp)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           }
)

create_quantile_segment_frame <- function(data, draw_quantiles, split = FALSE, grp = NULL) {
  dens <- cumsum(data$density) / sum(data$density)
  ecdf <- stats::approxfun(dens, data$y)
  ys <- ecdf(draw_quantiles)
  violin.xminvs <- (stats::approxfun(data$y, data$xminv))(ys)
  violin.xmaxvs <- (stats::approxfun(data$y, data$xmaxv))(ys)
  violin.xs <- (stats::approxfun(data$y, data$x))(ys)
  if (grp %% 2 == 0) {
    data.frame(
      x = ggplot2:::interleave(violin.xs, violin.xmaxvs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  } else {
    data.frame(
      x = ggplot2:::interleave(violin.xminvs, violin.xs),
      y = rep(ys, each = 2), group = rep(ys, each = 2)
    )
  }
}

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, position = position, 
        show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

###############################################################################
# Annotation with relative position to plot
###############################################################################
annotate_textp <- function(label, x, y, facets = NULL, hjust = 0, vjust = 0, color = 'black', alpha = NA,
                           family = thm$text$family, size = thm$text$size, fontface = 1, lineheight = 1.0,
                           box_just = ifelse(c(x,y) < 0.5,0,1), margin = unit(size/2, 'pt'), thm = theme_get()) {
  x <- scales::squish_infinite(x)
  y <- scales::squish_infinite(y)

  tg <- grid::textGrob(
    label, x = 0, y = 0, hjust = hjust, vjust = vjust,
    gp = grid::gpar(col = alpha(color, alpha), fontsize = size, fontfamily = family, fontface = fontface, lineheight = lineheight)
  )
  ts <- grid::unit.c(grid::grobWidth(tg), grid::grobHeight(tg))
  vp <- grid::viewport(x = x, y = y, width = ts[1], height = ts[2], just = box_just)
  tg <- grid::editGrob(tg, x = ts[1]*hjust, y = ts[2]*vjust, vp = vp)
  inner <- grid::grobTree(tg, vp = grid::viewport(width = unit(1, 'npc') - margin*2, height = unit(1, 'npc') - margin*2))
  
  layer(
    data = NULL,
    stat = StatIdentity,
    position = PositionIdentity,
    geom = GeomCustomAnn,
    inherit.aes = TRUE,
    params = list(
      grob = grid::grobTree(inner), 
      xmin = -Inf, 
      xmax = Inf, 
      ymin = -Inf, 
      ymax = Inf
    )
  )
}
###############################################################################


# For test
# file_input = "D:\\_MATLAB_2PI\\amplitude.csv"
# file_output = "D:\\_MATLAB_2PI\\test.pdf"
# data_type = "amplitude"
# data_limits = "0.0, 10.0"
# verbose = TRUE

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
colnames(data) <- c("group", "day_1", "day_2")

# Set colors
green = c(0, 158, 115) / 255
magenta = c(204, 121, 167) / 255
palette = c(rgb(green[1], green[2], green[3]), rgb(magenta[1], magenta[2], magenta[3]))

# Assign datapoints in each group to a quantile based on day 1
data.wide <- data
data.wide$qnt_idx <- 0
for (group in group_names) {
  group_idx = which(data.wide$group == group)
  x <- data.wide$day_1[group_idx]
  data.wide$qnt_idx[group_idx] <- as.numeric(cut(x, breaks = quantile(x, probs = c(0, .25, .50, .75, 1.0), type = 8), include.lowest = TRUE))
}
# Convert data to long (i.e., tall) format
day_1 <- data$day_1
day_2 <- data$day_2
data.long <- tidyr::gather(data = data.wide, key = "day", value = "value", day_1, day_2)
data.long$group <- as.factor(data.long$group)
data.long$day  <- as.factor(data.long$day)

# Perform basic tests
if (verbose) {
  wt = wilcox.test(data.wide$day_1, data.wide$day_2)
  print(paste("Wilcoxon signed-rank test: p-value", wt$p.value))
}

# Compute shift function for dependent groups
if (verbose) {
  print("Computing shift function")
}
data.sf <- mkt2(day_1, day_2, gr_names = "day", obs_names = "value", group_labels = c("day_1", "day_2"))
# sf <- shiftdhd(data = data.sf, formula = value ~ day, nboot = 1000, todo = list(c("day_2", "day_1")))

# 1D scatterplots + superimposed deciles
if (verbose) {
  print("Plotting distributions with deciles")
}
if (data_type == "amplitude") {
  ylabel = "Event amplitude (a.u.)"
} else {
  ylabel = "Event frequency (Hz)"
}
# Remove 0s from plot
data.toplot <- data.sf[data.sf$value != 0,]
n.active = c(length(which(data.toplot$day == "day_1")), length(which(data.toplot$day == "day_2")))
n.total = c(length(which(data.sf$day == "day_1")), length(which(data.sf$day == "day_2")))
n.inactive = n.total - n.active
data.toplot$x = 1
title_offset = .2
zero.position = 10^-5
n_offset = .001

graph1 <- ggplot(data.toplot, aes(x = .data$x, y = .data$value, fill = .data$day)) +
  # Add split violin
  geom_split_violin(adjust = .75, draw_quantiles = c(.25, .50, .75), alpha = .75, lwd = 1, scale = "count", trim = FALSE) +
  # Add text
  annotate_textp(day_names_str[1], x = .25, y = 1, hjust = .5, vjust = 1, size = 18, box_just = c(.5, 1)) + 
  annotate_textp(day_names_str[2], x = .75, y = 1, hjust = .5, vjust = 1, size = 18, box_just = c(.5, 1)) + 
  annotate_textp(paste("n active", n.active[1], "\nn inactive", n.inactive[1]), x = .25, y = 0, hjust = .5, vjust = 0, size = 14, box_just = c(.5, 0)) + 
  annotate_textp(paste("n active", n.active[2], "\nn inactive", n.inactive[2]), x = .75, y = 0, hjust = .5, vjust = 0, size = 14, box_just = c(.5, 0)) + 
  # Fix axes appearance  
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(limits = data_limits, expand = expand_scale(mult = c(0, 0))) +
  theme_bw() + 
  theme(panel.border = element_blank(), 
        axis.line = element_line(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        legend.title = element_blank(),
        legend.position = "none") + 
  labs(x = "", y = ylabel) + 
  scale_fill_manual(labels = day_names_str, values = palette)

# Apply logarithmic scale
if (data_type == "amplitude") {
  if (data_limits[1] == 0){
    data_limits[1] = 1e-2
  }
  graph1 <- graph1 + 
    scale_y_log10(breaks = exp(log(10) * seq(-5, 1)), labels = function(x) sprintf("%.1f", x), limits = data_limits, expand = expand_scale(mult = c(0, 0))) + 
    annotation_logticks(base = 10, sides = "l", short = unit(-0.1, "cm"), mid = unit(-0.2, "cm"), long = unit(-0.3, "cm")) + 
    theme(axis.text.y = element_text(margin = margin(r = 10)))
  # Rotate logticks to outside the plot area
  graph1.build.gtable <- ggplot_gtable(ggplot_build(graph1))
  graph1.build.gtable$layout$clip
  graph1.build.gtable$layout$clip <- c("on", "off", "off", "off", "off", "off", "off", "off", "off")
}

# # Plot shift function
# if (verbose) {
#   print("Plotting shift function")
# }
# graph2 <- plot_sf(sf, plot_theme = 1, diffint_col = "black", diffint_size = 1, q_line_col = "black", q_line_size = 0.5)
# graph2[[1]] <- graph2[[1]] + theme(panel.border = element_blank(), 
#                              axis.line = element_line(),
#                              panel.background = element_blank()) +
#   xlab(paste("Deciles for", day_names_str[2])) + 
#   ylab(paste("Decile difference", day_names_str[2], "-", day_names_str[1]))


# Linked stripcharts
if (verbose) {
  print("Plotting distributions per group")
}
# Set position for dodging distributions
pjd <- position_jitterdodge(seed = 17)
pd <- position_dodge(width = pjd$dodge.width)
# Draw plot
graph3 <- ggplot(data.long, aes(x = .data$group, y = .data$value, group = .data$day)) +
  geom_quasirandom(aes(colour = .data$qnt_idx), method = "pseudorandom", size = .5, alpha = .75, varwidth = TRUE, dodge.width = .75) +
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
  ylab(ylabel) + 
  scale_y_continuous(limits = data_limits, expand = expand_scale(mult = c(0, 0)))
# Apply logarithmic scale
if (data_type == "amplitude") {
  graph3 <- graph3 + 
    scale_y_log10(breaks = exp(log(10) * seq(-5, 1)), labels = function(x) sprintf("%.1f", x), limits = data_limits, expand = expand_scale(mult = c(0, 0))) + 
    annotation_logticks(base = 10, sides = "l", short = unit(-0.1, "cm"), mid = unit(-0.2, "cm"), long = unit(-0.3, "cm")) + 
    theme(axis.text.y = element_text(margin = margin(r = 10)))
}


# # Distribution of differences
# if (verbose) {
#   print("Plotting difference distributions by group")
# }
# data.wide$difference <- data.wide$day_2 / data.wide$day_1
# data.wide$difference[is.na(data.wide$difference)] <- 0
# # Because differences will also be negative, if the scale has to be trasnformed in log, many data points will be hidden.
# # Therefore, we directly transform positive values to log scale
# data.wide$difference_raw <- data.wide$difference
# data.wide$difference[data.wide$difference_raw > 1] = log2(data.wide$difference[data.wide$difference_raw > 1])
# data.wide$difference[data.wide$difference_raw <= 1] = log2(data.wide$difference[data.wide$difference_raw <= 1])
# 
# graph4 <- ggplot(data.wide, aes(x = .data$group, y = .data$difference)) + 
#   geom_abline(intercept = 0, slope = 0, linetype = 2) +
#   geom_quasirandom(aes(color = .data$group), method = "pseudorandom", size = .75, alpha = .75, varwidth = TRUE, dodge.width = .75) +
#   theme_bw() +
#   theme(legend.position = "none",
#         axis.text.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14),
#         axis.title.y = element_text(size = 16, face = "bold"),
#         plot.title = element_text(colour = "black", size = 20),
#         panel.border = element_blank(), 
#         axis.line = element_line(),
#         panel.background = element_blank()) +
#   # scale_colour_gradient2(midpoint = 0, low = "blue", mid = "gray80", high = "red", breaks = c(1, 0, -1), labels = c("+", "~", "-")) +
#   scale_colour_manual(values = palette) +
#   ylab(paste(day_names_str[2], "-", day_names_str[1], "difference in", paste0(tolower(substr(ylabel, 1, 1)), substr(ylabel, 2, nchar(ylabel))))) +
#   xlab(""); graph4 # + scale_y_log10()
#   # scale_y_log10(breaks = exp(log(10) * seq(-5, 1)), labels = function(x) sprintf("%.2f", x)) 

# # Calculate confidence of interval
# if (verbose) {
#   out <- hdpbci(data.wide$difference)
#   print(paste("Confidence interval of the difference", paste0("[", out$ci[1], ", ", out$ci[2], "]")))
# }


# Print plots on multiple pages
if (verbose) {
  print("Printing figures")
}
pdf(file = file_output, width = 10, height = 10)
if (data_type == "amplitude") {
  grid.draw(graph1.build.gtable)
} else{
  print(graph1)
}
# print(graph2)
print(graph3)
dev.off()

# Log end of script
if (verbose) {
  print('done')
}
