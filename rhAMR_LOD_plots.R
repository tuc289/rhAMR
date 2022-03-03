library(reshape2)
library(ggplot2)
library(ggthemes)
library(cowplot)


regression_ra <- read.table("~/regression_ra_rhamr.csv", sep=",", header=T, row.names=1)
regression_abs <- read.table("~/Desktop/regression_rhAMR.csv", sep =",", header=T, row.names=1)
plot_data <- regression_ra
data_for_sorting <- cbind(plot_data[2:4], regression_abs[,3])

plot <- plot_data[,2:4]
plot <- cbind("Gene" = row.names(plot), plot)
colnames(plot) <- c("Gene", "100", "1000", "10000")
threshold <- rep(NA, nrow(plot))

for (i in 1:nrow(plot)){
  if (plot[i,2]> 1e-03) {
    threshold[i] = "Pass"
  } else {
      threshold[i] = "Fail"}
}

plot <- cbind(plot, threshold)

plot_melt <- melt(plot, id.vars=c("Gene", "threshold"))
plot_melt$variable <- as.numeric(as.character(plot_melt$variable))
plot_all <- ggplot(plot_melt, aes(x = variable, y = value, color = threshold)) + 
  geom_line(aes(group = Gene)) + 
  scale_y_log10() +
  scale_x_log10() + 
  geom_point(aes(group = Gene)) + 
  theme_classic() + 
  labs(x = "Dilution factor", y = "Log relative abundance") + 
  scale_color_manual(name = "cut-off", values = c("Red", "Blue"))+
  theme(legend.position = "none")
plot_all
ggsave("lineplot.tiff", plot = plot_all, device= "tiff", width = 7.5, height = 5, units = "in", 
       dpi = 600)

plot_reg <- ggplot(data = plot_melt, aes(x = variable, y = value, Color = Gene)) +
  geom_point(aes(color = threshold)) + scale_y_log10() + scale_x_log10() + 
  geom_smooth(method = "glm", se = FALSE, aes(color = threshold), linetype = 1, size = 0.3) +
  scale_color_manual(name = "cut-off", values = c("Red", "Blue")) + 
  theme_classic() + 
  labs(x = "Dilution factor", y = "Log relative abundance") + 
  scale_color_manual(name = "cut-off", values = c("Red", "Blue"))+
  theme(legend.position = "none")
ggsave("lineplot_regression.tiff", plot = plot_reg, device= "tiff", width = 7.5, height = 5, units = "in")  



good_amplification <- subset.data.frame(data_for_sorting, data_for_sorting$p0_1.100 > 1e-03)
good_amplification_plot <- good_amplification[,1:4]
good_amplification_plot <- cbind("Gene" = row.names(good_amplification_plot), good_amplification_plot)
colnames(good_amplification_plot) <- c("Gene", "100", "1000", "10000")
good_amplification_plot <- good_amplification_plot[, 1:4]
good_melt <- melt(good_amplification_plot, id.vars = "Gene")
good_melt$variable <- as.numeric(as.character(good_melt$variable))
good <- ggplot(good_melt, aes(x = variable, y = value)) + 
  geom_line(aes(group = Gene), color = "Blue") + 
  scale_y_log10() +
  scale_x_log10() + 
  geom_point(aes(group = Gene), color = "Blue") + 
  scale_color_manual(values = c(rep("red", 34), rep("blue", 3)))+
  theme_classic() + 
  labs(x = "Dilution factor", y = "Log relative abundance") + 
  scale_color_manual(name = "cut-off", values = c("Red", "Blue"))+
  theme(legend.position = "none")

bad_amplification <- subset.data.frame(data_for_sorting, data_for_sorting$p0_1.100 < 1e-03)
bad_amplification_plot <- bad_amplification[,1:3]
bad_amplification_plot <- cbind("Gene" = row.names(bad_amplification_plot), bad_amplification_plot)
colnames(bad_amplification_plot) <- c("Gene", "100", "1000", "10000")
bad_melt <- melt(bad_amplification_plot)
bad_melt$variable <- as.numeric(as.character(bad_melt$variable))
bad <- ggplot(bad_melt, aes(x = variable, y = value)) + 
  geom_line(aes(group = Gene), color = "Red") + 
  scale_y_log10() +
  scale_x_log10() + 
  geom_point(aes(group = Gene), color = "Red") + 
  scale_color_manual(values = c(rep("red", 34), rep("blue", 3)))+
  theme_classic() + 
  labs(x = "Dilution factor", y = "Log relative abundance") + 
  scale_color_manual(name = "cut-off", values = c("Red", "Blue"))+
  theme(legend.position = "none")

panel <- plot_grid(good, bad, labels = c("A", "B"))
panel
ggsave("seperate_line_plot.tiff", panel, device = "tiff", width = 10, height = 5, units = "in", dpi = 600)
