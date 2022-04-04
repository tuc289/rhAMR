library(reshape2)
library(ggplot2)
library(ggthemes)
library(cowplot)
library(ggpubr)


regression_ra <- read.table("~/regression_ra_rhamr.csv", sep=",", header=T, row.names=1)
regression_abs <- read.table("~/Desktop/regression_rhAMR.csv", sep =",", header=T, row.names=1)
plot_data <- regression_ra

plot <- regression_abs[,3:5]
plot <- cbind("Gene" = row.names(plot), plot)
colnames(plot) <- c("Gene", "100", "1000", "10000")
threshold <- rep(NA, nrow(plot))

for (i in 1:nrow(plot)){
  if (plot[i,2]> 1e-03) {
    threshold[i] = "Pass"
  } else {
      threshold[i] = "Fail"}
}


plot_melt <- melt(plot, id.vars=c("Gene"))
plot_all <- ggplot(plot_melt, aes(x = variable, y = value, fill = Gene)) + 
  geom_bar(stat = "Identity") +
  theme_classic() + 
  labs(x = "Dilution factor", y = "Gene Counts") + 
  theme(legend.position = "none") + facet_wrap(~ Gene, ncol = 12) + 
  theme(axis.title = element_text(size=18)) +
  scale_y_log10()

plot_all
ggsave("barplot.tiff", plot = plot_all, device= "tiff", width = 10, height = 7.5, units = "in", 
       dpi = 600)

good_amplification <- subset.data.frame(plot_data, plot_data$p0_1.100 > 1e-03)
good_amplification_plot <- good_amplification[,2:4]
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
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(face="bold", size =16)) +
  theme(axis.text.y = element_text(face="bold", size =16)) + 
  theme(axis.title = element_text(size=18))

bad_amplification <- subset.data.frame(plot_data, plot_data$p0_1.100 < 1e-03)
bad_amplification_plot <- bad_amplification[,2:4]
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
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(face="bold", size =16)) +
  theme(axis.text.y = element_text(face="bold", size =16)) + 
  theme(axis.title = element_text(size=18))

panel <- plot_grid(good, bad, labels = c("A", "B"), label_size = 25)
panel 
ggsave("seperate_line_plot.tiff", panel, device = "tiff", width = 10, height = 5, units = "in", dpi = 600)


regression_ra

p<-ggplot(regression_ra, aes(x=box_plot, y=r.squared, color=box_plot, fill=box_plot)) +
  geom_boxplot(color = "black")
p_fin <- p  + theme_classic() + scale_fill_manual(values = c("Blue", "Red")) +   theme(legend.position = "none") + 
  xlab("") + 
  ylab(bquote('Regression Coefficient ' (R^2))) + 
  scale_x_discrete(labels= c("b" = "ARGs \n counts <25", "a" = " ARGs \n counts >25")) +
  theme(axis.text.x = element_text(face="bold", size =12)) +
  theme(axis.text.y = element_text(face="bold", size =12)) + 
  theme(axis.title = element_text(size=15))
ggsave("box_plot_LOD.tiff", plot = p_fin, device = "tiff", dpi= 600, width =3 , height =5 , units = "in")


