library(reshape2)
library(ggplot2)
library(dplyr)

for_R <- read.table("~/Downloads/cdc_panel.csv", sep =",", header = T, row.names=1)
str(for_R)
colnames(for_R) <- c("Type", "Class", "Mechanism", "Gene", "p0_1:100", "p0_1:1000", 
                     "p0_1:10000", "p2_1:1", "p2_1:100", "p2_1:1000", "p2_1:10000")
for_plot <- for_R[,4:11]


p0_gene_list <- c("ACRA", "ACRB", "ACRD", "ACRE", "ACRF", "ACRS", "AMPH", "ASMA", "BACA", 
                  "BAER", "BAES", "BLAEC", "CATA", "CPXAR", "CRP", "DFRA", "EMRA", "EMRB", "EMRD",
                  "EMRK", "EMRR", "EMRY", "EPTA", "EVGS", "GADW", "HNS", "KDPE", "KPNO", "MARA", "MARR", 
                  "MDFA", "MDTA", "MDTB", "MDTC", "MDTE", "MDTF", "MDTG", "MDTH", "MDTI", "MDTJ",
                  "MDTK", "MDTN", "MDTO", "MDTP", "MPHA", "MPHB", "MSBA", "MVRC", "PBP4B", "PMRF", 
                  "ROBA", "SOXS", "TETB")
p2_gene_list <- c("ACRB", "DFRA", "GADW", "MDTK")

p0_plot <- filter(for_plot, Gene %in% p0_gene_list)
p2_plot <- filter(for_plot, Gene %in% p2_gene_list)
p2_plot_genes <- p2_plot[1]

p0_plot <- p0_plot[1:4]
p2_plot <- p2_plot[5:8]
p2_plot <- cbind(p2_plot_genes, p2_plot)

p0_plot

p0_plot <- aggregate(.~ Gene, data = p0_plot, FUN=sum)

melted_p0 <- melt(p0_plot)
melted_p2 <- melt(p2_plot)

head(melted_p0)
melted_p0
head(melted_p2)
melted_p2

results_df <- data.frame(matrix(ncol=3, nrow = length(p0_gene_list)))
colnames(results_df) <- c("Gene", "r-squared", "p-value")
par(mfrow=c(8,8))
for (x in 1:length(p0_gene_list)){
  tryCatch({
    new_df <- filter(melted_p0, Gene %in% p0_gene_list[x])
    results_df[x,1] <- p0_gene_list[x]
    new_df$variable <- c("1", "2", "3")
    model <- lm(variable ~ log(value), data = new_df)
    summary <- summary(model)
    results_df[x,2] <- summary$r.squared
    results_df[x,3] <- summary$coefficients[8]
  }, error=function(e){cat("ERROR","from",p0_gene_list[x],"\n")})
  }

melted_p0_ABRB <- filter(melted_p0, Gene %in% "ACRB")
melted_p0_ABRB$variable <- c("1", "2", "3")
model <- lm(variable ~ log(value), data = melted_p0_ABRB)
summary(model)
plot(melted_p0_ABRB$variable, log(melted_p0_ABRB$value))
lines(melted_p0_ABRB$variable, log(melted_p0_ABRB$value))


p0_plot = ggplot(melted_p0, aes(x = variable, fill = Gene, y = value)) + 
  #facet_grid(vars(Gene), scales = "free") +
  geom_bar(stat = "identity") + 
  facet_wrap(~ Gene, nrow=10, ncol=10, scales = "free") + 
  #theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
  #      axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
  #      legend.text = element_text(size = 12, face = "bold", colour = "black"), 
  #      axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(trans='log10') +
  guides(fill="none") +
  xlab("Ratio")+
  ylab("Count")
p0_plot

ggsave("p0_plots.tiff", plot=p0_plot, device=tiff, width=30, height=15, units=c("in"))


p2_plot = ggplot(melted_p2, aes(x = variable, fill = Gene, y = value)) + 
  #facet_grid(vars(Gene), scales = "free") +
  geom_bar(stat = "identity") + 
  facet_wrap(~ Gene, nrow=10, ncol=10, scales = "free") + 
  #theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
  #      axis.title.y = element_text(size = 16, face = "bold"), legend.title = element_text(size = 16, face = "bold"), 
  #      legend.text = element_text(size = 12, face = "bold", colour = "black"), 
  #      axis.text.y = element_text(colour = "black", size = 12, face = "bold")) + 
  scale_y_continuous(trans='log10') +
  guides(fill="none") +
  xlab("Ratio")+
  ylab("Count")
p2_plot

ggsave("p2_plots.tiff", plot=p2_plot, device=tiff, width=30, height=15, units=c("in"))



