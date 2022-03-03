#Clustering ABRICATE results on AMRg detection using WGS data from 29 strains (PSU strains)
#This script is modified from original tutorial - http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/114-mca-multiple-correspondence-analysis-in-r-essentials/

install.packages(c("FactoMineR", "factoextra"))

library("FactoMineR")
library("factoextra")

data <- read.table("~/Desktop/cluster.csv", header =T, sep =",",row.names=1)
for(i in 2:150){
  data[,i] <- as.factor(data[,i])
}
str(data)
res.mca <- MCA(data, quali.sup = 1, graph = FALSE)

plot <- fviz_mca_ind(res.mca, repel=TRUE, ggtheme= theme_minimal(), max.overlaps=10)
plot_data <- plot$data
plot_data <- cbind(data[,1:2], plot_data)
p <- ggplot(plot_data, aes(x=x, y=y, color=Species)) + 
  geom_point() +
  geom_text_repel(label= rownames(plot_data), nudge_x= 0.1, size = 3) + 
  theme_classic() + labs(x = "PC1", y = "PC2")
p
ggsave("MCA.tiff", plot = p, device = "tiff", width = 7.5, height = 5, units = "in")

#Human fecal samples clustering
data[] <- lapply(data, function(x) gsub("P", "1", x))
data[] <- lapply(data, function(x) gsub("A", "0", x))
for(i in 2:150){
  data[,i] <- as.numeric(data[,i])
}
str(data)

df <- data.frame(matrix(ncol=8, nrow = 29))
colnames(df) <- row.names(data)[30:37]
rownames(df) <- row.names(data)[1:29]
for (x in 1:29){
  strain <- data[x,]
  for (y in 30:37) {
    human <- data[y,]
    combined <- rbind(strain, human)
    combined <- combined[,-1]
    for (z in 1:ncol(combined)){
      if(combined[1,z] == 1 & combined[2,z] == 0){tmp = 1}
      else{tmp = 0}
      a[z] = tmp
      df[x, y-29] = sum(a)
    } 
  }
}
df

pair <- as.data.frame(sapply(df, function(x) head(row.names(df)[order(x, decreasing = TRUE)],1)))
pair
df
write.csv(df, "table.csv")


##Strain pool##
cluster <- plot_data[,4:5]
fit <- kmeans(cluster, 4)
cluster_result <- fviz_cluster(fit, data = cluster)
fit$cluster
with_group <- cbind(data, group = fit$cluster)
with_group <- with_group[1:29,]

df_strain <- data.frame(matrix(ncol=3, nrow = ))
colnames(df_strain) <- row.names(data)[1:29]
rownames(df_strain) <- row.names(data)[1:29]
group1 <- with_group[with_group$group ==1,]
group2 <- with_group[with_group$group ==2,]
group3 <- with_group[with_group$group ==3,]

group1_optimal = NULL
for (i in 1:8){
  a = rbind(group2[i,-151],group1[,-151])
  for (z in 2:ncol(a)){
    if(a[1,z] == 1 & sum(a[-1,z]) == 0) {group1_optimal[z-1]= 1}
    else{group1_optimal[z-1] = 0}
  }
  print(sum(b))
}

  
  


