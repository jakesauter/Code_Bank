# library
library(ggplot2)
library(dplyr)
library(viridis)
library(hrbrthemes)

# import seurat_object
seurat_object <-
  readRDS('data/umap_seurat_object_with_10_by_10_SOM_metaclusters.Rds')

# Plot the breakdown of disease per cluster 
metadata <- 
  seurat_object@meta.data %>% 
  mutate(metacluster = leiden_clusters)

total_cells_per_disease <- 
  metadata$cancer_type %>% 
  table()

metadata <- 
  metadata %>% 
  count(metacluster, cancer_type) %>% 
  group_by(metacluster) %>% 
  mutate(n = n / sum(n)) %>% 
  ungroup()
  
cell_count_per_cluster <- 
  seurat_object@meta.data %>% 
  mutate(metacluster = as.integer(leiden_clusters)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


max_metacluster <- 
          max(as.integer(metadata$metacluster))

metadata %>% 
  ggplot() + 
  geom_bar(aes(fill=cancer_type, y=n, x=metacluster), 
           position="stack", stat="identity") +
  # scale_fill_manual(values = c("#77d1e5",
  #                              "#e8b8a1",
  #                              "#c6b6e5",
  #                              "#b0d9b1")) +
  scale_fill_brewer(palette = 'Paired') + 
  geom_text(data = cell_count_per_cluster, 
            aes(x = metacluster - 0.3, y = 1.05, label=n), 
            size = 3, 
            angle = 45, 
            hjust = 0) + 
  xlab('Metacluster') + 
  ylab("Proportion of Cluster From Disease") + 
  guides(fill = guide_legend(title = '')) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = 'bottom') + 
  coord_cartesian(xlim=c(1, max_metacluster+1), ylim=c(0, 1.2)) + 
  annotate(x=0, xend=max_metacluster+0.5, y=-.05, yend=-.05, 
           lwd=0.75, geom="segment") + 
  annotate(x=0.45, xend=0.45, y=-0.05, yend=1.03, 
           lwd=0.85, geom="segment") 

# Lets make a heatmap of cancer type / cancer type detailed by
# metacluster
metadata <- 
  seurat_object@meta.data %>% 
  mutate(metacluster = leiden_clusters)

total_cells_per_disease <- 
  metadata$cancer_type_detailed %>% 
  table()

metadata <- 
  metadata %>% 
  count(metacluster, cancer_type_detailed) %>% 
  group_by(metacluster) %>% 
  mutate(n = n / sum(n)) %>% 
  ungroup()

max_metacluster <- 
  metadata$metacluster %>% 
  as.integer() %>% 
  max()

cell_count_per_cluster <- 
  seurat_object@meta.data %>% 
  mutate(metacluster = as.integer(leiden_clusters)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


metadata %>% 
  ggplot() + 
  geom_bar(aes(fill=cancer_type_detailed, y=n, x=metacluster), 
           position="stack", stat="identity") +
  scale_fill_manual(values = c("#63a11b",
                               "#6260d7",
                               "#a1cd44",
                               "#a44cbe",
                               "#b3e36d",
                               "#003089",
                               "#8fe987",
                               "#bf0c78",
                               "#01933e",
                               "#ce3ea3",
                               "#59ecbd",
                               "#c91d40",
                               "#006a27",
                               "#ff96f0",
                               "#4d6a00",
                               "#3596ff",
                               "#908a00",
                               "#0150a3",
                               "#ff9d45",
                               "#780058",
                               "#6b9e5f",
                               "#ff619d",
                               "#8c7400",
                               "#d4a6eb",
                               "#983900",
                               "#7d3261",
                               "#fa7742",
                               "#630027",
                               "#ff97af",
                               "#930012")) +
  geom_text(data = cell_count_per_cluster, 
            aes(x = metacluster - 0.3, y = 1.05, label=n), 
            size = 3, 
            angle = 45, 
            hjust = 0) + 
  xlab('Metacluster') + 
  ylab("Proportion of Cluster From Disease") + 
  guides(fill = guide_legend(title = '')) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "bottom", 
        legend.text = element_text(size = 6.5)) + 
  coord_cartesian(xlim=c(1, max_metacluster+1), ylim=c(0, 1.2)) + 
  annotate(x=0, xend=max_metacluster+0.5, y=-.05, yend=-.05, 
           lwd=0.75, geom="segment") + 
  annotate(x=0.45, xend=0.45, y=-0.05, yend=1.03, 
           lwd=0.85, geom="segment") 





