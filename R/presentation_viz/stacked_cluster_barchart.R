# library
library(ggplot2)
library(viridis)
library(hrbrthemes)

# import seurat_object
seurat_object <-
  readRDS('data/umap_seurat_object_with_10_by_10_SOM_metaclusters.Rds')

# Plot the breakdown of disease per cluster 
metadata <- 
  seurat_object@meta.data

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
  mutate(metacluster = as.integer(metacluster)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


max_metacluster <- 
          max(as.integer(metadata$metacluster))

metadata %>% 
  ggplot() + 
  geom_bar(aes(fill=cancer_type, y=n, x=metacluster), 
           position="stack", stat="identity") +
  scale_fill_manual(values = c("#77d1e5",
                               "#e8b8a1",
                               "#c6b6e5",
                               "#b0d9b1")) +
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
  seurat_object@meta.data

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
  mutate(metacluster = as.integer(metacluster)) %>% 
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



## Lets also make the presentation quality heatmap here
metadata <- 
  seurat_object@meta.data

total_cells_per_disease <- 
  metadata$cancer_type %>% 
  table()

metadata %>% 
  count(metacluster, cancer_type) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_disease[cancer_type])




# Lets make a heatmap of cancer type / cancer type detailed by
# metacluster
total_cells_per_disease <- 
  metadata$cancer_type_detailed %>% 
  table()

patient_count_per_disease <- 
  seurat_object@meta.data %>% 
  count(cancer_type_detailed) %>% 
  mutate(n = n / 100e3) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))

max_metacluster <- 
  metadata$metacluster %>% 
  as.integer() %>% 
  max()

num_diseases <- 
  metadata$cancer_type_detailed %>% 
  unique() %>% 
  length()

levels(metadata$metacluster) <- 
  rev(levels(metadata$metacluster))

metadata %>% 
  count(metacluster, cancer_type_detailed) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_disease[cancer_type_detailed]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = cancer_type_detailed, 
                y = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = 'Prop Disease Represented') + 
  geom_text(data = patient_count_per_disease, 
            aes(x = cancer_type_detailed, 
                y = max_metacluster + 1, label = n), 
            angle = 0, size = 3, 
            fontface = 'italic') +
  geom_text(data = cell_count_per_cluster, 
            aes(x = num_diseases-1.5, 
                y = metacluster, 
                label = n), 
            size = 3, 
            hjust = -1,
            angle = 0, 
            fontface = 'italic') + 
  xlab('') + ylab('Metcluster') + 
  theme(axis.text.x = element_text(angle = -35, 
                                   hjust = 0), 
        axis.text.y.left = element_text(hjust = 1),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(0, max_metacluster+2), 
                  xlim=c(0.8, num_diseases+4)) + 
  annotate(y=0.45, yend=0.45, 
           x=0.5, xend=num_diseases+0.5, 
           lwd=0.85, geom="segment") 

