# Lets make a heatmap of cancer type / cancer type detailed by
# metacluster

metadata <- 
  seurat_object@meta.data %>% 
  filter(metacluster != 5)

metadata$metacluster <- 
  factor(metadata$metacluster, 
         levels = sort(unique(metadata$metacluster)))

total_cells_per_disease <- 
  metadata$cancer_type_detailed %>% 
  table()

patient_count_per_disease <- 
  seurat_object@meta.data %>% 
  count(cancer_type_detailed) %>% 
  mutate(n = n / 100e3) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))

cell_count_per_cluster <- 
  metadata %>% 
  mutate(metacluster = as.integer(metacluster)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


max_metacluster <- 
  metadata$metacluster %>% 
  as.integer() %>% 
  max()

num_diseases <- 
  metadata$cancer_type_detailed %>% 
  unique() %>% 
  length()

p1 <- metadata %>% 
  count(metacluster, cancer_type_detailed) %>% 
  tidyr::complete(metacluster, cancer_type_detailed) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_disease[cancer_type_detailed]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = cancer_type_detailed,
                y = metacluster,
                fill = n)) +
  scale_fill_viridis_c(name = 'Prop. of Disease Represented', 
                       na.value = '#4b0057') +
  geom_text(data = patient_count_per_disease,
            aes(x = cancer_type_detailed,
                y = max_metacluster + 1, label = n),
            angle = 0, size = 3,
            fontface = 'italic') +
  geom_text(data = cell_count_per_cluster,
            aes(x = num_diseases-.5,
                y = metacluster,
                label = n),
            size = 3,
            hjust = -1,
            angle = 0,
            fontface = 'italic') +
  xlab('') + ylab('SOM Metacluster') +
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

# Lets make a heatmap of cancer type / cancer type detailed by
# metacluster
metadata <- seurat_object@meta.data

metadata <- 
  seurat_object@meta.data %>% 
  mutate(metacluster = louvain_clusters)

metadata$metacluster <- 
  factor(metadata$metacluster, 
         levels = sort(unique(metadata$metacluster)))

total_cells_per_disease <- 
  metadata$cancer_type_detailed %>% 
  table()

patient_count_per_disease <- 
  seurat_object@meta.data %>% 
  count(cancer_type_detailed) %>% 
  mutate(n = n / 100e3) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))

cell_count_per_cluster <- 
  metadata %>% 
  mutate(metacluster = as.integer(metacluster)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


max_metacluster <- 
  metadata$metacluster %>% 
  as.integer() %>% 
  max()

num_diseases <- 
  metadata$cancer_type_detailed %>% 
  unique() %>% 
  length()

p2 <- metadata %>% 
  count(metacluster, cancer_type_detailed) %>% 
  tidyr::complete(metacluster, cancer_type_detailed) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_disease[cancer_type_detailed]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = cancer_type_detailed,
                y = metacluster,
                fill = n)) +
  scale_fill_viridis_c(name = 'Prop. of Disease Represented', 
                       na.value = '#4b0057') +
  geom_text(data = patient_count_per_disease,
            aes(x = cancer_type_detailed,
                y = max_metacluster + 1, label = n),
            angle = 0, size = 3,
            fontface = 'italic') +
  geom_text(data = cell_count_per_cluster,
            aes(x = num_diseases-.5,
                y = metacluster,
                label = n),
            size = 3,
            hjust = -1,
            angle = 0,
            fontface = 'italic') +
  xlab('') + ylab('Louvain Community') +
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





p1 

p2
