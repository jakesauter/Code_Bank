metadata <- 
  seurat_object@meta.data %>% 
  filter(metacluster != 5) %>%
  filter(flow_cytometer %in% 
           c("R65822R1018",
             "R658222R1012", 
             "R66093700082", 
             "R66093700081")) 

metadata$metacluster <- 
  factor(metadata$metacluster, 
         levels = sort(unique(metadata$metacluster)))

total_cells_per_machine <- 
  metadata$flow_cytometer %>% 
  table()

cell_count_per_cluster <- 
  metadata %>% 
  mutate(metacluster = as.integer(metacluster)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


patient_count_per_machine <- 
  metadata %>% 
  count(flow_cytometer) %>% 
  mutate(n = n / 100e3) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))

max_metacluster <- 
  metadata$metacluster %>% 
  as.integer() %>% 
  max()

num_machines <- 
  metadata$flow_cytometer %>% 
  unique() %>% 
  length()


p1 <- metadata %>% 
  count(metacluster, flow_cytometer) %>% 
  tidyr::complete(metacluster, flow_cytometer) %>%
  rowwise() %>% 
  mutate(n = n / total_cells_per_machine[flow_cytometer]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = flow_cytometer, 
                y = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = 'Prop. of Cells Represented', 
                       na.value = '#4b0057') + 
  geom_text(data = cell_count_per_cluster, 
            aes(x = num_machines, 
                y = metacluster, 
                label = n), 
            size = 3, 
            hjust = -1,
            angle = 0, 
            fontface = 'italic') + 
  xlab('Flow Cytometer') + ylab('SOM Metacluster') + 
  theme(axis.text.x = element_text(angle = -25, 
                                   hjust = 0), 
        axis.text.y.left = element_text(hjust = 1),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(0, max_metacluster+2), 
                  xlim=c(1, num_machines+1)) + 
  annotate(y=0.45, yend=0.45, 
           x=0.5, xend=num_machines+0.5, 
           lwd=0.85, geom="segment") 


metadata <- 
  seurat_object@meta.data %>% 
  mutate(metacluster = louvain_clusters) %>%
  filter(flow_cytometer %in% 
           c("R65822R1018",
             "R658222R1012", 
             "R66093700082", 
             "R66093700081")) 

metadata$metacluster <- 
  factor(metadata$metacluster, 
         levels = sort(unique(metadata$metacluster)))

total_cells_per_machine <- 
  metadata$flow_cytometer %>% 
  table()

cell_count_per_cluster <- 
  metadata %>% 
  mutate(metacluster = as.integer(metacluster)) %>% 
  count(metacluster) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))


patient_count_per_machine <- 
  metadata %>% 
  count(flow_cytometer) %>% 
  mutate(n = n / 100e3) %>% 
  mutate(n = format(round(as.numeric(n), 1), big.mark=","))

max_metacluster <- 
  metadata$metacluster %>% 
  as.integer() %>% 
  max()

num_machines <- 
  metadata$flow_cytometer %>% 
  unique() %>% 
  length()


p2 <- metadata %>% 
  count(metacluster, flow_cytometer) %>% 
  tidyr::complete(metacluster, flow_cytometer) %>%
  rowwise() %>% 
  mutate(n = n / total_cells_per_machine[flow_cytometer]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = flow_cytometer, 
                y = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = 'Prop. of Cells Represented', 
                       na.value = '#4b0057') + 
  geom_text(data = cell_count_per_cluster, 
            aes(x = num_machines, 
                y = metacluster, 
                label = n), 
            size = 3, 
            hjust = -1,
            angle = 0, 
            fontface = 'italic') + 
  xlab('Flow Cytometer') + ylab('Louvain Community') + 
  theme(axis.text.x = element_text(angle = -25, 
                                   hjust = 0), 
        axis.text.y.left = element_text(hjust = 1),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  coord_cartesian(ylim=c(0, max_metacluster+2), 
                  xlim=c(1, num_machines+1)) + 
  annotate(y=0.45, yend=0.45, 
           x=0.5, xend=num_machines+0.5, 
           lwd=0.85, geom="segment") 

p1 + p2