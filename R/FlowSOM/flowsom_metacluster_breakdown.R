library(dplyr)
library(Seurat)
library(ggplot2)
library(FlowSOM)
library(flowCore)
library(magrittr)
library(patchwork)

filter <- dplyr::filter

seurat_object <- 
  readRDS('data/umap_seurat_object_with_metaclusters.Rds')


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


metadata %>% 
  count(metacluster, cancer_type_detailed) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_disease[cancer_type_detailed]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = cancer_type_detailed, 
                y = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = '% Disease Represented') + 
  theme(axis.text.x = element_text(angle = -35, 
                                   hjust = 0), 
        panel.background=element_rect(fill="#450d54", 
                                      color = '#450d54'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab('') + ylab('Metcluster')


# Lets make a heatmap of cancer type / cancer type detailed by
# metacluster
total_cells_per_disease <- 
  metadata$cancer_type %>% 
  table()


metadata %>% 
  count(metacluster, cancer_type) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_disease[cancer_type]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = cancer_type, 
                y = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = '% Disease Represented') + 
  theme(axis.text.x = element_text(angle = -10, 
                                   hjust = 0), 
        panel.background=element_rect(fill="#450d54", 
                                      color = '#450d54'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab('') + ylab('Metcluster')


# Lets make a heatmap of flow cytometer by
# metacluster
flow_cyto <- 
  metadata$flow_cytometer %>% 
  {dplyr::case_when(. == "R658222R1012" ~ "Fortessa I", 
                   . == "R65822R1018" ~  "Fortessa II",
                   . == "R66093700072" ~ "Symphony I",
                   . == "R66093700081" ~ "Symphony II",
                   . == "R66093700082" ~ "Symphony III",
                   . == "V657338000098" ~ 'Canto 5')}


total_cells_per_flow_cytometer <- 
  flow_cyto %>% 
  table()


metadata %>% 
  mutate(flow_cytometer = flow_cyto) %>% 
  count(metacluster, flow_cytometer) %>% 
  rowwise() %>% 
  mutate(n = n / total_cells_per_flow_cytometer[flow_cytometer]) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(x = flow_cytometer, 
                y = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = '% Disease Represented') + 
  theme(axis.text.x = element_text(angle = -10, 
                                   hjust = 0), 
        panel.background=element_rect(fill="#450d54", 
                                      color = '#450d54'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab('') + ylab('Metcluster')




metadata %>% 
  mutate(date = lubridate::floor_date(experiment_date, 'month')) %>%
  filter(date > lubridate::ymd('2020-01-01')) %>% 
  count(flow_cytometer, date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, 
                y = flow_cytometer, 
                fill = n)) 




### Patient contribution to each cluster



metadata %>% 
  count(metacluster, flow_dir) %>% 
  rowwise() %>% 
  mutate(n = n / 100e3) %>%
  ungroup() %>% 
  ggplot() + 
  geom_tile(aes(y = flow_dir, 
                x = metacluster, 
                fill = n)) + 
  scale_fill_viridis_c(name = "% Patient's Cells Represented") + 
  theme(axis.text.y = element_text(size = 6), 
        panel.background=element_rect(fill="#450d54", 
                                      color = '#450d54'), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  ylab('') + xlab('Metcluster')




