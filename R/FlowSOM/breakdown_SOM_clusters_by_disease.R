library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(flowCore)
library(stringr)
library(RColorBrewer)

filter <- dplyr::filter


m1_som <- 
  readRDS('data/m1_som_100_cells_per_patient.Rds')

flow_ids <- 
  read.table('data/filename_order_100k_per_patient.csv', 
             sep = ',') %>% 
  unlist() %>% 
  str_extract('F[0-9\\-]+$')

patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE)

patient_df <- 
  patient_df %>% 
  filter(Exclusion.Reason == "Not Excluded", 
         !str_detect(Flow.Data.Folder, ','), 
         !str_detect(M1_path, ','), 
         !str_detect(M2_path, ',')) 


## Getting patient df in the same order
# as cells generated in the SOM 
patient_df <- 
  patient_df %>% 
  filter(Accession.Number %in% flow_ids)


rownames(patient_df) <- 
  patient_df$Accession.Number

patient_df <- 
  patient_df[flow_ids, ]

# Now get the cells to metacluster
# mapping
cells_to_cluster_mapping <- 
  FlowSOM::GetClusters(m1_som)

cluster_to_metacluster_mapping <- 
  m1_som$metaclustering

cells_to_metacluster_mapping <- 
  cells_to_cluster_mapping %>% 
  cluster_to_metacluster_mapping[.]


# Now associdate each cell to the disease the 
# patient it came from had 

cells_to_disease_mapping <- 
  patient_df$CbioPortal.Cancer.Type  %>% 
  lapply(function(x) rep(x, 100e3)) %>% 
  unlist()

cells_to_specific_disease_mapping <- 
  patient_df$CbioPortal.Cancer.Type.Detailed  %>% 
  lapply(function(x) rep(x, 100e3)) %>% 
  unlist()


detailed_cancer_types <- 
  patient_df$CbioPortal.Cancer.Type.Detailed %>% 
  unique()

high_level_cancer_types <- 
  vector('character', length(detailed_cancer_types))

for (i in seq_along(detailed_cancer_types)) {
 
  detailed_cancer_type <- 
    detailed_cancer_types[i]
  
  high_level_cancer_types[i] <- 
    patient_df %>% 
    filter(CbioPortal.Cancer.Type.Detailed == detailed_cancer_type) %>% 
    select(CbioPortal.Cancer.Type) %>% 
    head(1) %>% unlist()
}

names(high_level_cancer_types) <- 
  detailed_cancer_types

# Now for each cluster, breakdown inclusion 
# of diseases

unique_diseases <- unique(cells_to_specific_disease_mapping)

df <- matrix(ncol = length(unique_diseases), 
             nrow = length(levels(m1_som$metaclustering)), 
             0) %>% 
  set_colnames(unique_diseases)

for (i in seq_along(levels(m1_som$metaclustering))) {

  print(i)
  
  df[i, ] <-
    cells_to_specific_disease_mapping[
      cells_to_metacluster_mapping == i] %>% 
    table() %>% .[unique_diseases] 
}

df <- 
  as.data.frame(df) 

df[is.na(df)] = 0

# Plot this dataframe as a heatmap

plot_df <- df %>% 
  apply(1, function(x) x / sum(x)) %>% 
  t() %>% as.data.frame()
  
plot_df$cluster <- seq(1, nrow(df))


plot_df <- plot_df %>% 
  tidyr::pivot_longer(-c(cluster))

plot_df$cluster <- 
  factor(plot_df$cluster, 
         levels = rev(levels(m1_som$metaclustering))) 

cancer_type_levels <- 
  high_level_cancer_types %>%
  sort(index.return = TRUE) %>% 
  .$ix %>% detailed_cancer_types[.]

plot_df$name <- 
  factor(plot_df$name, 
         levels = cancer_type_levels)



# Heatmap 
plot_df %>% 
  ggplot() + 
  geom_tile(aes(name, cluster, fill= value)) + 
  ylab('Cluster Label') + xlab('') + 
  # scale_x_discrete(labels = high_level_cancer_types) + 
  theme(axis.text.x = element_text(angle = -25, hjust = 0))
  # theme(axis.text.x = element_blank())

