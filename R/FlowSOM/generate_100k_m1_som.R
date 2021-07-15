
library(arrow)
library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(FlowSOM)
library(RColorBrewer)

filter <- dplyr::filter

# Gather the order of the flow ids in the UMAP results
flow_ids <- 
  read.table('data/filename_order_100k_per_patient.csv', 
             sep = ',') %>% 
  unlist() %>% 
  str_extract('F[0-9\\-]+$')


flow_dirs_root <-
  "/Users/sauterj1/Documents/Woodlist/Anon_Flow_Folders_Server_Version"

flow_directories <-
  list.files(flow_dirs_root,
             full.names = TRUE)

som_data_flow_ids <- basename(flow_directories)
current_flow_id_order <- som_data_flow_ids
som_data_flow_ids <- seq_along(som_data_flow_ids)
names(som_data_flow_ids) <- current_flow_id_order
som_data_flow_id_order <- som_data_flow_ids[flow_ids]

flow_directories <- flow_directories[som_data_flow_id_order]

m1_fcs_data_list <-
  vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  
  cat('Reading fcs file:', i, 'of',
      length(flow_directories), '\n')
  
  fcs_file <- file.path(flow_directory,
                        'processed_M1_subsampled_100k.fcs')
  
  m1_fcs_data_list[[i]] <-
    read.FCS(fcs_file,
             transformation = FALSE) %>%
    exprs()
  
}

m1_fcs_data_mat <-
  m1_fcs_data_list %>%
  do.call(rbind, .)

m1_fcs_ff <-
  flowFrame(m1_fcs_data_mat)

m1_som <- FlowSOM::FlowSOM(m1_fcs_ff,
                           xdim = 10,
                           ydim = 10,
                           nClus = 20,
                           colsToUse = seq_len(ncol(m1_fcs_data_list[[1]])))


saveRDS(m1_som, 'data/m1_som_100_cells_per_patient_10_by_10_dim_20_clus.Rds')