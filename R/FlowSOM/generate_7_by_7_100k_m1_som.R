
library(dplyr)
library(magrittr)
library(flowCore)
library(FlowSOM)

filter <- dplyr::filter

flow_dirs_root <-
  "/home/sauterj1/shared_data_folder/sauterj1/Anon_Flow_Folders"

flow_directories <-
  list.files(flow_dirs_root,
             full.names = TRUE)

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
                           xdim = 7,
                           ydim = 7,
                           maxMeta = 40,
                           colsToUse = seq_len(ncol(m1_fcs_data_list[[1]])))

bmp('7_by_7_metaclustering_elbow_curve.bmp', 
    width = 640, 
    height = 400)

meta_clustering <- 
  FlowSOM::MetaClustering(
    data = m1_som$map$codes, 
    method = "metaClustering_consensus", 
    max = 40, 
    seed = 497, 
    plot = TRUE)


dev.off()

saveRDS(meta_clustering, 
        '7_by_7_metaclustering_results.Rds')


cat('\n\nDetermined number of metaclusters is: ', FlowSOM::NMetaclusters(m1_som))

saveRDS(m1_som, '7_by_7_SOM_empirical_num_clusters.Rds')