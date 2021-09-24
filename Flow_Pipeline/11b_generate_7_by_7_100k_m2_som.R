
library(dplyr)
library(magrittr)
library(Seurat)
library(flowCore)
library(FlowSOM)

filter <- dplyr::filter

seurat_object <- 
  readRDS('/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_centralized_Seurat_object.Rds')

m2_fcs_data_mat <-
  seurat_object@assays$MultiParamFlowCyto@data %>% 
  as.matrix()

m2_fcs_ff <-
  flowFrame(m2_fcs_data_mat)

m2_som <- FlowSOM::FlowSOM(m2_fcs_ff,
                           xdim = 7,
                           ydim = 7,
                           maxMeta = 40,
                           colsToUse = seq_len(ncol(m2_fcs_data_list[[1]])))

bmp('M2_7_by_7_metaclustering_elbow_curve.bmp', 
    width = 640, 
    height = 400)

meta_clustering <- 
  FlowSOM::MetaClustering(
    data = m2_som$map$codes, 
    method = "metaClustering_consensus", 
    max = 40, 
    seed = 497, 
    plot = TRUE)


dev.off()

saveRDS(meta_clustering, 
        'M2_7_by_7_metaclustering_results.Rds')


cat('\n\nDetermined number of metaclusters is: ', FlowSOM::NMetaclusters(m2_som))

saveRDS(m2_som, 'M2_7_by_7_SOM_empirical_num_clusters.Rds')