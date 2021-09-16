
library(dplyr)
library(magrittr)
library(Seurat)
library(flowCore)
library(FlowSOM)

filter <- dplyr::filter

seurat_object <- 
  readRDS('/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_centralized_Seurat_object.Rds')

m1_fcs_data_mat <-
  seurat_object@assays$MultiParamFlowCyto@data %>% 
  as.matrix()

m1_fcs_ff <-
  flowFrame(m1_fcs_data_mat)

m1_som <- FlowSOM::FlowSOM(m1_fcs_ff,
                           xdim = 7,
                           ydim = 7,
                           maxMeta = 40,
                           colsToUse = seq_len(ncol(m1_fcs_data_list[[1]])))

bmp('M1_7_by_7_metaclustering_elbow_curve.bmp', 
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
        'M1_7_by_7_metaclustering_results.Rds')


cat('\n\nDetermined number of metaclusters is: ', FlowSOM::NMetaclusters(m1_som))

saveRDS(m1_som, 'M1_7_by_7_SOM_empirical_num_clusters.Rds')