
library(magrittr)
library(flowCore)
library(FlowSOM)

filter <- dplyr::filter

m1_som <- 
  readRDS('m1_som_100_cells_per_patient_Flowsom_chosen_n_clusters.Rds')


bmp('Metaclustering_results.bmp', 
    width = 1280, 
    height = 800)

m1_som <- 
  FlowSOM::MetaClustering(
              data = m1_som$map$codes, 
              method = "metaClustering_consensus", 
              max = 50, 
              seed = 497, 
              plot = TRUE)

saveRDS(m1_som, '15_by_15_SOM_chosen_clusters.Rds')

dev.off()