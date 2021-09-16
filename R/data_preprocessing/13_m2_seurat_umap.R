
library(dplyr)
library(magrittr)
library(Seurat)
library(flowCore)
library(FlowSOM)
library(future)

filter <- dplyr::filter

seurat_object <- 
  readRDS('/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_centralized_Seurat_object.Rds')

# Make a multi-core plan
plan("multiprocess", workers = 30)

cat('\n\nRunning UMAP\n\n')

#  Run umap
seurat_object <-
  RunUMAP(seurat_object,
          features = rownames(seurat_object))

saveRDS(seurat_object, 
        '/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_seurat_object_with_umap.R')


