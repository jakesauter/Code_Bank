library(Seurat)
library(magrittr)
library(patchwork)

seurat_object <- 
  readRDS('data/umap_seurat_object.Rds')

# Read in patient metadata, and correlate all
# interested variables with metadata of each
# cell based on the flow_id / flow_dir metadata



Seurat::DimPlot(
  seurat_object, 
  reduction = 'umap')


seurat_object[[]]$flow_dir %>% unique()

