library(Seurat)
library(future) 

plan("multiprocess", 
     workers = 10)

seurat_object <- 
  readRDS('umap_seurat_object_with_metaclusters.Rds')


seurat_object <- 
  FindNeighbors(seurat_object, 
                features = rownames(seurat_object), 
                dims = NULL)

saveRDS(seurat_object, 'seurat_object_after_find_neighbors.RDS')

seurat_object <- 
  FindClusters(seurat_object, 
               resolution = 0.8)

saveRDS(seurat_object, 'seurat_object_after_find_clusters.RDS')

