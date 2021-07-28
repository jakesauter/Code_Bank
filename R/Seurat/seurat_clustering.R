library(Seurat)
library(future) 


# TODO: scanpy: https://scanpy.readthedocs.io/en/stable/
# https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
# https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb

plan("multiprocess", 
     workers = 10)

seurat_object <- 
  readRDS('data/umap_seurat_object_with_7_by_7_SOM_metaclusters.Rds')

# seurat_expression_data <- 
#   t(as.matrix(GetAssayData(seurat_object))) %>% 
#   as.data.frame()
# 
# arrow::write_feather(seurat_expression_data, 
#                      'seurat_expression_data.feather')

Seurat:::FindNeighbors.Seurat
Seurat:::FindNeighbors.default

AnnoySearch <- Seurat:::AnnoySearch
debug(AnnoySearch)

undebug(Seurat:::NNHelper) 

undebug(Seurat:::AnnoyNN)
debug(Seurat:::AnnoyBuildIndex)

seurat_object <- 
  FindNeighbors(seurat_object, 
                features = rownames(seurat_object), 
                dims = NULL)

saveRDS(seurat_object, 'seurat_object_after_find_neighbors.RDS')

seurat_object <- 
  FindClusters(seurat_object, 
               resolution = 0.8)

saveRDS(seurat_object, 'seurat_object_after_find_clusters.RDS')

