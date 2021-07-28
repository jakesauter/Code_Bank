library(Seurat)
library(reticulate)

reticulate::use_python(
  '/gpfs/mskmind_ess/sauterj1/miniconda3/envs/Seurat/bin/python3', 
  required = TRUE)

# scanpy: https://scanpy.readthedocs.io/en/stable/
# https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html
# https://github.com/LuckyMD/Code_snippets/blob/master/Seurat_to_anndata.ipynb
scanpy <- reticulate::import('scanpy')
ad <- import("anndata", convert = FALSE)

seurat_object <- 
  readRDS('/gpfs/mskmind_ess/sauterj1/seurat/umap_seurat_object_with_7_by_7_SOM_metaclusters.Rds')

seurat_ad <- scanpy$AnnData(
  X   = t(as.matrix(GetAssayData(seurat_object)))
  # ,
  # obs = seurat_object[[]],
  # var = t(as.matrix(GetAssay(seurat_object)[[]]))
)

# make scanpy object? 
neighbors_result <- 
  scanpy$pp$neighbors(
  seurat_ad, 
  use_rep=reticulate::py_none(), 
  n_pcs=0, 
  method='umap', 
  n_neighbors = 30
)

saveRDS(seurat_object, 'seurat_object_after_find_clusters.RDS')

