library(dplyr)
library(Seurat)
library(future)
library(flowCore)
library(patchwork)

# First Need to get the cells x markers matrix
flow_directories <- 
  list.files('/home/sauterj1/shared_data_folder/sauterj1/Anon_Flow_Folders', 
             full.names = TRUE)

fcs_dat_list <-
  vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  cat('Reading fcs file:', i, 'of', 
      length(flow_directories), '\n')
  
  m1_fcs_path <- 
    file.path(flow_directories[[i]], 
              'processed_M1_subsampled_100k.fcs')
  
  fcs_dat <- 
    exprs(read.FCS(m1_fcs_path,
                      transformation = FALSE))
  
  rownames(fcs_dat) <- 
    rep(basename(flow_directories[i]), 
        nrow(fcs_dat))
  
  fcs_dat_list[[i]] <- fcs_dat
  
}

cat('\n\nCombining data into single matrix\n\n')
fcs_data <- 
  do.call(rbind, fcs_dat_list)

# Seurat expects cells as columns and
# features as rows
fcs_data <- t(fcs_data)

# Create the metadata matrix
meta.data <- 
  data.frame(flow_dir = colnames(fcs_data))

# Each cell needs a uniqe name (lets make it the index)
colnames(fcs_data) <- seq(1, ncol(fcs_data))

cat('\n\nCreating Seurat object\n\n')
# Create seurat object
seurat_obj <- 
  CreateSeuratObject(
    counts = fcs_data,
    project = 'FCS Example',
    assay = 'MultiParamFlowCyto', 
    meta.data = meta.data)

# Make a multi-core plan
plan("multiprocess", workers = 20)

cat('\n\nRunning UMAP\n\n')
#  Run umap
seurat_obj <- 
  RunUMAP(seurat_obj, 
          features = rownames(seurat_obj))


cat('\n\nSaving Seurat object\n\n')
# Save Seurat object
saveRDS(seurat_obj, 'umap_seurat_object.Rds')



