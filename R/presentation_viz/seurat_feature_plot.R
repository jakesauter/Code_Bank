
library(dplyr)
library(Seurat)
library(magrittr)
library(patchwork)

filter <- dplyr::filter

seurat_object <- 
  readRDS('data/umap_seurat_object_with_7_by_7_SOM_metaclusters.Rds')


bmp('images/seurat_feature_plot.bmp',
    width = 400,
    height = 1000)

p <- 
  FeaturePlot(seurat_object, 
            features = rownames(seurat_object))

p + 
  scale_colour_gradientn(
    colours = rev(brewer.pal(n = 11, name = "RdBu")))

p

dev.off()