library(dplyr)
library(Seurat)
library(FlowSOM)
library(flowCore)
library(magrittr)
library(patchwork)

filter <- dplyr::filter


# Need to run seurat/build_seurat_object_with_SOM_clusters.R
seurat_object <- 
  readRDS('data/umap_seurat_object_with_7_by_7_SOM_metaclusters.Rds')

umap_points <- 
  seurat_object@reductions$umap@cell.embeddings

louvain_clusters <- 
  arrow::read_feather('data/louvain_clustering_results.feather')

seurat_object@meta.data$louvain_clusters <- 
  louvain_clusters %>% 
  unlist()

# Define function to plot UMAP results
plot_umap <- function(points,
                      labels,
                      colors, 
                      alpha = 0.85, 
                      cex=0.01) {
  
  pad=0.2
  pch=19
  
  xylim = range(points)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  
  par(mar=c(0.2,0.7,1.2,0.7), ps=10)
  plot(xylim, xylim, type="n", axes=F, frame=F)
  rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  
  
  trans_colors <- scales::alpha(colors, alpha = alpha)
  
  points(points[,1], points[,2], col=trans_colors[as.integer(labels)],
         cex=cex, pch=pch)
}

labels <- seurat_object@meta.data$louvain_clusters
labels <- factor(labels)


colors <-
  c("#45286c",
    "#443248",
    "#592aae",
    "#64432e",
    "#3a5144",
    "#973a53",
    "#5f4faa",
    "#54692b",
    "#b73995",
    "#c7423f",
    "#606be5",
    "#d9472b",
    "#ad4be3",
    "#d5437f",
    "#bb6152",
    "#b17232",
    "#d84ac5",
    "#6186d5",
    "#ea5a17",
    "#bb74a5",
    "#ad74db",
    "#9a8fc1",
    "#57a0cf",
    "#c88c9b",
    "#77a3b6",
    "#cf85d0",
    "#dd9628",
    "#60b88c",
    "#bba6a4",
    "#75bd58",
    "#cfa971",
    "#b0b28a",
    "#bdb549",
    "#7acd49",
    "#59d087",
    "#53ccc5",
    "#a0c2b0",
    "#bacd84",
    "#86e63a",
    "#dad833")


# I think for the UMAP I might stick to 
# my manual way of plotting
bmp('images/seurat_umap_with_louvain_clusters.bmp', 
    width = 800, 
    height = 600)

plot_umap(umap_points, 
          labels = labels, 
          colors = colors, 
          alpha = 0.35)

dev.off()


# saveRDS(seurat_object, 'data/seurat_object_with_louvain_clusters.Rds')










