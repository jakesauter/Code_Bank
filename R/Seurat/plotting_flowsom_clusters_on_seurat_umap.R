library(dplyr)
library(Seurat)
library(FlowSOM)
library(flowCore)
library(magrittr)
library(patchwork)

filter <- dplyr::filter


# Need to run seurat/build_seurat_object_with_SOM_clusters.R
seurat_object <- 
  readRDS('data/umap_seurat_object_with_10_by_10_SOM_metaclusters.Rds')

umap_points <- 
  seurat_object@reductions$umap@cell.embeddings

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
  
  labels.u = levels(labels)
  legend.pos = "topleft"
  legend.text = as.character(labels.u)
  
  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[seq_along(labels.u)],
         bty="n", pch=pch, cex = 1.5, pt.cex = 3.2)
  
}

labels <- seurat_object@meta.data$metacluster
labels <- factor(labels)


# 7x7 color palette
# colors <- 
#   c("#e07200",
#      "#4a80ff",
#      "#ae9500",
#      "#002e9f",
#      "#00bf67",
#      "#820093",
#      "#4c7400",
#      "#ff6fe7",
#      "#b20500",
#      "#004584",
#      "#fa8d69",
#      "#e293bd",
#      "#6e2107")

colors <- 
  c("#a60051",
     "#eac302",
     "#0080fc",
     "#ff8c0d",
     "#003e9b",
     "#d74800",
     "#a585ff",
     "#3a6100",
     "#ff69ff",
     "#e3c376",
     "#ea019b",
     "#00726c",
     "#ffa5fb",
     "#670006",
     "#015e8c",
     "#efb7c7")


# I think for the UMAP I might stick to 
# my manual way of plotting
bmp('images/seurat_umap_with_10_by_10_SOM_metaclusters.bmp', 
    width = 800, 
    height = 500)

plot_umap(umap_points, 
          labels = labels, 
          colors = colors, 
          alpha = 0.35)

dev.off()

