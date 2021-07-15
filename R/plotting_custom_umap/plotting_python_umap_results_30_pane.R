
library(arrow)
library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(FlowSOM)
library(RColorBrewer)

filter <- dplyr::filter

# Read in UMAP results
results <- 
  arrow::read_feather('data/umap_results_100k_per_patient.feather')


# Gather the order of the flow ids in the UMAP results
flow_ids <- 
  read.table('data/filename_order_100k_per_patient.csv', 
             sep = ',') %>% 
  unlist() %>% 
  str_extract('F[0-9\\-]+$')

m1_som <- readRDS('data/m1_som_100_cells_per_patient_30_clus.Rds')



# Define function to plot UMAP results
plot_umap <- function(points,
                      labels,
                      colors, 
                      cex=0.01) {
  
  pad=0.1
  pch=19
  
  xylim = range(points)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  
  par(mar=c(0.2,0.7,1.2,0.7), ps=10)
  plot(xylim, xylim, type="n", axes=F, frame=F)
  rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  
  
  points(points[,1], points[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
}



cluster_to_metacluster_mapping <- 
  m1_som$metaclustering

cells_to_cluster_mapping <- 
  FlowSOM::GetClusters(m1_som)

cells_to_metacluster_mapping <- 
  cells_to_cluster_mapping %>% 
  cluster_to_metacluster_mapping[.]


n <- length(levels(labels))
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
colors <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colors <- c(scales::alpha('grey', 0.1), 'black')




plot_results <- 
  as.matrix(results)

bmp('images/umap_som_clusters/umap_96_patients_30_meta_clusters.bmp', 
    width = 1280*2, 
    height = 800*2)

par(mfrow = c(5, 6))

for (label in levels(cells_to_metacluster_mapping)) {
  
  cat('Plotting UMAP for SOM cluster: ', label, '\n')
  
  
  labels <- cells_to_metacluster_mapping
  labels <- ifelse(labels == label, "Included", "Not Included")
  labels <- factor(labels, levels = c("Not Included", "Included"))
  
  plot_umap(plot_results,
            labels,
            colors)
  
}

dev.off()
