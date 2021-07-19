library(dplyr)
library(Seurat)
library(magrittr)
library(patchwork)

filter <- dplyr::filter

seurat_object <- 
  readRDS('data/umap_seurat_object.Rds')

flow_ids <- seurat_object[[]]$flow_dir

# Read in patient metadata, and correlate all
# interested variables with metadata of each
# cell based on the flow_id / flow_dir metadata
patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE) %>%
  filter(Accession.Number %in% flow_ids)

m1_meta_info_df <- 
  readRDS('/Users/sauterj1/Documents/Patient_Folder_Analysis/data/m1_meta_info_df.Rds') %>% 
  filter(flow_id %in% flow_ids)

patient_and_meta_info_df <- 
  dplyr::left_join(patient_df, m1_meta_info_df, 
                   by = c('Accession.Number' = 'flow_id'))

rownames(patient_and_meta_info_df) <- 
  patient_and_meta_info_df$Accession.Number


# Add Date, Machine, Course Disease and Fine Disease to
# Seurat metadata
seurat_object@meta.data$cancer_type <-
  seurat_object@meta.data$flow_dir %>% 
  patient_and_meta_info_df[., "CbioPortal.Cancer.Type"]

seurat_object@meta.data$cancer_type_detailed <- 
  seurat_object@meta.data$flow_dir %>% 
  patient_and_meta_info_df[., "CbioPortal.Cancer.Type.Detailed"] %>% 
  unlist()

seurat_object@meta.data$flow_cytometer <-
  seurat_object@meta.data$flow_dir %>% 
  patient_and_meta_info_df[., "cyt_num"] %>% 
  unlist()

seurat_object@meta.data$experiment_date <- 
  seurat_object@meta.data$flow_dir %>% 
  patient_and_meta_info_df[., "date"] %>% 
  unlist()


bmp('images/seurat_umap.bmp', 
    width = 1280, 
    height = 800)

# Plot 
Seurat::DimPlot(
  seurat_object, 
  reduction = 'umap', 
  pt.size = 0.1, 
  group.by = 'flow_cytometer')

dev.off()

# I think for the UMAP I might stick to 
# my manual way of plotting

# Define function to plot UMAP results
plot_umap <- function(points,
                      labels,
                      colors, 
                      alpha = 0.85, 
                      cex=0.01) {
  
  pad=0.1
  pch=19
  
  xylim = range(points)
  xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  
  xlim <- xylim
  # xlim[1] <- xlim[1] - .5 * abs(xlim[1])
  
  par(mar=c(0.2,0.7,1.2,0.7), ps=10)
  plot(xlim, xylim, type="n", axes=F, frame=F)
  rect(xlim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  
  
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

umap_points <- 
  seurat_object@reductions$umap@cell.embeddings

labels <- seurat_object@meta.data$flow_cytometer

labels <- 
  dplyr::case_when(labels == "R658222R1012" ~ "Fortessa I", 
                   labels == "R65822R1018" ~  "Fortessa II",
                   labels == "R66093700072" ~ "Symphony I",
                   labels == "R66093700081" ~ "Symphony II",
                   labels == "R66093700082" ~ "Symphony III",
                   labels == "V657338000098" ~ 'Canto 5')

include_point <- 
  labels != "Symphony I" & labels != 'Canto 5'

umap_points <- umap_points[include_point, ]

labels <- labels[include_point]
labels <- factor(labels)

colors <-
  c("#fe00b1",
   "#01f2a5",
   "#006fa6",
   "#ffe6fc")

bmp('images/seurat_umap_colored_by_cytometer_no_symph_1_or_canto_5.bmp', 
    width = 1280, 
    height = 800)


plot_umap(umap_points, 
          labels = labels, 
          colors = colors)

dev.off()


bmp('images/seurat_feature_plot.bmp', 
    width = 1280, 
    height = 800)


# Plotting by patient disease
labels <- seurat_object@meta.data$cancer_type_detailed

labels <- factor(labels)

colors <-
  c("#296600",
   "#6738e4",
   "#80f934",
   "#64009b",
   "#6eff7c",
   "#d10083",
   "#01b81f",
   "#5f0068",
   "#f7e800",
   "#ac90ff",
   "#9bd100",
   "#00396a",
   "#01e378",
   "#ff4e6c",
   "#57ffa1",
   "#eb5200",
   "#02dcfb",
   "#ffb057",
   "#0088c0",
   "#fffe93",
   "#290400",
   "#a2ffda",
   "#003332",
   "#ffaee3",
   "#1f3f00",
   "#a7b4ff",
   "#785b00",
   "#aaf3ff",
   "#ff968c",
   "#005746")

bmp('images/seurat_umap_colored_by_minor_disease.bmp', 
    width = 1280, 
    height = 800)


plot_umap(umap_points, 
          labels = labels, 
          colors = colors)

dev.off()


bmp('images/seurat_feature_plot.bmp', 
    width = 1280, 
    height = 800)

FeaturePlot(seurat_object, 
            features = rownames(seurat_object))

dev.off()
