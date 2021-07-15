

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
    

# Gather patient info and metadata for plot coloring
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

# Define point labels

labels <-
  sapply(flow_ids, function(x) rep(x, 100e3))

patient_and_meta_info_df <- 
  patient_and_meta_info_df %>% 
  mutate(cytometer = 
           case_when(
             cyt_num == "V657338000021" ~ "Canto 4",
             cyt_num == "V657338000098" ~ "Canto 5",
             cyt_num == "V657338000099" ~ "Canto 6",
             cyt_num == "R658222R1012" ~ "Fortessa 1", 
             cyt_num == "R65822R1018"  ~ "Fortessa 2", 
             cyt_num == "R66093700072" ~ "Symphony 1", 
             cyt_num == "R66093700081" ~ "Symphony 2", 
             cyt_num == "R66093700082" ~ "Symphony 3", 
             
           ))


# Get patient meta info df in same order as points
patient_and_meta_info_df <- 
  patient_and_meta_info_df %>% 
  set_rownames(
    patient_and_meta_info_df$Accession.Number
  ) %>%  
  .[flow_ids, ]

labels <- 
  patient_and_meta_info_df %>% 
  .$cytometer %>% 
  unlist() %>% 
  lapply(function(x) rep(x, 40e3)) %>% 
  unlist() 


levels <- 
  patient_and_meta_info_df$cytometer %>% 
  unlist() %>% 
  unique() %>% 
  sort() 

labels <- factor(labels, levels)


colors <-
  RColorBrewer::brewer.pal(length(unique(labels)), 'Set1')

results <- 
  as.matrix(results)

bmp('images/umap_96_patients_100k_cells.bmp', 
    width = 1280*1.2, 
    height = 800*1.2)

plot_umap(results,
          labels,
          colors)

dev.off()


#############################################################
# Now lets color them by which cluster they fell into in the 
# Self organzing map 

# I can't ensure that the order of the flow directories in the 
# UMAP was the same as the order of the flow directories in the 
# SOM matrix, so the saftest way to do this would be to read in the 
# data from the flow directories in the order that they were 
# included in the UMAP data matrix, and freshly classify the points
# on the SOM. 

# Actually, since the subsampled data on the server is not the same as 
# the subsampled data on my local machine, I will first have to download
# that subsampled data and just build the SOM again

# flow_dirs_root <- 
#   "/Users/sauterj1/Documents/Woodlist/Anon_Flow_Folders_Server_Version"
# 
# flow_directories <- 
#   list.files(flow_dirs_root, 
#              full.names = TRUE)
# 
# 
# m1_fcs_data_list <- 
#   vector('list', length(flow_directories))
# 
# for (i in seq_along(flow_directories)) {
#   
#   flow_directory <- flow_directories[[i]]
#   
#   cat('Reading fcs file:', i, 'of', 
#       length(flow_directories), '\n')
#                                     
#   fcs_file <- file.path(flow_directory, 
#                         'processed_M1_subsampled_100k.fcs')
#   
#   m1_fcs_data_list[[i]] <- 
#     read.FCS(fcs_file, 
#              transformation = FALSE) %>% 
#     exprs()
#   
# }
# 
# m1_fcs_data_mat <- 
#   m1_fcs_data_list %>% 
#   do.call(rbind, .)
# 
# m1_fcs_ff <- 
#   flowFrame(m1_fcs_data_mat)
# 
# m1_som <- FlowSOM::FlowSOM(m1_fcs_ff,
#                            xdim = 15, 
#                            ydim = 15, 
#                            nClus = 40, 
#                            colsToUse = seq_len(ncol(m1_fcs_data_list[[1]])))


m1_som <- readRDS('data/m1_som_100_cells_per_patient.Rds')

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

colors <- c(scales::alpha('grey', 0.05), 'blue')


bmp('images/umap_96_patients_colored_by_SOM_select_clusters.bmp', 
    width = 1280*1.2, 
    height = 800*1.2)

par(mfrow = c(8, 5))


plot_results <- 
  as.matrix(results)

for (label in levels(cells_to_metacluster_mapping)) {
  cat('Plotting UMAP for SOM cluster: ', label, '\n')
  
  labels <- cells_to_metacluster_mapping
  labels <- ifelse(labels == label, "Included", "Not Included")
  labels <- factor(labels, levels = c("Not Included", "Included"))
  
  plot_umap(plot_results[1:100, ],
            labels,
            colors)
}

dev.off()

