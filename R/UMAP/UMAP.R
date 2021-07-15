#!/bin/Rscript


### Input data
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Anon_Flow_Folders/"

library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(flowCore)
library(flowStats)
library(patchwork)
library(lubridate)
library(RColorBrewer)
library(knitr)
library(umap)

filter <- dplyr::filter


#########################################################
# Next Step: 
# In order to get the most value of a UMAP projection, 
# we should choose our samples wisely. Firstly, 
# all of our samples should have the same diagnosis.
# then we should choose our samples evenly over the
# machines that they were recorded on, ideally selecting
# 2-3 samples per machine. We should also keep in mind 
# the batches / change in configurations that we have
# pointed out before when selecting by machine / batch
#
# First lets see what distribution of patients / machines
# we have that are already currently parsed, then if 
# we need to expand the set to woodlist 2 patient set
# we can convert some files manually if needed
#########################################################
patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE)

flow_ids_processed <- 
  list.files('~/Documents/Woodlist/Anon_Flow_Folders/')


patient_df <- 
  patient_df %>% 
  filter(Exclusion.Reason == "Not Excluded", 
         !str_detect(Flow.Data.Folder, ','), 
         !str_detect(M1_path, ','), 
         !str_detect(M2_path, ',')) 


patient_df <- 
  patient_df %>% 
  filter(Accession.Number %in% flow_ids_processed)


# Need to pick 2 patients per machine within the same
# back for the machine with all patients having the
# same disease

m1_meta_info_df <- 
  readRDS('/Users/sauterj1/Documents/Patient_Folder_Analysis/data/m1_meta_info_df.Rds')


patient_and_meta_df <- 
  dplyr::left_join(patient_df, m1_meta_info_df, 
                    by = c('Accession.Number' = 'flow_id'))


df <- 
  patient_and_meta_df %>% 
  filter(CbioPortal.Cancer.Type.Detailed == 
           "Acute Myeloid Leukemia") %>% 
  group_by(cyt_num) %>% 
  slice_head(n = 2) %>% 
  ungroup()




flow_directories <- 
  paste0('/Users/sauterj1/Documents/Woodlist/Anon_Flow_Folders/', 
         df$Accession.Number)

flow_directories <- 
  list.files('/Users/sauterj1/Documents/Woodlist/Anon_Flow_Folders/', 
             full.names = TRUE)


fcs_dat_list <-
  vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  cat('Reading fcs file:', i, 'of', 
        length(flow_directories), '\n')
  
  m1_fcs_path <- 
    file.path(flow_directories[[i]], 
              'processed_M1_subsampled.fcs')
  
  fcs_dat_list[[i]] <- read.FCS(m1_fcs_path,
                                transformation = FALSE)
  
}


data_mat <- 
  lapply(fcs_dat_list, 
       function(x) exprs(x)) %>% 
  do.call(rbind, .)

tic <- Sys.time()


library(reticulate)
reticulate::use_python('/Users/sauterj1/anaconda3/bin/python3.8', 
                       required = TRUE)

fcs_umap <- umap::umap(data_mat, method = 'umap-learn')

saveRDS(fcs_umap, 
        '/Users/sauterj1/Documents/Patient_Folder_Analysis/data/fcs_umap.Rds')

toc <- Sys.time()

timediff <- toc - tic

cat('UMAP took: ')

print(timediff)




if (FALSE) {
  
  fcs_umap <- readRDS('/Users/sauterj1/Documents/Patient_Folder_Analysis/data/fcs_umap.Rds')
  
  
  plot_umap <- function(points,
                        labels,
                        colors) {
  
    pad=0.1
    cex=0.6
    pch=19
  
    xylim = range(points)
    xylim = xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  
    
    trans_colors <- scales::alpha(colors, alpha = 0.85)
    
    points(points[,1], points[,2], col=trans_colors[as.integer(labels)],
           cex=0.01, pch=pch)
  
    labels.u = levels(labels)
    legend.pos = "topleft"
    legend.text = as.character(labels.u)
  
    legend(legend.pos, legend=legend.text, inset=0.03,
           col=colors[seq_along(labels.u)],
           bty="n", pch=pch, cex = 1.5, pt.cex = 3.2)
  
  }
  
  
  labels <-
    sapply(1:9, function(x) rep(x, 40e3))
  
  df <- 
    df %>% 
    mutate(cytometer = 
      case_when(
        cyt_num == "R658222R1012" ~ "Fortessa 1", 
        cyt_num == "R65822R1018"  ~ "Fortessa 2", 
        cyt_num == "R66093700072" ~ "Symphony 1", 
        cyt_num == "R66093700081" ~ "Symphony 2", 
        cyt_num == "R66093700082" ~ "Symphony 3", 
        
      ))
  
  
  df <- 
    df %>% mutate(cyt_and_date = paste0(cytometer, '_', date)) 
  
  labels <- 
    df %>% 
    .$cyt_and_date %>% 
    unlist() %>% 
    sapply(function(x) rep(x, 40e3))
  
  levels <- 
    df$cyt_and_date %>% 
    unlist() %>% 
    unique() %>% 
    sort() 
  
  labels <- factor(labels, levels)
  
  
  colors <-
    RColorBrewer::brewer.pal(length(unique(labels)), 'Set1')
  
  
  bmp('images/umap.bmp', 
      width = 1280*1.2, 
      height = 800*1.2)
  
  plot_umap(fcs_umap$layout,
            labels,
            colors)
  
  dev.off()
    
  
}



if (FALSE) {
  
  library(FlowSOM)
  
  ff_data_mat <- 
    flowFrame(data_mat)
  
  # Next George's idea was to color each 
  m1_som <- readRDS('data/m1_som.Rds')
  
  new_data_som <- 
    FlowSOM::NewData(m1_som, ff_data_mat)
  
  meta_clusters <- 
    new_data_som$metaclustering
  
  data_clustering <- 
    GetClusters(new_data_som)
  
  data_meta_clustering <- 
    data_clustering %>% 
    meta_clusters[.]
  
  labels <- 
    df %>% 
    .$cyt_and_date %>% 
    unlist() %>% 
    sapply(function(x) rep(x, 40e3))
  
  levels <- 
    df$cyt_and_date %>% 
    unlist() %>% 
    unique() %>% 
    sort() 
  
  labels <- factor(labels, levels)
  
  
  colors <-
    RColorBrewer::brewer.pal(length(unique(labels)), 'Set1')
  
  
  bmp('images/umap.bmp', 
      width = 1280*1.2, 
      height = 800*1.2)
  
  plot_umap(fcs_umap$layout,
            labels,
            colors)
  
  dev.off()

  
  
}



bmp('images/fsom.bmp', 
    width = 1280*1.2, 
    height = 800*1.2)



