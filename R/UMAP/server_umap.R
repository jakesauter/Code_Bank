#!/bin/Rscript

library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
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

flow_directories <- 
  list.files("/home/sauterj1/shared_data_folder/sauterj1/Anon_Flow_Folders", 
             full.names = TRUE)


fcs_dat_list <-
  vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  cat('Processing directory:', i, 'of:', length(flow_directories), '\n')
  
  flow_directory <- flow_directories[i]
  
  cat('Flow Directory: ', flow_directory, '\n\n')
  
  
  m1_fcs_file <- file.path(flow_directory, 'processed_M1.fcs')
  m2_fcs_file <- file.path(flow_directory, 'processed_M2.fcs')
  
  m1_subsampled_fcs_file <-
    m1_fcs_file %>% 
    strsplit('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_subsampled.fcs')
  
  m2_subsampled_fcs_file <-
    m2_fcs_file %>% 
    strsplit('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_subsampled.fcs')
  
  if (!file.exists(m1_subsampled_fcs_file)) {
    m1_fcs_dat <- 
      read.FCS(m1_fcs_file, 
               transformation = FALSE)
    
    
    exprs(m1_fcs_dat) <- 
      exprs(m1_fcs_dat)[sample(1:nrow(m1_fcs_dat), 
                               40e3), ]
    
    write.FCS(m1_fcs_dat, 
              m1_subsampled_fcs_file)
    
  } 
  
  if (!file.exists(m2_subsampled_fcs_file)) {
    
    m2_fcs_dat <- 
      read.FCS(m2_fcs_file, 
               transformation = FALSE)
    
    exprs(m2_fcs_dat) <- 
      exprs(m2_fcs_dat)[sample(1:nrow(m2_fcs_dat), 
                               40e3), ]
    
    write.FCS(m2_fcs_dat, 
              m2_subsampled_fcs_file)
  } 
  
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

fcs_umap <- umap::umap(data_mat)

saveRDS(fcs_umap, 'fcs_umap_subsampled_40k.Rds')

toc <- Sys.time()

timediff <- toc - tic

cat('UMAP took: ', timediff, ' to generate\n')









