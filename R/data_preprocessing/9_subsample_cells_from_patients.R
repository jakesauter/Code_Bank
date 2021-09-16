#!/bin/Rscript

library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(umap)

filter <- dplyr::filter

flow_directories <- 
  list.files("/gpfs/mskmind_ess/sauterj1/Full_Cohort_Anon_Flow_Folders", 
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
    str_c('_subsampled_100k.fcs')
  
  m2_subsampled_fcs_file <-
    m2_fcs_file %>% 
    strsplit('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_subsampled_100k.fcs')
  
  if (!file.exists(m1_subsampled_fcs_file)) {
    m1_fcs_dat <- 
      read.FCS(m1_fcs_file, 
               transformation = FALSE)
    
    
    exprs(m1_fcs_dat) <- 
      exprs(m1_fcs_dat)[sample(1:nrow(m1_fcs_dat), 
                               100e3), ]
    
    write.FCS(m1_fcs_dat, 
              m1_subsampled_fcs_file)
    
  } 
  
  if (!file.exists(m2_subsampled_fcs_file)) {
    
    m2_fcs_dat <- 
      read.FCS(m2_fcs_file, 
               transformation = FALSE)
    
    exprs(m2_fcs_dat) <- 
      exprs(m2_fcs_dat)[sample(1:nrow(m2_fcs_dat), 
                               100e3), ]
    
    write.FCS(m2_fcs_dat, 
              m2_subsampled_fcs_file)
  } 
}








