#!/bin/Rscript

flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Flow_Folders"


library(dplyr)
library(flowAI)
library(stringr)
library(magrittr)
library(flowCore)
library(flowStats)

filter <- dplyr::filter


#====================================================
# Gathering all QC Statistics
#----------------------------------------------------

flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE)

flow_directories <- sample(flow_directories, 10)

qc_results_list <- vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  
  # Parse and better save QC Results
  qc_results <- 
    read.table(file.path(flow_directory, 
                         'flowAI_QC_Results',
                         'QCmini.txt'), 
               skip = 1, 
               sep = '\t') %>% 
    .[, c(1, 3, 5, 6, 7)] %>% 
    set_colnames(c('filename', 
                   'total_filtered_perc',
                   'FR_filtered_perc', 
                   'FS_filtered_perc', 
                   'FM_filtered_perc (not used)'))
  
  qc_results$filename <-
    file.path(flow_directory, qc_results$filename)
  
  
  qc_results_list[[i]] <- qc_results
  
}

full_qc_results <- do.call(rbind, qc_results_list)

full_qc_results <- 
  full_qc_results %>% 
  dplyr::distinct() %>%
  arrange(total_filtered_perc)
#====================================================



#====================================================
# Reading in High Quality M1 files 
#---------------------------------------------------

flow_dirs_passed_qc <- 
  full_qc_results %>% 
  filter(total_filtered_perc < 80) %>% 
  mutate(dirname = dirname(filename)) %>% 
  count(dirname) %>% 
  filter(n > 1) %>% 
  .$dirname

flow_directories <- 
  flow_dirs_passed_qc

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')]
  
  m1_fcs_file <- 
    str_detect(fcs_files, 
               'M1|m1') %>% 
    fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  m2_fcs_file <- 
    str_detect(fcs_files, 
               'M2|m2') %>% 
    fcs_files[.] %>% 
    file.path(flow_directory, .)
  

  m1_fcs_dat <- 
    flowCore::read.FCS(m1_fcs_file, 
                       transformation = FALSE)
  
  m2_fcs_dat <- 
    flowCore::read.FCS(m2_fcs_file, 
                       transformation = FALSE)
  
  # Transform
  m1_meta <- pData(parameters(m1_fcs_dat))
  m1_param_names <- as.character(m1_meta$name)
  m1_cols_to_transform <- m1_param_names[c(1:4, which(!is.na(m1_meta$desc)))]
  m1_cols_to_transform <- unname(m1_cols_to_transform)
  
  m2_meta <- pData(parameters(m2_fcs_dat))
  m2_param_names <- as.character(m2_meta$name)
  m2_cols_to_transform <- m2_param_names[c(1:4, which(!is.na(m2_meta$desc)))]
  m2_cols_to_transform <- unname(m2_cols_to_transform)
  
  
  exprs(m1_fcs_dat) <- exprs(m1_fcs_dat)[, c(m1_cols_to_transform, 'Time')]
  exprs(m2_fcs_dat) <- exprs(m2_fcs_dat)[, c(m2_cols_to_transform, 'Time')]
  
  
  m1_fcs_dat <- 
    flowTrans::flowTrans(m1_fcs_dat,
                         'mclMultivArcSinh',
                         m1_cols_to_transform, 
                         FALSE,
                         FALSE)$result
  
  m2_fcs_dat <- 
    flowTrans::flowTrans(m2_fcs_dat,
                         'mclMultivArcSinh',
                         m2_cols_to_transform, 
                         FALSE,
                         FALSE)$result
  
  # Save
  m1_transformed_fcs_file <- 
    m1_fcs_file %>% 
    str_split('.fcs') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_and_transformed.fcs')  
  
  m2_transformed_fcs_file <- 
    m2_fcs_file %>% 
    str_split('.fcs') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_and_transformed.fcs')  
  
  flowCore::write.FCS(m1_fcs_dat, 
                      m1_transformed_fcs_file)
  
  flowCore::write.FCS(m2_fcs_dat, 
                      m2_transformed_fcs_file)
}










#====================================================