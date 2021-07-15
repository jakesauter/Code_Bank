#!/bin/Rscript

flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Correct_QC_Flow_Folders/"


library(dplyr)
library(flowAI)
library(stringr)
library(magrittr)
library(flowCore)
library(flowStats)
library(flowTrans)

filter <- dplyr::filter

#====================================================
# Reading in High Quality M1 files 
#---------------------------------------------------

flow_directories <- 
  list.dirs(flow_dir_save_loc, recursive = FALSE)

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")

  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')] %>% 
    .[!str_detect(., 'transformed')]
  
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
  
  cat('Reading FCS File: ', m1_fcs_file, '\n')
  
  
  m1_fcs_dat <- 
    flowCore::read.FCS(m1_fcs_file, 
                       transformation = FALSE)
  
  cat('Reading FCS File: ', m2_fcs_file, '\n')
  
  
  m2_fcs_dat <- 
    flowCore::read.FCS(m2_fcs_file, 
                       transformation = FALSE)
  
  # Transform
  m1_meta <- pData(parameters(m1_fcs_dat))
  m1_param_names <- as.character(m1_meta$name)
  m1_cols_to_transform <- m1_param_names[which(!is.na(m1_meta$desc))]
  m1_cols_to_transform <- unname(m1_cols_to_transform)
  m1_cols_to_transform <- c('SSC-A', 'SSC-H', m1_cols_to_transform)
  
  m2_meta <- pData(parameters(m2_fcs_dat))
  m2_param_names <- as.character(m2_meta$name)
  m2_cols_to_transform <- m2_param_names[which(!is.na(m2_meta$desc))]
  m2_cols_to_transform <- unname(m2_cols_to_transform)
  m2_cols_to_transform <- c('SSC-A', 'SSC-H', m2_cols_to_transform)
  
  exprs(m1_fcs_dat) <- exprs(m1_fcs_dat)[, c('FSC-A', 'FSC-H', 
                                             m1_cols_to_transform)]
  rownames(parameters(m1_fcs_dat)) <- paste0('$P', seq_len(ncol(exprs(m1_fcs_dat))))
  
  exprs(m2_fcs_dat) <- exprs(m2_fcs_dat)[, c('FSC-A', 'FSC-H', 
                                             m2_cols_to_transform)]
  rownames(parameters(m2_fcs_dat)) <- paste0('$P', seq_len(ncol(exprs(m2_fcs_dat))))
  

  tl <-
    flowCore::transformList(from = m1_cols_to_transform,
                            tfun = function(x) asinh(x/150),
                            transformationId = 'M1')
  m1_fcs_dat <-
    flowCore::transform(m1_fcs_dat, tl)
  
  exprs(m1_fcs_dat)[, 'FSC-A'] <- 
    exprs(m1_fcs_dat)[, 'FSC-A'] %>% 
    {. / max(.) * 8.159}
  
  exprs(m1_fcs_dat)[, 'FSC-H'] <- 
    exprs(m1_fcs_dat)[, 'FSC-H'] %>% 
    {. / max(.) * 8.159}
  
  
  tl <-
    flowCore::transformList(from = m2_cols_to_transform,
                            tfun = function(x) asinh(x/150),
                            transformationId = 'M2')
  m2_fcs_dat <-
    flowCore::transform(m2_fcs_dat, tl)
  
  exprs(m2_fcs_dat)[, 'FSC-A'] <- 
    exprs(m2_fcs_dat)[, 'FSC-A'] %>% 
    {. / max(.) * 8.159}
  
  exprs(m2_fcs_dat)[, 'FSC-H'] <- 
    exprs(m2_fcs_dat)[, 'FSC-H'] %>% 
    {. / max(.) * 8.159}
  
  
  # Save
  m1_transformed_fcs_file <-
    m1_fcs_file %>%
    str_split('.fcs') %>%
    .[[1]] %>% .[1] %>%
    str_c('_lin_and_asinh_transformed.fcs')

  m2_transformed_fcs_file <-
    m2_fcs_file %>%
    str_split('.fcs') %>%
    .[[1]] %>% .[1] %>%
    str_c('_lin_and_asinh_transformed.fcs')
  
  flowCore::write.FCS(m1_fcs_dat, 
                      m1_transformed_fcs_file)
  
  flowCore::write.FCS(m2_fcs_dat, 
                      m2_transformed_fcs_file)

}








#====================================================