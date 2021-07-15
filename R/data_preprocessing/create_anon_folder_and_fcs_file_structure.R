#!/bin/Rscript

### Input data
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Correct_QC_Flow_Folders/"

library(fs)
library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)

filter <- dplyr::filter

# Gather all QC Statistics
flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE)

anon_flow_dir_loc <- 
  path(flow_dir_save_loc, 
        '../', 
      'Anon_Flow_Folders')
  
if (!dir_exists(anon_flow_dir_loc)) {
  dir_create(anon_flow_dir_loc) 
}

anon_flow_dir_loc <- 
  path_real(anon_flow_dir_loc)

anon_flow_directories <- 
  flow_directories %>% 
  basename() %>% 
  str_extract('F[0-9]+-[0-9]+') %>% 
  path(anon_flow_dir_loc, .)

anonymize_fcs_header <- function(fcs) {
  description(fcs)$`$FIL` <- ""
  description(fcs)$`$SRC` <- ""
  description(fcs)$`EXPERIMENT NAME` <- ""
  description(fcs)$`FILENAME` <- ""
  
  return(fcs)
}

coalesce_fcs_colnames <- function(fcs_data, tube = 'M1') {
  
  if (tube == 'M1') {
    names_order <- c('FSC-A', 'FSC-H', 'SSC-A', 'SSC-H', 
                     'CD13', 'CD15', 'CD19', 'CD33', 'CD34', 
                     'CD38', 'CD45', 'CD71', 'CD117', 'HLA-DR')
  } else if (tube == 'M2') {
    names_order <- c('FSC-A', 'FSC-H', 'SSC-A', 'SSC-H', 
                     'CD11b', 'CD13', 'CD14', 'CD16', 'CD34', 
                     'CD38', 'CD45', 'CD64', 'CD123', 'HLA-DR') 
  } else {
    stop(paste0('TUBE NAME: ', tube, ' NOT DEFINED')) 
  }
    
    
  fcs_expr <- exprs(fcs_data)
  metadata <- pData(parameters(fcs_data))
  
  to_keep <- ifelse(is.na(metadata$desc), FALSE, TRUE)
  to_keep[1:4] <- TRUE
  
  names <- metadata$desc[to_keep]
  names[1:4] <- colnames(fcs_expr)[1:4]
  names <- str_extract(names, 
                       'FSC-A|FSC-H|SSC-A|SSC-H|CD[0-9]+b?|HLA')
  names[names == 'HLA'] = 'HLA-DR'
  
  new_expr <- fcs_expr[, to_keep] 
  colnames(new_expr) <- names
  
  # To ensure same order of column names.
  new_expr <- new_expr[, names_order] 

  if (ncol(fcs_expr) != length(to_keep)) {
    stop('DROPPED COLUMNS IN NAME ORDERING') 
  }
  
  return(flowFrame(new_expr))
}

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')] %>% 
    .[str_detect(., '_lin_and_asinh_transformed')] %>% 
    .[str_detect(., 'singlets')]
  
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
  
  cat('Reading FCS file: ', m1_fcs_file, '\n')
  
  m1_fcs_dat <- 
    flowCore::read.FCS(m1_fcs_file, 
                       transformation = FALSE)
  
  cat('Reading FCS file: ', m2_fcs_file, '\n')
  
  m2_fcs_dat <- 
    flowCore::read.FCS(m2_fcs_file, 
                       transformation = FALSE)
  
  
  # Anonymyze fcs header
  m1_anon_fcs_data <- 
    anonymize_fcs_header(m1_fcs_dat) %>% 
    coalesce_fcs_colnames(tube = 'M1')
  
  m2_anon_fcs_data <- 
    anonymize_fcs_header(m2_fcs_dat) %>% 
    coalesce_fcs_colnames(tube = 'M2')
  
  
  ## Saving
  if (!dir_exists(anon_flow_directories[[i]])) {
    dir_create(anon_flow_directories[[i]])
  }
  
  m1_anon_fcs_file <-
    path(anon_flow_directories[[i]], 
         "processed_M1.fcs")
  
  m2_anon_fcs_file <-
    path(anon_flow_directories[[i]], 
         "processed_M2.fcs")
  
  
  flowCore::write.FCS(m1_anon_fcs_data, 
                      m1_anon_fcs_file)
  
  flowCore::write.FCS(m2_anon_fcs_data, 
                      m2_anon_fcs_file)  
  
}
