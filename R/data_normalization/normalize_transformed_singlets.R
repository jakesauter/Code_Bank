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

# flow_directories <- sample(flow_directories, 25)


m1_fcs_list <- vector('list', length(flow_directories))
m2_fcs_list <- vector('list', length(flow_directories))

names(m1_fcs_list) <- flow_directories
names(m2_fcs_list) <- flow_directories

m1_fcs_input_files <- vector('character', length(flow_directories))
m2_fcs_input_files <- vector('character', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")

  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')] %>% 
    .[str_detect(., 'asinh_transformed')] %>% 
    .[str_detect(., 'singlets')] %>% 
    .[!str_detect(., 'normalized')]
  
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
  
  
  m1_fcs_list[[i]] <- 
    flowCore::read.FCS(m1_fcs_file, 
                       transformation = FALSE)
  
  cat('Reading FCS File: ', m2_fcs_file, '\n')
  
  m2_fcs_list[[i]] <- 
    flowCore::read.FCS(m2_fcs_file, 
                       transformation = FALSE)
  
  m1_fcs_input_files[i] <- m1_fcs_file
  m2_fcs_input_files[i] <- m2_fcs_file
}


## Coalesce expression marker names in order to make flowset
#===============================================================================
coalesce_fcs_list_to_flowset <- function(fcs_list, cols_to_keep = NULL) {
  for(i in seq_along(fcs_list)) {
    fcs_data <- fcs_list[[i]]
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
    
    # To ensure same order of column names. Needed
    # to convert to flowset
    if (i > 1) {
      names <- colnames(exprs(fcs_list[[1]]))
      new_expr <- new_expr[, names] 
    }
    
    if (!is.null(cols_to_keep)) {
      new_expr <- new_expr[, cols_to_keep]
    }
    
    
    fcs_list[[i]] <- flowFrame(new_expr)
  }
  
  fcs_flowset <- as(fcs_list, 'flowSet')
}


m1_fcs_flowset <- 
  coalesce_fcs_list_to_flowset(m1_fcs_list, 
                               cols_to_keep = c('FSC-H', 'FSC-A', 'SSC-H', 
                                                'SSC-A', 'CD34', 'CD45', 'CD117'))

m2_fcs_flowset <- 
  coalesce_fcs_list_to_flowset(m2_fcs_list, 
                               cols_to_keep = c('FSC-H', 'FSC-A', 'SSC-H', 
                                                'SSC-A', 'CD38', 'CD45', 'CD64'))



#===============================================================================


# Normalization
#===============================================================================
## M1 

cat('Normalizing all M1 files ... \n')

channel_names <- colnames(exprs(m1_fcs_flowset[[1]]))

m1_normalized_fcs_flowset <- 
  flowStats::gaussNorm(m1_fcs_flowset, 
                       channel.names = channel_names, 
                       max.lms = c(3, 2, 3, 2, 2, 3, 2))

m1_normalization_confidence <- 
  m1_normalized_fcs_flowset$confidence

cat('M1 confidence: ', m1_normalization_confidence, '\n\n')

m1_normalized_fcs_flowset <- 
  m1_normalized_fcs_flowset$flowset



## M2

cat('Normalizing all M2 files ... \n')

channel_names <- colnames(exprs(m2_fcs_flowset[[1]]))

m2_normalized_fcs_flowset <- 
  flowStats::gaussNorm(m2_fcs_flowset, 
                       channel.names = channel_names, 
                       max.lms = c(3, 2, 3, 3, 3, 3, 3))

m2_normalization_confidence <- 
  m2_normalized_fcs_flowset$confidence

cat('M2 confidence: ', m2_normalization_confidence, '\n\n')


m2_normalized_fcs_flowset <- 
  m2_normalized_fcs_flowset$flowset
#===============================================================================



# Save
#===============================================================================
for (i in seq_along(m1_fcs_list)) {
  m1_transformed_fcs_file <- 
    m1_fcs_input_files[i] %>% 
    str_split('.fcs') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_gaussNorm_normalized.fcs')  
  
  m2_transformed_fcs_file <- 
    m2_fcs_input_files[i] %>% 
    str_split('.fcs') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_gaussNorm_normalized.fcs')  
  
  cat('\nSaving: ', m1_transformed_fcs_file, '\n')
  
  flowCore::write.FCS(m1_normalized_fcs_flowset[[i]], 
                      m1_transformed_fcs_file)
  
  cat('\nSaving: ', m2_transformed_fcs_file, '\n')
  
  flowCore::write.FCS(m2_normalized_fcs_flowset[[i]], 
                      m2_transformed_fcs_file)
  
}
#===============================================================================
