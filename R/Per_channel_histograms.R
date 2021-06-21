#!/bin/Rscript

### Input data

flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(flowCore)
filter <- dplyr::filter

# Gather all QC Statistics
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

flow_dirs_passed_qc <- 
  full_qc_results %>% 
  filter(total_filtered_perc < 80) %>% 
  mutate(dirname = dirname(filename)) %>% 
  count(dirname) %>% 
  filter(n > 1) %>% 
  .$dirname

flow_directories <- 
  flow_dirs_passed_qc

# flow_directories <- 
#   sample(flow_directories, 10)

fcs_list <- vector('list', length(flow_directories))
names(fcs_list) <- flow_directories

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')]
  
  fcs_file <- 
    str_detect(fcs_files, 
               'M1|m1') %>% 
    fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  print(fcs_file)
  
  fcs_list[[i]] <- 
    flowCore::read.FCS(fcs_file, 
                       transformation = FALSE)
  
}



### Merging all data into flowFrame object

# First big problem, need to modify the metadata so that we can merge together into a flowFrame


library(flowStats)
library(flowCore)
library(stringr)

# Need to compose each flowFrame to have the
# same parameter names
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
  
  fcs_list[[i]] <- flowFrame(new_expr)
}

fcs_flowset <- as(fcs_list, 'flowSet')


### Plotting Per-Channel KDEs


## First, gather metadata on patients
  
  
patient_df <- read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
                       stringsAsFactors = FALSE)




new_fcs_list <- fcs_list
names(new_fcs_list) <- basename(names(fcs_list))

# arrange by cancer type to make insights easier 
patient_df <- 
  patient_df %>% 
  arrange(CbioPortal.Cancer.Type, CbioPortal.Cancer.Type.Detailed) %>% 
  filter(basename(Flow.Data.Folder) %in% names(new_fcs_list))

new_fcs_list <- 
  new_fcs_list[basename(patient_df$Flow.Data.Folder)]


for (marker in colnames(exprs(new_fcs_list[[1]]))) {
  cat('generating plots for marker: ', marker, '\n')
  pdf(paste0('~/Documents/Patient_Folder_Analysis/KDE_pdfs/', marker, '.pdf'))
  
  for (i in seq_along(new_fcs_list)) {
    patient_directory <- names(new_fcs_list)[i]
    patient_info <- 
      patient_df %>% 
      filter(patient_directory == basename(Flow.Data.Folder))
    
    p <- ggplot(data.frame(x = exprs(new_fcs_list[[i]])[, marker])) + 
      geom_density(aes(x=x)) + 
      ggtitle(paste0('Cancer_Type= ', patient_info$CbioPortal.Cancer.Type, 
                     '\nDetailed_Cancer_Type= ', patient_info$CbioPortal.Cancer.Type.Detailed)) 
    
    print(p)
  }
  dev.off()
}


