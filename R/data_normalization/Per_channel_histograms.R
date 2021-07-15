#!/bin/Rscript

### Output Location 
pdf_dir <- '~/Documents/Patient_Folder_Analysis/optim_transformed_KDE_pdfs/'

### Input data
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Correct_QC_Flow_Folders/"

patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
            stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(flowCore)
library(flowStats)

filter <- dplyr::filter

# Gather all QC Statistics
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


m1_fcs_list <- vector('list', length(flow_directories))
m2_fcs_list <- vector('list', length(flow_directories))

names(m1_fcs_list) <- flow_directories
names(m2_fcs_list) <- flow_directories

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')] %>% 
    .[str_detect(., 'transformed')] %>% 
    .[!str_detect(., 'asinh')] %>% 
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
  
  m1_fcs_list[[i]] <- 
    flowCore::read.FCS(m1_fcs_file, 
                       transformation = FALSE)
  
  cat('Reading FCS file: ', m2_fcs_file, '\n')
  
  m2_fcs_list[[i]] <- 
    flowCore::read.FCS(m2_fcs_file, 
                       transformation = FALSE)
  
}



### Merging all data into flowFrame object
#
# Need to modify the metadata so that we can 
# merge together into a flowFrame
coalesce_fcs_list <- function(fcs_list) {
  for(i in seq_along(fcs_list)) {
    fcs_expr <- exprs(fcs_list[[i]])
    metadata <- pData(parameters(fcs_list[[i]]))
    
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
  return(fcs_list)
}

cat('\n\nExtracting M1 marker columns ... \n')
m1_fcs_list <- coalesce_fcs_list(m1_fcs_list)

cat('\nExtracting M2 marker columns ... \n\n')
m2_fcs_list <- coalesce_fcs_list(m2_fcs_list)

### Plotting Per-Channel KDEs
plot_kde_for_fcs_list <- function(fcs_list, pdf_dir) {
  
  if (!dir.exists(pdf_dir)) {
    dir.create(pdf_dir, recursive = TRUE) 
  }
  
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
    pdf(file.path(pdf_dir, paste0(marker, '.pdf')), 
        width = 7, 
        height = 4)
    
    max_marker_value <- max(sapply(fcs_list, function(x) max(exprs(x)[, marker])))
    min_marker_value <- min(sapply(fcs_list, function(x) min(exprs(x)[, marker])))
    
    for (i in seq_along(new_fcs_list)) {
      patient_directory <- names(new_fcs_list)[i]
      patient_info <- 
        patient_df %>% 
        filter(patient_directory == basename(Flow.Data.Folder))
      
      p <- ggplot(data.frame(x = exprs(new_fcs_list[[i]])[, marker])) + 
        geom_density(aes(x=x, y = ..scaled..)) + 
        ggtitle(paste0('Cancer_Type= ', patient_info$CbioPortal.Cancer.Type, 
                       '\nDetailed_Cancer_Type= ', patient_info$CbioPortal.Cancer.Type.Detailed)) + 
        xlim(c(min_marker_value, max_marker_value)) + 
        xlab(paste0(marker, ' (Arcsinh Transformed)')) + 
        ylab('Density')
      
      print(p)
    }
    dev.off()
  }
  
  invisible('Completed Sucessfully')
}

cat('Plotting M1 KDEs ... \n')
plot_kde_for_fcs_list(m1_fcs_list, file.path(pdf_dir, 'M1'))

cat('Plotting M2 KDEs ... \n')
plot_kde_for_fcs_list(m2_fcs_list, file.path(pdf_dir, 'M2'))

