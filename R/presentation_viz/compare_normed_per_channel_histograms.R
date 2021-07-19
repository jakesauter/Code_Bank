#!/bin/Rscript

### Output Location 
pdf_dir <- '~/Documents/Patient_Folder_Analysis/normalized_KDE_smaller_kernel_pdfs/'

### Input data
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE)

library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(flowCore)
library(flowStats)
library(patchwork)

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


m1_normed_fcs_list <- vector('list', length(flow_directories))
m1_non_normed_fcs_list <- vector('list', length(flow_directories))

m2_normed_fcs_list <- vector('list', length(flow_directories))
m2_non_normed_fcs_list <- vector('list', length(flow_directories))

names(m1_normed_fcs_list) <- flow_directories
names(m2_normed_fcs_list) <- flow_directories
names(m1_non_normed_fcs_list) <- flow_directories
names(m2_non_normed_fcs_list) <- flow_directories

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')] %>% 
    .[!str_detect(., 'transformed')] %>% 
    .[!str_detect(., 'singlets')] 
  
  normed_fcs_files <- 
    fcs_files %>% 
    .[str_detect(., 'normalized')]
  
  non_normed_fcs_files <- 
    fcs_files %>% 
    .[!str_detect(., 'normalized')]
  
  m1_normed_fcs_file <- 
    str_detect(normed_fcs_files, 
               'M1|m1') %>% 
    normed_fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  m1_non_normed_fcs_file <- 
    str_detect(non_normed_fcs_files, 
               'M1|m1') %>% 
    non_normed_fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  
  m2_normed_fcs_file <- 
    str_detect(normed_fcs_files, 
               'M2|m2') %>% 
    normed_fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  m2_non_normed_fcs_file <- 
    str_detect(non_normed_fcs_files, 
               'M2|m2') %>% 
    non_normed_fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  cat('Reading FCS file: ', m1_normed_fcs_file, '\n')
  
  m1_normed_fcs_list[[i]] <- 
    flowCore::read.FCS(m1_normed_fcs_file, 
                       transformation = FALSE)
  
  
  cat('Reading FCS file: ', m1_non_normed_fcs_file, '\n')
  
  m1_non_normed_fcs_list[[i]] <- 
    flowCore::read.FCS(m1_non_normed_fcs_file, 
                       transformation = FALSE)
  
  
  cat('Reading FCS file: ', m2_normed_fcs_file, '\n')
  
  m2_normed_fcs_list[[i]] <- 
    flowCore::read.FCS(m2_normed_fcs_file, 
                       transformation = FALSE)
  
  cat('Reading FCS file: ', m2_non_normed_fcs_file, '\n')
  
  m2_non_normed_fcs_list[[i]] <- 
    flowCore::read.FCS(m2_non_normed_fcs_file, 
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
m1_normed_fcs_list <- coalesce_fcs_list(m1_normed_fcs_list)
m1_non_normed_fcs_list <- coalesce_fcs_list(m1_non_normed_fcs_list)

cat('\nExtracting M2 marker columns ... \n\n')
m2_normed_fcs_list <- coalesce_fcs_list(m2_normed_fcs_list)
m2_non_normed_fcs_list <- coalesce_fcs_list(m2_non_normed_fcs_list)

### Plotting Per-Channel KDEs
plot_kde_for_fcs_list <- function(normed_fcs_list, 
                                  non_normed_fcs_list, 
                                  pdf_dir) {
  
  if (!dir.exists(pdf_dir)) {
    dir.create(pdf_dir, recursive = TRUE) 
  }
  
  names(normed_fcs_list) <- basename(names(normed_fcs_list))
  names(non_normed_fcs_list) <- basename(names(non_normed_fcs_list))
  
  # arrange by cancer type to make insights easier 
  patient_df <- 
    patient_df %>% 
    arrange(CbioPortal.Cancer.Type, CbioPortal.Cancer.Type.Detailed) %>% 
    filter(basename(Flow.Data.Folder) %in% names(normed_fcs_list))
  
  normed_fcs_list <- 
    normed_fcs_list[basename(patient_df$Flow.Data.Folder)]
  
  non_normed_fcs_list <- 
    non_normed_fcs_list[basename(patient_df$Flow.Data.Folder)]
  
  flow_ids <- patient_df$Accession.Number
  
  for (i in seq_len(ncol(exprs(normed_fcs_list[[1]])))) {
    marker <- colnames(exprs(normed_fcs_list[[1]]))[i]
    cat('generating plots for marker: ', marker, '\n')
    pdf(file.path(pdf_dir, paste0(marker, '.pdf')), 
        width = 7, 
        height = 6)
    
    max_marker_value <- max(sapply(normed_fcs_list, function(x) max(exprs(x)[, marker])))
    max_marker_value <- max(c(max_marker_value, 
                              sapply(non_normed_fcs_list, function(x) max(exprs(x)[, marker]))))
    
    
    min_marker_value <- min(sapply(normed_fcs_list, function(x) min(exprs(x)[, marker])))
    min_marker_value <- min(c(min_marker_value, 
                              sapply(normed_fcs_list, function(x) min(exprs(x)[, marker]))))
    
    for (i in seq_along(normed_fcs_list)) {
      patient_directory <- names(normed_fcs_list)[i]
      patient_info <- 
        patient_df %>% 
        filter(patient_directory == basename(Flow.Data.Folder))
      
      p1 <- ggplot(data.frame(x = exprs(non_normed_fcs_list[[i]])[, marker])) + 
        geom_density(aes(x=x), n=200) + 
        ggtitle(paste0(flow_ids[[i]], " Not Normalized")) + 
        xlim(c(min_marker_value, max_marker_value)) + 
        xlab(paste0(marker, ' (Arcsinh Transformed)')) + 
        ylab('Density') 
      
      p2 <- ggplot(data.frame(x = exprs(normed_fcs_list[[i]])[, marker])) + 
        geom_density(aes(x=x), n=200) + 
        ggtitle(paste0(flow_ids[[i]], ' Normalized')) +
        xlim(c(min_marker_value, max_marker_value)) + 
        xlab(paste0(marker, ' Normalized (Arcsinh Transformed)')) + 
        ylab('Density') 
      
      p3 <- ggplot(data.frame(x1 = exprs(non_normed_fcs_list[[i]])[, marker], 
                        x2 = exprs(normed_fcs_list[[i]])[, marker])) + 
        geom_density(aes(x=x1, 
                         fill = 'non-normed'), 
                     alpha = 0.5, 
                     lwd = .7, 
                     n = 200) + 
        geom_density(aes(x=x2, 
                         fill = 'normed'), 
                     alpha = 0.5, 
                     lwd = .7) + 
        xlim(c(min_marker_value, max_marker_value)) + 
        xlab(paste0(marker, ' (Arcsinh Transformed)')) + 
        ylab('Density') + 
        scale_fill_manual(values = c('non-normed' = 'blue', 
                                      'normed' = 'red'), 
                          labels = c('non-normed' = 'Not Normalized', 
                                     'normed' = 'Normalized')) + 
        theme(legend.title = element_blank())
      
      p <- (p1 + p2) / p3
      p <- p + plot_annotation(title = paste0('Cancer_Type= ', patient_info$CbioPortal.Cancer.Type, 
                                              '\nDetailed_Cancer_Type= ', patient_info$CbioPortal.Cancer.Type.Detailed))
      
      print(p)
    }
    dev.off()
  }
  
  invisible('Completed Sucessfully')
}

cat('Plotting M1 KDEs ... \n')
plot_kde_for_fcs_list(m1_normed_fcs_list,
                      m1_non_normed_fcs_list, 
                      file.path(pdf_dir, 'M1'))

cat('Plotting M2 KDEs ... \n')
plot_kde_for_fcs_list(m2_normed_fcs_list, 
                      m2_non_normed_fcs_list, 
                      file.path(pdf_dir, 'M2'))

