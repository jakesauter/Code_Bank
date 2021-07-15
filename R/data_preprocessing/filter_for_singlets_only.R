#!/bin/Rscript


### Input data
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Correct_QC_Flow_Folders/"


library(dplyr)
library(ggcyto)
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

is_single <- function(ff, plot = FALSE, ...) {
  fsc_a <- flowCore::exprs(ff)[,"FSC-A"]
  fsc_h <- flowCore::exprs(ff)[,"FSC-H"]
  
  bins <- cut(fsc_a, 10)
  
  ratios <- fsc_h / fsc_a
  slope_per_bin <- tapply(ratios, bins, mean)
  expected_values <- fsc_a * slope_per_bin[bins]
  deviations <- abs(fsc_h - expected_values)
  
  x <- tapply(fsc_a, bins, mean)
  e <- tapply(expected_values, bins, mean)
  d <- tapply(deviations, bins, function(x){mean(x) + 2*sd(x)})
  y <- e - d
  
  spl <- splinefun(x, y)
  
  if (plot) {
    flowDensity::plotDens(ff, c("FSC-A", "FSC-H"), ...)
    points(x, e, col = "red", pch = 19)
    points(x, y, col = "red", pch = 19)
    lines(seq(1, 300000, by = 1000), 
          spl(seq(1, 300000, by = 1000)),
          col = "red",
          lwd = 2)
  }
  
  selection <- fsc_h > spl(fsc_a)
  return(selection)
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
    .[!str_detect(., 'singlets')]
  
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
  
  
  m1_singlet_events <- is_single(m1_fcs_dat)
  m2_singlet_events <- is_single(m2_fcs_dat)
  
  exprs(m1_fcs_dat) <- exprs(m1_fcs_dat)[m1_singlet_events, ]
  exprs(m2_fcs_dat) <- exprs(m2_fcs_dat)[m2_singlet_events, ]
  
  cat(length(which(!m1_singlet_events))/length(m1_singlet_events), '% doublet events detected in M1\n')
  cat(length(which(!m2_singlet_events))/length(m2_singlet_events), '% doublet events detected in M2\n')
  
  
  ## Saving
  m1_transformed_fcs_file <- 
    m1_fcs_file %>% 
    str_split('.fcs') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_singlets_only.fcs')  
  
  m2_transformed_fcs_file <- 
    m2_fcs_file %>% 
    str_split('.fcs') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_singlets_only.fcs')  
  
  flowCore::write.FCS(m1_fcs_dat, 
                      m1_transformed_fcs_file)
  
  flowCore::write.FCS(m2_fcs_dat, 
                      m2_transformed_fcs_file)  
  
}


