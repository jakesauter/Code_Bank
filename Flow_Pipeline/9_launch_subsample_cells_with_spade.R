#!/bin/Rscript

set.seed(497)

options(warn = -1)

flow_directories <-
  list.files("/gpfs/mskmind_ess/sauterj1/Full_Cohort_Anon_Flow_Folders",
             full.names = TRUE)

# flow_directories <- 
#   list.files("~/Documents/Woodlist/Full_Cohort_Anon_Flow_Folders/", 
#              full.names = TRUE)

for (i in seq_along(flow_directories)) {
  
  cat('Processing directory:', i, 'of:', length(flow_directories), '\n')
  
  flow_directory <- flow_directories[i]
  
  cat('Flow Directory: ', flow_directory, '\n\n')
  
  tic <- Sys.time()
  
  system(paste0(
    'Rscript --no-save --no-restore launch_spade_for_dir.R ', 
    flow_directories[[i]]
  ), intern = TRUE)
  
  toc <- Sys.time()
  
  cat('\n')
  print(toc-tic)
  cat('\n')
  
  if (file.exists(file.path(flow_directory, 'error_file.txt'))) {
    cat('\n\n ERROR OCCURRED WHEN SAMPLING DIRECTORY: ', flow_directory, '\n\n') 
  }
  
}



if (FALSE) {
  
  
  library(flowCore)
  
  flow_directories <-
    list.files("/gpfs/mskmind_ess/sauterj1/Full_Cohort_Anon_Flow_Folders",
               full.names = TRUE)
  
  for (i in seq_along(flow_directories)) {
    
    cat('Processing directory:', i, 'of:', length(flow_directories), '\n')
    
    flow_directory <- flow_directories[i]
    
    cat('Flow Directory: ', flow_directory, '\n\n')

    
    m1_rand_fcs   <- read.FCS(file.path(flow_directory, "processed_M1_randomly_subsampled_seed_497.fcs"))
    m2_rand_fcs   <- read.FCS(file.path(flow_directory, "processed_M2_randomly_subsampled_seed_497.fcs"))
    m1_spade_fcs  <- read.FCS(file.path(flow_directory, "processed_M1_subsampled_100k.fcs"))
    m2_spade_fcs  <- read.FCS(file.path(flow_directory, "processed_M2_subsampled_100k.fcs"))
    
    
    cat('m1_spade_fcs: ', nrow(exprs(m1_spade_fcs)), '\n')
    cat('m2_spade_fcs: ', nrow(exprs(m2_spade_fcs)), '\n')
    cat('m1_rand_fcs: ', nrow(exprs(m1_rand_fcs)), '\n')
    cat('m2_rand_fcs: ', nrow(exprs(m2_rand_fcs)), '\n')
  
  }
  
  
  
  
}






