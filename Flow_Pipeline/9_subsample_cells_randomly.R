  
library(flowCore)
library(magrittr)
library(stringr)

set.seed(497)

flow_directories <-
  list.files("/gpfs/mskmind_ess/sauterj1/Full_Cohort_Anon_Flow_Folders",
             full.names = TRUE)

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
    str_c('_randomly_subsampled_seed_497.fcs')
  
  m2_subsampled_fcs_file <-
    m2_fcs_file %>% 
    strsplit('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_randomly_subsampled_seed_497.fcs')
  
  
  m1_fcs_dat <-
    read.FCS(m1_fcs_file,
             transformation = FALSE)


  exprs(m1_fcs_dat) <-
    exprs(m1_fcs_dat)[sample(1:nrow(exprs(m1_fcs_dat)),
                             100e3), ]

  write.FCS(m1_fcs_dat,
            m1_subsampled_fcs_file)


  m2_fcs_dat <-
    read.FCS(m2_fcs_file,
             transformation = FALSE)

  exprs(m2_fcs_dat) <-
    exprs(m2_fcs_dat)[sample(1:nrow(exprs(m2_fcs_dat)),
                             100e3), ]

  write.FCS(m2_fcs_dat,
            m2_subsampled_fcs_file)
  
  spade::SPADE.removeExistingDensityAndClusterColumns(m1_subsampled_fcs_file)
  spade::SPADE.removeExistingDensityAndClusterColumns(m2_subsampled_fcs_file)
  
  system(paste0('rm ', file.path(flow_directory, 'processed_M1_randomly_subsampled_seed_497.fcs.orig1')))
  system(paste0('rm ', file.path(flow_directory, 'processed_M2_randomly_subsampled_seed_497.fcs.orig1')))
  
  
}



if (FALSE) {
  
   
  library(flowCore)
  
  flow_directories <-
    list.files("/gpfs/mskmind_ess/sauterj1/Full_Cohort_Anon_Flow_Folders",
               full.names = TRUE)
  
  
  for (i in seq_along(flow_directories)) {
    
    flow_directory <- flow_directories[i]
    
    cat('Flow Directory: ', i, ' of: ', length(flow_directories),  '\n')
    cat('Flow Directory: ', flow_directory, '\n\n')
    
    m1_file <- 
      file.path(flow_directory, 
                'processed_M1_randomly_subsampled_seed_497.fcs')
    
    m2_file <- 
      file.path(flow_directory, 
                'processed_M2_randomly_subsampled_seed_497.fcs')
    
    spade::SPADE.removeExistingDensityAndClusterColumns(m1_file)
    spade::SPADE.removeExistingDensityAndClusterColumns(m2_file)

    system(paste0('rm ', file.path(flow_directory, 'processed_M1_randomly_subsampled_seed_497.fcs.orig1')))
    system(paste0('rm ', file.path(flow_directory, 'processed_M2_randomly_subsampled_seed_497.fcs.orig1')))
  }
  
  
  
  
}









