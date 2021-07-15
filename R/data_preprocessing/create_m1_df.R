
library(dplyr)
library(stringr)
library(ggplot2)
library(parallel)
library(flowCore)
library(lubridate)

get_fcs_markers <- function(path) {
  sample <- flowCore::read.FCS(path, transformation = FALSE)
  flowCore::markernames(sample)
}

patient_df <- 
  read.csv('~/Documents/Woodlist/jake_processed_flow_df_ids_051321.csv', 
                       stringsAsFactors = FALSE)

m1_df <- 
  patient_df %>% 
  dplyr::filter(str_starts(M1_path, '/')) %>% 
  mutate(m1_markers = 
           mclapply(M1_path, function(m1_path) {
             
             if (!file.exists(m1_path)) {
               return('Invalid Filepath')
             }
             
             print(m1_path)
             
             # currently there is the possiblity that multiple paths
             # are in M1 path or M1 path columns
             m1_path <- 
               m1_path %>% 
               str_split(',') %>% 
               .[[1]] %>% .[[1]] %>% 
               str_trim()
             
             get_fcs_markers(m1_path)
           }, mc.cores = 10))

saveRDS(m1_df, '~/Documents/Patient_Folder_Analysis/data/m1_df.rds')





