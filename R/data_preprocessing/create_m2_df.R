
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

patient_df <- read.csv('~/Documents/Woodlist/jake_processed_flow_df_ids_051321.csv', 
                       stringsAsFactors = FALSE)

m2_df <- 
  patient_df %>% 
  dplyr::filter(str_starts(M2_path, '/')) %>% 
  mutate(m2_markers = 
           mclapply(M2_path, function(m2_path) {
             
             if (!file.exists(m2_path)) {
               return('Invalid Filepath')
             }
             
             print(m2_path)
             
             # currently there is the possiblity that multiple paths
             # are in M1 path or M2 path columns
             m2_path <- 
               m2_path %>% 
               str_split(',') %>% 
               .[[1]] %>% .[[1]] %>% 
               str_trim()
             
             get_fcs_markers(m2_path)
           }, mc.cores = 10))

saveRDS(m2_df, '~/Documents/Patient_Folder_Analysis/data/m2_df.rds')