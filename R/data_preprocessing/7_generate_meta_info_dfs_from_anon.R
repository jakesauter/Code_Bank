library(dplyr)
library(ggplot2)
library(stringr)
library(magrittr)
library(flowCore)
library(flowStats)
library(patchwork)
library(lubridate)
library(stringr)
library(knitr)

filter <- dplyr::filter


ANON_COHORT_PATH <- 
  "/Users/sauterj1/Documents/Woodlist/Full_Cohort_Anon_Flow_Folders/"

### Input data
patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE) 

patient_df <- 
  patient_df %>% 
  filter(Exclusion.Reason == "Not Excluded")




m1_meta_info_list <- 
  vector('list', 
         length(flow_directories))

m2_meta_info_list <- 
  vector('list', 
         length(flow_directories))

parse_meta_info <- function(parsed_meta) {
  
  cyt_num   <- parsed_meta$CYTNUM
  date      <- parsed_meta$`$DATE`
  
  param_names <- 
    colnames(parsed_meta$SPILL)

  meta_info_list <- 
    list(
      cyt_num = cyt_num, 
      date = date
    ) 
  
  meta_info_list
}


dirs <- list.dirs(ANON_COHORT_PATH)
counter <- 0

for (dir in dirs) {
  
  flow_directory <- basename(dir)
  counter <- counter + 1
  cat('Processing flow directory:', counter, 'out of:', length(dirs))
  cat('\nFlow directory: ', flow_directory, '\n\n')
  
  # Need to resolve which M1/M2 here? 
  
  m1_fcs_file <- 
    file.path('~/Documents/Woodlist/Full_Cohort_Anon_Flow_Folders/', 
              patient_df$Accession.Number[[i]], 
              'processed_M1.fcs')
  
  m2_fcs_file <- 
    file.path('~/Documents/Woodlist/Full_Cohort_Anon_Flow_Folders/', 
              patient_df$Accession.Number[[i]], 
              'processed_M2.fcs')
  
  parsed_meta <- description(read.FCS(m1_fcs_file, 
                                      transformation = FALSE, 
                                      which.lines = 1))
  
  
  m1_meta_info_list[[counter]] <- parse_meta_info(parsed_meta)
  
  parsed_meta <- description(read.FCS(m2_fcs_file, 
                                      transformation = FALSE, 
                                      which.lines = 1))
  
  # m2_meta_info_list[[counter]] <- 
    parsed_meta %>% 
    data.frame(cyt_num = .$CYTNUM, 
               date = .$`$DATE`)
    
  
}

m1_meta_info_df <- 
  do.call(rbind, m1_meta_info_list) %>% 
  as.data.frame() %>% 
  mutate(
    date = lubridate::dmy(date)) %>% 
  mutate(flow_id = patient_df$Accession.Number) %>% 
  arrange(cyt_num, date)

m2_meta_info_df <- 
  do.call(rbind, m2_meta_info_list) %>% 
  as.data.frame() %>% 
  mutate(
    date = lubridate::dmy(date)) %>% 
  mutate(flow_id = patient_df$Accession.Number) %>% 
  arrange(cyt_num, date)


saveRDS(m1_meta_info_df, 'data/m1_meta_info_df.Rds')
saveRDS(m2_meta_info_df, 'data/m2_meta_info_df.Rds')
