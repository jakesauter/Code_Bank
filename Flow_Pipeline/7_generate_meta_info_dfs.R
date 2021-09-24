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


choose_file_to_use <- function(filenames) {
  
  if (length(filenames) == 1)
    return(filenames)
  
  # Remove any filenames that have "DO NOT USE"
  contains_do_not_use <-
    str_detect(basename(filenames),
               'DO NOT USE')
  
  filenames <-
    filenames[!contains_do_not_use]
  
  if (length(filenames) == 1)
    return(filenames[[1]])
  
  # # Use repeat
  contains_repeat <-
    str_detect(basename(filenames), 'RPT') |
    str_detect(basename(filenames), 'REPEAT')
  
  if (length(which((contains_repeat))) == 1)
    return(filenames[contains_repeat])
  # 
  # # Use wash 
  contains_wash <-
    str_detect(basename(filenames), 'WASH') |
    str_detect(basename(filenames), 'XW')
  
  if(length(which(contains_wash)) == 1)
    return(filenames[contains_wash])
  
  
  # If not choose the latest file
  file <-
    filenames %>% 
    file.mtime() %>% 
    lubridate::as_datetime() %>% 
    which.max() %>% 
    filenames[.]
  
  return(file)
}


### Input data
patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE) 

patient_df <- 
  patient_df %>% 
  filter(Exclusion.Reason == "Not Excluded")




m1_meta_info_list <- 
  vector('list', 
         nrow(patient_df))

m2_meta_info_list <- 
  vector('list', 
         nrow(patient_df))

parse_meta_info <- function(parsed_meta) {
  baseline_date      <- parsed_meta$`CST BASELINE DATE` %>% 
    str_split('T') %>% .[[1]] %>% .[1]
  
  setup_date         <- parsed_meta$`CST SETUP DATE` %>% 
    str_split('T') %>% .[[1]] %>% .[1]
  
  config_create_date <- parsed_meta$`CYTOMETER CONFIG CREATE DATE` %>% 
    str_split('T') %>% .[[1]] %>% .[1]
  
  beads_lot_id       <- parsed_meta$`CST BEADS LOT ID`
  config_name        <- parsed_meta$`CYTOMETER CONFIG NAME`
  
  cytometer <- parsed_meta$`$CYT`
  cyt_num   <- parsed_meta$CYTNUM
  date      <- parsed_meta$`$DATE`
  
  param_names <- 
    colnames(parsed_meta$SPILL)
  
  # marker names
  
  marker_names <- 
    names(parsed_meta) %>% 
    str_detect('P[0-9]+S') %>% 
    parsed_meta[.] %>% 
    unname() %>% unlist() %>% 
    sort()
  
  meta_info_list <- 
    list(
      cytometer = cytometer, 
      cyt_num = cyt_num, 
      date = date, 
      baseline_date = baseline_date, 
      setup_date = setup_date, 
      config_create_date = config_create_date, 
      beads_lot_id = beads_lot_id, 
      config_name = config_name, 
      param_names = param_names, 
      marker_names = marker_names
    ) 
  
  meta_info_list
}

for (i in seq_len(nrow(patient_df))) {
  
  flow_directory <- patient_df[i, "Flow.Data.Folder"]
  
  cat('Processing flow directory:', i, 'out of:', nrow(patient_df))
  cat('\nFlow directory: ', flow_directory, '\n')
  
  # Need to resolve which M1/M2 here? 
  
  m1_fcs_file <-
    flow_directory %>% 
    list.files(full.names = TRUE, 
               pattern = '.*\\.fcs') %>% 
    .[str_detect(., 'm1|M1')] %>% 
    choose_file_to_use()
  
  m2_fcs_file <-
    flow_directory %>%     
    list.files(full.names = TRUE, 
               pattern = '.*\\.fcs') %>% 
    .[str_detect(., 'm2|M2')] %>% 
    choose_file_to_use()

  
  
  cat('\n\nReading: ', basename(m1_fcs_file), ' \n')
  parsed_meta <- description(read.FCS(m1_fcs_file,
                                      transformation = FALSE,
                                      which.lines = 1))


  m1_meta_info_list[[i]] <- parse_meta_info(parsed_meta)

  cat('Reading: ', basename(m2_fcs_file), ' \n\n')
  parsed_meta <- description(read.FCS(m2_fcs_file,
                                      transformation = FALSE,
                                      which.lines = 1))

  m2_meta_info_list[[i]] <- parse_meta_info(parsed_meta)
  
}

m1_meta_info_df <- 
  do.call(rbind, m1_meta_info_list) %>% 
  as.data.frame() %>% 
  mutate(
    date = lubridate::dmy(date)) %>% 
  mutate(baseline_date = lubridate::ymd(baseline_date), 
         setup_date = lubridate::ymd(setup_date), 
         flow_id = patient_df$Accession.Number) %>% 
  arrange(cyt_num, date, baseline_date, setup_date)

m2_meta_info_df <- 
  do.call(rbind, m2_meta_info_list) %>% 
  as.data.frame() %>% 
  mutate(
    date = lubridate::dmy(date)) %>% 
  mutate(baseline_date = lubridate::ymd(baseline_date), 
         setup_date = lubridate::ymd(setup_date), 
         flow_id = patient_df$Accession.Number) %>% 
  arrange(cyt_num, date, baseline_date, setup_date)


save_loc <- "/Users/sauterj1/Documents/Patient_Folder_Analysis/data"

saveRDS(m1_meta_info_df, 
        file.path(save_loc, 'm1_meta_info_df.Rds'))

saveRDS(m2_meta_info_df, 
        file.path(save_loc, 'm2_meta_info_df.Rds'))
