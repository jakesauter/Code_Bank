#!/bin/Rscript


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
library(lubridate)
library(knitr)

filter <- dplyr::filter

flow_directories <- 
  patient_df %>% 
  filter(Exclusion.Reason == "Not Excluded", 
         !str_detect(Flow.Data.Folder, ','), 
         !str_detect(M1_path, ','), 
         !str_detect(M2_path, ',')) %>% 
  .$Flow.Data.Folder
  


m1_meta_info_list <- vector('list', length(flow_directories))
m2_meta_info_list <- vector('list', length(flow_directories))

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

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[i]
  
  cat('Processing flow directory:', i, 'out of:', length(flow_directories))
  cat('\nFlow directory: ', flow_directory, '\n\n')
  
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[!str_detect(., 'compensated')] %>% 
    .[!str_detect(., 'high_quality_events')] %>% 
    .[!str_detect(., 'transformed')] %>% 
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
  
  
  parsed_meta <- description(read.FCS(m1_fcs_file, 
                                      transformation = FALSE, 
                                      which.lines = 1))
  
  m1_meta_info_list[[i]] <- parse_meta_info(parsed_meta)
  
  parsed_meta <- description(read.FCS(m2_fcs_file, 
                                      transformation = FALSE, 
                                      which.lines = 1))
  m2_meta_info_list[[i]] <- parse_meta_info(parsed_meta)
 
}

m1_meta_info_df <- 
  do.call(rbind, m1_meta_info_list) %>% 
  as.data.frame() %>% 
  mutate(
    date = date %>% 
      str_replace_all('JAN', '1') %>% 
      str_replace_all('FEB', '2') %>% 
      str_replace_all('MAR', '3') %>% 
      str_replace_all('APR', '4') %>% 
      str_replace_all('MAY', '5') %>% 
      str_replace_all('JUN', '6') %>% 
      str_replace_all('JUL', '7') %>% 
      str_replace_all('AUG', '8') %>% 
      str_replace_all('SEP', '9') %>% 
      str_replace_all('OCT', '10') %>% 
      str_replace_all('NOV', '11') %>% 
      str_replace_all('DEC', '12') %>% 
      lubridate::dmy()
  ) %>% 
  mutate(baseline_date = lubridate::ymd(baseline_date), 
         setup_date = lubridate::ymd(setup_date)) %>% 
  arrange(cyt_num, date, baseline_date, setup_date)

m2_meta_info_df <- 
  do.call(rbind, m2_meta_info_list) %>% 
  as.data.frame() %>% 
  mutate(
    date = date %>% 
      str_replace_all('JAN', '1') %>% 
      str_replace_all('FEB', '2') %>% 
      str_replace_all('MAR', '3') %>% 
      str_replace_all('APR', '4') %>% 
      str_replace_all('MAY', '5') %>% 
      str_replace_all('JUN', '6') %>% 
      str_replace_all('JUL', '7') %>% 
      str_replace_all('AUG', '8') %>% 
      str_replace_all('SEP', '9') %>% 
      str_replace_all('OCT', '10') %>% 
      str_replace_all('NOV', '11') %>% 
      str_replace_all('DEC', '12') %>% 
      lubridate::dmy()
  ) %>% 
  mutate(baseline_date = lubridate::ymd(baseline_date), 
         setup_date = lubridate::ymd(setup_date)) %>% 
  arrange(cyt_num, date, baseline_date, setup_date)



cytometer_range <- function(df, cytom_num) {
  df <- 
    df %>% 
    arrange(date) %>% 
    filter(cyt_num == cytom_num)
  
  first <- 
    df %>% 
    head(1) %>% 
    .$date
  
  last <- 
    df %>% 
    tail(1) %>% 
    .$date

  c(first, last)
}

cytometer_ids <- 
  m1_meta_info_df %>% 
  .$cyt_num %>% 
  unique() %>% 
  unlist()

cytometer_ids %>% 
  lapply(function(x) cytometer_range(m1_meta_info_df, x)) %>% 
  set_names(cytometer_ids)

# Looks like machines for M1 were changed on 2020-7-11 / 2020-10-13

# Also looks like multiple cytometers were used at the same time

m1_meta_info_df %>%
  
  mutate(month = month(date), 
         year = year(date)) %>% 
  group_by(year, month, cyt_num, cytometer) %>% 
  count() %>% 
  mutate(cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  arrange(year, month, cyt_num) %>% 
  print(n = Inf)

# Starts with V -- FACSCanto
# Starts with R -- LSRFortessa

unique_cytometers_used <- 
  sorted_meta_info_df$cyt_num %>% 
  unique()

# 8 different cytometers used for m1 over the years

sorted_meta_info_df %>% 
  mutate(year = year(date), 
         cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  count(cyt_num, year) %>% 
  arrange(year)

# same three used from 2016-2019 (21, 98, 99)
# new machines introduced in 2020 and 2021


# make plot -- time on x axis, machine on y axis, 
# heatmap value of how many uses

p1 <- 
  m1_meta_info_df %>% 
  mutate(date = floor_date(date, 'quarter'),  
         cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  count(cyt_num, date) %>% 
  arrange(date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y = cyt_num, fill = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('') + ylab('') + 
  ggtitle("M1 Cytometers Used Over Time") + 
  theme(legend.position = 'None')

p2 <- 
  m2_meta_info_df %>% 
  mutate(date = floor_date(date, 'quarter'), 
         cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  count(cyt_num, date) %>% 
  arrange(date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y = cyt_num, fill = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('') + ylab('') + 
  ggtitle("M2 Cytometers Used Over Time") + 
  theme(legend.position = 'None')


p1 / p2


# Same flow cytometers used for M1 / M2? 
m1_meta_info_df %>% 
  filter(cyt_num == "V657338000099") %>% 
  select(date, ) %>% 
  arrange(date) %>% 
  head()


m2_meta_info_df %>% 
  filter(cyt_num == "V657338000099") %>% 
  select(date) %>% 
  arrange(date) %>% 
  head()


# What about parameters? 
#
# For each machine determine if the paramters have changed 

param_names <- 
  m1_meta_info_df %>% 
  filter(cyt_num == "R66093700082") %>% 
  select(param_names) %>% 
  unlist()



p3 <- 
  m1_meta_info_df %>% 
  mutate(date = floor_date(date, 'quarter'), 
         cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  rowwise() %>% 
  mutate(param_names = paste0(param_names, collapse = ',')) %>%
  tidyr::separate_rows(param_names, sep = ',') %>% 
  count(param_names, date) %>% 
  arrange(date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y = param_names, fill = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('') + ylab('') + 
  ggtitle("Cytometer Channels Over Time")

p4 <- 
  m2_meta_info_df %>% 
  mutate(date = floor_date(date, 'quarter'), 
         cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  rowwise() %>% 
  mutate(param_names = paste0(param_names, collapse = ',')) %>%
  tidyr::separate_rows(param_names, sep = ',') %>% 
  count(param_names, date) %>% 
  arrange(date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y = param_names, fill = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('') + ylab('') + 
  ggtitle("Channels Over Time")


p4 / p1



# marker names over time 


p5 <- 
  m1_meta_info_df %>% 
  mutate(date = floor_date(date, unit = 'quarter')) %>% 
  rowwise() %>% 
  mutate(marker_names = paste0(marker_names, collapse = ',')) %>%
  tidyr::separate_rows(marker_names, sep = ',') %>% 
  count(marker_names, date) %>% 
  arrange(date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y = marker_names, fill = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  # guides(fill = guide_legend(title = '# Flow Ids')) +
  theme(legend.position = 'None') + 
  xlab('') + ylab('') + 
  ggtitle("M1 Markers Over Time")

p6 <- 
  m2_meta_info_df %>% 
  mutate(date = floor_date(date, 'quarter')) %>% 
  rowwise() %>% 
  mutate(marker_names = paste0(marker_names, collapse = ',')) %>%
  tidyr::separate_rows(marker_names, sep = ',') %>% 
  count(marker_names, date) %>% 
  arrange(date) %>% 
  ggplot() + 
  geom_tile(aes(x = date, y = marker_names, fill = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('') + ylab('') + 
  ggtitle("M2 Markers Over Time") + 
  theme(legend.position = 'None')


p1 / p5

p3

p2 / p6

###################
# Baseline dates
##################

m1_batches <- 
  m1_meta_info_df %>% 
  count(baseline_date, cyt_num) %>% 
  arrange(cyt_num, baseline_date) 

m1_batches$cyt_num <- unlist(m1_batches$cyt_num)

p1 + 
  geom_point(data = batches, 
             aes(x = baseline_date, 
                 y = cyt_num, 
                 shape = 'baseline date'), 
             size = 4) + 
  theme(legend.position = 'right') + 
  scale_shape_manual(values = c('baseline date' = 4)) + 
  theme(legend.title = element_blank()) + 
  ggtitle('Cytometers Used Over Time')


p3 <- m1_meta_info_df %>% 
  mutate(cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  count(cyt_num, baseline_date) %>% 
  arrange(baseline_date) %>%
  mutate(baseline_date = as.character(baseline_date)) %>% 
  ggplot() + 
  geom_tile(aes(x = baseline_date, y = cyt_num, fill = n)) + 
  geom_text(aes(x = baseline_date, y = cyt_num, label = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('Baseline Date') + ylab('') + 
  ggtitle("Cytometers and Baseline Dates") + 
  theme(legend.position = 'None')



p1 <- m1_meta_info_df %>% 
  mutate(cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  count(cyt_num, config_create_date) %>% 
  arrange(config_create_date) %>%
  mutate(config_create_date = as.character(config_create_date)) %>% 
  ggplot() + 
  geom_tile(aes(x = config_create_date, y = cyt_num, fill = n)) + 
  geom_text(aes(x = config_create_date, y = cyt_num, label = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('Config Create Date') + ylab('') + 
  ggtitle("Cytometers and Config Create Dates") + 
  theme(legend.position = 'None')




p2 <- m1_meta_info_df %>% 
  mutate(cyt_num = as.character(cyt_num), 
         cytometer = as.character(cytometer)) %>% 
  count(cyt_num, setup_date) %>% 
  arrange(setup_date) %>%
  mutate(setup_date = as.character(setup_date)) %>% 
  ggplot() + 
  geom_tile(aes(x = setup_date, y = cyt_num, fill = n)) + 
  geom_text(aes(x = setup_date, y = cyt_num, label = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('Setup Date') + ylab('') + 
  ggtitle("Cytometers and Setup Dates") + 
  theme(legend.position = 'None')

p1 / p3 / p2


cyt_nums <- 
  c("R66093700082", "V657338000098", "V657338000021", "V657338000099", 
    "R65822R1018", "R658222R1012", "R66093700081", "R66093700072")

for (spec_cyt_num in cyt_nums) {
  p <- 
    m1_meta_info_df %>% 
    filter(cyt_num == spec_cyt_num) %>% 
    mutate(cyt_num = as.character(cyt_num), 
           cytometer = as.character(cytometer)) %>% 
    count(date, baseline_date) %>% 
    arrange(date) %>%
    mutate(baseline_date = as.character(baseline_date)) %>% 
    ggplot() + 
    geom_point(aes(x = date, y = baseline_date)) + 
    scale_fill_distiller(palette = 'Reds', 
                         direction = 1) + 
    guides(fill = guide_legend(title = '# Flow Ids')) +
    xlab('Experiment Date') + ylab('Baseline Date') + 
    ggtitle(paste0("Cyt: ", spec_cyt_num, " Experiment Date vs Baseline Dates")) + 
    theme(legend.position = 'None') 
  
  print(p)
  
  readline(prompt = 'Press [Enter] to contiue')
}


# This problem is only with cytometer R81


# So what are the date ranges for the batches for each cytometer? 
p1
# Though have to deal with cytometer R81




## Making a table of the batches

# Make sure to account for strange issue in R81 
# This means that we cannot just group by unique 
# baseline dates, but have to make batches by ordering
# by date, then splitting whenever we see the baseline
# date change 

# This makes it an interesting case, also a question for Misha then, 
# should those cases be in the same batch, or different batches? 

# With all cases that have the same baseline date in the same batch 
m1_meta_info_df %>% 
  group_by(cyt_num, baseline_date) %>% 
  mutate(batch_start = min(date), 
         batch_end = max(date)) %>% 
  count(cyt_num, baseline_date, batch_start, batch_end) %>% 
  mutate(cyt_num = unlist(cyt_num)) %>% 
  arrange(batch_start) %>% 
  select(cyt_num, batch_start, batch_end, baseline_date, n) %>% 
  rename('Cytometer_ID' = cyt_num,
         'Batch_Start_Date' = batch_start, 
         'Batch_End_Date' = batch_end, 
         'Baseline_Date' = baseline_date, 
         'Num_Cases' = n)


# With batches defined as contiguous chunks


batch_df_list <- 
  m1_meta_info_df %>% 
  select(date, cyt_num, baseline_date) %>% 
  arrange(date) %>%
  group_by(cyt_num) %>% 
  group_split()

final_batch_df_list <- list()
batch_num <- 1

for (i in seq_along(batch_df_list)) {
  df <- batch_df_list[[i]]
  
  # find location of non-contiguous baseline dates
  break_locations <- c()
  for (j in seq_len(nrow(df)-1)) {
     break_found <- 
       df[j+1, 'baseline_date'] != df[j, 'baseline_date']
     
     if (break_found) {
       break_locations <- c(break_locations, j+1)
     }
  }
  
  if (length(break_locations) != 0) {
    break_locations <- c(1, break_locations, nrow(df))
    batches <- 
      cut(1:nrow(df), 
          break_locations, 
          include.lowest = TRUE, 
          right = FALSE) %>% 
      as.integer()
    
    for (j in seq_along(unique(batches))) {
      batch_start_idx <- min(which(batches == j))
      batch_end_idx   <- max(which(batches == j))
      
      final_batch_df_list[[batch_num]] <- 
          df[batch_start_idx:batch_end_idx, ] %>% 
          mutate(batch_number = batch_num)
      
      batch_num <- batch_num + 1
    }
  } else {
    final_batch_df_list[[batch_num]] <- 
      df %>% 
      mutate(batch_number = batch_num)
    
    batch_num <- batch_num + 1
  }
}

final_batch_df_list

final_batch_df_list %>% 
  do.call(rbind, .) %>% 
  group_by(batch_number) %>% 
  mutate(batch_start = min(date), 
         batch_end = max(date)) %>% 
  ungroup(batch_number) %>% 
  count(cyt_num, baseline_date, batch_start, batch_end, batch_number) %>% 
  mutate(cyt_num = unlist(cyt_num)) %>% 
  arrange(cyt_num, batch_start) %>% 
  select(cyt_num, batch_start, batch_end, baseline_date, n) %>% 
  rename('Cytometer_ID' = cyt_num,
         'Batch_Start_Date' = batch_start, 
         'Batch_End_Date' = batch_end, 
         'Baseline_Date' = baseline_date, 
         'Num_Cases' = n) %>% 
  kable()









