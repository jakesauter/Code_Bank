#!/bin/Rscript

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

### Input data
patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv', 
           stringsAsFactors = FALSE) 

patient_df <- 
  patient_df %>% 
  filter(Exclusion.Reason == "Not Excluded", 
         !str_detect(M1_path, ','), 
         !str_detect(M2_path, ','))
         


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

for (i in seq_len(nrow(patient_df))) {
  
  flow_directory <- patient_df[i, "Flow.Data.Folder"]
  
  cat('Processing flow directory:', i, 'out of:', nrow(patient_df))
  cat('\nFlow directory: ', flow_directory, '\n\n')
  
  m1_fcs_file <- 
    patient_df[i, "M1_path"]
  
  m2_fcs_file <- 
    patient_df[i, "M2_path"]
  
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


saveRDS(m1_meta_info_df, 'data/m1_meta_info_df.Rds')
saveRDS(m2_meta_info_df, 'data/m2_meta_info_df.Rds')



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


p3 / p4


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


p5 / p6

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



factored_config_name <-
  m1_meta_info_df %>% 
  select(config_create_date,
         config_name) %>% 
  mutate(config_create_date = lubridate::ymd(config_create_date)) %>% 
  arrange(config_create_date) %>% 
  .$config_name %>% 
  str_replace_all('-', ' ') %>% 
  unlist() %>% 
  factor(levels = unique(.))


m1_meta_info_df %>% 
  mutate(cytometer = as.character(cytometer), 
         config_name = factored_config_name, 
         cyt_num = unlist(cyt_num)) %>% 
  count(cyt_num, config_name) %>% 
  mutate(config_name = as.integer(config_name)) %>% 
  ggplot() + 
  geom_tile(aes(x = config_name, y = cyt_num, fill = n)) + 
  scale_x_continuous(breaks = 1:length(levels(factored_config_name)), 
                     labels = levels(factored_config_name)) + 
  geom_text(aes(x = config_name, y = cyt_num, label = n)) + 
  scale_fill_distiller(palette = 'Reds', 
                       direction = 1) + 
  guides(fill = guide_legend(title = '# Flow Ids')) +
  xlab('') + ylab('') + 
  ggtitle("Cytometers and Config Names") + 
  theme(legend.position = 'None', 
        axis.text.x = element_text(angle = -35, hjust = 0), 
        plot.margin=unit(c(.5,4.3,.5,.5),"cm"))





# Need to confirm that the date of the experiments
# that used varying configs for the Fortessa's
# (12, 18) changed between 7/13/2020 and 10/5/2020


m1_meta_info_df %>% 
  select(config_create_date,
         config_name) %>% 
  mutate(config_create_date = lubridate::ymd(config_create_date)) %>% 
  arrange(config_create_date) %>% 
  count(config_create_date, config_name) %>% 
  select(config_name, config_create_date) %>% 
  kable()


m1_meta_info_df <- 
  m1_meta_info_df %>% 
  rowwise() %>% 
  mutate(num_params = length(param_names))

m1_meta_info_df %>% 
  mutate(cyt_num = unlist(cyt_num), 
         config_name = unlist(config_name)) %>% 
  count(cyt_num, num_params, config_name) %>% 
  select(config_name, cyt_num, num_params) %>% 
  arrange(num_params)


# Config change by experiment date

df <- 
  m1_meta_info_df %>% 
  group_by(cyt_num, config_name) %>% 
  arrange(date) %>% 
  slice_head(n = 1) %>% 
  rename(earliest_exp_date = date) %>% 
  ungroup() 


exp_count <- 
  m1_meta_info_df %>% 
  count(cyt_num, config_name) %>% 
  .$n

latest_exp_date <- 
  m1_meta_info_df %>% 
  group_by(cyt_num, config_name) %>% 
  arrange(date) %>% 
  slice_tail(n = 1) %>% 
  ungroup() %>% 
  .$date

df$latest_exp_date <- latest_exp_date
df$n <- exp_count
  

df %>% 
  arrange(earliest_exp_date) %>% 
  mutate(cyt_num = unlist(cyt_num), 
         config_name = unlist(config_name)) %>% 
  select(cyt_num, earliest_exp_date, latest_exp_date, config_name, n) %>% 
  arrange(cyt_num, earliest_exp_date) %>% 
  kable()


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









