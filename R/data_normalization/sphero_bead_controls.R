

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

sphero_fcs_path <- 
  "/Volumes/Patients FCS files/2021/F1/12032021/Fortessa I SpherotechLOT#AM01 WEEK OF 030821/Spherotech Lot AL01_RPT_031221.fcs"

exp_fcs_path <- 
  "/Volumes/Patients FCS files/2021/F1/12032021/F21-2446 KAPPEN 031221/F1F21-2446KAPP_M1.fcs"

sphero_fcs_data <- read.FCS(sphero_fcs_path, transformation = FALSE)
exp_fcs_data    <- read.FCS(exp_fcs_path, transformation = FALSE)

meta <- description(sphero_fcs_data)


sphero_fcs_data %>% 
  parameters() %>% 
  pData()

exp_fcs_data %>% 
  parameters() %>% 
  pData()

desc <- 
  description(sphero_fcs_data)

desc$CYTNUM
desc$`$DATE`



cd_45_sphero <- 
  sphero_fcs_data %>% 
  exprs() %>% 
  .[, 'BUV805-A'] 


sphero_fcs_data %>% 
  exprs() %>% 
  .[, 'BUV805-A'] %>% 
  .[. > 0] %>% 
  log() %>% 
  density() %>% 
  plot()

abline(v = log(124), 
       col = 'red', 
       lty = 2)

abline(v = log(150), 
       col = 'red', 
       lty = 2)










# Manually going to find the spherotech runs for the 10
# samples I used for UMAP

m1_paths <- 
  c("/Volumes/PATIENTS/S/SAUL_JULIE_35536305/F21-766 SAUL 012521/S3F21-766SAUL_M1.fcs", 
    "/Volumes/PATIENTS/K/KAPOOR_GOURI ARVIND_38014929/F21-2241 KAPOOR 030821/S3F21-2241KAPO_M1.fcs", 
    "/Volumes/PATIENTS/S/SAUL_JULIE_35536305/F21-767 SAUL 012521/S2F21-767SAUL_M1.fcs", 
    "/Volumes/PATIENTS/B/BEARD_LORIE_38130682/F21-1865 BEARD 022521/S2F21-1865BEAR_M1.fcs", 
    "/Volumes/PATIENTS/E/EISNER_ANNIE_38098201/F21-1288 EISNER 020921/F1F21-1288EISN_M1.fcs", 
    "/Volumes/PATIENTS/R/RAMIREZ_GRACE L_38106336/F20-10179 RAMIREZ 120320/F1F20-10179RAMI_M1.fcs", 
    "/Volumes/PATIENTS/M/MORIN_HOPE CATHERINE_38127184/F20-8371 MORIN 101220/F2F20-8371MORI M1 001.fcs", 
    "/Volumes/PATIENTS/S/SMILEY_ZORIEAN_38132309/F20-10111 SMILEY 120120/F2F20-10111SMIL M1 001.fcs", 
    "/Volumes/PATIENTS/S/SIEGEL_JUDITH_00370285/F21-2131 SIEGEL 030321/S1F21-2131SIEG_M1.fcs"
  )

meta_info <- 
  m1_paths %>% 
  lapply(function(x) {
    x %>% 
      read.FCS(transformation = FALSE, 
               which.lines = 1) %>% 
      description()
  })

useful_info <- 
  meta_info %>% 
  lapply(function(x) {
    year <- 
      x$`$DATE` %>% 
      lubridate::dmy() %>% 
      year
    
    date <- 
      x$`$DATE` %>% 
      lubridate::dmy() %>% 
      as.character()
    
    cyt_num <- 
      x$CYTNUM
    
    list(cyt_num, year, date)
  }) %>% 
 do.call(rbind, .)

useful_info %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  set_colnames(c('cyt_num', 'year', 'date')) %>% 
  mutate(
    cytometer_dir =
      case_when(
        cyt_num == "V657338000021" ~ "CANTO IV",
        cyt_num == "V657338000098" ~ "CANTO V",
        cyt_num == "V657338000099" ~ "CANTO VI",
        year >= 2018 & year <= 2020 & cyt_num == "R658222R1012" ~ "Fortessa I",
        year >= 2018 & year <= 2020 & cyt_num == "R65822R1018"  ~ "Fortessa II",
        year >= 2021 & cyt_num == "R658222R1012"  ~ "F1",
        year >= 2021 & cyt_num == "R65822R1018"   ~ "F2",
        cyt_num == "R66093700072" ~ "S1",
        cyt_num == "R66093700081" ~ "S2",
        cyt_num == "R66093700082" ~ "S3"
  )) %>% 
  select(cyt_num, year, cytometer_dir, date)



sphero_fcs <- 
  c("/Volumes/Patients FCS files/2021/S3/25012021/DAILY QC SYMPHONY III LOT AM01 WEEK OF 012521/LOT AL01_ADJUST_012521_25012021094612.fcs",
    "/Volumes/Patients FCS files/2021/S3/08032021/DAILY QC SYMPHONY III LOTAM01 WEEK OF 030821/LOT AM01_030821.fcs", 
    "/Volumes/Patients FCS files/2021/S2/25012021/DAILY QC SYMPHONY II AM01 WEEK OF 012521/LOT AL01_adjusted_012521_25012021104441.fcs", 
    "", 
    "", 
    "/Volumes/Patients FCS files/2020/Fortessa I/03122020/Fortessa I Spherotech LOT AM01 WEEK OF 113020/Spherotech Lot AM01_ADJUST_120320.fcs", 
    "13102020 Missing"
    "1201202 Also Missing"
  )




# have to map machine name by year and cytometer number
# m1_meta_info_df <- readRDS('data/m1_meta_info_df.Rds')
# 
# m1_meta_info_df <- 
#   m1_meta_info_df %>% 
#   mutate(year = year(date), 
         # cytometer_dir =
         #   case_when(
         #     cyt_num == "V657338000021" ~ "CANTO IV",
         #     cyt_num == "V657338000098" ~ "CANTO V",
         #     cyt_num == "V657338000099" ~ "CANTO VI",
         #     year >= 2018 & year <= 2020 & cyt_num == "R658222R1012" ~ "Fortessa I",
         #     year >= 2018 & year <= 2020 & cyt_num == "R65822R1018"  ~ "Fortessa II",
         #     year >= 2021 & cyt_num == "R658222R1012"  ~ "F1",
         #     year >= 2021 & cyt_num == "R65822R1018"   ~ "F2",
         #     cyt_num == "R66093700072" ~ "S1",
         #     cyt_num == "R66093700081" ~ "S2",
         #     cyt_num == "R66093700082" ~ "S3"
         #   ))
# 
# 
# root_dir <- "/Volumes/Patients FCS files/"
# 
# 
# year <- m1_meta_info_df$year[1]
# cyt_dir <- m1_meta_info_df$cytometer_dir[1]
# flow_id <- m1_meta_info_df$flow_id[1]
# 
# 
# 
# dir <- 
#   file.path(root_dir, year, cyt_dir) %>% 
#   fs::path_real()





