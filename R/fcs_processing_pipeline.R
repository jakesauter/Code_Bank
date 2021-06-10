library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(flowAI)
filter <- dplyr::filter

##################################################################################
# Parameter Variables
#=================================================================================
INPUT_FLOW_DF_CSV <- "/Users/sauterj1/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv"

parse_comp_mat_python_path <-
  "/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/parse_flow_folder_compensation.py"

flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Flow_Folders"
##################################################################################




##################################################################################
# Reading in initial Patients DF and filtering for only Woodlist 3 directories
#  that have no ambiguity in nlys files or m1/m2 fcs files
#=================================================================================
cat('\n\nDetermining Woodlist 3 analysis directories ...\n\n')
patient_df <- 
  read.csv(INPUT_FLOW_DF_CSV, 
           stringsAsFactors = FALSE)

WL3_patients <-
  patient_df %>%
  rowwise() %>%
  mutate(WL3_nlys_file_count =
           length(
             list.files(Flow.Data.Folder,
                        pattern = '\\.nlys$')
           )) 

WL3_patients <- 
  WL3_patients %>% 
  filter(WL3_nlys_file_count == 1) %>% 
  filter(!str_detect(M1_path, ','), 
         !str_detect(M2_path, ',')) 
##################################################################################


##################################################################################
# Saving Flow Folders locally (only with fcs and nlys files that we need)
# as well as parsing and saving compensation matrices
#=================================================================================
cat('\n\nDownloading flow directories and parsing compensation ...\n\n')

flow_directories <- WL3_patients$Flow.Data.Folder

for (flow_directory in flow_directories) {
  
  cat("\nDownloading Flow Directory: ", flow_directory, "\n\n")
  
  # New path we will save to
  local_flow_directory <- 
    file.path(flow_dir_save_loc, 
              basename(flow_directory))
  
  # Need to copy flow folders to local 
  # in order to modify and save compensation matrices
  if (!dir.exists(local_flow_directory)) {
    
    dir.create(local_flow_directory)
    
    files_to_copy <- 
      list.files(flow_directory, 
                 pattern = 'M1.*\\.fcs$|M2.*\\.fcs$|.*\\.nlys$', 
                 full.names = TRUE)
    
    file.copy(files_to_copy, 
              local_flow_directory)
  }
  
  if (!dir.exists(
    file.path(local_flow_directory, 
              'compensation_matrices')
  )) {
    command <- 
      paste0('/Users/sauterj1/anaconda3/bin/python3.8 ', 
             parse_comp_mat_python_path, 
             " \"", local_flow_directory, "\"")
    
    system(command)
  }
}
##################################################################################


##################################################################################
# QCing all m1 and m2 files in flow directories with flowAI and saving high 
# quality fcs file
#=================================================================================
cat('\n\nPerforming quality checks and saving ...\n\n')

library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(flowAI)
filter <- dplyr::filter

INPUT_FLOW_DF_CSV <- "/Users/sauterj1/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv"

parse_comp_mat_python_path <-
  "/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/parse_flow_folder_compensation.py"

flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Flow_Folders"


flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE) 

counter <- 0


for (flow_directory in flow_directories) {
  
  counter <- counter + 1
  cat('\n\nFlow Directory: ', counter, " out of: ", length(flow_directories), '\n\n')
  
  flowAI_QC_output_folder <- 
    file.path(flow_directory, 'flowAI_QC_Results')
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[!str_detect(., 'compensated')] %>% 
    .[!str_detect(., 'high_quality_events')]
  
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
  
  
  m1_fcs <- flowCore::read.FCS(m1_fcs_file, 
                               transformation = FALSE)
  m1_metadata <- pData(parameters(m1_fcs))
  
  m2_fcs <- flowCore::read.FCS(m2_fcs_file, 
                               transformation = FALSE)
  m2_metadata <- pData(parameters(m2_fcs))
  
  # Only use channel names with markers to check for 
  # signal quality
  exclude_channels <- 
    m1_metadata$name[is.na(m1_metadata$desc)]
  
  exclude_channels <- 
    exclude_channels[exclude_channels != "Time"]
  
  m1_fcs_flowai_qc <- 
    flowAI::flow_auto_qc(m1_fcs, 
                         remove_from = 'FR_FS', 
                         ChExcludeFS = exclude_channels,
                         ChExcludeFM = exclude_channels,
                         folder_results = flowAI_QC_output_folder, 
                         html_report = "_QC",
                         max_cptFS = 30,
                         fcs_QC = FALSE)
  
  
  # Only use channel names with markers to check for 
  # signal quality
  exclude_channels <- 
    m2_metadata$name[is.na(m2_metadata$desc)]
  
  exclude_channels <- 
    exclude_channels[exclude_channels != "Time"]
  
  m2_fcs_flowai_qc <- 
    flowAI::flow_auto_qc(m2_fcs, 
                         remove_from = 'FR_FS', 
                         ChExcludeFS = exclude_channels,
                         ChExcludeFM = exclude_channels,
                         folder_results = flowAI_QC_output_folder, 
                         html_report = "_QC",
                         max_cptFS = 30,
                         fcs_QC = FALSE)
  
  # Save high quality m1 and m2 fcs 
  m1_high_qual_fcs_filename <- 
    paste0(str_split(m1_fcs_file, '\\.fcs')[[1]][1], 
           '_high_quality_events.fcs')
  
  flowCore::write.FCS(m1_fcs_flowai_qc, 
                      m1_high_qual_fcs_filename)
  
  m2_high_qual_fcs_filename <- 
    paste0(str_split(m2_fcs_file, '\\.fcs')[[1]][1], 
           '_high_quality_events.fcs')
  
  flowCore::write.FCS(m2_fcs_flowai_qc, 
                      m2_high_qual_fcs_filename)
}
##################################################################################


##################################################################################
# Applying compensation to high quality fcs files
#=================================================================================
library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(flowAI)
filter <- dplyr::filter

INPUT_FLOW_DF_CSV <- "/Users/sauterj1/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv"

parse_comp_mat_python_path <-
  "/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/parse_flow_folder_compensation.py"

flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

cat('\n\nApplying compensation and saving ...\n\n')

flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE)

for (flow_directory in flow_directories) {
  
  cat('\nWorking on flow directory: ', flow_directory, '\n\n')
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'high_quality_events.fcs')]
  
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
  
  m1_fcs <- flowCore::read.FCS(m1_fcs_file, 
                               transformation = FALSE)
  m2_fcs <- flowCore::read.FCS(m2_fcs_file, 
                               transformation = FALSE)
  
  m1_comp_mat <- 
    m1_fcs_file %>% 
    basename() %>% 
    str_split('_high_quality_events\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('.csv') %>% 
    file.path(flow_directory, 
              'compensation_matrices',
              .) %>% 
    read.csv(header = FALSE) %>% 
    set_colnames(pData(parameters(m1_fcs))$name)
  
  m2_comp_mat <- 
    m2_fcs_file %>% 
    basename() %>% 
    str_split('_high_quality_events\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('.csv') %>% 
    file.path(flow_directory, 
              'compensation_matrices',
              .) %>% 
    read.csv(header = FALSE) %>% 
    set_colnames(pData(parameters(m2_fcs))$name)
  
  # compensate 
  m1_fcs <- compensate(m1_fcs, m1_comp_mat)
  m2_fcs <- compensate(m2_fcs, m2_comp_mat)
  
  compensated_m1_file <- 
    str_c(
      str_split(m1_fcs_file, '\\.')[[1]][1], 
      '_compensated.fcs'
    )
  
  compensated_m2_file <- 
    str_c(
      str_split(m2_fcs_file, '\\.')[[1]][1], 
      '_compensated.fcs'
    )
  
  
  flowCore::write.FCS(m1_fcs, 
                      compensated_m1_file)
  
  flowCore::write.FCS(m2_fcs, 
                      compensated_m2_file)
}
##################################################################################


# Gather all QC Statistics
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE)


qc_results_list <- vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  
  # Parse and better save QC Results
  qc_results <- 
    read.table(file.path(flow_directory, 
                         'flowAI_QC_Results',
                         'QCmini.txt'), 
               skip = 1, 
               sep = '\t') %>% 
    .[, c(1, 3, 5, 6, 7)] %>% 
    set_colnames(c('filename', 
                   'total_filtered_perc',
                   'FR_filtered_perc', 
                   'FS_filtered_perc', 
                   'FM_filtered_perc (not used)'))
  
  qc_results$filename <-
    file.path(flow_directory, qc_results$filename)
  
  
  qc_results_list[[i]] <- qc_results
  
}

full_qc_results <- do.call(rbind, qc_results_list)

full_qc_results <- 
  full_qc_results %>% 
  dplyr::distinct() %>%
  arrange(total_filtered_perc)

print(full_qc_results)

##################################################################################


##################################################################################
# QCing all m1 and m2 files WITH FLOWCLEAN in flow directories with flowAI and saving high 
# quality fcs file
#=================================================================================
# cat('\n\nPerforming quality checks and saving ...\n\n')
# 
# library(dplyr)
# library(stringr)
# library(magrittr)
# library(flowCore)
# library(flowClean)
# filter <- dplyr::filter
# 
# INPUT_FLOW_DF_CSV <- "/Users/sauterj1/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv"
# 
# parse_comp_mat_python_path <-
#   "/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/parse_flow_folder_compensation.py"
# 
# flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Flow_Folders"
# 
# 
# flow_directories <- 
#   list.files(flow_dir_save_loc, 
#              full.names = TRUE) 
# 
# counter <- 0
# 
# 
# for (flow_directory in flow_directories) {
#   
#   counter <- counter + 1
#   cat('\n\nFlow Directory: ', counter, " out of: ", length(flow_directories), '\n\n')
#   
#   cat('Directory: ', flow_directory, '\n')
#   
#   fcs_files <- 
#     list.files(flow_directory, 
#                pattern = "\\.fcs") %>% 
#     .[!str_detect(., 'compensated')] %>% 
#     .[!str_detect(., 'high_quality_events')]
#   
#   m1_fcs_file <- 
#     str_detect(fcs_files, 
#                'M1|m1') %>% 
#     fcs_files[.] %>% 
#     file.path(flow_directory, .)
#   
#   m2_fcs_file <- 
#     str_detect(fcs_files, 
#                'M2|m2') %>% 
#     fcs_files[.] %>% 
#     file.path(flow_directory, .)
#   
#   
#   m1_fcs <- flowCore::read.FCS(m1_fcs_file)
#   m1_metadata <- pData(parameters(m1_fcs))
#   
#   m2_fcs <- flowCore::read.FCS(m2_fcs_file)
#   m2_metadata <- pData(parameters(m2_fcs))
#   
#   # Only use channel names with markers to check for 
#   # signal quality
#   include_channels <- which(!is.na(m1_metadata$desc))
#   
#   flowclean_output_prefix <- 
#     file.path(flow_directory, 
#               'flowClean_QC_Results') %>% 
#     paste0('/', basename(m1_fcs_file))
#   
#   if (!dir.exists(dirname(flowclean_output_prefix))) {
#     dir.create(dirname(flowclean_output_prefix))
#   }
#   
#   
#   m1_good_vs_bad <- clean(m1_fcs, 
#                           vectMarkers=include_channels,
#                           filePrefixWithDir=flowclean_output_prefix,
#                           returnVector = TRUE, 
#                           ext = 'fcs', 
#                           diagnostic=TRUE)
#   
#   exprs(m1_fcs) <- 
#     exprs(m1_fcs)[m1_good_vs_bad <= 10e3 ,]
#   
#   ## Save QC Results
#   m1_perc_removed <- 
#     ifelse(m1_good_vs_bad <= 10e3, 1, 10e3) %>% 
#     {length(which(. > 10e3)) / 
#         length(.)}
#   
#   
#   # Only use channel names with markers to check for 
#   # signal quality
#   include_channels <- which(!is.na(m2_metadata$desc))
#   
#   flowclean_output_prefix <- 
#     file.path(flow_directory, 
#               'flowClean_QC_Results') %>% 
#     paste0('/', basename(m2_fcs_file))
#   
#   m2_good_vs_bad <- clean(m2_fcs, 
#                           vectMarkers=include_channels,
#                           filePrefixWithDir=flowclean_output_prefix,
#                           returnVector = TRUE, 
#                           ext = 'fcs'
#                           diagnostic=TRUE)
#   
#   exprs(m2_fcs) <- 
#     exprs(m2_fcs)[m2_good_vs_bad <= 10e3 ,]
#   
#   m2_perc_removed <- 
#     ifelse(m2_good_vs_bad <= 10e3, 1, 10e3) %>% 
#     {length(which(. > 10e3)) / 
#         length(.)}
#   
#   # Save QC Data
#   
#   qc_file <- file.path(flow_directory, 
#                        'flowClean_QC_Results', 
#                        'flowClean_removal_rate.csv')
#   
#   data.frame(filename = c(m1_fcs_file, m2_fcs_file), 
#              perc_removed = c(m1_perc_removed, m2_perc_removed)) %>% 
#     write.csv(qc_file, row.names = FALSE)
#   
#   # Save high quality m1 and m2 fcs 
#   m1_high_qual_fcs_filename <- 
#     paste0(str_split(m1_fcs_file, '\\.fcs')[[1]][1], 
#            '_high_quality_events.fcs')
#   
#   flowCore::write.FCS(m1_fcs, 
#                       m1_high_qual_fcs_filename)
#   
#   m2_high_qual_fcs_filename <- 
#     paste0(str_split(m2_fcs_file, '\\.fcs')[[1]][1], 
#            '_high_quality_events.fcs')
#   
#   flowCore::write.FCS(m2_fcs, 
#                       m2_high_qual_fcs_filename)
# }
# ##################################################################################
# 
# 
# 
# # Gather all QC Statistics
# 
# qc_results_list <- vector('list', length(flow_directories))
# 
# for (i in seq_along(flow_directories)) {
#   
#   flow_directory <- flow_directories[[i]]
#   
#   # Parse and better save QC Results
#   qc_results <- 
#     read.csv(file.path(flow_directory, 
#                          'flowClean_QC_Results',
#                          'flowClean_removal_rate.csv'), 
#              header = TRUE)
#   
# 
#   qc_results_list[[i]] <- qc_results
#   
# }
# 
# full_qc_results <- do.call(rbind, qc_results_list)
# 
# full_qc_results <- 
#   full_qc_results %>% 
#   dplyr::distinct() %>%
#   arrange(perc_removed)
# 
# print(full_qc_results)

