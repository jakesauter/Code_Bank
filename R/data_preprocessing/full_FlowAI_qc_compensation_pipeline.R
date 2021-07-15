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

flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Correct_QC_Flow_Folders"
##################################################################################

if (!dir.exists(flow_dir_save_loc)) 
  dir.create(flow_dir_save_loc)



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

# New: also need to filter by m1 and m2 fcs file size
# lets say needs to be larger than 10MB, though will check this
# threshold when I have this data for each of the cases we have
# downloaded for now
WL3_patients <-
  WL3_patients %>%
  mutate(M1_file_size = file.info(M1_path)$size / (1024^2),
         M2_file_size = file.info(M2_path)$size / (1024^2)) %>%
  filter(M1_file_size > 10,
         M2_file_size > 10)
##################################################################################


##################################################################################
# Saving Flow Folders locally (only with fcs and nlys files that we need)
# as well as parsing and saving compensation matrices
#=================================================================================
cat('\n\nDownloading flow directories and parsing compensation ...\n\n')

flow_directories <- WL3_patients$Flow.Data.Folder

counter <- 0

for (flow_directory in flow_directories) {
  
  counter <- counter + 1
  cat('\n\nFlow Directory: ', counter, " out of: ", length(flow_directories), '\n\n')

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
    
    copied_files <- 
      list.files(local_flow_directory,
                 pattern = 'M1.*\\.fcs$|M2.*\\.fcs$|.*\\.nlys$',
                 full.names = TRUE)    
    
    # Check for error in file transfer 
    if (length(copied_files) != length(files_to_copy)) {
      error(paste0("Not all files copied for flow dir: ", flow_directory, 
                   " to location: ", local_flow_directory))
    } else {
      files_to_copy <- sort(files_to_copy)
      copied_files <- sort(copied_files)
      
      for (i in seq_along(files_to_copy)) {
        file_to_copy_size <- file.info(files_to_copy[i])$size / (1024^2)
        copied_file_size <- file.info(copied_files[i])$size / (1024^2)
        if (!near(file_to_copy_size, copied_file_size, 1))
          error(paste0('Error in copying file: ', files_to_copy[i], 
                       ' to file: ', copied_files[i]))
      }
    }
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
  
  Sys.sleep(time = 1)
}
##################################################################################


##################################################################################
# QCing all m1 and m2 files in flow directories with flowAI and saving high
# quality fcs file
#=================================================================================
cat('\n\nPerforming quality checks and saving ...\n\n')

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
                         remove_from = 'FS',
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
                         remove_from = 'FS',
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
# Remove entire flow directories from cohort that had too many 
# events removed from quality control, as I am getting errors in 
# applying compensation to files with no events due to 100% filtering
#=================================================================================
# 
# flow_directories <- 
#   list.files(flow_dir_save_loc, 
#              full.names = TRUE)
# 
# 
# qc_results_list <- vector('list', length(flow_directories))
# 
# for (i in seq_along(flow_directories)) {
#   
#   flow_directory <- flow_directories[[i]]
#   
#   # Parse and better save QC Results
#   qc_results <- 
#     read.table(file.path(flow_directory, 
#                          'flowAI_QC_Results',
#                          'QCmini.txt'), 
#                skip = 1, 
#                sep = '\t') %>% 
#     .[, c(1, 3, 5, 6, 7)] %>% 
#     set_colnames(c('filename', 
#                    'total_filtered_perc',
#                    'FR_filtered_perc', 
#                    'FS_filtered_perc', 
#                    'FM_filtered_perc (not used)'))
#   
#   qc_results$filename <-
#     file.path(flow_directory, qc_results$filename)
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
#   arrange(total_filtered_perc)
# 
# flow_dirs_passed_qc <- 
#   full_qc_results %>% 
#   filter(total_filtered_perc < 80) %>% 
#   mutate(dirname = dirname(filename)) %>% 
#   count(dirname) %>% 
#   filter(n > 1) %>% 
#   .$dirname
##################################################################################


##################################################################################
# Applying compensation to high quality fcs files
#=================================================================================
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

