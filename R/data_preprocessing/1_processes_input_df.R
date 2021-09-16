#!/bin/Rscript 


################################################
# Input: dat_file -- location of cohort csv
#===============================================
# Output: `jake_processed` appended filename in same 
#  location as original csv file, with flow folder
#  and M1/M2 filepaths resolved for each patient, as 
#  well as reasons that the files could not be found
#  or accessed if that is the case.
###############################################

#########################################################################
# Modify input data_file parameter and run for new cohort definitions
#========================================================================
data_file = '~/Documents/Woodlist/Input_Flow_DFs/flow_df_ids_052621.csv'

output_csv_path <- 
  file.path(dirname(data_file), 
            paste0('jake_processed_', basename(data_file)))
#########################################################################


library(fs)
library(DT)
library(dplyr)
library(knitr)
library(ggplot2)
library(stringr)
library(magrittr)
library(patchwork)
library(lubridate)
library(reticulate)


find_patient_toplevel_folders <- function(MRNs, 
                                         server_path = '/Volumes/PATIENTS/') {
  
  input_MRNs <- MRNs
  
  MRNs <- MRNs %>% 
    as.character() %>% 
    unique()
  
  paths <- rep(NA, length(MRNs))
  names(paths) <- MRNs
  
  for (letter in letters) {
    letter <- toupper(letter)    
    path <- file.path(server_path, letter)
    dir_names <- list.files(path = path)
    
    for (MRN in MRNs) {
      # if path has already been found, skip
      if(!is.na(paths[MRN])) next
      
      grep_res <- grep(MRN, dir_names, value = TRUE)
      
      if (length(grep_res) == 1) {
        # if only one folder was found (success)
        paths[MRN] <- 
          fs::path(server_path, letter, grep_res) %>% 
          stringr::str_replace_all('//', '/')
      } else if(length(grep_res) > 1) {
          # if more than one folder was found
          paths[MRN] <- 
            fs::path(server_path, letter, grep_res) %>% 
            paste(collapse = ', ')
      }
    }
    
    Sys.sleep(0.01)
  }
  
  output_paths <- 
    sapply(input_MRNs, 
           function(MRN) paths[as.character(MRN)])
  
  invisible(output_paths)
}

find_patient_flow_directory <- function(paths, flow_ids) {

  flow_id_folders <- rep('', length(paths))
  
  for (i in seq_along(paths)) {
  
    # Will work for patients with one or multiple files, 
    # as if they do not have a comma, only the filepath will
    # be returned
    split_paths <- str_split(paths[i], ',') 
    
    all_patient_files <- 
      sapply(split_paths, function(path) {
        # need to trim the path as there is a space
        # at the start of any path other than the first
        # so that printing the dataframe is more user friendly
        path <- str_trim(path)
        list.files(path = path, 
                  full.names = TRUE)
        })
  
    flow_id_folder <- 
      grep(flow_ids[i], 
           all_patient_files, 
           value = TRUE) 
    
    if (length(flow_id_folder) == 0) {
      flow_id_folder <- 'No Flow ID Folder Found'
    }
    
    # Sometimes the flow id folder is copied to the 
    # renamed paitient folder, though the patient flow
    # files are always exactly the same, so just picking 
    # the first one here
    if(length(flow_id_folder) > 1) {
      flow_id_folder <- flow_id_folder[1]
    }
    
    flow_id_folders[i] <- flow_id_folder
    Sys.sleep(0.001)
  }
  
  flow_id_folders
}

find_patient_M1_and_M2 <- 
  function(flow_id_folders) {
        
    df <-
        data.frame(flow_id_folder = flow_id_folders, 
                   M1_path = rep('FOLDER DID NOT EXIST', length(flow_id_folders)), 
                   M2_path = rep('FOLDER DID NOT EXIST', length(flow_id_folders)), 
                   stringsAsFactors = FALSE)
    
    for (i in seq_along(flow_id_folders)) {
     
      if (!dir.exists(flow_id_folders[i])) next
        
      m1_path <-
          grep('M1.*\\.fcs', 
               list.files(flow_id_folders[i], 
                          full.names = TRUE), 
               ignore.case = TRUE, 
               value = TRUE)
      
      m2_path <-
            grep('M2.*\\.fcs', 
                 list.files(flow_id_folders[i], 
                            full.names = TRUE), 
                 ignore.case = TRUE, 
                 value = TRUE)
      
      if (length(m1_path) > 1) {
        m1_path <- 
          m1_path %>% 
          paste(collapse = ', ')
      } else if (length(m1_path) == 0) {
        m1_path <- 'NO M1 FOUND' 
      }      
      
      if (length(m2_path) > 1) {
        m2_path <- 
          m2_path %>% 
          paste(collapse = ', ')
      } else if (length(m2_path) == 0) {
        m2_path <- 'NO M2 FOUND' 
      }
      
      df[i, 'M1_path'] <- m1_path
      df[i, 'M2_path'] <- m2_path
      Sys.sleep(0.001)
    }
    
  df[, c('M1_path', 'M2_path')]                 
}

patient_df <- read.csv(data_file, 
                       stringsAsFactors = FALSE) 

patient_df <- 
  patient_df %>% 
  mutate(Path.Procedure.Date = as_date(Path.Procedure.Date)) %>% 
  mutate(Top.Level.Folder = 
           find_patient_toplevel_folders(MRN)) %>% 
  mutate(Flow.Data.Folder = 
           find_patient_flow_directory(Top.Level.Folder, 
                                      Accession.Number)) %>% 
  mutate(Permission.To.Read = 
           as.logical(file.access(Flow.Data.Folder, mode = 4) + 1))

paths_df <- 
  find_patient_M1_and_M2(patient_df$Flow.Data.Folder)

# Bind M1 and M2 paths to patient data frame
patient_df <- 
  patient_df %>% 
  cbind(paths_df)

# Exclusion Reasons
patient_df <- 
  patient_df %>% 
  mutate(Exclusion.Reason = 
           if_else(is.na(Top.Level.Folder), 
                   'Top Level Folder Not Found', 'Not Excluded')) %>% 
   mutate(Exclusion.Reason = 
           if_else(Flow.Data.Folder == 'No Flow ID Folder Found' & 
                   Exclusion.Reason == 'Not Excluded',  
                    Flow.Data.Folder, Exclusion.Reason)) %>% 
   mutate(Exclusion.Reason = 
           if_else(Permission.To.Read == FALSE & 
                    Exclusion.Reason == 'Not Excluded', 
                   'Flow Directory Read Permissions', 
                   Exclusion.Reason)) %>% 
  mutate(Exclusion.Reason = 
           if_else(Permission.To.Read == TRUE & 
                      M1_path == 'NO M1 FOUND' & 
                      M2_path == 'NO M2 FOUND' & 
                      Exclusion.Reason == 'Not Excluded', 
                   'No M1.fcs or M2.fcs found', 
                    Exclusion.Reason)) %>% 
  mutate(Exclusion.Reason = 
           if_else(Permission.To.Read == TRUE &
                     M1_path == 'NO M1 FOUND' &
                     M2_path != 'NO M2 FOUND' & 
                     Exclusion.Reason == 'Not Excluded', 
                    'No M1.fcs found, but M2.fcs found', 
                    Exclusion.Reason)) %>% 
  mutate(Exclusion.Reason = 
           if_else(Permission.To.Read == TRUE &
                     M1_path != 'NO M1 FOUND' &
                     M2_path == 'NO M2 FOUND' & 
                     Exclusion.Reason == 'Not Excluded', 
                    'M1.fcs found, but no M2.fcs found', 
                    Exclusion.Reason))  

patient_df <- 
  patient_df %>% 
  arrange(Exclusion.Reason)


# Write out final df
write.csv(patient_df, 
          output_csv_path, 
          row.names = FALSE)


