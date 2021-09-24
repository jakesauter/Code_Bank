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

flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Full_Correct_QC_Flow_Folders"
##################################################################################

if (!dir.exists(flow_dir_save_loc)) {
  dir.create(flow_dir_save_loc)
}



##################################################################################
# Reading in initial Patients DF and filtering for only Woodlist 3 directories
#  that have no ambiguity in nlys files or m1/m2 fcs files
#=================================================================================
cat('\n\nDetermining patients with non-ambiguous M1 and M2 files ...\n\n')

patient_df <-
  read.csv(INPUT_FLOW_DF_CSV,
           stringsAsFactors = FALSE) %>%
  filter(Exclusion.Reason == "Not Excluded")

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
  latest_file_idx <-
    filenames %>% 
    file.mtime() %>% 
    lubridate::as_date() %>% 
    which.max()
  
  return(filenames[latest_file_idx])
}


# New: also need to filter by m1 and m2 fcs file size
# lets say needs to be larger than 10MB, though will check this
# threshold when I have this data for each of the cases we have
# downloaded for now
filtered_patient_df <-
  patient_df %>%
    rowwise() %>% 
    mutate(
      M1_path = {
        paths <- 
          M1_path %>% str_split(',') %>% 
            .[[1]] %>% str_trim()
        
        choose_file_to_use(paths) 
      },
      M2_path = {
        paths <- 
          M2_path %>% str_split(',') %>% 
          .[[1]] %>% str_trim()        
        choose_file_to_use(paths) 
      }, 
    )


filtered_patient_df <- 
  filtered_patient_df %>% 
  mutate(M1_file_size = file.info(M1_path)$size / (1024^2),
         M2_file_size = file.info(M2_path)$size / (1024^2)) %>% 
  filter(!is.na(M1_file_size) & M1_file_size > 10,
         !is.na(M1_file_size) & M2_file_size > 10)

    
  

  
orig_patient_df <- patient_df
patient_df <- filtered_patient_df


##################################################################################

##################################################################################
# Saving Flow Folders locally (only with fcs and nlys files that we need)
# as well as parsing and saving compensation matrices
#=================================================================================
cat('\n\nDownloading flow directories ...\n\n')

flow_directories <- patient_df$Flow.Data.Folder

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
                 # pattern = 'M1.*\\.fcs$|M2.*\\.fcs$|.*\\.nlys$|^[^.]*$',
                 full.names = TRUE)
    
    file.copy(files_to_copy,
              local_flow_directory)
    
    copied_files <- 
      list.files(local_flow_directory,
                 # pattern = 'M1.*\\.fcs$|M2.*\\.fcs$|.*\\.nlys$|^[^.]*$',
                 full.names = TRUE)    
    
    # Check for error in file transfer 
    if (length(copied_files) != length(files_to_copy)) {
      message(paste0("Not all files copied for flow dir: ", flow_directory, 
                   " to location: ", local_flow_directory))
    } else {
      files_to_copy <- sort(files_to_copy)
      copied_files <- sort(copied_files)
      
      for (i in seq_along(files_to_copy)) {
        file_to_copy_size <- file.info(files_to_copy[i])$size / (1024^2)
        copied_file_size <- file.info(copied_files[i])$size / (1024^2)
        if (!near(file_to_copy_size, copied_file_size, 1))
          message(paste0('Error in copying file: ', files_to_copy[i], 
                       ' to file: ', copied_files[i]))
      }
    }
  }

  Sys.sleep(time = 1)
}



