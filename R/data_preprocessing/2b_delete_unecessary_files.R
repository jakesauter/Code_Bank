dirs <- 
  list.files('~/Documents/Woodlist/Full_Correct_QC_Flow_Folders/', 
             full.names = TRUE) %>% 
  .[fs::is_dir(.)]


choose_file_to_use <- function(filenames) {
  
  if (length(filenames) == 1)
    return(filenames)
  
  # remove any pre-processed files for now
  contains_high_qual_events <- 
    str_detect(basename(filenames), 
               'high_quality_events.fcs')
  
  filenames <- 
    filenames[contains_high_qual_events]
  
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




for (dir in dirs) {
  
  files <- 
    list.files(dir, 
               full.names = TRUE, 
               all.files = TRUE)
  
  m1_fcs_file <-
    (
      str_detect(files, '\\.fcs') &
      str_detect(files, 'm1|M1') 
    ) %>% 
    files[.] %>%
    .[str_detect(., 'singlets')]
    # choose_file_to_use() %>% 
    fs::path_real()
  
  m2_fcs_file <-
    (
      str_detect(files, '\\.fcs') &
        str_detect(files, 'm2|M2') 
    ) %>% 
    files[.] %>%
    choose_file_to_use() %>% 
    fs::path_real()
  
  
  to_delete <-
    files[fs::path_real(files) != m1_fcs_file &
          fs::path_real(files) != m2_fcs_file &
          !str_detect(files, '\\.nlys$') &
          !str_detect(files, '^[^.]*$') &
          !fs::is_dir(files)]
  
  sapply(to_delete, 
         function(file) {
            res <- try(file.remove(file), 
                        silent = FALSE)
            
            if ('try-error' %in% class(res)) {
              message(paste0('File: ', file, 
                             ' could not be deleted'))
            }
  })
  
}




# Found out how many and which directories don't
# have the correct compensation matrix for their
# version of the m1 and m2 panel

bad_dirs <- c()

for (flow_directory in flow_directories) {

  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs", 
               full.names = TRUE) %>% 
    .[str_detect(., 'high_quality_events.fcs')]
  
  m1_fcs_filename <-
    str_detect(fcs_files,
               'M1|m1') %>%
    fcs_files[.] %>%
    choose_file_to_use() %>%
    fs::path_real()
  
  
  m2_fcs_filename <-
    str_detect(fcs_files,
               'M2|m2') %>%
    fcs_files[.] %>%
    choose_file_to_use() %>%
    fs::path_real()
  
  
  m1_comp_mat_filename <-
    file.path(flow_directory, 
              'compensation_matrices') %>% 
    list.files(pattern = 'M1.*\\.csv')
  
  m2_comp_mat_filename <-
    file.path(flow_directory, 
              'compensation_matrices') %>% 
    list.files(pattern = 'M2.*\\.csv')
  

  file.path(flow_directory, 
            'compensation_matrices') %>% 
    list.files(pattern = 'M1.*\\.csv') %>% 
    print()
  
  file.path(flow_directory, 
            'compensation_matrices') %>% 
    list.files(pattern = 'M2.*\\.csv') %>% 
    print()
  
  print(m1_comp_mat_filename)
  print(m2_comp_mat_filename)
  
  if (
      (
        length(m1_fcs_filename) &
        length(m1_comp_mat_filename) &
        length(m2_fcs_filename) &
        length(m2_comp_mat_filename) 
      ) &&
      ({
        !str_detect(
        str_split(basename(m1_fcs_filename), '\\.')[[1]][1], 
        str_split(basename(m1_comp_mat_filename), '\\.')[[1]][1])
       } | 
      {
        !str_detect(
          str_split(basename(m2_fcs_filename), '\\.')[[1]][1], 
          str_split(basename(m2_comp_mat_filename), '\\.')[[1]][1])
      })) {
    
    bad_dirs <- 
      c(bad_dirs,
        dirname(m1_fcs_filename)
      )
    
  }
  
}















