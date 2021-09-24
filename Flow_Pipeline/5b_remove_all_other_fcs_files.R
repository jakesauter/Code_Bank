dirs <- 
  list.files('~/Documents/Woodlist/Full_Correct_QC_Flow_Folders/', 
             full.names = TRUE) %>% 
  .[fs::is_dir(.)]




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
    .[str_detect(., 'singlets')] %>% 
    fs::path_real()
  
  m2_fcs_file <-
    (
      str_detect(files, '\\.fcs') &
        str_detect(files, 'm2|M2') 
    ) %>% 
    files[.] %>%
    .[str_detect(., 'singlets')] %>% 
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



