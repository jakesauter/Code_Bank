
# This cohort will be defined by the previous cohort plus
# additional patients to reach 200-250 patients with a decently
# fair distribution between diseases

previous_flow_ids
  
patient_df <- 
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv')

previous_patient_df <- 
  patient_df %>% 
  filter(Accession.Number %in% previous_flow_ids)

previous_patient_df %>% 
  count(CbioPortal.Cancer.Type)

sample_df <- 
  patient_df %>% 
  filter(!(Accession.Number %in% previous_flow_ids)) %>% 
  filter(CbioPortal.Cancer.Type != "Leukemia") %>% 
  filter(CbioPortal.Cancer.Type != 
           "Myeloid Neoplasms with Germ Line Predisposition", 
         Exclusion.Reason == "Not Excluded") 

new_patient_df <- 
  previous_patient_df %>% 
  rbind({sample_df %>% 
            filter(CbioPortal.Cancer.Type == "Myelodysplastic Syndromes") %>% 
            sample_n(49)
        }) %>% 
  rbind({sample_df %>% 
            filter(CbioPortal.Cancer.Type == "Myelodysplastic/Myeloproliferative Neoplasms") %>% 
            sample_n(39)
        }) %>%
  rbind({sample_df %>% 
            filter(CbioPortal.Cancer.Type == "Myeloproliferative Neoplasms") %>% 
            sample_n(55)
        }) 

new_patient_df %>% 
  count(CbioPortal.Cancer.Type)
  

cypher_df <- 
  read.csv('data/flow_df_ids_052621.csv')

new_patient_df <- 
  new_patient_df %>% 
  rowwise() %>% 
  mutate(anon_flow_folder_id = 
           sapply(Accession.Number, 
                  function(x) {
                    cypher_df[cypher_df$Accession.Number == x, ]$PRPT_PATH_RPT_ID
                  }) %>% unname()
         ) %>% 
  ungroup() 

write.csv(new_patient_df, 
          'data/239_cohort.csv')

# TODO: maybe want to incorporate flow cytometer
# the sample was measured on and the date that 
# the sample was measured. Take your time in 
# selecting a balanced cohort as the job 
# will take a while to run and this could yield
# robust interpretable results if done right


m1_meta_info_df <- 
  readRDS(file.path('/Users/sauterj1/', 
                    'Documents/Patient_Folder_Analysis/', 
                    'data/m1_meta_info_df.Rds'))



patient_df <- 
  read.csv(file.path('/Users/sauterj1/Documents/', 
                     'Woodlist/Input_Flow_DFs/',
                     'jake_processed_flow_df_ids_052621.csv')) %>% 
  filter(CbioPortal.Cancer.Type != 
         "Myeloid Neoplasms with Germ Line Predisposition", 
         Exclusion.Reason == "Not Excluded")

  
patient_and_meta_df <-
  dplyr::left_join(patient_df, m1_meta_info_df, 
                   by = c('Accession.Number' = 'flow_id')) %>% 
  mutate(cyt_num = unlist(cyt_num)) %>% 
  filter(Accession.Number %in% basename(list.dirs('~/Documents/Woodlist/Full_Cohort_Anon_Flow_Folders/')))


patient_and_meta_df %>% 
  group_by(CbioPortal.Cancer.Type) %>% 
  sample_n(46) %>% 
  ungroup() %>% 
  count(CbioPortal.Cancer.Type, cyt_num) %>% 
  print(n = Inf)



