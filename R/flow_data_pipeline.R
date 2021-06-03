parse_comp_mat_python_path <- "/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/parse_flow_folder_compensation.py"

flow_dir_save_loc <- "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

WL3_patients <- readRDS('~/Documents/Patient_Folder_Analysis/WL3_patients.Rds')

WL3_patients <- 
  WL3_patients %>% 
  filter(WL3_nlys_file_count == 1)

flow_directories <- WL3_patients$Flow.Data.Folder

# Testing on 10 for now
flow_directories <- flow_directories[1:10]

compensated_fcs_data <- vector('list', length(flow_directories))
names(compensated_fcs_data) <- basename(flow_directories)

# Now load fcs, load compensation, and apply
for (flow_directory in flow_directories) {
  
  cat("Working on Flow Directory: ", flow_directory, "\n\n")
  
  # New path we will save to
  local_flow_directory <- 
    file.path(flow_dir_save_loc, 
              basename(flow_directory))
  
  # Need to copy flow folders to local 
  #in order to modify and save compensation matrices
  if (!dir.exists(local_flow_directory)) {
    file.copy(flow_directory, 
              flow_dir_save_loc, 
              recursive = TRUE, 
              overwrite = FALSE)
  }
  
  # Parse the compensation matrices and save
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
  
  
  m1_file <- 
    WL3_patients %>% 
    filter(basename(Flow.Data.Folder) == 
             basename(local_flow_directory)) %>% 
    .$M1_path %>% basename() %>% 
    file.path(local_flow_directory, .)
  
  m2_file <- 
    WL3_patients %>% 
    filter(basename(Flow.Data.Folder) == 
             basename(local_flow_directory)) %>% 
    .$M2_path %>% basename() %>% 
    file.path(local_flow_directory, .)
  
  m1_fcs <- flowCore::read.FCS(m1_file)
  m2_fcs <- flowCore::read.FCS(m2_file)
  
  m1_comp_mat <- 
    m1_file %>% 
    basename() %>% 
    str_split('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('.csv') %>% 
    file.path(local_flow_directory, 
              'compensation_matrices',
              .) %>% 
    read.csv(header = FALSE) %>% 
    set_colnames(pData(parameters(m1_fcs))$name)
  
  m2_comp_mat <- 
    m2_file %>% 
    basename() %>% 
    str_split('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('.csv') %>% 
    file.path(local_flow_directory, 
              'compensation_matrices',
              .) %>% 
    read.csv(header = FALSE) %>% 
    set_colnames(pData(parameters(m2_fcs))$name)
  
  # compensate 
  m1_fcs <- compensate(m1_fcs, m1_comp_mat)
  m2_fcs <- compensate(m2_fcs, m1_comp_mat)
  
  # save
  compensated_fcs_data[[basename(flow_directory)]] <- 
    list('m1_fcs' = m1_fcs, 'm2_fcs' = m2_fcs)
  
  # later probably going to want to do the analysis here as not going to want to save 1000 flowFrames. The ncdfFlow package supposidely can perform all the same operations as flowCore without storing all flowSets in Ram, though we can take a look at this when needed
  
}

library(flowAI)

for (comp_fcs_data in compensated_fcs_data) {
  
  flowAI_QC_output_folder <- "/Users/sauterj1/Documents/Patient_Folder_Analysis/flowAI_QC_Results"
  
  m1_fcs <- comp_fcs_data$m1_fcs
  m1_metadata <- pData(parameters(m1_fcs))
  m2_fcs <- comp_fcs_data$m2_fcs
  m2_metadata <- pData(parameters(m2_fcs))
  
  # Only use channel names with markers to check for 
  # signal quality
  exclude_channels <- 
    m1_metadata$name[is.na(m1_metadata$desc)]
  
  exclude_channels <- 
    exclude_channels[exclude_channels != "Time"]
  
  m1_fcs_flowai_qc <- 
    flowAI::flow_auto_qc(m1_fcs, 
                         ChExcludeFS = exclude_channels,
                         ChExcludeFM = exclude_channels,
                         folder_results = flowAI_QC_output_folder, 
                         mini_report = FALSE, 
                         html_report = "_QC",
                         fcs_QC = FALSE)
  
  
  # Only use channel names with markers to check for 
  # signal quality
  exclude_channels <- 
    m2_metadata$name[is.na(m2_metadata$desc)]
  
  exclude_channels <- 
    exclude_channels[exclude_channels != "Time"]
  
  m2_fcs_flowai_qc <- 
    flowAI::flow_auto_qc(m2_fcs, 
                         ChExcludeFS = exclude_channels,
                         ChExcludeFM = exclude_channels,
                         folder_results = flowAI_QC_output_folder, 
                         mini_report = FALSE, 
                         html_report = "_QC",
                         fcs_QC = FALSE)
}