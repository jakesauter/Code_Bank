library(flowSOM)
library(stringr)
library(flowCore)
library(magrittr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)


### Input data
flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Anon_Flow_Folders/"

flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE)


m1_fcs_list <- vector('list', length(flow_directories))
m2_fcs_list <- vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  cat('Processing directory:', i, 'of:', length(flow_directories), '\n')
  cat('Flow Directory: ', flow_directory, '\n\n')

  flow_directory <- flow_directories[i]
  
  m1_fcs_file <- file.path(flow_directory, 'processed_M1.fcs')
  m2_fcs_file <- file.path(flow_directory, 'processed_M2.fcs')
  
  m1_subsampled_fcs_file <-
    m1_fcs_file %>% 
    strsplit('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_subsampled.fcs')
  
  m2_subsampled_fcs_file <-
    m2_fcs_file %>% 
    strsplit('\\.') %>% 
    .[[1]] %>% .[1] %>% 
    str_c('_subsampled.fcs')
  
  if (!file.exists(m1_subsampled_fcs_file)) {
    m1_fcs_dat <- 
      read.FCS(m1_fcs_file, 
               transformation = FALSE)
    
    
    exprs(m1_fcs_dat) <- 
      exprs(m1_fcs_dat)[sample(1:nrow(m1_fcs_dat), 
                               40e3), ]
    
    write.FCS(m1_fcs_dat, 
              m1_subsampled_fcs_file)
    
  } else {
    m1_fcs_dat <- 
      read.FCS(m1_subsampled_fcs_file, 
               transformation = FALSE)
  }
  
 if (!file.exists(m2_subsampled_fcs_file)) {

    m2_fcs_dat <- 
      read.FCS(m2_fcs_file, 
               transformation = FALSE)
    
    exprs(m2_fcs_dat) <- 
      exprs(m2_fcs_dat)[sample(1:nrow(m2_fcs_dat), 
                               40e3), ]
    
    write.FCS(m2_fcs_dat, 
              m2_subsampled_fcs_file)
 } else {
   m2_fcs_dat <- read.FCS(m2_subsampled_fcs_file, 
                          transformation = FALSE)
 }
  
  m1_fcs_list[[i]] <- m1_fcs_dat
  m2_fcs_list[[i]] <- m2_fcs_dat
}

m1_flowset <- as(m1_fcs_list, 'flowSet')
m2_flowset <- as(m2_fcs_list, 'flowSet')


m1_som <- FlowSOM::FlowSOM(m1_flowset, 
                        xdim = 15, 
                        ydim = 15, 
                        nClus = 40, 
                        colsToUse = seq_len(ncol(m1_fcs_list[[1]])))


m2_som <- FlowSOM::FlowSOM(m2_flowset, 
                        xdim = 15, 
                        ydim = 15, 
                        nClus = 40, 
                        colsToUse = seq_len(ncol(m2_fcs_list[[1]])))



## Next thing that might be interesting: 

# What is the distribution of clusters
# per patient

# probably will have to import patient df
# to correlate clinincal information with 
# cluster info

# can also try to quickly pump these features
# through a random classifier, though this
# will probably have to wait till next week

# metacluster assignments for each of
# the 225 nodes
som$metaclustering

## weight values for each of the
# 225 nodes
som$FlowSOM$map$codes


# best matching unit for every
# individual cell used to build
# the map
FlowSOM::GetClusters(som)

# How to classify new points? 
new_data_som <- 
  FlowSOM::NewData(som, m1_fcs_list[[1]])



# Goal: for each patient determine the distribution
# of events over the SOM for all events in the file

m1_event_classifications <- 
  vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  cat('Processing directory:', i, 'of:', length(flow_directories), '\n')
  cat('Flow Directory: ', flow_directory, '\n\n')
  
  flow_directory <- flow_directories[i]
  
  m1_fcs_file <- file.path(flow_directory, 'processed_M1.fcs')
  m2_fcs_file <- file.path(flow_directory, 'processed_M2.fcs')
  
  m1_fcs_dat <- 
    read.FCS(m1_fcs_file, 
             transformation = FALSE)
  
  m1_event_classifications[[i]] <- 
    FlowSOM::NewData(m1_som, m1_fcs_dat)
  
  # m2_fcs_dat <- 
  #   read.FCS(m2_fcs_file, 
  #            transformation = FALSE)
  
}

saveRDS(m1_event_classifications, 'data/m1_event_classifications.Rds')
# Question: When FlowSOM is "scaling the data" in the NewData()
# function, is it using the same scaling min/max values defined in 
# training (should be yes though make sure)



m1_event_classifications <- readRDS('data/m1_event_classifications.Rds')


# get the cluster and metacluster classifications from 
# each file

m1_cluster_classifications <- 
  vector('list', length(flow_directories))

m1_meta_cluster_classifications <- 
  vector('list', length(flow_directories))

m1_metaclustering <- 
  m1_event_classifications[[1]]$metaclustering

for (i in seq_along(flow_directories)) {
  
  m1_cluster_classifications[[i]] <- 
    m1_event_classifications[[i]] %>% 
    FlowSOM::GetClusters()
  
  
  m1_meta_cluster_classifications[[i]] <-
    m1_cluster_classifications[[i]] %>% 
    m1_metaclustering[.]
}


meta_clustering_mat <- matrix(nrow = length(flow_directories), 
                             ncol = length(unique(m1_metaclustering)))


for (i in seq_along(flow_directories)) {
  meta_clustering_mat[i, ] <- 
    table(m1_meta_cluster_classifications[[i]]) / 
    length(m1_meta_cluster_classifications[[i]])
  
}

colnames(meta_clustering_mat) <- 
  seq(1, ncol(meta_clustering_mat))


meta_clustering_df <- 
  as.data.frame(meta_clustering_mat) %>% 
  mutate(patient_id = dplyr::row_number()) %>% 
  tidyr::pivot_longer(-c(patient_id)) %>% 
  mutate(name = factor(as.integer(name)))

for (i in seq(1, length(unique(m1_metaclustering)))) {
  p <- 
    meta_clustering_df %>% 
    filter(name == paste0('V', i)) %>% 
    ggplot() + 
    geom_histogram(aes(x = value)) + 
    xlab('') + xlim(c(0, 1)) + 
    ggtitle(paste0('Cluster Inclusion: ', i))
  
  print(p)
  
  readline(prompt = 'Press [Enter] to continue')
  
}

plot_list <- vector('list', length(flow_directories))

pdf('images/patient_meta_cluster_distributions')

for (i in seq(1, length(flow_directories))) {
  plot_list[[i]] <- 
    meta_clustering_df %>% 
    filter(patient_id == i) %>% 
    ggplot(aes(x="", y=value, fill=name)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0, direction = -1) + 
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank()) +
    xlab('') + ylab('') + 
    ggtitle(paste0('Patient: ', i))
  
  print(plot_list[[i]])
}

dev.off()


meta_clustering_df %>% 
  filter(patient_id == i) %>% 
  ggplot(aes(x="", y=value, fill=name)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0, direction = -1) + 
  xlab('') + ylab('') + 
  ggtitle(paste0('Patient: ', i))
