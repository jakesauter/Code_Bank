
library(dplyr)
library(stringr)
library(magrittr)
library(flowCore)
library(FlowSOM)
library(parallel)

filter <- dplyr::filter

flow_dirs_root <-
  "/home/sauterj1/shared_data_folder/sauterj1/Anon_Flow_Folders"

flow_directories <-
  list.files(flow_dirs_root,
             full.names = TRUE) %>% 
  sort()

m1_fcs_data_list <-
  vector('list', length(flow_directories))

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  
  cat('Reading fcs file:', i, 'of',
      length(flow_directories), '\n')
  
  fcs_file <- file.path(flow_directory,
                        'processed_M1_subsampled_100k.fcs')
  
  m1_fcs_data_list[[i]] <-
    read.FCS(fcs_file,
             transformation = FALSE) %>%
    exprs()
  
}

m1_fcs_data_mat <-
  m1_fcs_data_list %>%
  do.call(rbind, .)

m1_fcs_ff <-
  flowFrame(m1_fcs_data_mat)

# Coefficient of variation is sdev(y) / mean(y), though this 
# does not make sense for data in which 0 is not the absolute 
# minimum value
# FlowSOM::GetMetaclusterCVs


# Thus I will implement calculating the **variance** for every 
# marker in every metacluster myself to use this as a metric
# to plot for the elbow method

GetMetaclusterVariance <- 
  function(fsom) {
    variance <- 
      sapply(seq_along(levels(fsom$metaclustering)), 
          function(i) {
            apply(
              subset(
                    fsom$data, 
                    fsom$metaclustering[GetClusters(fsom)] == i), 2, 
                      function(y) ifelse (length(y) > 0, stats::sd(y), NA)
            )
          }) 
    
    t(variance)
  }

cat('\nBuilding FlowSOM object ...\n')
m1_som <- FlowSOM::FlowSOM(m1_fcs_ff,
                           xdim = 15,
                           ydim = 15,
                           nClus = 50,
                           colsToUse = seq_len(ncol(m1_fcs_data_mat)))

cat('\nSaving FlowSOM object ...\n')
saveRDS(m1_som, 'm1_som_elbow_method.Rds')


nclust_to_try <- seq(50, 10, -1)
variance_scores <- rep(0, length(nclust_to_try))

cat('\nBegining metacluster attempts ...\n')

variance_scores <- 
  parallel::mclapply(
    seq_along(nclust_to_try), 
    function(i) {
    
      nclust <- nclust_to_try[i]
      
      cl <- as.factor(
        metaClustering_consensus(
          data = m1_som$map$codes, 
          k = nclust, 
          seed = 497
        )
      )
      
      m1_som <- 
        FlowSOM::ReassignMetaclusters(m1_som, cl)
      
      sum(GetMetaclusterVariance(m1_som))
  }, mc.cores = 60)


cat('\nSaving metacluster variance scores ...\n')
saveRDS(variance_scores, 'FlowSOM_variance_scores.Rds')






