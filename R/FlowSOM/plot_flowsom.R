library(FlowSOM)
library(stringr)
library(magrittr)

m1_som <- readRDS('data/10_by_10_SOM_empirical_num_clusters.Rds')

FlowSOM::NMetaclusters(m1_som)
fsom <- m1_som

end_slice_idxs <- 
  m1_som$prettyColnames %>% 
  str_locate("<") %>% .[, 1] 
  
names_colnames <- 
  names(m1_som$prettyColnames)

fsom$prettyColnames <-
  sapply(seq_along(end_slice_idxs), 
         function(idx) m1_som$prettyColnames[idx] %>% 
           str_sub(1, end_slice_idxs[idx]-2))

names(fsom$prettyColnames) <- 
  names_colnames

PlotStars(fsom,
          backgroundValues = m1_som$metaclustering)



FlowSOM::FlowSOMmary(fsom, 
                     plotFile = '10_by_10_SOMmary.pdf')


PlotFlowSOM(fsom, view = 'MST', equalNodeSize = TRUE, 
            title = "FlowSOM Metaclusters") %>% 
    AddNodes(values = fsom$metaclustering, 
             showLegend = TRUE, label = "Metaclusters") %>% 
    AddLabels(labels = fsom$metaclustering)

