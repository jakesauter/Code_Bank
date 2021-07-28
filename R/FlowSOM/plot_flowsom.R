library(FlowSOM)
library(stringr)
library(magrittr)

m1_som <- readRDS('data/7_by_7_SOM_empirical_num_clusters.Rds')

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


p <- PlotFlowSOM(fsom, view ='MST', 
            equalNodeSize = TRUE, 
            title = "FlowSOM Metaclusters") %>%
    AddNodes(values = fsom$metaclustering, 
             showLegend = TRUE, label = "Metaclusters", 
             colorPalette = c("#e07200",
                                   "#4a80ff",
                                   "#ae9500",
                                   "#002e9f",
                                   "#00bf67",
                                   "#820093",
                                   "#4c7400",
                                   "#ff6fe7",
                                   "#b20500",
                                   "#004584",
                                   "#fa8d69",
                                   "#e293bd",
                                   "#6e2107")) %>%
    AddLabels(labels = fsom$metaclustering)


p
