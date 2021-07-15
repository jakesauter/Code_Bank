library(FlowSOM)
library(stringr)
library(magrittr)

m1_som <- readRDS('data/m1_som_100_cells_per_patient_30_clus.Rds')

m1_som$FlowSOM$map$xdim
m1_som$FlowSOM$map$ydim

end_slice_idxs <- 
  m1_som$FlowSOM$prettyColnames %>% 
  str_locate("<") %>% .[, 1] 
  
m1_som$FlowSOM$prettyColnames <-
  sapply(seq_along(end_slice_idxs), 
         function(idx) m1_som$FlowSOM$prettyColnames[idx] %>% 
           str_sub(1, end_slice_idxs[idx]-2))

fsom <- m1_som

fsom$FlowSOM$MST$size <- m1_som$FlowSOM$MST$size / 2
fsom$FlowSOM$MST$size <-  2

PlotStars(fsom$FlowSOM, 
          backgroundValues = fsom$metaclustering)

PlotStars(m1_som$FlowSOM, 
          backgroundValues = m1_som$metaclustering, 
          view = 'grid')

fsom$FlowSOM$MST$size <- 3

PlotStars(fsom$FlowSOM, 
          backgroundValues = m1_som$metaclustering, 
          view = 'tSNE')

