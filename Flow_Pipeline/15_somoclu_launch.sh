#!/bin/bash



mpirun -np 60 somoclu -x 10 -y 10 -e 1000 -v 2 \
  /gpfs/mskmind_ess/sauterj1/pipeline_objects/M1_random_sample_seurat_expression_data.table \
  M1_random_sample_somoclust

mpirun -np 60 somoclu -x 10 -y 10 -e 1000 -v 2 \
  /gpfs/mskmind_ess/sauterj1/pipeline_objects/M1_SPADE_seurat_expression_data.table \
  M1_spade_sample_somoclust

mpirun -np 60 somoclu -x 10 -y 10 -e 1000 -v 2 \
  /gpfs/mskmind_ess/sauterj1/pipeline_objects/M2_random_sample_seurat_expression_data.table \
  M2_random_sample_somoclust

mpirun -np 60 somoclu -x 10 -y 10 -e 1000 -v 2 \
  /gpfs/mskmind_ess/sauterj1/pipeline_objects/M2_SPADE_seurat_expression_data.table \
  M2_spade_sample_somoclust


# library('Rsomoclu')
# library('kohonen')
# library(magrittr)

# seurat_object <- 
#   readRDS('data/M2_centralized_Seurat_object.Rds')

# input_data <-
#   seurat_object@assays$MultiParamFlowCyto@data %>%
#   as.matrix() %>% t()

# write.table(input_data,
#           '~/Documents/Patient_Folder_Analysis/data/somoclu_input_data.txt', 
#           row.names = FALSE)

# input_data <- 
  # arrow::read_feather("/gpfs/.../.feather")

# 
# 
# input_data <- data.matrix(input_data)
# 
# input_data <- t(input_data)
# 
# nSomX <- 10
# nSomY <- 10
# nEpoch <- 1
# radius0 <- 0
# radiusN <- 0
# radiusCooling <- "linear"
# scale0 <- 0
# scaleN <- 0.01
# scaleCooling <- "linear"
# kernelType <- 0
# mapType <- "planar"
# gridType <- "rectangular"
# compactSupport <- FALSE
# codebook <- NULL
# neighborhood <- "gaussian"
# stdCoeff <- 0.5
# 
# res <- Rsomoclu.train(input_data, nEpoch, nSomX, nSomY,
#                       radius0, radiusN,
#                       radiusCooling, scale0, scaleN,
#                       scaleCooling,
#                       kernelType, mapType, gridType, compactSupport,
#                       neighborhood, stdCoeff, codebook)
# 
# Rsomoclu:::Rsomoclu.kohonen
# 
# ## Convert to kohonen object for plotting
# sommap = Rsomoclu.kohonen(input_data, res)
# ## Show 'codebook'
# plot(sommap, type="codes", main = "Codes")
# ## Show 'component planes'
# plot(sommap, type = "property", property = sommap$codes[[1]][,1],
#      main = colnames(sommap$codes)[1])
# ## Show 'U-Matrix'
# plot(sommap, type="dist.neighbours")
