

set.seed(497)

suppressPackageStartupMessages(library(spade))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(flowCore))

options(warn = -1)

filter <- dplyr::filter

args <- commandArgs(trailingOnly=TRUE)


flow_directory <- args[[1]]

owd <- setwd(flow_directory)

cat('Flow Directory: ', flow_directory, '\n\n')

system("echo 'If this file exists, there was an error in subsampling' > error_file.txt")

tic <- Sys.time()

m1_fcs_file <- file.path(flow_directory, 'processed_M1.fcs')
m2_fcs_file <- file.path(flow_directory, 'processed_M2.fcs')

SPADE.removeExistingDensityAndClusterColumns(m1_fcs_file)
SPADE.removeExistingDensityAndClusterColumns(m2_fcs_file)


m1_subsampled_fcs_file <-
  m1_fcs_file %>% 
  strsplit('\\.') %>% 
  .[[1]] %>% .[1] %>% 
  str_c('_subsampled_100k.fcs')

m2_subsampled_fcs_file <-
  m2_fcs_file %>% 
  strsplit('\\.') %>% 
  .[[1]] %>% .[1] %>% 
  str_c('_subsampled_100k.fcs')

if (file.exists(m1_subsampled_fcs_file)) {
  file.remove(m1_subsampled_fcs_file)
}

if (file.exists(m2_subsampled_fcs_file)) {
  file.remove(m2_subsampled_fcs_file)
}

SPADE.addDensityToFCS(m1_fcs_file, 
                      m1_fcs_file, 
                      cols = NULL, 
                      transforms=NULL,
                      kernel_mult = 5, 
                      apprx_mult = 1.5,
                      med_samples = 10e3,
                      comp=FALSE)

SPADE.addDensityToFCS(m2_fcs_file, 
                      m2_fcs_file, 
                      cols = NULL, 
                      transforms=NULL,
                      kernel_mult = 5, 
                      apprx_mult = 1.5,
                      med_samples = 10e3,
                      comp=FALSE)


SPADE.downsampleFCS(m1_fcs_file, 
                    m1_subsampled_fcs_file, 
                    target_percent = NULL, 
                    target_number = 100e3)

SPADE.downsampleFCS(m2_fcs_file, 
                    m2_subsampled_fcs_file, 
                    target_percent = NULL, 
                    target_number = 100e3)

SPADE.removeExistingDensityAndClusterColumns(m1_subsampled_fcs_file)
SPADE.removeExistingDensityAndClusterColumns(m2_subsampled_fcs_file)


system('rm *.orig1')


system('rm error_file.txt')

setwd(owd)
