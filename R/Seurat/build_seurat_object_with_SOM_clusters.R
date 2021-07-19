library(dplyr)
library(Seurat)
library(FlowSOM)
library(flowCore)
library(magrittr)
library(patchwork)

filter <- dplyr::filter

cat('Reading in FlowSOM object ... \n')

m1_som <-
  readRDS('data/10_by_10_SOM_empirical_num_clusters.Rds')

cat('Reading in Seurat object ... \n')


seurat_object <-
  readRDS('data/umap_seurat_object.Rds')

flow_ids <- seurat_object[[]]$flow_dir


cat('Adding metadata to Seurat object ... \n')


# Read in patient metadata, and correlate all
# interested variables with metadata of each
# cell based on the flow_id / flow_dir metadata
patient_df <-
  read.csv('~/Documents/Woodlist/Input_Flow_DFs/jake_processed_flow_df_ids_052621.csv',
           stringsAsFactors = FALSE) %>%
  filter(Accession.Number %in% flow_ids)

m1_meta_info_df <-
  readRDS('/Users/sauterj1/Documents/Patient_Folder_Analysis/data/m1_meta_info_df.Rds') %>%
  filter(flow_id %in% flow_ids)

patient_and_meta_info_df <-
  dplyr::left_join(patient_df, m1_meta_info_df,
                   by = c('Accession.Number' = 'flow_id'))

rownames(patient_and_meta_info_df) <-
  patient_and_meta_info_df$Accession.Number


# Add Date, Machine, Course Disease and Fine Disease to
# Seurat metadata
seurat_object@meta.data$cancer_type <-
  seurat_object@meta.data$flow_dir %>%
  patient_and_meta_info_df[., "CbioPortal.Cancer.Type"]

seurat_object@meta.data$cancer_type_detailed <-
  seurat_object@meta.data$flow_dir %>%
  patient_and_meta_info_df[., "CbioPortal.Cancer.Type.Detailed"] %>%
  unlist()

seurat_object@meta.data$flow_cytometer <-
  seurat_object@meta.data$flow_dir %>%
  patient_and_meta_info_df[., "cyt_num"] %>%
  unlist()

seurat_object@meta.data$experiment_date <-
  seurat_object@meta.data$flow_dir %>%
  patient_and_meta_info_df[., "date"] %>%
  unlist()



cat('Computing metaclusters of Seurat object ... \n')


# Add FlowSOM metacluster as a meta.data
# variable
data <-
  seurat_object@assays$MultiParamFlowCyto@data %>%
  as.matrix() %>%
  t()

ex_data <-
  read.FCS('~/Documents/Woodlist/Anon_Flow_Folders_Server_Version/F18-10168/processed_M1_subsampled_100k.fcs')

data_ff <-
  flowCore::flowFrame(data,
                      parameters = parameters(ex_data))

new_data_som <-
  FlowSOM::NewData(m1_som,
                   data_ff)

data_clusters <-
  FlowSOM::GetClusters(new_data_som)

data_metaclusters <-
  data_clusters %>%
  m1_som$metaclustering[.]

cat('Adding metacluster info to Seurat object ... \n')

seurat_object@meta.data$metacluster <-
  data_metaclusters

cat('Saving Seurat object ... \n')


saveRDS(seurat_object,
        'data/umap_seurat_object_with_10_by_10_SOM_metaclusters.Rds')

