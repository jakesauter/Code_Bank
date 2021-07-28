library(dplyr)
library(Seurat)
library(ggplot2)
library(FlowSOM)
library(flowCore)
library(magrittr)
library(patchwork)

filter <- dplyr::filter


Sys.setenv('R_MAX_VSIZE'=100e9)

cat('\nReading in Seurat object ...\n')

seurat_object <- 
  readRDS('seurat_object_with_louvain_clusters.Rds')


cat('\nConverting Seurat object to matrix ...\n')

data_mat <- 
  seurat_object@assays$MultiParamFlowCyto@data %>% 
  as.matrix() %>% 
  t()

# Found this issue when trying to plot whole object, 
# should remove this in quality checking 
rows_to_keep <- 
  data_mat[, 'FSC-A'] > -5

clean_data_mat <- 
  data_mat[rows_to_keep, ]

cat('\nConverting matrix to dataframe ...\n')

df <- data.frame(clean_data_mat)

df[, 'metacluster'] <- 
  seurat_object$louvain_clusters[rows_to_keep]

cat('\nPivoting dataframe ...\n')

longer_df <- 
  tidyr::pivot_longer(df, cols = -c(metacluster))

longer_df$name <- 
  longer_df$name %>% 
  factor(levels = 
           c("FSC.A", "FSC.H", "SSC.A", "SSC.H", "CD13", "CD15", "CD19", 
             "CD33", "CD34", "CD38", "CD45", "CD71", "CD117", "HLA.DR"))

metaclusters <- 
  unique(df$metacluster) %>% 
  sort()

longer_df$metacluster <- 
  longer_df$metacluster %>% 
  factor(levels = metaclusters, 
         labels = paste0('Community_', metaclusters))

cat('\nPlotting first plot ...\n')

bmp('louvain_marker_dist_across_communities.bmp', 
    height = 600, 
    width = 800)

longer_df %>% 
  ggplot() + 
  geom_violin(aes(x = metacluster, 
                  y = value, 
                  fill = metacluster)) + 
  xlab('') + ylab('') + 
  facet_wrap(vars(name)) + 
  guides(fill=guide_legend(title="Metacluster")) + 
  theme(axis.text.x = element_text(angle = -45,
                                   hjust =  0)) 

dev.off()

cat('\nPlotting second array of plots ...\n')

dir.create('marker_expression_per_louvain_community')
setwd('marker_expression_per_louvain_community')

for (community in unique(df$metacluster)) {

  bmp(paste0('community_', community, '_marker_distributions.bmp'), 
      height = 600, 
      width = 800)
  
  longer_df %>% 
    ggplot() + 
    geom_violin(aes(x = name, 
                    y = value, 
                    fill = name)) + 
    xlab('') + ylab('') + 
    facet_wrap(vars(metacluster)) + 
    guides(fill=guide_legend(title="Marker")) + 
    theme(axis.text.x = element_text(angle = -45,
                                     hjust =  0)) 
  
  dev.off()
}



