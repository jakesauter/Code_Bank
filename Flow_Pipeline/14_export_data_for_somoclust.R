



#  OMP_NUM_THREADS=50 somoclu -x 10 -y 10 -e 1000 -v 2 /gpfs/mskmind_ess/sauterj1/seurat/objects/M1_SPADE_seurat_expression_data.table m1_spade_somoclust
# OMP_NUM_THREADS=50 somoclu -x 10 -y 10 -e 1000 -v 2 /gpfs/mskmind_ess/sauterj1/seurat/objects/M1_random_sample_seurat_expression_data.table m1_random_somoclust


data <- 
  arrow::read_feather('/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_SPADE_seurat_expression_data.feather')

write.table(data, 
            '/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_SPADE_seurat_expression_data.table', 
            row.names = FALSE, 
            col.names = FALSE)

data <- 
  arrow::read_feather('/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_SPADE_seurat_expression_data.feather')

write.table(data, 
            '/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_SPADE_seurat_expression_data.table', 
            row.names = FALSE, 
            col.names = FALSE)

data <- 
  arrow::read_feather('/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_random_sample_seurat_expression_data.feather')

write.table(data, 
            '/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_random_sample_seurat_expression_data.table', 
            row.names = FALSE, 
            col.names = FALSE)

data <- 
  arrow::read_feather('/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_random_sample_seurat_expression_data.feather')

write.table(data, 
            '/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_random_sample_seurat_expression_data.table', 
            row.names = FALSE, 
            col.names = FALSE)