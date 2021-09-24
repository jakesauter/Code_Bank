import time
import scipy
import scanpy
import pandas as pd
import pyarrow

seurat_expression_data = \
pd.read_feather('/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_random_sample_seurat_expression_data.feather')

seurat_ad = \
scanpy.AnnData(X = seurat_expression_data)

start = time.time()

scanpy.pp.neighbors(
  seurat_ad, 
  use_rep=None, 
  n_pcs=0, 
  method='umap', 
  n_neighbors = 30)

end = time.time()
print(end-start)

datetime = time.strftime("%Y%m%d-%H%M%S")

start2 = time.time()

scanpy.tl.leiden(seurat_ad, 
                 directed = False, 
                 use_weights = True,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
                 n_iterations = 3)

end2 = time.time()
print(end2-start2)


pyarrow.feather.write_feather(seurat_ad.obs['leiden'].to_frame(), 
                              'M1_RANDOM_leiden_clustering_results_' + datetime + '.feather')
