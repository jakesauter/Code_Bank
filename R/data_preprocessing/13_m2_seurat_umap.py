import os
import csv
import umap
# import fcsparser
import numpy as np
import pandas as pd
import pyarrow.feather as feather

# Setting number of threads available for Numna to use
# os.environ['NUMBA_NUM_THREADS'] = '20' 

print('Reading in Seurat Data ...')


m2_data = pd.read_feather('/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_seurat_expression_data.feather')

# flow_dir_save_loc = "/home/sauterj1/shared_data_folder/sauterj1/Anon_Flow_Folders"
# 
# print('Reading data from flow directories ...')
# 
# flow_dirs = os.listdir(flow_dir_save_loc)
# flow_dirs = [os.path.join(flow_dir_save_loc, flow_dir) for flow_dir in flow_dirs if
#              os.path.isdir(os.path.join(flow_dir_save_loc, flow_dir))]
# 
# m1_files = [os.path.join(x, 'processed_M1_subsampled_100k.fcs') for x in flow_dirs]
# 
# m1_data_list = [fcsparser.parse(x)[1] for x in m1_files]
# 
# print('Concatenating data ...')
# 
# m1_data_mat = np.concatenate(m1_data_list)

print('Creating embedding ...')
reducer = umap.UMAP(n_jobs = 40)
embedding = reducer.fit_transform(m2_data)

print('Saving Embedding ...')

# Save embedding and list of directories in order 
feather.write_feather(pd.DataFrame(embedding), 'M2_umap_results_100k_per_patient.feather')
