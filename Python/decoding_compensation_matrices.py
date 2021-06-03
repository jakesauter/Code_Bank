#!/bin/python3.8

# In order to have access to ccl_bplist module
import sys
sys.path.append('/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/')

import re
import os
import json
import struct
import pprint
import ccl_bplist
import numpy as np

def decode_comp_matrix(comp_matrix_string, num_params): 
  comp_matrix = struct.unpack('>' + 'd' * (len(list(comp_matrix_string)) // 8), comp_matrix_string)
  comp_matrix = np.reshape(comp_matrix, (num_params, num_params))
  return comp_matrix
  
# In order to have access to ccl_bplist module
import sys
sys.path.append('/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/')

import re
import os
import json
import struct
import pprint
import ccl_bplist
import numpy as np

def decode_comp_matrix(comp_matrix_string, num_params): 
  comp_matrix = struct.unpack('>' + 'd' * (len(list(comp_matrix_string)) // 8), comp_matrix_string)
  comp_matrix = np.reshape(comp_matrix, (num_params, num_params), order = 'F')
  return comp_matrix
  
def parse_comp_mats_for_flow_dir(flow_directory): 
  owd = os.getcwd()
  os.chdir(flow_directory)
  flow_directory_files = os.listdir(flow_directory)
  nlys_files = [file for file in flow_directory_files if file.endswith(".nlys")]

  if len(nlys_files) > 1: 
    print("MORE THAN ONE NLYS FILE FOUND, PROCEDING WITH: " + nlys_files[0])
    
  nlys_filename = nlys_files[0]
  xml_filename = nlys_filename.split('.')[0] + '.xml'
  
  os.system('plutil -convert xml1 -o ' + xml_filename + ' ' + nlys_filename)
  
  xml_file = open(xml_filename, 'r')
  xml_lines = xml_file.readlines()

  cmmatrix_in_lines = ['CMMatrix' in line for line in xml_lines]
  cmmatrix_line_idxs = [i for i, x in enumerate(cmmatrix_in_lines) if x]
  
  filename_in_lines = [re.search('<string>.*\.fcs</string>', line) is not None for line in xml_lines]
  filename_line_idxs = [i for i, x in enumerate(filename_in_lines) if x]
  filename_lines = [xml_lines[i] for i in filename_line_idxs]
  filenames = [re.search('>.*<', filename_line).group(0)[1:-1] for filename_line in filename_lines]

  # Associate CMMatrix Start Line with Filename

  filename_cmmatrix_loc_dict = {}

  for cmmatrix_line_idx in cmmatrix_line_idxs: 
    assoc_filename_line = max(filter(lambda filename_line_idx: filename_line_idx < cmmatrix_line_idx, filename_line_idxs))
    filename = filenames[filename_line_idxs.index(assoc_filename_line)]
    filename_cmmatrix_loc_dict[filename] = cmmatrix_line_idx
  
  filename_comp_mat_dict = {}

  for filename, start_idx in filename_cmmatrix_loc_dict.items(): 

    lines = xml_lines[start_idx:start_idx + 8]

    # Find which line has UID
    # parse uid out of <integer> </integer> tags from next line
    uid_line_detected = ['UID' in line for line in lines]
    uid_line_idx = [i for i, x in enumerate(uid_line_detected) if x][0]
    uid_int_line = lines[uid_line_idx + 1]
    cmmatrix_uid = int(re.search('>.*<', uid_int_line).group(0)[1:-1])

    # Find which line has CMNumberParams 
    # Parse number params out of next line 
    params_line_detected = ['CMNumberParams' in line for line in lines]
    params_line_idx = [i for i, x in enumerate(params_line_detected) if x][0]
    params_int_line = lines[params_line_idx + 1]
    num_params = int(re.search('>.*<', params_int_line).group(0)[1:-1])

    filename_comp_mat_dict[filename] = {'cmmatrix_uid': cmmatrix_uid, 'num_params': num_params}

  f = open(nlys_filename, "rb")
  plist = ccl_bplist.load(f)
  object_table = plist['$objects']

  comp_mats_dict = {}

  for filename, comp_dict in filename_comp_mat_dict.items(): 
    comp_mat_string = object_table[comp_dict['cmmatrix_uid']]
    comp_mats_dict[filename] = decode_comp_matrix(comp_mat_string, comp_dict['num_params'])
    
  # Save comp mats
  keys = [str(x) for x in comp_mats_dict.keys()]
  m1_or_m2_keys = [re.search('m1|m2', key, flags = re.IGNORECASE) for key in keys]
  m1_or_m2_keys = [keys[i] for i in range(len(keys)) if m1_or_m2_keys[i] is not None]

  comp_mat_dir = os.path.join(flow_directory, 'compensation_matrices')
  if not os.path.isdir(comp_mat_dir): os.mkdir(comp_mat_dir)

  for key in m1_or_m2_keys: 
    comp_mat = comp_mats_dict[key]
    filepath = os.path.join(comp_mat_dir, key.split('.')[0] + '.csv')
    np.savetxt(filepath, comp_mat, delimiter=",")

  # change back to orignal working directory
  os.chdir(owd)
            
  return(comp_mats_dict)


if __name__ == "__main__":
  parse_comp_mats_for_flow_dir(sys.argv[1])
  print("Saved M1 and M2 compensation matrices to " + os.path.join(sys.argv[1], 'compensation_matrices') + '\n\n')
