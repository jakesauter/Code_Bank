#!/bin/python3.8

import sys
sys.path.append('/Users/sauterj1/Documents/Patient_Folder_Analysis/Python/')

import os
import re
import sys
import csv
import pprint
import ccl_bplist

# Wanted to break up logic in calculating the proportion info as 
# this might get messy with accounting for possible spelling mistakes

def compute_pop_proportions(pop_info_dict): 
  
  population_proportion_dict = {}

  wbc_key_pairs = {
    "CD34+ Blasts" : 'CD34+',
    "CD117+ Cells" : 'CD117',
    "B Cells" : 'B cells',
    "Monocytes": 'Monocytes',
    "Lymphocytes": 'Lymphocytes',
    "Basophils": 'Basophils',
    "Eosinophils": 'Eosynophils',
    "PCDC":'PCDC',
    "Plasma Cells": 'Plasma cells'
  }
  
  for key, value in wbc_key_pairs.items(): 
    if value in pop_info_dict['M1']: 
      population_proportion_dict[key] = \
        pop_info_dict['M1'][value] / pop_info_dict['M1']['WBC']
    elif value in pop_info_dict['M2']: 
       population_proportion_dict[key] = \
        pop_info_dict['M2'][value] / pop_info_dict['M2']['WBC']
    else: 
      population_proportion_dict[key] = -1
  
  if "Eryhthroids" in pop_info_dict['M1']: 
    population_proportion_dict["Eryhthroids"] = \
      pop_info_dict['M1']['Eryhthroids'] / pop_info_dict['M1']['Singlets']
  elif "Eryhthroids" in pop_info_dict['M2']: 
    population_proportion_dict["Eryhthroids"] = \
      pop_info_dict['M2']['Eryhthroids'] / pop_info_dict['M2']['Singlets']
  else: 
    population_proportion_dict["Eryhthroids"] = -1
    
  if "WBC" in pop_info_dict['M1']: 
    population_proportion_dict["WBC"] = \
      pop_info_dict['M1']['WBC'] / pop_info_dict['M1']['Singlets']
  elif "WBC" in pop_info_dict['M2']: 
    population_proportion_dict["WBC"] = \
      pop_info_dict['M2']['WBC'] / pop_info_dict['M2']['Singlets']
  else: 
    population_proportion_dict["Eryhthroids"] = -1
    
  for k,v in population_proportion_dict.items(): 
    # population_proportion_dict[k] = round(population_proportion_dict[k]*100, 2)
    population_proportion_dict[k] = population_proportion_dict[k]*100

  return(population_proportion_dict)

def write_pop_info_to_csv(pop_prop_dict, csv_filepath): 
  with open(csv_filepath, 'w') as f: 
      w = csv.DictWriter(f, pop_prop_dict.keys())
      w.writeheader()
      w.writerow(pop_prop_dict)
      
  return True
  

def parse_pop_info_for_flow_dir(flow_directory): 
  
  populations = ['CD34+', 'WBC', 'CD117', 'B cells', 'Lymphocytes', 
               'Eryhthroids', 'Singlets', 'Plasma cells', 'Granulocytes', 
               'WBC', 'Monocytes', 'Basophils', 'Eosynophils', 'PCDC']
  
  owd = os.getcwd()
  os.chdir(flow_directory)
  flow_directory_files = os.listdir(flow_directory)
  nlys_files = [file for file in flow_directory_files if file.endswith(".nlys")]
  
  if len(nlys_files) > 1: 
    print("MORE THAN ONE NLYS FILE FOUND, PROCEDING WITH: " + nlys_files[0])
    
  nlys_filename = nlys_files[0]
  xml_filename = nlys_filename.split('.')[0] + '.xml'
  
  os.system('plutil -convert xml1 -o ' + xml_filename + ' ' + nlys_filename)
  
  plist = ccl_bplist.load(open(nlys_filename, 'rb'))
  object_table = plist['$objects']

  xml_file = open(xml_filename, 'r')
  xml_lines = xml_file.readlines() 

  # Now filter these lines for lines that the next line contains "<integer>"
  number_lines = [":Number" in line for line in xml_lines]
  number_line_idxs =  [i for i, x in enumerate(number_lines) if x]

  # Considering whether the line after OR second line after are integer lines
  pop_number_lines = [x for x in number_line_idxs if '<integer>' in xml_lines[x+1] or '<integer>' in xml_lines[x+2]]
  pop_number_strings = [re.search('>.*:', xml_lines[i]).group(0)[1:-1] for i in pop_number_lines]

  # Have to choose one line or the other
  pop_nums = []
  for i in pop_number_lines: 
    pop_num_1 = re.search('integer>.*<', xml_lines[i+1])
    pop_num_2 = re.search('integer>.*<', xml_lines[i+2])

    if pop_num_1 is None: 
      pop_num = pop_num_2.group(0)[8:-1]
    else: 
      pop_num = pop_num_1.group(0)[8:-1]

    pop_nums.append(pop_num)

  filename_pop_number_loc_dict = {}
  
  filename_in_lines = [re.search('<string>.*\.fcs</string>', line) is not None for line in xml_lines]
  filename_line_idxs = [i for i, x in enumerate(filename_in_lines) if x]
  filename_lines = [xml_lines[i] for i in filename_line_idxs]
  filenames = [re.search('>.*<', filename_line).group(0)[1:-1] for filename_line in filename_lines]

  for pop_number_line_idx in pop_number_lines: 
    assoc_filename_line = max(filter(lambda filename_line_idx: filename_line_idx < pop_number_line_idx, filename_line_idxs))
    filename = filenames[filename_line_idxs.index(assoc_filename_line)]
    filename_pop_number_loc_dict[pop_number_line_idx] = filename

  pop_number_dict = {}

  for i in range(len(pop_nums)):
    line_idx = pop_number_lines[i]
    filename = filename_pop_number_loc_dict[line_idx]
    if filename not in pop_number_dict: 
      pop_number_dict[filename] = {pop_number_strings[i]: int(pop_nums[i])}
    else: 
      pop_number_dict[filename][pop_number_strings[i]] = int(pop_nums[i])
      
  # Resolve unfound populations
  pop_num_dict_keys = list(pop_number_dict.keys())
  found_keys = list(pop_number_dict[pop_num_dict_keys[0]].keys()) + \
               list(pop_number_dict[pop_num_dict_keys[1]].keys())

  for population in populations: 
    if population not in found_keys: 
      number_lines = [population + ":Number" in line for line in xml_lines]
      number_line_idx =  [i for i, x in enumerate(number_lines) if x][0]
      ## TODO: Could make this more robust by searching for last ANValue tag
      pop_number_string =  re.search('integer>.*<', xml_lines[number_line_idx-3])
      try: 
        pop_uid = int(pop_number_string.group(0)[8:-1])
        pop_num = object_table[pop_uid]
        # Figure out which M1/M2 it belongs to 
        assoc_filename_line = max(filter(lambda filename_line_idx: filename_line_idx < number_line_idx, filename_line_idxs))
        filename = filenames[filename_line_idxs.index(assoc_filename_line)]
        pop_number_dict[filename][population] = int(pop_num)
      except: 
        pass
      
  simplified_pop_number_dict = {}

  for k,v in pop_number_dict.items():
    key = re.search("m1|m2", k, re.IGNORECASE)
    if key is None: continue
    simplified_pop_number_dict[key.group(0).upper()] = pop_number_dict[k]
    
  # change back to orignal working directory
  os.chdir(owd)
            
  return(simplified_pop_number_dict)


if __name__ == "__main__":
  pp = pprint.PrettyPrinter(indent=4)
  output = parse_pop_info_for_flow_dir(sys.argv[1])
  output = compute_pop_proportions(output)
  output_csv_filepath =  os.path.join(sys.argv[1], 'population_proportion_info.csv')
  write_pop_info_to_csv(output, output_csv_filepath)
  print('\n\nSaved population proportion information to: ' + output_csv_filepath + '\n')
  pp.pprint(output)
  print('\n\n')