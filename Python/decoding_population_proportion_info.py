#!/bin/python3.8

import os
import re
import sys
import pprint

def compute_pop_proportions(pop_info_dict): 
  population_proportion_dict = {
    "CD34+ Blasts" : pop_info_dict['M1']['CD34+'] / pop_info_dict['M1']['WBC'],
    "CD117+ Cells" : pop_info_dict['M1']['CD117'] / pop_info_dict['M1']['WBC'],
    "B Cells" : pop_info_dict['M1']['B cells'] / pop_info_dict['M1']['WBC'],
    "Granulocytes" : pop_info_dict['M2']['Granulocytes'] / pop_info_dict['M2']['WBC'],
    "Monocytes": pop_info_dict['M2']['Monocytes'] / pop_info_dict['M2']['WBC'],
    "Lymphocytes": pop_info_dict['M1']['Lymphocytes'] / pop_info_dict['M1']['WBC'],
    "Basophils": pop_info_dict['M2']['Basophils'] / pop_info_dict['M2']['WBC'],
    "Eosinophils": pop_info_dict['M2']['Eosynophils'] / pop_info_dict['M2']['WBC'],
    "PCDC": pop_info_dict['M2']['PCDC'] / pop_info_dict['M2']['WBC'],
    "Erythroids": pop_info_dict['M1']['Eryhthroids'] / pop_info_dict['M1']['Singlets'],
    "Plasma Cells": pop_info_dict['M1']['Plasma cells'] / pop_info_dict['M1']['WBC']
  }


  for k,v in population_proportion_dict.items(): 
    population_proportion_dict[k] = round(population_proportion_dict[k]*100, 2)

  return(population_proportion_dict)

def parse_pop_info_for_flow_dir(flow_directory): 
  
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
      
  simplified_pop_number_dict = {}

  for k,v in pop_number_dict.items():
    key = re.search("m1|m2", k, re.IGNORECASE)
    if key is None: continue
    simplified_pop_number_dict[key.group(0).upper()] = pop_number_dict[k]
    
  # change back to orignal working directory
  os.chdir(owd)
            
  return(simplified_pop_number_dict)


# TODO: denote M1 vs M2 when adding to the dictionary to denote M1_WBC and M2_WBC, etc
# check if all events is relative to tube or all events (probably tube specific)
if __name__ == "__main__":
  pp = pprint.PrettyPrinter(indent=4)
  output = parse_pop_info_for_flow_dir(sys.argv[1])
  output = compute_pop_proportions(output)
  pp.pprint(output)
