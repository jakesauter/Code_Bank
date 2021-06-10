---
  title: "FCS Transformation and Normalization"
author: "Jake Sauter"
date: "6/10/2021"
output: html_document
editor_options: 
  chunk_output_type: console
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE,
                      message = FALSE, 
                      comment = NA)
```

### Input data

```{r}
library(flowCore)

flow_dir_save_loc <- 
  "/Users/sauterj1/Documents/Woodlist/Flow_Folders"

flow_directories <- 
  list.files(flow_dir_save_loc, 
             full.names = TRUE)

flow_directories <- 
  sample(flow_directories, 10)

fcs_list <- vector('list', length(flow_directories))
names(fcs_list) <- flow_directories

for (i in seq_along(flow_directories)) {
  
  flow_directory <- flow_directories[[i]]
  cat("\n\nFlow Directory: ", flow_directory, "\n\n")
  
  fcs_files <- 
    list.files(flow_directory, 
               pattern = "\\.fcs") %>% 
    .[str_detect(., 'compensated')] %>% 
    .[str_detect(., 'high_quality_events')]
  
  fcs_file <- 
    str_detect(fcs_files, 
               'M2|m2') %>% 
    fcs_files[.] %>% 
    file.path(flow_directory, .)
  
  print(fcs_file)
  
  fcs_list[[i]] <- 
    flowCore::read.FCS(fcs_file, 
                       transformation = FALSE)
  
}
```

### Normalization 

First big problem, need to modify the metadata so that we can merge together into a flowFrame

```{r}
library(flowStats)
library(flowCore)
library(stringr)

# Need to compose each flowFrame to have the
# same parameter names
for(i in seq_along(fcs_list)) {
  fcs_data <- fcs_list[[i]]
  fcs_expr <- exprs(fcs_data)
  metadata <- pData(parameters(fcs_data))
  
  to_keep <- ifelse(is.na(metadata$desc), FALSE, TRUE)
  to_keep[1:4] <- TRUE
  
  names <- metadata$desc[to_keep]
  names[1:4] <- colnames(fcs_expr)[1:4]
  names <- str_extract(names, 
                       'FSC-A|FSC-H|SSC-A|SSC-H|CD[0-9]+b?|HLA')
  names[names == 'HLA'] = 'HLA-DR'
  
  new_expr <- fcs_expr[, to_keep] 
  colnames(new_expr) <- names
  
  # To ensure same order of column names. Needed
  # to convert to flowset
  if (i > 1) {
    names <- colnames(exprs(fcs_list[[1]]))
    new_expr <- new_expr[, names] 
  }
  
  fcs_list[[i]] <- flowFrame(new_expr)
}

fcs_flowset <- as(fcs_list, 'flowSet')
```

### Normalizing using gaussNorm

```{r}
library(flowStats)

channel_names <- colnames(exprs(fcs_flowset[[1]]))
channel_names <- channel_names[5:length(channel_names)]

normalized_fcs_flowset <- 
  flowStats::warpSet(fcs_flowset, channel_names)
```








