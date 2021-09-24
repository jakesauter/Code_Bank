library(dplyr)
library(Seurat)
library(future)
library(flowCore)
library(patchwork)

SAVE_DIR <- "/home/sauterj1/mskmind_ess/seurat/objects"

COHORT <- c("3517371", "3471488", "3509885", "3460637", "3438445", "3497635", 
            "3497636", "3522361", "3512289", "3451481", "3433416", "3446317", 
            "3522435", "3459765", "3441110", "3059919", "3521647", "3407728", 
            "3509114", "3442027", "3435058", "3475561", "3452175", "3509175", 
            "3409341", "3506501", "3436846", "3493527", "3498381", "3469698", 
            "3448770", "3525781", "3447324", "3444351", "3447212", "3501884", 
            "3407722", "3446250", "3439300", "3478154", "3519922", "3446217", 
            "3509123", "3495775", "3446313", "3454373", "3470647", "3446251", 
            "3489962", "3440076", "3486706", "3520651", "3486030", "3516901", 
            "3459713", "3435057", "3438439", "3489948", "3468337", "3444378", 
            "3445395", "3407815", "3453122", "3462324", "3454298", "3463338", 
            "3511375", "3466022", "3498352", "3488263", "3516314", "3479927", 
            "3486634", "3505610", "3513838", "3490013", "3494186", "3494232", 
            "3519924", "3495780", "3496636", "3490886", "3501823", "3500059", 
            "3513804", "3511446", "3507200", "3509041", "3516375", "3516326", 
            "3519098", "3519909", "3520648", "3520759", "3523058", "3523934", 
            "3017935", "3161169", "3165874", "3107852", "3083434", "3222180", 
            "3188921", "2864151", "2860466", "3181538", "3398024", "3480817", 
            "3143119", "3188256", "3003151", "3220501", "2948833", "3373798", 
            "3109450", "3122225", "3218713", "3140541", "2941774", "3235190", 
            "2793301", "2761521", "3175622", "2973311", "2967895", "3216611", 
            "3118187", "3186522", "2967897", "2953538", "3102788", "3415201", 
            "3195222", "3224351", "3396954", "3394653", "3185533", "3071309", 
            "3244064", "3202287", "3107845", "2930920", "3490873", "2944746", 
            "3107036", "3376751", "3138131", "3212524", "3218698", "3224361", 
            "3006135", "3108648", "3200074", "3180617", "3426278", "3186573", 
            "3167829", "3235249", "3072949", "3283545", "3194445", "3193480", 
            "3400206", "3181574", "3059226", "3407813", "3150884", "3213259", 
            "3147164", "3210890", "3107875", "3110429", "3183214", "3374515", 
            "2930140", "3385615", "3118186", "3102046", "3232798", "3124259", 
            "2982912", "2793992", "3112061", "3107869", "3266538", "3382034", 
            "3151768", "3120040", "3120953", "3100260", "3086637", "3211620", 
            "3225939", "3375680", "3423620", "3185517", "3102047", "3127982", 
            "2861257", "3060881", "2956575", "3010680", "3429125", "3093541", 
            "3030356", "3002221", "3389418", "3204115", "3430703", "3161190", 
            "3108652", "3136771", "3190031", "3416735", "3180540", "3080723", 
            "3224358", "3069565", "3231002", "3403408", "3164106", "2997907", 
            "3479058", "3232797", "2980676", "3052003", "2844061", "3420660", 
            "3156134", "3234838", "2989108", "3096847", "3162676", "2968642", 
            "3150907", "3073800", "3094355", "2993499", "3196751")

# First Need to get the cells x markers matrix
flow_directories <- 
  list.files('/gpfs/mskmind_ess/sauterj1/Full_Cohort_Anon_Flow_Folders', 
             full.names = TRUE) %>% 
  .[basename(.) %in% COHORT]

for (assay in c('M1', 'M2')) {
  
  fcs_dat_list <-
    vector('list', length(flow_directories))
  
  for (i in seq_along(flow_directories)) {
    
    
    cat('Reading fcs file:', i, 'of', 
        length(flow_directories), '\n')
    
    fcs_path <- 
      file.path(flow_directories[[i]], 
                paste0('processed_', assay, '_randomly_subsampled_seed_497.fcs'))
    
    fcs_dat <- 
      exprs(read.FCS(fcs_path,
                     transformation = FALSE))
    
    rownames(fcs_dat) <- 
      rep(basename(flow_directories[i]), 
          nrow(fcs_dat))
    
    fcs_dat_list[[i]] <- fcs_dat
    
  }
  
  cat('\n\nCombining data into single matrix\n\n')
  fcs_data <- 
    do.call(rbind, fcs_dat_list)
  
  # Seurat expects cells as columns and
  # features as rows
  fcs_data <- t(fcs_data)
  
  # Create the metadata matrix
  meta.data <- 
    data.frame(flow_dir = colnames(fcs_data))
  
  # Each cell needs a uniqe name (lets make it the index)
  colnames(fcs_data) <- seq(1, ncol(fcs_data))
  
  cat('\n\nCreating Seurat object\n\n')
  # Create seurat object
  seurat_obj <- 
    CreateSeuratObject(
      counts = fcs_data,
      project = 'FCS Example',
      assay = 'MultiParamFlowCyto', 
      meta.data = meta.data)
  
  
  cat('\n\nSaving Seurat object\n\n')
  # Save Seurat object
  saveRDS(seurat_obj, 
          file.path(SAVE_DIR, paste0(assay, '_random_sample_centralized_Seurat_object.Rds')))
  
  seurat_expression_data <-
    t(as.matrix(GetAssayData(seurat_obj))) %>%
    as.data.frame()
  
  arrow::write_feather(seurat_expression_data,
                       file.path(SAVE_DIR, paste0(assay, '_random_sample_seurat_expression_data.feather')))
}



# 
# library(dplyr)
# library(Seurat)
# library(future)
# library(flowCore)
# library(patchwork)
# 
# 
# seurat_obj <- 
#   readRDS('/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_centralized_Seurat_object.Rds')
# 
# 
# seurat_expression_data <-
#   t(as.matrix(GetAssayData(seurat_obj))) %>%
#   as.data.frame()
# 
# arrow::write_feather(seurat_expression_data,
#                      '/gpfs/mskmind_ess/sauterj1/seurat/objects/M1_seurat_expression_data.feather')
# 
# 
# 
# seurat_obj <- 
#   readRDS('/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_centralized_Seurat_object.Rds')
# 
# 
# seurat_expression_data <-
#   t(as.matrix(GetAssayData(seurat_obj))) %>%
#   as.data.frame()
# 
# arrow::write_feather(seurat_expression_data,
#                      '/gpfs/mskmind_ess/sauterj1/seurat/objects/M2_seurat_expression_data.feather')
# 
# 
# 
# 

