library(WGCNA)
library(tibble)
library(readr)
library(dplyr)
library(circlize)
library(MotrpacRatTraining6moData)
library(tidyr)


RTHG <- RAT_TO_HUMAN_GENE
pass1b_pheno <- PHENO
bat_t <- TRNSCRPT_BAT_NORM_DATA
hip_t <- TRNSCRPT_HIPPOC_NORM_DATA
adr_t <- TRNSCRPT_ADRNL_NORM_DATA
col_t <- TRNSCRPT_COLON_NORM_DATA
cor_t <- TRNSCRPT_CORTEX_NORM_DATA
hea_t <- TRNSCRPT_HEART_NORM_DATA
hyp_t <- TRNSCRPT_HYPOTH_NORM_DATA
kid_t <- TRNSCRPT_KIDNEY_NORM_DATA
liv_t <- TRNSCRPT_LIVER_NORM_DATA
lun_t <- TRNSCRPT_LUNG_NORM_DATA
ova_t <- TRNSCRPT_OVARY_NORM_DATA
sgn_t <- TRNSCRPT_SKMGN_NORM_DATA
svl_t <- TRNSCRPT_SKMVL_NORM_DATA
sml_t <- TRNSCRPT_SMLINT_NORM_DATA
spl_t <- TRNSCRPT_SPLEEN_NORM_DATA
tes_t <- TRNSCRPT_TESTES_NORM_DATA
ven_t <- TRNSCRPT_VENACV_NORM_DATA
wat_t <- TRNSCRPT_WATSC_NORM_DATA

control_8w_sample <- rownames(pass1b_pheno %>%
                                filter(study_group_timepoint == "8 week control group (pass1b)"))

training_8w_sample <- rownames(pass1b_pheno %>%
                                 filter(study_group_timepoint == "8 weeks training group (pass1b)"))
training_1w_sample <- rownames(pass1b_pheno %>%
                                 filter(study_group_timepoint == "1 week training group (pass1b)"))
training_2w_sample <- rownames(pass1b_pheno %>%
                                 filter(study_group_timepoint == "2 weeks training group (pass1b)"))
training_4w_sample <- rownames(pass1b_pheno %>%
                                 filter(study_group_timepoint == "4 weeks training group (pass1b)"))

subset_seq_data5 <- function(seq_data, control_samples, training_1w_samples, training_2w_samples, training_4w_samples, training_8w_samples, pheno_data) {
  matching_control_pids <- intersect(control_samples, colnames(seq_data))
  matching_training_1w_pids <- intersect(training_1w_samples, colnames(seq_data))
  matching_training_2w_pids <- intersect(training_2w_samples, colnames(seq_data))
  matching_training_4w_pids <- intersect(training_4w_samples, colnames(seq_data))
  matching_training_8w_pids <- intersect(training_8w_samples, colnames(seq_data))
  
  control_data <- seq_data %>% dplyr::select(feature_ID, tissue, assay, all_of(matching_control_pids))
  training_1w_data <- seq_data %>% dplyr::select(feature_ID, tissue, assay, all_of(matching_training_1w_pids))
  training_2w_data <- seq_data %>% dplyr::select(feature_ID, tissue, assay, all_of(matching_training_2w_pids))
  training_4w_data <- seq_data %>% dplyr::select(feature_ID, tissue, assay, all_of(matching_training_4w_pids))
  training_8w_data <- seq_data %>% dplyr::select(feature_ID, tissue, assay, all_of(matching_training_8w_pids))
  
  return(list(
    control = control_data, training_1w = training_1w_data,
    training_2w = training_2w_data, training_4w = training_4w_data,
    training_8w = training_8w_data
  ))
}
setwd("../data")
tissue_list <- list(
  "BAT" = bat_t,
  "HIPPO" = hip_t,
  "ADREN" = adr_t,
  "COLON" = col_t,
  "CORTEX" = cor_t,
  "HEART" = hea_t,
  "HYPOTH" = hyp_t,
  "KIDNEY" = kid_t,
  "LIVER" = liv_t,
  "LUNG" = lun_t,
  "SKM-GN" = sgn_t,
  "SKM-VL" = svl_t,
  "SMLINT" = sml_t,
  "SPLEEN" = spl_t,
  "VENACV" = ven_t,
  "WATSC" = wat_t
)
run_bicorrelation_analysis5 <- function(
    origin_tissue,
    origin_tissue_name,
    tissue_list, 
    control_samples, 
    training_1w_samples, 
    training_2w_samples, 
    training_4w_samples, 
    training_8w_samples, 
    pheno_data
) {
  results_dir <- paste0("bicor_results/", origin_tissue_name, "/")
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  for (target_tissue_name in names(tissue_list)) {
    if (target_tissue_name == origin_tissue_name) next
    
    cat("\nProcessing:", target_tissue_name, "...\n")
    
    target_tissue <- tissue_list[[target_tissue_name]]
    
    tissue_groups <- subset_seq_data5(
      target_tissue, control_samples, training_1w_samples,
      training_2w_samples, training_4w_samples, training_8w_samples, pheno_data
    )
    
    origin_groups <- subset_seq_data5(
      origin_tissue, control_samples, training_1w_samples,
      training_2w_samples, training_4w_samples, training_8w_samples, pheno_data
    )
    
    cor_results <- list()
    for (grp in names(tissue_groups)) {
      # Fix x:
      x_df <- tissue_groups[[grp]] %>%
        dplyr::select(-tissue, -assay)
      
      x <- x_df %>%
        column_to_rownames("feature_ID") %>%
        as.matrix() %>%
        t()
      
      # Fix y:
      y_df <- origin_groups[[grp]] %>%
        dplyr::select(-tissue, -assay)
      
      y <- y_df %>%
        column_to_rownames("feature_ID") %>%
        as.matrix() %>%
        t()
      
      cor_result <- bicorAndPvalue(x, y)
      
      # Add row/col names for traceability
      rownames(cor_result$bicor) <- colnames(x)
      colnames(cor_result$bicor) <- colnames(y)
      rownames(cor_result$p) <- colnames(x)
      colnames(cor_result$p) <- colnames(y)
      
      cor_results[[grp]] <- cor_result
    }
    
    saveRDS(cor_results, file = paste0(results_dir, "bicor_", target_tissue_name, ".rds"))
    cor_results <- NULL
    gc()
    
    cat(" - Correlation completed and saved for", target_tissue_name, "\n")
  }
  
  cat("\n✅ All correlations completed for origin tissue:", origin_tissue_name, "\n")
}

run_bicorrelation_analysis5(
  origin_tissue = wat_t,    # Change this to any origin tissue
  origin_tissue_name = "WATSC",  # Make sure to match tissue name format
  tissue_list = tissue_list, # The list of all tissues
  control_samples = control_8w_sample, 
  training_1w_samples = training_1w_sample,
  training_2w_samples = training_2w_sample,
  training_4w_samples = training_4w_sample,
  training_8w_samples = training_8w_sample,
  pheno_data = pass1b_pheno
) # This generates one RDS file per origin -> target tissue pair. Each correlation can take hours to complete. With 92GB RAM, average time to complete each correlation is ~20min


