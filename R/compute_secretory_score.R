library(WGCNA)
library(tibble)
library(readr)
library(dplyr)
library(circlize)
library(MotrpacRatTraining6moData)
library(tidyr)


### Load Uniprot DB for secretomes
sec_db <- read.delim("uniprotkb_proteome_UP000002494_AND_keyw_2025_02_17.tsv",
                     header = TRUE, sep = "\t", quote = "")
sec_db_expanded <- sec_db %>%
  separate_rows(Gene.Names, sep = " ")

# Merge sec_db with RAT_TO_HUMAN_GENE using Gene.Names (sec_db) and RAT_SYMBOL (RAT_TO_HUMAN_GENE)
sec_db_annotated <- sec_db_expanded %>%
  left_join(RAT_TO_HUMAN_GENE %>% select(RAT_SYMBOL, RAT_ENSEMBL_ID),
            by = c("Gene.Names" = "RAT_SYMBOL"))
sum(!is.na(sec_db_annotated$RAT_ENSEMBL_ID)) # 2897 non-NA RAT ENSEMBL

# Function to compute Ssec for each training group, only on filtered secretory features. 
compute_Ssec <- function(data, group_type, sec_db_annotated) {
  # Step 1: Extract p-value matrix
  group_data <- data[[group_type]]$p
  
  # Step 2: Valid Ensembl IDs
  valid_ensembl_ids <- sec_db_annotated$RAT_ENSEMBL_ID
  
  # Step 3: Subset matrix to secreted proteins
  group_data_filtered <- group_data[, colnames(group_data) %in% valid_ensembl_ids]
  
  # Step 4: Transform p-values
  group_data_transformed <- -log(group_data_filtered)
  group_data_transformed[is.infinite(group_data_transformed)] <- NA
  
  # Step 5: Compute Ssec
  Ssec <- colSums(group_data_transformed, na.rm = TRUE) / nrow(group_data_transformed)
  
  # Step 6: Build base data frame
  Ssec_df <- data.frame(
    RAT_ENSEMBL_ID = names(Ssec),
    Ssec = Ssec,
    row.names = NULL
  )
  
  # Step 7: Join gene symbol (may produce duplicates)
  Ssec_annotated <- dplyr::left_join(
    Ssec_df,
    sec_db_annotated[, c("RAT_ENSEMBL_ID", "Gene.Names")],
    by = "RAT_ENSEMBL_ID"
  ) %>%
    dplyr::rename(RAT_GENE_SYMBOL = Gene.Names)
  
  # ✅ Step 8: Keep only distinct (Ensembl ID + Symbol + Ssec)
  Ssec_unique <- Ssec_annotated %>%
    dplyr::distinct(RAT_ENSEMBL_ID, RAT_GENE_SYMBOL, Ssec, .keep_all = TRUE)
  
  # Step 9: Save
  file_name <- paste0(deparse(substitute(data)), "_", group_type, "_Ssec.csv")
  write.csv(Ssec_unique, file = file_name, row.names = FALSE)
  
  return(Ssec_unique)
}

# Load RDS files that were generated from compute_correlation.R
#1 WAT and BAT
wat_bat<- readRDS("./bicor_results/WATSC/bicor_BAT.rds")

# Example usage for both 'control' and 'training' groups - WATSC to BAT
control_Ssec_df <- compute_Ssec(wat_bat, "control", sec_db_annotated)
training1w_Ssec_df <- compute_Ssec(wat_bat, "training_1w", sec_db_annotated)
training2w_Ssec_df <- compute_Ssec(wat_bat, "training_2w", sec_db_annotated)
training4w_Ssec_df <- compute_Ssec(wat_bat, "training_4w", sec_db_annotated)
training8w_Ssec_df <- compute_Ssec(wat_bat, "training_8w", sec_db_annotated)

#load off wat_bat
rm(wat_bat)
gc()