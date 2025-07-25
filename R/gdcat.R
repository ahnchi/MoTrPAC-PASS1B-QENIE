## Gene Derived Correlation Across Tissues (GD-CAT)
## This runs fgsea in target tissue genes that correlate with the given origin tissue gene. 
## Correlation files need to be generated from compute_correlation.R in order to process


remotes::install_github("igordot/msigdbr@v7.5.1")

library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)
library(rlang)


gene_sets <- msigdbr(species = "Rattus norvegicus", category = "C5") 
contrast_order <- c("training_1w", "training_2w", "training_4w", "training_8w", "control")



check_fgsea <- function(Direction, target_gene) {
  
  # Step 1: Extract bicor vectors for the target gene across conditions
  extract_bicor_values <- function(Direction, target_gene) {
    bicor_values <- list()
    for (condition in names(Direction)) {
      bicor_matrix <- Direction[[condition]]$bicor
      if (target_gene %in% colnames(bicor_matrix)) {
        bicor_values[[condition]] <- bicor_matrix[, target_gene]
      } else {
        warning(paste("Gene", target_gene, "not found in condition", condition))
        bicor_values[[condition]] <- NULL
      }
    }
    return(bicor_values)
  }
  
  bicor_values <- extract_bicor_values(Direction, target_gene)
  
  # Step 2: Remove any NULL or empty entries
  bicor_values <- lapply(bicor_values, function(vec) {
    if (is.null(vec)) return(NULL)
    vec <- vec[!is.na(names(vec)) & names(vec) != ""]
    return(vec)
  })
  bicor_values <- Filter(Negate(is.null), bicor_values)  # remove NULLs
  
  # Step 3: Prepare FGSEA gene sets (Ensembl-based)
  gene_sets_formatted <- gene_sets %>%
    dplyr::select(gs_name, ensembl_gene) %>%
    dplyr::group_by(gs_name) %>%
    dplyr::summarise(genes = list(unique(ensembl_gene)), .groups = "drop")
  
  pathways <- setNames(gene_sets_formatted$genes, gene_sets_formatted$gs_name)
  
  # Step 4: Run FGSEA for each contrast
  fgsea_results_list <- list()
  
  for (contrast in names(bicor_values)) {
    vec <- bicor_values[[contrast]]
    
    ranked_df <- data.frame(
      gene = names(vec),
      score = vec,
      stringsAsFactors = FALSE
    ) %>%
      dplyr::filter(!is.na(gene) & gene != "unknown") %>%
      dplyr::group_by(gene) %>%
      dplyr::summarise(score = ifelse(all(score > 0), max(score),
                                      ifelse(all(score < 0), min(score),
                                             ifelse(abs(max(score)) > abs(min(score)), max(score), min(score)))),
                       .groups = "drop") %>%
      dplyr::arrange(desc(score))
    
    stats <- ranked_df$score
    names(stats) <- ranked_df$gene
    stats <- stats[is.finite(stats)]
    
    # Safety check for sufficient overlap
    overlap <- length(intersect(names(stats), unlist(pathways)))
    if (overlap < 15) {
      message("Skipping contrast ", contrast, ": insufficient overlap with gene sets.")
      next
    }
    
    fgsea_res <- fgsea::fgseaMultilevel(
      pathways = pathways,
      stats = stats,
      minSize = 15,
      maxSize = 500
    )
    
    fgsea_res$contrast <- contrast
    fgsea_results_list[[contrast]] <- fgsea_res
  }
  
  if (length(fgsea_results_list) == 0) {
    warning("No FGSEA results computed.")
    return(NULL)
  }
  
  all_fgsea_results <- do.call(rbind, lapply(fgsea_results_list, as.data.frame))
  assign("FGSEA_RESULTS", all_fgsea_results, envir = .GlobalEnv)
  return(all_fgsea_results)
}

#Example use
fgsea_tgfb2 <- check_fgsea(wat_bat, "ENSRNOG00000002418") # ENSRNOG00000002418 = TGFb2. Testing BAT genes that correlate with WATSC-TGFb2

