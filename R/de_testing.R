library(reticulate)

# - CPM means -----------------------------------------------------------------
cpm_means <- function(cpm, info, group) {
  unique_groups <- unique(info[[group]])
  mean_vals <- function(x) {
    df <- cpm[, filter(info, !!sym(group) == x)$Sample_ID]
    if (is.null(ncol(df))) {
      warning(paste(x, "has only 1 sample. Cannot calculate mean."),
              call. = FALSE)
      return(df)
    } else {
      return(2^rowMeans(log2(df + 1))-1)
    }
  }
  df <- map(unique_groups, mean_vals) %>% bind_cols()
  df <- set_names(df, unique_groups)
  df$gene_id <- cpm$gene_id
  df <- select(df, gene_id, unique_groups)
  return(df)
}

# - Intermine -----------------------------------------------------------------
intermine <- function(genes) {
  source_python("~/Programs/rnaseq/python/mousemine.py")
  result <- clean_result(query_intermine(genes))
  result <- as.data.frame(result)
  result <- result[!duplicated(result), ]
  return(result)
}

# - GO analysis ---------------------------------------------------------------
run_go <- function(deseq_results, go) {
  deseq_results <- filter(deseq_results, !is.na(padj))
  deseq_results <- filter(deseq_results, gene_name %in% unique(go$gene))
  go <- filter(go, gene %in% deseq_results$gene_name)
  
  # set up population values for hypergeometric test
  sig_genes <- sum(deseq_results$padj < 0.05)
  all_genes <- nrow(deseq_results)
  go_test <- function(go_term) {
    all_go_genes <- sum(go$GO == go_term)
    sig_go_genes <- sum(
      go$GO == go_term & 
        go$gene %in% filter(deseq_results, padj < 0.05)$gene_name
    )
    phyper(sig_go_genes, sig_genes, all_genes - sig_genes, all_go_genes)
  }
  
  # run hypergeometric test on all GO terms and adjust with Bonferroni
  go_terms <- unique(go$GO)
  results <- map_dbl(go_terms, go_test)
  results <- p.adjust(results, method = "bonferroni")
  names(results) <- go_terms
  return(results)
}
