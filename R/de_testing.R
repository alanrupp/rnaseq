library(reticulate)

# - CPM means -----------------------------------------------------------------
cpm_means <- function(cpm, info, group) {
  unique_groups <- unique(info[[group]])
  mean_vals <- function(x) {
    df <- cpm[, dplyr::filter(info, !!sym(group) == x)$Sample_ID]
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
  df <- dplyr::select(df, gene_id, unique_groups)
  return(df)
}

# - Intermine -----------------------------------------------------------------
intermine <- function(genes) {
  source_python("~/Programs/rnaseq/python/mousemine.py")
  result <- clean_result(genes)
  result <- result[!duplicated(result), ]
  return(result)
}

# - GO analysis ---------------------------------------------------------------
run_go <- function(deseq_results, ontology = NULL) {
  go <- read_csv("~/Programs/rnaseq/data/mgi.csv.gz", col_types = "cc")
  terms <- read_tsv("~/Programs/rnaseq/data/go_terms.mgi", col_names = FALSE,
                    col_types = "ccc")
  if (!is.null(ontology)) {
    terms <- filter(terms, X1 == ontology)
    go <- filter(go, GO %in% terms$X2)
  }
  deseq_results <- filter(deseq_results, !is.na(padj))
  deseq_results <- filter(deseq_results, gene_name %in% unique(go$gene))
  go <- filter(go, gene %in% deseq_results$gene_name)
  
  # set up population values for hypergeometric test
  sig_genes <- sum(deseq_results$padj < 0.05)
  all_genes <- nrow(deseq_results)
  
  # remove GO terms that are underpowered (don't even test)
  max_go_genes <- seq(max(table(go$GO)))
  power_test <- phyper(max_go_genes-1, sig_genes, all_genes-sig_genes, 
                       max_go_genes, lower.tail = FALSE)
  if (sum(power_test >= 0.05) > 0) {
    underpowered <- which(power_test >= 0.5)
    print(paste("Removing", max(underpowered), "GO categories because they are",
                "underpowered."))
    underpowered_categories <- go %>% group_by(GO) %>% count() %>%
      filter(n <= max(underpowered)) %>% .$GO
    go <- filter(go, !GO %in% underpowered_categories)
    deseq_results <- filter(deseq_results, gene_name %in% unique(go$gene))
    go <- filter(go, gene %in% deseq_results$gene_name)
  } else {
    print("Running on all GO terms because there is enough power to detect.")
  }
  
  # run hypergeometric test on all GO terms and adjust with FDR
  go_terms <- unique(go$GO)
  print(paste("Running GO for", length(go_terms), "terms"))
  go_test <- function(go_term) {
    all_go_genes <- sum(go$GO == go_term)
    sig_go_genes <- sum(
      go$GO == go_term & 
        go$gene %in% filter(deseq_results, padj < 0.05)$gene_name
    )
    phyper(sig_go_genes-1, sig_genes, all_genes-sig_genes, all_go_genes,
           lower.tail = FALSE)
  }
  results <- pbapply::pbsapply(go_terms, go_test)
  results <- p.adjust(results, method = "BH")
  
  # export as data.frame
  results <- data.frame("GO" = go_terms, "P" = results)
  results <- left_join(results, terms, by = c("GO" = "X2"))
  results <- dplyr::rename(results, "Ontology" = X1, "Name" = X3)
  results <- arrange(results, P)
  return(dplyr::select(results, GO, Name, Ontology, P))
}