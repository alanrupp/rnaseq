library(tidyverse)

# - Read counts ---------------------------------------------------------------
read_counts <- function(files, column = "X2") {
  if (str_detect(files, "Sample")) {
    samples <- str_extract(files, "Sample_[0-9]+")
  } else {
    samples <- paste0("Sample_", str_extract(files, "[:alnum:]+(?=Reads)"))
  }
  
  # read data
  if (length(files) > 1) {
    map(files, ~ read_tsv(.x, skip = 4, col_names = FALSE)) %>%
      bind_cols() %>%
      select(`X1...1`, starts_with(column)) %>%
      set_names("gene_id", samples)
  } else {
    read_tsv(files, skip = 4, col_names = FALSE) %>%
      select(`X1`, !!sym(column)) %>%
      set_names("gene_id", samples)
  }
}

# - Make counts to CPM --------------------------------------------------------
make_cpm <- function(counts, log2 = FALSE) {
  if (log2) {
    df <- mutate_if(counts, is.numeric, ~ log2(.x / sum(.x) * 10^6 + 1))
  } else {
    df <- mutate_if(counts, is.numeric, ~ .x / sum(.x) * 10^6)
  }
  return(df)
}

# - Find genes expressed above given threshold --------------------------------
expressed_genes <- function(counts, cpm_threshold = 1, min_samples = 3,
			    gene_col = "gene_id") {
  counts <- make_cpm(counts)
  above_threshold <- mutate_if(counts, is.numeric, ~ .x > cpm_threshold)
  above_threshold <- rowSums(select_if(counts, is.numeric)) >= min_samples
  return(counts[above_threshold, ] %>% pull(gene_col))
}

# - Choose PCs from PCA by bootstrapping --------------------------------------
choose_pcs <- function(pca, matrix, reps = 100, sig_level = 0.05,
		       progress_bar = FALSE) {
  mtrx <- as.matrix(matrix)
  randomize <- function(mtrx) {
    result <- sapply(seq(nrow(mtrx)), function(x) sample(mtrx[x, ]))
    return(result)
  }
  find_eigenvalues <- function(matrix) {
    pca <- prcomp(randomize(matrix), scale. = TRUE)
    return(pca$sdev^2)
  }
  if (progress_bar) replicate <- pbapply::pbreplicate
  eigenvalues <- replicate(reps, find_eigenvalues(mtrx))
  # find the frequency that real eigenvalue is above bootstrap eigenvalue
  prob <- 1 - rowSums(pca$sdev^2 > eigenvalues) / reps
  sig_pcs <- which(prob < sig_level)
  return(sig_pcs)
}

# - Correlation ---------------------------------------------------------------
correlation <- function(counts, genes = NULL) {
  if (!is.null(genes)) {
    counts <- filter(counts, gene_id %in% genes) 
  }
  corr <- cor(select(counts, -gene_id))
  corr <- as.data.frame(corr) %>%
    rownames_to_column("Sample_A") %>%
    gather(-Sample_A, key = "Sample_B", value = "corr")
  return(corr)
}
