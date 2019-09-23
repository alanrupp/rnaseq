library(tidyverse)

# - Read counts ---------------------------------------------------------------
read_counts <- function(files, column = "X2") {
  samples <- str_extract(files, "Sample_[0-9]+")
  map(files, ~ read_tsv(.x, skip = 4, col_names = FALSE)) %>%
    bind_cols() %>%
    select(X1, starts_with(column)) %>%
    set_names("gene_id", samples)
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
expressed_genes <- function(counts, cpm_threshold = 1, min_samples = 3) {
  counts <- make_cpm(counts)
  above_threshold <- mutate_if(counts, is.numeric, ~ .x > cpm_threshold)
  above_threshold <- rowSums(select_if(counts, is.numeric)) >= min_samples
  return(counts$gene_id[above_threshold])
}

