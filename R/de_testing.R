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

unique_groups <- unique(info$Treatment)
map(unique_groups,
    ~ apply(., 2, function(x) 
      mean(log2(cpm[, filter(bead_info, Treatment == x)$Sample_ID] + 1)))
    )


apply(cpm[, filter(bead_info, Treatment == "PBS")$Sample_ID], 1, mean)
