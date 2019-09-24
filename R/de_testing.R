library(reticulate)

# - CPM means -----------------------------------------------------------------
cpm_means <- function(cpm, info, group) {
  map(unique(as.data.frame(info)[, group]),
        ~ rowMeans(log2(cpm[, filter(info, !!sym(group) == .x)$Sample_ID] + 1))
        ) %>%
    as.data.frame() %>%
    set_names(unique(as.data.frame(info)[, group])) %>%
    mutate_all(~ 2^.x - 1) %>%
    as.data.frame() %>%
    mutate("gene_id" = cpm$gene_id) %>%
    select(gene_id, unique(as.data.frame(info)[, group]))
}

# - Intermine -----------------------------------------------------------------
intermine <- function(genes) {
  source_python("~/Programs/rnaseq/python/mousemine.py")
  result <- clean_result(query_intermine(genes))
  result <- as.data.frame(result)
  result <- result[!duplicated(result), ]
  return(result)
}


