library(reticulate)
# - Intermine -----------------------------------------------------------------
intermine <- function(genes) {
  source_python("~/Programs/rnaseq/python/mousemine.py")
  result <- clean_result(query_intermine(genes))
  result <- as.data.frame(result)
  result <- result[!duplicated(result), ]
  return(result)
}
