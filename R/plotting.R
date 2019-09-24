library(tidyverse)
library(ggrepel)
library(wesanderson)

# - Read depth ----------------------------------------------------------------
# ENCODE suggests read depth > 30M
read_depth <- function(counts) {
  depth <- colSums(select_if(counts, is.numeric))
  return(depth)
}

# - Plot count data -----------------------------------------------------------
counts_plot <- function(counts, gene_ids, metadata, group = NULL,
                        color = NULL, shape = NULL, dodge = NULL,
                        pair = FALSE, cpm = TRUE) {
  # set up data.frame with counts and metadata
  metadata <- select(metadata, Sample_ID, pair, cells, quo(color), quo(shape))
  
  # make CPM
  if (cpm) {
    df <- make_cpm(counts)
  } else {
    df <- counts
  }
  
  # add metadata
  df <- df %>%
    filter(gene_id %in% gene_ids) %>%
    gather(key = "Sample_ID", value = "expr") %>%
    left_join(., metadata, by = "Sample_ID") %>%
    left_join(., genes, by = )
  
  # plot
  p <- ggplot(df, aes(x = !!dodge, y = expr, group = !!group)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(trans = "log2") +
    theme_classic() +
    ylab("Expression (CPM)")
  if (pair) {
    
  }
  # custom color aesthetic
  if (!is.null(color)) {
    p <- p + 
      geom_line(aes(color = !!color)) + 
      geom_point(aes(color = !!color))
  }
  # custom shape aesthetic
  if (!is.null(shape)) {
    p <- p + 
      geom_line(aes(shape = shape)) + 
      geom_point(aes(shape = shape))
  }
  # facet for multiple genes
  if (length(gene_ids) > 1) {
    p <- p + facet_wrap(~gene_ids)
  }
  return(p)
}

# - Counts heatmap ------------------------------------------------------------
counts_heatmap <- function(counts, metadata = NULL, annotation = NULL) {
  
}


# - Volcano plot --------------------------------------------------------------
volcano_plot <- function(results, method = "DESeq2", label = NULL) {
  # set up y-axis scale
  padj_min <- min(results$padj, na.rm = TRUE)
  if (-log10(padj_min) < 10) {
    
  } else {
    y_scale <- 10^-(0:as.integer(-log10(padj_min)+1))
    if (length(y_scale) > 10) {
      y_scale <- seq(from = -log10(y_scale[1]), 
                     to = -log10(y_scale[length(y_scale)]), 
                     length.out = 10)
    }
  }
  # set up x axis scale
  min_fc <- abs(max(filter(results, padj < 0.05)))
  

  # make plot
  p <- ggplot(results, 
              aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)) +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
    geom_vline(aes(xintercept = 0)) +
    geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
    geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") +
    geom_point(show.legend = FALSE, stroke = 0, alpha = ) +
    scale_color_manual(values = c("gray", "firebrick3")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = expression("Fold change (log"[2]*")"), 
         y = expression(italic("P")*" value")) +
    scale_x_continuous(limits = c(-min_fc, min_fc)) #+
    #scale_y_continuous(breaks = y_scale,
    #                   labels = 10^-y_scale,
    #                   expand = c(0, 0)) 
  
  if (label) {
    df <- left_join(results, genes, by = "gene_id") %>%
      filter(gene_name %in% label)
    p <- p + 
      geom_text_repel(data = df, aes(label = gene_name), color = "black")
  }
}


# - Heatmap plot --------------------------------------------------------------
heatmap_plot <- function(counts, genes, 
                         info = NULL, annotation = NULL) {
  counts <- make_cpm(counts, log2 = TRUE)
  counts <- filter(counts, gene_id %in% genes)
  
  # cluster genes and samples
  gene_clust <- counts %>% as.data.frame() %>%
    column_to_rownames("gene_id") %>% dist() %>% hclust()
  sample_clust <- hclust(dist(t(counts[,2:ncol(counts)])))
  
  # make tidy
  df <- gather(counts, -gene_id, key = "Sample_ID", value = "counts") %>%
    mutate(gene_id = factor(gene_id, 
                            levels = gene_clust$labels[gene_clust$order]),
           Sample_ID = factor(Sample_ID,
                              levels = sample_clust$labels[sample_clust$order]))
  
  # plot
  plt <- ggplot(df, aes(x = Sample_ID, y = as.numeric(gene_id), fill = counts)) +
    geom_tile() +
    scale_fill_gradient(low = "#273046", high = "#FAD510",
                        name = expression(underline("log"[2]*"CPM"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     color = "black"),
          axis.text.y = element_blank()) +
    ylab(NULL) +
    xlab(NULL) +
    scale_y_continuous(expand = c(0, 0))
  
  # add annotation
  if (!is.null(info)) {
    anno <- left_join(df, info, by = "Sample_ID")
    anno <- anno %>% mutate(
      Sample_ID = factor(Sample_ID, levels = sample_clust$labels[sample_clust$order])
    ) %>% as.data.frame()
  }
  
  # add annotation (optional)
  add_annotation_rect <- function(i, j) {
    annotate("rect", 
             xmin = which(levels(anno$Sample_ID) == samples[j]) - 0.5, 
             xmax = which(levels(anno$Sample_ID) == samples[j]) + 0.5, 
             ymin = max_yval + block_size * (i - 0.5), 
             ymax = max_yval + block_size * (i + 0.5),
             fill = colors[anno[anno$Sample_ID == samples[j], annotation[i]]]
    )
  }
  add_annotation_text <- function(i) {
    annotate("label", x = length(samples)/2 + 0.5, 
             y = max_yval + block_size * i, 
             label = annotation[i], hjust = 0.5, color = "gray90", 
             fill = "gray10", alpha = 0.5, label.size = NA)
  }
  if (!is.null(info) & !is.null(annotation)) {
    max_yval <- length(unique(df$gene_id))
    block_size <- max_yval * 0.05
    samples <- unique(df$Sample_ID)
    palettes <- base::sample(1:length(wes_palettes), length(annotation),
                             replace = FALSE)
    for (i in 1:length(annotation)) {
      colors <- wes_palette(palettes[i])
      classes <- unique(anno[, annotation[i]])
      colors <- base::sample(colors, length(classes), replace = FALSE)
      names(colors) <- classes
      for (j in 1:length(samples)) {
        plt <- plt + add_annotation_rect(i, j)
      }
      plt <- plt + add_annotation_text(i)
    }
  } else if (!is.null(annotation) & is.null(info)) {
    message("Cannot add annotation to a plot without an info object")
  }
  
  # return plot
  return(plt)
}


# - TSNE ----------------------------------------------------------------------
tsne_plot <- function(counts, info = NULL, color = NULL, shape = NULL,
                      text = NULL) {
  df <- counts %>%
    filter(gene_id %in% expressed_genes(counts)) %>%
    make_cpm(., log2 = TRUE)
  
  # make df with rownames
  df <- df %>% as.data.frame() %>% column_to_rownames("gene_id")
  
  # run PCA
  pca <- prcomp(t(df), scale. = TRUE)
  pcs <- choose_pcs(pca, df)
  
  # run TSNE from PCA
  set.seed(32)
  tsne <- Rtsne::Rtsne(pca$x[, pcs], perplexity = (nrow(pca$x) - 1) / 3)
  coord <- 
    as.data.frame(tsne$Y) %>%
    set_names("TSNE1", "TSNE2") %>%
    mutate("Sample_ID" = rownames(pca$x))
  
  if (!is.null(info)) {
    coord <- left_join(coord, info, by = "Sample_ID")
    color_var <- sym(color)
    shape_var <- sym(shape)
    text_var <- sym(text)
  }
  
  # plot
  plt <- ggplot(coord, aes(x = TSNE1, y = TSNE2)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (!is.null(info)) {
    plt <- plt + 
      geom_point(aes(color = !!color_var, shape = !!shape_var)) +
      geom_text_repel(aes(label = Sample_ID, color = !!text_var), 
                      show.legend = FALSE) +
      scale_color_manual(values = wes_palette("Cavalcanti1"))
  }
  return(plt)
}


# - Correlation ---------------------------------------------------------------
# ENCODE suggests correlation of replicates should have Spearman > 0.9
correlation_plot <- function(counts, genes = NULL, info = NULL,
                             annotation = NULL) {
  # grab data
  corr <- correlation(counts, genes)
  
  # cluster genes and samples
  sample_clust <- corr %>%
    spread(key = "Sample_B", value = "corr") %>%
    as.data.frame() %>%
    column_to_rownames("Sample_A") %>%
    dist() %>%
    hclust()
  
  corr <- corr %>% mutate(
    Sample_A = factor(Sample_A, levels = sample_clust$labels[sample_clust$order]),
    Sample_B = factor(Sample_B, levels = sample_clust$labels[sample_clust$order])
    )
  
  # plot
  plt <- ggplot(corr, aes(x = Sample_A, y = Sample_B, fill = corr)) +
    geom_tile() +
    scale_fill_gradient(low = "#273046", high = "#FAD510",
                         name = expression(underline("Correlation"))) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1),
          axis.title = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank())
  
  
  # add annotation
  if (!is.null(info)) {
    anno <- left_join(corr, info, by = c("Sample_A" = "Sample_ID"))
    anno <- anno %>% mutate(
      Sample_A = factor(Sample_A, levels = sample_clust$labels[sample_clust$order])
    )
  }

  # add annotation (optional)
  add_annotation_rect <- function(i, j) {
    annotate("rect", 
             xmin = which(levels(anno$Sample_A) == samples[j]) - 0.5, 
             xmax = which(levels(anno$Sample_A) == samples[j]) + 0.5, 
             ymin = max_yval + i - 0.5, 
             ymax = max_yval + i + 0.5,
             fill = colors[anno[anno$Sample_A == samples[j], annotation[i]]]
    )
  }
  add_annotation_text <- function(i) {
    annotate("label", x = length(samples)/2 + 0.5, y = max_yval + i, 
             label = annotation[i], hjust = 0.5, color = "gray90", 
             fill = "gray10", alpha = 0.5, label.size = NA)
  }
  if (!is.null(info) & !is.null(annotation)) {
    max_yval <- length(unique(corr$Sample_B))
    samples <- unique(corr$Sample_A)
    palettes <- base::sample(seq(length(wes_palettes)), length(annotation),
                             replace = FALSE)
    for (i in 1:length(annotation)) {
      colors <- wes_palette(palettes[i])
      classes <- unique(anno[,annotation[i]])
      colors <- base::sample(colors, length(classes), replace = FALSE)
      names(colors) <- classes
      for (j in 1:length(samples)) {
        plt <- plt + add_annotation_rect(i, j)
      }
      plt <- plt + add_annotation_text(i)
    }
  } else if (!is.null(annotation) & is.null(info)) {
    message("Cannot add annotation to a plot without an info object")
  }
  
  return(plt)
}


