library(tidyverse)
library(ggrepel)
library(wesanderson)
library(ggdendro)

# - Read depth ----------------------------------------------------------------
# ENCODE suggests read depth > 30M
read_depth <- function(counts) {
  depth <- colSums(select_if(counts, is.numeric))
  return(depth)
}

# - Plot count data -----------------------------------------------------------
counts_plot <- function(counts, gene_names, info, samples = NULL,
                        group = NULL, color = NULL, shape = NULL, dodge = NULL,
                        pair = FALSE, cpm = TRUE, n_col = NULL) {
  if (!is.null(samples)) {
    counts <- select(counts, gene_id, samples)
  }
  
  # make CPM
  if (cpm) {
    df <- make_cpm(counts)
  } else {
    df <- counts
  }
  
  # add metadata
  df <- left_join(df, genes, by = "gene_id") %>%
    filter(gene_name %in% gene_names) %>%
    select(gene_name, starts_with("Sample")) %>%
    gather(-gene_name, key = "Sample_ID", value = "expr") %>%
    left_join(., info, by = "Sample_ID")
  
  # plot
  p <- ggplot(df, aes(y = expr)) +
    scale_y_continuous(trans = "log2") +
    theme_minimal() +
    ylab("Expression (CPM)") +
    xlab(NULL)
  
  # custom color aesthetic
  if (!is.null(dodge)) {
    p <- p + geom_boxplot(aes(x = !!sym(dodge))) +
      geom_jitter(aes(x = !!sym(dodge)), height = 0, width = 0.2)
  } else {
    p <- p + geom_boxplot(aes(x = '')) +
      geom_jitter(aes(x = ''), height = 0, width = 0.2)
  }
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
  if (length(gene_names) > 1) {
    p <- p + facet_wrap(~gene_name, ncol = n_col, scales = "free_y")
  }
  return(p)
}


# - Volcano plot --------------------------------------------------------------
volcano_plot <- function(results, label = NULL, xmax = NULL, ymax = NULL) {
  # set up axis scales
  if (is.null(xmax)) {
    if (sum(results$padj < 0.05, na.rm = TRUE) > 0) {
      xmax <- max(abs(filter(results, padj < 0.05)$log2FoldChange), 
                  na.rm = TRUE)
    } else {
      xmax <- max(abs(results$log2FoldChange), na.rm = TRUE)
    }
  } 
  if (is.null(ymax)) {
    pvals <- -log10(results$padj)
    pvals <- pvals[is.finite(pvals)]
    ymax <- max(pvals)
  } 
  
  # remove NA values from plot
  results <- filter(results, !is.na(padj))
  
  # make plot
  p <- 
    ggplot(results, 
           aes(x = log2FoldChange, y = -log10(padj), color = padj < 0.05)) +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
    geom_vline(aes(xintercept = 0)) +
    geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
    geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") +
    geom_point(show.legend = FALSE, stroke = 0, alpha = 0.4) +
    scale_color_manual(values = c("gray", "firebrick3")) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(x = expression("Fold change (log"[2]*")"), 
         y = expression(italic("P")*" value (-log"[10]*")")) +
    scale_x_continuous(expand = c(0, 0), limits = c(-xmax, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax))
  
  # add gene labels
  if (!is.null(label)) {
    df <- left_join(results, genes, by = "gene_id") %>%
      mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
      filter(gene_name %in% label)
    p <- p + geom_text_repel(data = df, aes(label = gene_name),
                             show.legend = FALSE)
  }
  return(p)
}

# - Spectral color palette ----------------------------------------------------
# from colorlover python package
make_spectral <- function(n = 100) {
  colors <- c(rgb(158/255, 1/255, 66/255, 1), 
              rgb(213/255, 62/255, 79/255, 1),
              rgb(244/255, 109/255, 67/255, 1), 
              rgb(253/255, 174/255, 97/255, 1), 
              rgb(254/255, 224/255, 139/255, 1),
              rgb(255/255, 255/255, 191/255, 1), 
              rgb(230/255, 245/255, 152/255, 1), 
              rgb(171/255, 221/255, 164/255, 1), 
              rgb(102/255, 194/255, 165/255, 1), 
              rgb(50/255, 136/255, 189/255, 1), 
              rgb(94/255, 79/255, 162/255, 1))
  colorRampPalette(colors)(n)
} 

# - Heatmap plot --------------------------------------------------------------
heatmap_plot <- function(counts, gene_ids = NULL, 
                         info = NULL, annotation = NULL,
                         max_cpm = 10, make_cpm = TRUE,
                         label_genes = FALSE,
                         cluster_genes = TRUE, cluster_samples = TRUE,
                         draw_tree = FALSE,
                         color_palette = make_spectral(), 
                         tree_scaling = 1) {
  if (make_cpm) { counts <- make_cpm(counts, log2 = TRUE) }
  if (!is.null(genes)) { 
    counts <- filter(counts, gene_id %in% gene_ids) 
  } else {
    gene_ids <- unique(counts$gene_id)
  }
  
  # make tidy
  df <- gather(counts, -gene_id, key = "Sample_ID", value = "counts")
  
  # cluster genes and samples
  if (cluster_genes) {
    gene_clust <- counts %>% as.data.frame() %>%
      column_to_rownames("gene_id") %>% dist() %>% hclust()
    gene_levels <- gene_clust$labels[gene_clust$order]
  } else {
    gene_levels <- genes
  }
  if (cluster_samples) {
    sample_clust <- hclust(dist(t(counts[, 2:ncol(counts)])))
    sample_levels <- sample_clust$labels[sample_clust$order]
  } else {
    sample_levels <- colnames(counts)[colnames(counts) != "gene_id"]
  }
  df <- df %>%
    mutate(gene_id = factor(gene_id, levels = gene_levels),
           Sample_ID = factor(Sample_ID, levels = sample_levels))
    
  # clip max at given max value
  df <- mutate(df, counts = ifelse(counts > max_cpm, max_cpm, counts))
  
  # plot
  plt <- ggplot(df, aes(x = Sample_ID, y = as.numeric(gene_id), fill = counts)) +
    geom_tile() +
    xlab(NULL) + ylab(NULL) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     color = "black"),
          axis.text.y = element_blank())
  
  if (label_genes) {
    plt <- plt + theme(axis.text.y = element_text(color = "black"))
  }
  if (!is.null(color_palette)) {
    plt <- plt + 
      scale_fill_gradientn(colors = color_palette,
                           name = expression(underline("CPM\n(log"[2]*")")))
  }
  
  # add annotation
  if (!is.null(info)) {
    anno <- left_join(df, info, by = "Sample_ID")
    anno <- anno %>% mutate(
      Sample_ID = factor(Sample_ID, levels = sample_levels)
      ) %>% as.data.frame()
  }
  
  # add annotation (optional)
  add_annotation_rect <- function(i, j) {
    annotate("rect", 
             xmin = which(levels(anno$Sample_ID) == samples[j]) - 0.5, 
             xmax = which(levels(anno$Sample_ID) == samples[j]) + 0.5, 
             ymin = max_yval + block_size * (i - 1), 
             ymax = max_yval + block_size * i,
             fill = colors[anno[anno$Sample_ID == samples[j], annotation[i]]]
    )
  }
  add_annotation_text <- function(i) {
    annotate("label", x = length(samples)/2 + 0.5, 
             y = max_yval + block_size * (i - 0.5), 
             label = annotation[i], hjust = 0.5, color = "gray90", 
             fill = "gray10", alpha = 0.5, label.size = NA)
  }
  if (!is.null(info) & !is.null(annotation)) {
    max_yval <- length(levels(df$gene_id))
    block_size <- max_yval * 0.05
    samples <- unique(df$Sample_ID)
    brewer_palettes <- c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2",
                         "Set1", "Set2", "Set3")
    palettes <- base::sample(brewer_palettes, length(annotation))
    for (i in 1:length(annotation)) {
      classes <- unique(anno[, annotation[i]])
      colors <- RColorBrewer::brewer.pal(n = length(classes), 
                                         name = palettes[i])
      names(colors) <- classes
      for (j in 1:length(samples)) {
        plt <- plt + add_annotation_rect(i, j)
      }
      plt <- plt + add_annotation_text(i)
    }
  } else if (!is.null(annotation) & is.null(info)) {
    message("Cannot add annotation to a plot without an info object")
  } else {
    block_size <- 0
  }
  
  # return plot
  if (draw_tree) {
    # draw tree
    sticks <- dendro_data(sample_clust)$segments
    scaling <- length(gene_ids) * 0.1
    for (i in 1:nrow(sticks)) {
      plt <- plt +
        annotate("segment",
                 x = sticks[i, "x"], xend = sticks[i, "xend"],
                 y = (sticks[i, "y"] * tree_scaling) + 
                   (block_size * length(annotation)) +
                   length(gene_ids) + 0.5,
                 yend = (sticks[i, "yend"] * tree_scaling) +
                   (block_size * length(annotation)) + 
                   length(gene_ids) + 0.5)
    }
  }
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
                             annotation = NULL, threshold = NULL,
                             cluster_samples = TRUE,
                             draw_tree = FALSE) {
  # grab data
  corr <- correlation(counts, genes)
  
  # cluster genes and samples
  if (cluster_samples) {
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
    
  }
  
  # plot
  plt <- ggplot(corr, aes(x = Sample_A, y = Sample_B)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 1),
          axis.title = element_blank(),
          plot.background = element_blank(),
          panel.background = element_blank())
  if (is.null(threshold)) {
    plt <- plt + geom_tile(aes(fill = corr)) +
      scale_fill_gradient2(low = "#ca562c", mid = "#f6edbd", high = "#008080",
                           midpoint = 0.9,
                           name = expression(underline("Correlation")))
  } else {
    label <- paste("Correlation >", threshold)
    plt <- plt + geom_tile(aes(fill = corr > threshold)) +
      scale_fill_manual(values = c("#273046", "#FAD510"),
                        name = substitute(underline(label))
      )
  }
  
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
    brewer_palettes <- c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2",
                        "Set1", "Set2", "Set3")
    palettes <- base::sample(brewer_palettes, length(annotation))
    for (i in 1:length(annotation)) {
      classes <- unique(anno[,annotation[i]])
      colors <- RColorBrewer::brewer.pal(n = length(classes), 
                                         name = palettes[i])
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
  if (draw_tree) {
    n_samples <- length(sample_clust$labels)
    tree_plot <- ggdendrogram(dendro_data(sample_clust)) + 
      theme_void() +
      theme(plot.margin = unit(c(0, 1.1-(0.008*n_samples), 
                                 0, 0.2-(0.008*n_samples)), "in"))
    return(cowplot::plot_grid(plotlist = list(tree_plot, plt), ncol = 1,
                              rel_heights = c(0.2, 0.8))
    )
  } else {
    return(plt)
  }
}

# - Expression vs. enrichment plot --------------------------------------------
ma_plot <- function(results, cpm, label = NULL,
                    ymax = NULL, ylim = NULL) {
  # combine enrichment and expression data
  cpm <- data.frame("gene_id" = cpm$gene_id, 
                    "expr" = rowMeans(select(cpm, -gene_id))
                    )
  df <- left_join(results, cpm, by = "gene_id") %>%
    filter(!is.na(padj))
  
  # set axis limits
  if (is.null(ymax)) {
    if (sum(df$padj < 0.05) > 0) {
      ymax <- max(abs(filter(df, padj < 0.05)$log2FoldChange))
    } else {
      ymax <- max(abs(df$log2FoldChange))
    }
  }
  if (is.null(ylim)) {
    ylim <- c(1, 6000)
  }
  
  # plot
  p <- ggplot(df, aes(x = expr, y = log2FoldChange, color = padj < 0.05)) +
    geom_hline(aes(yintercept = -log2(1.5)), linetype = "dashed") +
    geom_hline(aes(yintercept = log2(1.5)), linetype = "dashed") +
    geom_point(alpha = 0.4, stroke = 0, show.legend = FALSE) +
    scale_x_continuous(breaks = c(1, 16, 256, 4096), limits = ylim,
                       expand = c(0, 0), trans = "log2") +
    scale_y_continuous(limits = c(-ymax, ymax), expand = c(0, 0)) +
    scale_color_manual(values = c("gray50", "firebrick3")) +
    geom_hline(aes(yintercept = 0)) +
    theme_bw() +
    labs(y = expression("Fold change (log"[2]*")"),
         x = "Expression (FPM)") +
    theme(panel.grid = element_blank())
  
  # label genes
  if (!is.null(label)) {
    df <- left_join(df, genes, by = "gene_id") %>% 
      mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
      filter(gene_name %in% label)
    p <- p + geom_text_repel(data = df, aes(label = gene_name),
                             show.legend = FALSE)
  }
  return(p)
}
