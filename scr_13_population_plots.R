#!/usr/bin/env Rscript
# Create PCA and ADMIXTURE plots for both VCF-derived datasets.
# The script expects outputs from scr_12_pca_admixture.sh.

# Configuration
ANALYSIS_DIR <- "./population_analysis"
PLOT_DIR <- file.path(ANALYSIS_DIR, "plots")
DATASETS <- c("modern", "with_all_samples")

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required. Install it with install.packages('ggplot2').")
}

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

for (dataset in DATASETS) {
  prefix <- file.path(ANALYSIS_DIR, dataset)
  eigenvec_file <- paste0(prefix, ".eigenvec")
  optimal_k_file <- paste0(prefix, "_optimal_K.txt")

  if (!file.exists(eigenvec_file)) {
    stop("PCA file not found: ", eigenvec_file)
  }
  if (!file.exists(optimal_k_file)) {
    stop("Optimal K file not found: ", optimal_k_file)
  }

  pca <- read.table(eigenvec_file, header = FALSE, stringsAsFactors = FALSE)
  if (ncol(pca) < 4) stop("Unexpected PCA file format: ", eigenvec_file)
  names(pca)[1:4] <- c("Family", "Sample", "PC1", "PC2")

  pca_plot <- ggplot2::ggplot(pca, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(size = 2, color = "steelblue") +
    ggplot2::labs(
      title = paste("PCA:", dataset),
      x = "PC1",
      y = "PC2"
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave(
    file.path(PLOT_DIR, paste0(dataset, "_PCA.png")),
    pca_plot,
    width = 8,
    height = 6,
    dpi = 300
  )

  best_k <- as.integer(readLines(optimal_k_file, n = 1))
  q_file <- paste0(prefix, ".", best_k, ".Q")
  fam_file <- paste0(prefix, ".fam")
  if (!file.exists(q_file) || !file.exists(fam_file)) {
    stop("ADMIXTURE output not found for ", dataset, " at K=", best_k)
  }

  q <- read.table(q_file, header = FALSE)
  fam <- read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  if (nrow(q) != nrow(fam)) stop("Q and FAM row counts differ for ", dataset)
  names(q) <- paste0("Cluster_", seq_len(ncol(q)))
  q$Sample <- fam$V2

  q_long <- reshape(
    q,
    varying = names(q)[seq_len(ncol(q) - 1)],
    v.names = "Ancestry",
    timevar = "Cluster",
    times = names(q)[seq_len(ncol(q) - 1)],
    idvar = "Sample",
    direction = "long"
  )
  q_long$Sample <- factor(q_long$Sample, levels = q$Sample)

  admixture_plot <- ggplot2::ggplot(
    q_long,
    ggplot2::aes(x = Sample, y = Ancestry, fill = Cluster)
  ) +
    ggplot2::geom_col(width = 1) +
    ggplot2::labs(
      title = paste("ADMIXTURE:", dataset, "(K =", best_k, ")"),
      x = "Sample",
      y = "Ancestry proportion"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(PLOT_DIR, paste0(dataset, "_ADMIXTURE_K", best_k, ".png")),
    admixture_plot,
    width = 14,
    height = 6,
    dpi = 300
  )
}

message("PCA and ADMIXTURE plots saved to: ", PLOT_DIR)
