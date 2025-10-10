# Differential Analysis
# differential -1
# two condition, 4 columns
# Post: Perform differential analysis between two conditions using limma-voom pipeline, analyzing each column cluster separately to identify significantly different genomic regions.
# Parameter: sample_names: Vector of sample names in order (first rep1 samples for condition1, then rep2 samples for condition2)
#            rep1: Number of replicates in condition 1 (minimum 2) --- case_num
#            rep2: Number of replicates in condition 2 (minimum 2) --- control_num
#            col_cluster_file: Path to column cluster assignment TSV file
#            wgc_file_path: Vector of paths to count matrix feather files for each sample
#            sig_result_dir: Directory to save differential analysis results
#            pseudocount: Pseudocount for normalization (default: 0.5)
#            normalization_factor: Scaling factor for CPM calculation (default: 1E6)
#            pseudocount_for_log: Pseudocount for log transformation (default: 1)
#            lowess_span: Span parameter for voom lowess fitting (default: 0.5)
#            l2fc_thres: Log2 fold change threshold for significance (default: 0.5)
#            mean_per_thres_list: Vector of mean expression percentile thresholds for filtering (default: c(0.25))
#            fdr_thres_list: Vector of FDR thresholds for significance calling (default: c(0.25))
# Output: Saves differential analysis results, filtered count matrices, and significant regions for each cluster and threshold combination
limma_column_cluster_differential_regions <- function(sample_names, rep1, rep2, col_cluster_file = NULL, wgc_file_path, sig_result_dir, normalization_factor = 1E6, lowess_span = 0.5, l2fc_thres = 0.5, mean_per_thres_list = c(0.25), fdr_thres_list = c(0.25)) { # cluster

  if(rep1 < 2 || rep2 < 2){
    stop("Each condition must have at least 2 replicates.")
  }
  if(length(sample_names) != (rep1+rep2)){
    stop("Length of sample_names must equal rep1+rep2.")
  }

  # load libraries
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(latex2exp)
    library(edgeR)
    library(matrixStats)
    library(limma)
  })

  group1 <- sample_names[1:rep1]
  group2 <- sample_names[(rep1+1):(rep1+rep2)]
  conditions <- c(rep("condition1", rep1), rep("condition2", rep2))
  coldata <- data.frame("condition" = conditions, row.names = sample_names)
  group <- factor(coldata$condition)
  condition_levels <- levels(group)
  mm <- model.matrix(~0 + group)

  wgc_list <- lapply(wgc_file_path, function(f) column_to_rownames(read_feather(f), var = "pos"))
  names(wgc_list) <- sample_names

  # Check if all WGC files have matching column names
  all_colnames <- lapply(wgc_list, colnames)
  all_colnames_sorted <- lapply(all_colnames, sort)
  colnames_strings <- sapply(all_colnames_sorted, paste, collapse = "|")

  if (length(unique(colnames_strings)) > 1) {
    mismatch_info <- c()
    for (i in 1:(length(all_colnames) - 1)) {
      for (j in (i + 1):length(all_colnames)) {
        if (!identical(all_colnames_sorted[[i]], all_colnames_sorted[[j]])) {
          mismatch_info <- c(mismatch_info,
                             paste0("  File ", i, " (", sample_names[i], ") vs File ", j, " (", sample_names[j], ")"))
        }
      }
    }
    stop("Column names do not match across WGC files:\n", paste(unique(mismatch_info), collapse = "\n"))
  }

  # Handle column cluster file
  if (is.null(col_cluster_file)) {
    # Get all unique column names from all WGC files
    all_features <- unique(unlist(lapply(wgc_list, colnames)))
    # Create a data frame where each feature gets its own label
    col_cluster <- data.frame(
      feature = all_features,
      label = seq_along(all_features),
      stringsAsFactors = FALSE
    )
  } else {
    # Load column cluster file
    col_cluster <- read.table(col_cluster_file, header = TRUE, sep = "\t", row.names = NULL)
  }

  # Check intersection between WGC columns and col_cluster features
  wgc_cols <- colnames(wgc_list[[1]])
  cluster_features <- col_cluster$feature

  wgc_not_in_cluster <- setdiff(wgc_cols, cluster_features)
  cluster_not_in_wgc <- setdiff(cluster_features, wgc_cols)

  if (length(wgc_not_in_cluster) > 0) {
    warning("WGC files contain ", length(wgc_not_in_cluster), " column(s) not found in col_cluster file")
  }

  if (length(cluster_not_in_wgc) > 0) {
    warning("col_cluster file contains ", length(cluster_not_in_wgc), " feature(s) not found in WGC files")
  }

  dir.create(sig_result_dir, recursive = TRUE, showWarnings = FALSE)

  # Loop over column clusters
  col_cluster <- col_cluster[col_cluster$feature %in% wgc_cols, ]
  col_label_list <- unique(col_cluster$label)
  for (col_label in col_label_list) {
    target_pair_select_list <- col_cluster[col_cluster$label == col_label, "feature"]

    # Combine counts for this cluster
    tmp_combine_orig <- do.call(
      cbind,
      lapply(sample_names, function(s) rowSums(wgc_list[[s]][, target_pair_select_list, drop=FALSE]))
    )
    colnames(tmp_combine_orig) <- sample_names

    # Keep rows with at least one group meeting both conditions:
    # Count non-zero replicates per group for each region
    group1_nonzero_count <- rowSums(tmp_combine_orig[, group1] != 0)
    group2_nonzero_count <- rowSums(tmp_combine_orig[, group2] != 0)

    # Condition a: more than 50% non-zero in at least one group
    group1_has_majority <- group1_nonzero_count > (rep1 * 0.5)
    group2_has_majority <- group2_nonzero_count > (rep2 * 0.5)

    # Condition b: at least 2 non-zero in at least one group
    group1_has_min2 <- group1_nonzero_count >= 2
    group2_has_min2 <- group2_nonzero_count >= 2

    keep_rows <- (group1_has_majority & group1_has_min2) | (group2_has_majority & group2_has_min2)
    tmp_combine <- tmp_combine_orig[keep_rows, ]


    # EdgeR object
    d0 <- DGEList(tmp_combine)
    d <- calcNormFactors(d0)
    norm_info <- d$samples
    lib.size <- norm_info$lib.size
    effect_libsize <- norm_info$lib.size * norm_info$norm.factors

    # voom (mean-variance plot)
    voom_plot_filename <- glue("{sig_result_dir}/voom_plot_col_cluster-{col_label}.pdf")
    pdf(voom_plot_filename)
    y <- voom(d, mm, span = lowess_span, plot = TRUE)
    dev.off()

    mean_log2_cpm <- rowMeans(y$E) + log2(exp(mean(log(lib.size + 1)))) - log2(normalization_factor)

    # # Fit full model
    # fit <- lmFit(y, mm)
    # fitted_values <- fit$coefficients %*% t(mm)
    # residuals <- y$E - fitted_values
    # sqrt_res_std <- sqrt(apply(residuals, 1, sd))

    # # Save "all regions" with p-value/log2FC
    # contr <- makeContrasts(grouptreated - groupuntreated, levels = colnames(coef(fit)))
    # tmp <- contrasts.fit(fit, contr)
    # tmp <- eBayes(tmp)
    # top.table <- topTable(tmp, sort.by = "P", n = Inf)
    # top.table_all <- rownames_to_column(top.table, var = "pos")
    # write_feather(top.table_all, glue("{sig_result_dir}/all_regions_col_cluster-{col_label}.feather"))

    # Loop over mean threshold filters
    for (mean_per_thres in mean_per_thres_list) {
      tmp_combine_filtered <- tmp_combine[mean_log2_cpm >= quantile(mean_log2_cpm, mean_per_thres), ]

      # Save filtered count matrix
      # double filter
      filtered_matrix_filename <- glue("{sig_result_dir}/filtered_counts_col_cluster-{col_label}_rowmean-{mean_per_thres}.feather")
      write_feather(rownames_to_column(as.data.frame(tmp_combine_filtered), var = "pos"), filtered_matrix_filename)

      # Re-fit with filtered
      d0_filtered <- DGEList(tmp_combine_filtered)
      d_filtered <- calcNormFactors(d0_filtered)
      y_filtered <- voom(d_filtered, mm, span = lowess_span, plot = FALSE)
      fit <- lmFit(y_filtered, mm)
      contr <- makeContrasts(groupcondition2 - groupcondition1, levels = colnames(coef(fit)))
      tmp <- contrasts.fit(fit, contr)
      tmp <- eBayes(tmp)
      top.table <- topTable(tmp, sort.by = "P", n = Inf)

      # Save full table for this threshold
      # significant count matrix
      top.table_all <- rownames_to_column(top.table, var = "pos")
      write_feather(top.table_all, glue("{sig_result_dir}/all_regions_rowmean-{mean_per_thres}_col_cluster-{col_label}.feather"))

      # Save significant regions for each FDR threshold
      for (fdr_thres in fdr_thres_list) {
        top.table_sig <- top.table[(top.table$adj.P.Val < fdr_thres) & (abs(top.table$logFC) > l2fc_thres), ]
        top.table_sig_filename <- glue("{sig_result_dir}/sig_regions_rowmean-{mean_per_thres}_FDR-{fdr_thres}_col_cluster-{col_label}.tsv")
        write.table(top.table_sig, top.table_sig_filename, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
      }
    }
  }
}
