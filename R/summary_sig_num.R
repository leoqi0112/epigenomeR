# differential -4
# sample_names should exactly in order (case_num then control_num)
# Post: Summarize the number of significant differentially accessible regions across column clusters using limma-voom or t-test methods, providing counts of total, upregulated, and downregulated regions for each cluster and threshold combination.
# Parameter: sample_names: Vector of sample names in exact order (first case_num samples for condition1, then control_num samples for condition2)
#            case_num: Number of replicates in condition 1 (minimum 2)
#            control_num: Number of replicates in condition 2 (minimum 2)
#            col_cluster_file: Path to TSV file containing column cluster assignments with 'label' and 'feature' columns
#            wgc_file_path: Vector of paths to feather files containing count matrices for each sample
#            sig_result_dir: Output directory where summary tables will be saved
#            l2fc_thres: Log2 fold change threshold for significance calling (default: 0.5)
#            mean_per_thres_list: Vector of mean expression percentile thresholds for filtering (default: c(0.25))
#            fdr_thres_list: Vector of FDR thresholds for significance testing (default: c(0.25))
#            normalization_factor: Scaling factor for CPM normalization (default: 1E6)
#            lowess_span: Span parameter for voom lowess fitting (default: 0.5)
#            method_name_list: Vector of statistical methods to use ("limma" or "ttest", default: c("limma"))
# Output: Saves TSV summary tables showing counts of significant, upregulated, and downregulated regions for each column cluster and parameter combination
summary_sig_num <- function(sample_names, case_num, control_num, wgc_file_path, sig_result_dir, col_cluster_file = NULL, l2fc_thres = 0.5, mean_per_thres_list = c(0.25), fdr_thres_list = c(0.25), normalization_factor = 1E6, lowess_span = 0.5, method_name_list = c("limma")) {
  if(case_num < 2 || control_num < 2){
    stop("Each condition must have at least 2 replicates.")
  }
  if(length(sample_names) != (case_num+control_num)){
    stop("Length of sample_names must equal case_num+control_num.")
  }

  # load libraries
  suppressPackageStartupMessages({
    library(arrow)
    library(tibble)
    library(glue)
    library(edgeR)
    library(matrixStats)
    library(matrixTests)
    library(limma)
    library(dplyr)
  })

  group1 <- sample_names[1:case_num]
  group2 <- sample_names[(case_num+1):(case_num+control_num)]
  conditions <- c(rep("condition1", case_num), rep("condition2", control_num))
  coldata <- data.frame("condition"=conditions, row.names = sample_names)
  group <- factor(coldata$condition)
  condition_levels <- levels(group)
  mm <- model.matrix(~0 + group)

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

  col_label_list <- unique(col_cluster$label)

  wgc_list <- lapply(wgc_file_path, function(f) column_to_rownames(read_feather(f), var="pos"))
  names(wgc_list) <- sample_names

  dir.create(sig_result_dir, recursive = TRUE, showWarnings = FALSE)

  for(fdr_thres in fdr_thres_list){
    for(mean_per_thres in mean_per_thres_list){

      fdr_sig_summary_table <- as.data.frame(matrix(0, nrow=length(col_label_list), ncol=3))
      colnames(fdr_sig_summary_table) <- c("sig","up","down")
      rownames(fdr_sig_summary_table) <- col_label_list

      for(method_name in method_name_list){
        for(col_label in col_label_list){

          target_pair_select_list <- col_cluster[col_cluster$label==col_label,"feature"]

          # replicates
          group1_mat <- sapply(group1, function(s) rowSums(wgc_list[[s]][, target_pair_select_list, drop=FALSE]))
          group2_mat <- sapply(group2, function(s) rowSums(wgc_list[[s]][, target_pair_select_list, drop=FALSE]))
          tmp_combine_orig <- data.frame(group1_mat, group2_mat)
          colnames(tmp_combine_orig) <- c(group1, group2)

          # >=1 group non zero
          group1_nonzero <- rowSums(tmp_combine_orig[, group1] == 0) == 0
          group2_nonzero <- rowSums(tmp_combine_orig[, group2] == 0) == 0
          tmp_combine <- tmp_combine_orig[group1_nonzero | group2_nonzero, ]

          # edgeR + voom
          d0 <- DGEList(tmp_combine)
          d <- calcNormFactors(d0)
          lib.size <- d$samples$lib.size
          y <- voom(d, mm, span=lowess_span, plot = FALSE)

          # filter
          mean_log2_cpm <- rowMeans(y$E) + log2(exp(mean(log(lib.size + 1)))) - log2(normalization_factor)
          tmp_combine_filtered <- tmp_combine[mean_log2_cpm >= quantile(mean_log2_cpm, mean_per_thres), ]

          d0_f <- DGEList(tmp_combine_filtered)
          d_f <- calcNormFactors(d0_f)
          y_f <- voom(d_f, mm, span=lowess_span, plot=FALSE)

          # differential analysis
          if(method_name=="limma"){
            fit <- lmFit(y_f, mm)
            contr <- makeContrasts(
              contrasts = paste0("group", condition_levels[2], " - group", condition_levels[1]),
              levels = colnames(coef(fit))
            )
            tmp <- contrasts.fit(fit, contr)
            tmp <- eBayes(tmp)
            top.table <- topTable(tmp, sort.by = "P", n = Inf)
          } else if(method_name=="ttest"){
            eff_libsize <- d_f$samples$lib.size * d_f$samples$norm.factors
            E_pre <- sweep((tmp_combine_filtered + 0.5), 2, eff_libsize, "/") * normalization_factor
            E_log <- log2(E_pre)
            g1_mean <- rowMeans(E_pre[, group1, drop=FALSE])
            g2_mean <- rowMeans(E_pre[, group2, drop=FALSE])
            log_fc <- log2((g2_mean + 1)/(g1_mean + 1))
            t_res <- row_t_welch(E_log[, group2, drop=FALSE], E_log[, group1, drop=FALSE])
            top.table <- data.frame(logFC=log_fc[rownames(t_res)], P.Value=t_res$pvalue)
            top.table <- top.table[complete.cases(top.table), ]
            top.table$adj.P.Val <- p.adjust(top.table$P.Value, method="BH")
          }

          top.table_sig <- top.table[(top.table$adj.P.Val<fdr_thres)&(abs(top.table$logFC)>l2fc_thres),]
          up_sig_num <- nrow(top.table[(top.table$adj.P.Val<fdr_thres)&(top.table$logFC>0),])
          down_sig_num <- nrow(top.table[(top.table$adj.P.Val<fdr_thres)&(top.table$logFC<=0),])

          fdr_sig_summary_table[col_label,"sig"] <- nrow(top.table_sig)
          fdr_sig_summary_table[col_label,"up"] <- up_sig_num
          fdr_sig_summary_table[col_label,"down"] <- down_sig_num

          message(glue("sig thres: {fdr_thres}, filter: {mean_per_thres}, {method_name}, column cluster: {col_label}, sig regions: {nrow(top.table_sig)}"))
        }
      }

      fdr_sig_summary_table <- tibble::rownames_to_column(fdr_sig_summary_table, var="column_cluster")
      out_file <- glue("summary_sig_num_{mean_per_thres}_FDR-{fdr_thres}_log2FC-{l2fc_thres}.tsv")
      write.table(fdr_sig_summary_table, file.path(sig_result_dir, out_file), sep="\t", quote=FALSE, row.names=FALSE)
    }
  }
}
