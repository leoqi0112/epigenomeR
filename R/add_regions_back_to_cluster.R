# Add Non-Informative Regions Back to Clusters Based on Correlation
# Post: Assign cluster labels to regions excluded from informative set by correlating them with 
#       existing cluster signatures. Filters non-informative regions by non-zero count, computes 
#       correlations with cluster-aggregated profiles, identifies best-matching clusters, and 
#       generates a comprehensive feature-label table with priority-based assignments.
# Parameter: nozero_file_path: Path to feather file containing all non-zero regions with 'pos' column.
#            transformed_file_path: Path to feather file with transformed count data (e.g., normalized, log-transformed).
#            informative_path: Path to feather file containing informative/significant regions used for clustering.
#            cluster_path: Path to TSV/TXT file with 'feature' and 'label' columns defining cluster assignments.
#            out_path: File path to save final feature-label table as feather format.
#            save_plot: Logical indicating whether to save correlation distribution histogram, default to FALSE.
#            plot_path: File path for saving histogram plot (PNG format). Required if save_plot = TRUE.
#            cutoff_non_zero: Minimum number of non-zero samples required per region, default to 10.
#            quantile_threshold: Quantile threshold (0-1) for filtering high-correlation regions, default to 0.75.
# Output: A dataframe with 'feature' and 'label' columns where labels follow priority hierarchy:
#         cluster_path assignments > correlation-based assignments > CRF_specific > Background.
#         Also saves result to out_path as feather file and optionally saves correlation histogram.
add_regions_back_to_cluster <- function(nozero_file_path, transformed_file_path, informative_path, cluster_path, out_path, save_plot = FALSE, plot_path = NULL, cutoff_non_zero = 10, quantile_threshold = 0.75) {
  # load libraries
  suppressPackageStartupMessages({
    library(arrow)
    library(dplyr)
    library(readr)
  })   
  
  informative_regions <- read_feather(informative_path)
  informative_regions <- informative_regions %>% tibble::column_to_rownames("pos")  
  informative_regions <- as.matrix(informative_regions) 
  transformed_regions <- read_feather(transformed_file_path)
  transformed_regions <- transformed_regions %>% tibble::column_to_rownames("pos")  
  transformed_regions <- as.matrix(transformed_regions)
  row_cluster_table <- read.table(cluster_path, header = TRUE)
  head(row_cluster_table)
  row_cluster_vector <- setNames(row_cluster_table$label, row_cluster_table$feature)
  row_cluster_vector
  
  df <- read_feather(nozero_file_path) 
  if ("pos" %in% colnames(df)) {
    df <- df %>% tibble::column_to_rownames("pos")
  }
  mat_regions <- rownames(informative_regions)
  df_regions <- rownames(df)
  reduced_df <- df[setdiff(rownames(df), rownames(informative_regions)), ]
  
  # filter <10
  filtered_df <- filter_greater_than_cutoff(reduced_df, cutoff_non_zero = cutoff_non_zero)
  filtered_transformed_df <- transformed_regions[rownames(transformed_regions) %in% rownames(filtered_df), ] # row: 387931
  
  cluster_sample_matrix <- reduce_cluster_group(informative_regions, row_cluster_vector) 
  
  correlation_result <- correlation_matrix(df = filtered_transformed_df, cluster_df = cluster_sample_matrix)
  
  correlation_result_max <- add_max_correlation(correlation_result)
  
  high_cor_regions <- filter_and_histogram_select_cutoff(correlation_result_max = correlation_result_max, quantile_threshold = quantile_threshold, save_plot = save_plot, plot_path = plot_path)
  
  output <- generate_result_df(full_df = transformed_regions, filtered_df = filtered_df, add_high_cor_regions = high_cor_regions, row_cluster_table = row_cluster_table, out_path = out_path)
  
  return(output)
}