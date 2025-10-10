# Generate Feature Label Assignment Table with Priority-Based Classification
# Post: Create a comprehensive feature labeling table by filtering zero-sum regions and assigning
#       labels based on a priority hierarchy: row_cluster_table (highest) > add_high_cor_regions >
#       filtered_df > Background (default). Labels are assigned sequentially, with later assignments
#       overwriting earlier ones.
# Parameter: full_df: Complete dataframe with features as rownames and samples as columns.
#                     Rows with all zero values are removed.
#            filtered_df: Subset dataframe containing CRF-specific regions (used for 2nd priority labeling).
#            add_high_cor_regions: Dataframe with rownames as features and 'best_cluster' column
#                                  indicating cluster assignments based on correlation (3rd priority).
#            row_cluster_table: Dataframe with 'feature' and 'label' columns for predefined cluster
#                               assignments (highest priority, overwrites all other labels).
#            out_path: Optional file path to save result as feather format. Creates directory if needed.
#                      If NULL, no file is saved.
# Output: A dataframe with two columns: 'feature' (region names) and 'label' (assigned cluster/category).
#         Labels follow priority: row_cluster_table > add_high_cor_regions > CRF_specific > Background.
#         Optionally saves result to feather file if out_path is specified.
generate_result_df <- function(full_df, filtered_df, add_high_cor_regions, row_cluster_table, out_path = NULL) {
  library(arrow)
  # filter row with all 0
  row_sums <- rowSums(full_df, na.rm = TRUE)
  full_df_filtered <- full_df[row_sums != 0, ]

  if (nrow(full_df_filtered) == 0) {
    warning("All regions are zero")
    return(data.frame())
  }

  result_df <- data.frame(
    feature = rownames(full_df_filtered),
    label = "Background",
    stringsAsFactors = FALSE
  )

  result_df$label[result_df$feature %in% rownames(filtered_df)] <- "CRF_specific"

  matched_features <- result_df$feature %in% rownames(add_high_cor_regions)
  result_df$label[matched_features] <- add_high_cor_regions[result_df$feature[matched_features], "best_cluster"]

  matched_features_info <- result_df$feature %in% row_cluster_table$feature
  matching_indices <- match(result_df$feature[matched_features_info], row_cluster_table$feature)
  result_df$label[matched_features_info] <- as.character(row_cluster_table$label[matching_indices])

  # result_df[result_df$feature == "chr22_20320001_20320800", ]
  if (!is.null(out_path)) {
    out_dir <- dirname(out_path)
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    # Check file extension and save accordingly
    if (grepl("\\.feather$", out_path, ignore.case = TRUE)) {
      write_feather(result_df, out_path)
    } else if (grepl("\\.tsv$", out_path, ignore.case = TRUE)) {
      write.table(result_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    } else if (grepl("\\.csv$", out_path, ignore.case = TRUE)) {
      write.csv(result_df, out_path, row.names = FALSE)
    } else {
      # Default to TSV
      write.table(result_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    }

    cat("Result saved to:", out_path, "\n")
  }
  return(result_df)
}
