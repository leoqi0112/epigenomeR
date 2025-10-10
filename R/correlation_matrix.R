# Compute Correlation Matrix Between Two Dataframes
# Post: Calculate pairwise correlations between rows of two dataframes across their common samples.
#       Filters both dataframes to retain only shared column names (samples), then computes 
#       correlation matrix where each element represents the correlation between a row from df 
#       and a row from cluster_df.
# Parameter: df: A dataframe or matrix where rows are features/regions and columns are samples.
#            cluster_df: A dataframe or matrix where rows are cluster features and columns are samples.
#                        Must share at least one column name with df.
# Output: A correlation matrix with dimensions [nrow(df) x nrow(cluster_df)].
#         Rows correspond to features from df, columns correspond to features from cluster_df.
#         Uses complete observations only (removes NA values pairwise).
correlation_matrix <- function(df, cluster_df) {
  common_samples <- intersect(colnames(df), colnames(cluster_df))
  
  if(length(common_samples) == 0) {
    stop("No common samples between df and cluster_df")
  }
  
  df_filtered <- df[, common_samples]
  cluster_filtered <- cluster_df[, common_samples]
  
  cor_matrix <- cor(t(df_filtered), t(cluster_filtered), use = "complete.obs")
  
  return(cor_matrix)
}