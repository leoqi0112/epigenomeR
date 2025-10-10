# Reduce Dataframe by Cluster Groups
# Post: Aggregate dataframe rows by cluster assignments using a specified summary function.
#       Matches regions between dataframe rownames and cluster assignments, groups rows by 
#       cluster labels, and applies summary function (e.g., mean, median) to all columns within each cluster.
# Parameter: dataframe: A dataframe or matrix where rownames represent genomic regions or features.
#            cluster_df: A named vector where names are region identifiers and values are cluster labels.
#                        Note: Code uses 'cluster_info' variable - ensure this is defined or rename parameter.
#            summary_func: Function to aggregate values within each cluster, default to mean.
#                          Can use median, sum, or any function accepting vectors.
# Output: A dataframe with clusters as rownames and summarized values for each column.
#         Only includes clusters/regions present in both input dataframe and cluster_df.
reduce_cluster_group <- function(dataframe, cluster_df, summary_func = mean) {    
  dataframe <- as.data.frame(dataframe)
  common_regions <- intersect(rownames(dataframe), names(cluster_info))
  
  if(length(common_regions) == 0) {
    stop("No matching regions between dataframe and cluster_info")
  }
  
  filtered_df <- dataframe[common_regions, ]
  filtered_cluster <- cluster_info[common_regions]
  
  cluster_summary <- filtered_df %>%
    mutate(cluster = filtered_cluster) %>%
    group_by(cluster) %>%
    summarise(across(everything(), ~ summary_func(.x, na.rm = TRUE)), .groups = 'drop') %>%
    tibble::column_to_rownames("cluster")
  
  return(cluster_summary)
}