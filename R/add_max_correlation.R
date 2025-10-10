# Add Maximum Correlation and Best Cluster to Correlation Matrix
# Post: Augment correlation matrix by adding two columns identifying the maximum correlation value
#       and its corresponding cluster for each row (feature/region).
#       Finds the highest correlation per row and the cluster name associated with that peak correlation.
# Parameter: correlation_result: A correlation matrix with features as rows and clusters as columns.
#                                Typically output from correlation_matrix() function.
# Output: A matrix with two additional columns appended:
#         - max_correlation: The maximum correlation value for each row.
#         - best_cluster: The column name (cluster identifier) corresponding to the maximum correlation.
#         Original correlation values are preserved in the first columns.
add_max_correlation <- function(correlation_result) {
  max_values <- apply(correlation_result, 1, max, na.rm = TRUE)
  
  max_clusters <- colnames(correlation_result)[apply(correlation_result, 1, which.max)]
  result <- cbind(correlation_result, 
                  max_correlation = max_values,
                  best_cluster = max_clusters)
  
  return(result)
}