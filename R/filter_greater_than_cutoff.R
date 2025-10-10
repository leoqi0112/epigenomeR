# Filter Rows by Non-Zero Count
# Post: Filter dataframe to retain only rows with non-zero values exceeding a specified threshold.
#       Counts the number of non-zero (positive) values per row and keeps rows where this count 
#       is greater than the cutoff.
# Parameter: df: A dataframe or matrix with numeric values.
#            cutoff_non_zero: Minimum number of non-zero values required per row, default to 10.
#                             Rows with MORE than this many non-zero values are kept.
# Output: Filtered dataframe containing only rows that exceed the non-zero count threshold.
filter_greater_than_cutoff <- function(df, cutoff_non_zero = 10) {
  zero_counts <- rowSums(df > 0, na.rm = TRUE)
  filtered_df <- df[zero_counts > cutoff_non_zero, ]
  return(filtered_df)
}