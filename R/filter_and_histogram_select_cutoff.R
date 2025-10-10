# Filter Correlation Results by Quantile Threshold with Histogram Visualization
# Post: Filter correlation matrix to retain only rows exceeding a specified quantile threshold
#       of maximum correlation values. Optionally generates and saves a histogram showing the 
#       distribution of correlations with the threshold marked. Prints filtering statistics.
# Parameter: correlation_result_max: Matrix with 'max_correlation' column, typically output from 
#                                     add_max_correlation() function.
#            quantile_threshold: Quantile value (0-1) for filtering, default to 0.75.
#                                Regions with max_correlation > this quantile are retained.
#            save_plot: Logical indicating whether to save histogram plot, default to FALSE.
#            plot_path: File path for saving plot (PNG format). Required if save_plot = TRUE.
# Output: Filtered matrix containing only rows where max_correlation exceeds the quantile threshold.
#         Prints quantile value, region counts, and percentage retained to console.
#         If save_plot = TRUE, saves histogram with threshold line as PNG file.
filter_and_histogram_select_cutoff <- function(correlation_result_max, quantile_threshold = 0.75, save_plot = FALSE, plot_path = NULL) {
  
  max_cor_values <- as.numeric(correlation_result_max[, "max_correlation"])
  threshold <- quantile(max_cor_values, quantile_threshold, na.rm = TRUE)
  cat("Quantile", quantile_threshold*100, "%:", threshold, "\n")
  
  if(save_plot) {
    png(plot_path, width = 800, height = 600)
    hist(max_cor_values, 
         main = paste("Distribution of Maximum Correlations"),
         xlab = "Max Correlation", 
         ylab = "Frequency",
         bins = 1000,
         #breaks = 30,
         col = "lightblue",
         border = "black")
    abline(v = threshold, col = "red", lty = 2, lwd = 2)
    legend("topright", paste(quantile_threshold*100, "% quantile =", round(threshold, 3)), 
           col = "red", lty = 2)
    dev.off()
    cat("Plot saved as:", plot_path, "\n") 
  }
  
  filtered_data <- correlation_result_max[max_cor_values > threshold, ]
  
  cat("Original regions:", nrow(correlation_result_max), "\n")
  cat("Filtered regions (>", quantile_threshold*100, "%):", nrow(filtered_data), "\n")
  cat("Percentage retained:", round(nrow(filtered_data)/nrow(correlation_result_max)*100, 2), "%\n")
  
  return(filtered_data)
}