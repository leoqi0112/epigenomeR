# R/handlers/count_matrix_handler.R
# Handler function for count-matrix subcommand

handle_count_matrix <- function(args, verbose) {
  
  if (verbose) {
    cat("=== Count Matrix Generation ===\n")
    cat("Parsing arguments...\n")
  }
  
  # Parse BAM paths
  bam_paths <- trimws(strsplit(args$bam_paths, ",")[[1]])
  if (!all(file.exists(bam_paths))) {
    stop("Error: Not all BAM files exist. Check paths:\n", 
         paste(bam_paths[!file.exists(bam_paths)], collapse = "\n"))
  }
  if (verbose) {
    cat("BAM files:", length(bam_paths), "\n")
    cat("  -", paste(basename(bam_paths), collapse = "\n  - "), "\n")
  }
  
  # Parse regions (either integer or file paths)
  regions <- tryCatch({
    as.integer(args$regions)
  }, warning = function(w) {
    trimws(strsplit(args$regions, ",")[[1]])
  }, error = function(e) {
    trimws(strsplit(args$regions, ",")[[1]])
  })
  
  if (is.numeric(regions)) {
    if (verbose) cat("Using bin size:", regions, "\n")
  } else {
    if (!all(file.exists(regions))) {
      stop("Error: Not all region files exist. Check paths:\n", 
           paste(regions[!file.exists(regions)], collapse = "\n"))
    }
    if (verbose) {
      cat("Region files:", length(regions), "\n")
      cat("  -", paste(basename(regions), collapse = "\n  - "), "\n")
    }
  }
  
  # Parse transformations
  transformations <- NULL
  if (!is.null(args$transformations)) {
    transformations <- trimws(strsplit(args$transformations, ",")[[1]])
    if (verbose) {
      cat("Transformations:", paste(transformations, collapse = ", "), "\n")
    }
  }
  
  # Parse boolean flags
  apply_transformation <- truthy(args$apply_transformation)
  save_each_step <- truthy(args$save_each_step)
  do_qc <- truthy(args$do_qc)
  
  if (verbose) {
    cat("Apply transformation:", apply_transformation, "\n")
    cat("Save each step:", save_each_step, "\n")
    cat("QC filtering:", do_qc, "\n")
    if (do_qc) {
      cat("QC percentile threshold:", args$qc_filtered_percentile, "\n")
    }
    cat("\n")
  }
  
  # Call the main function
  count_matrix_function_with_qc(
    bam_path = bam_paths,
    regions = regions,
    save_dir = args$save_dir,
    libnorm_type = args$libnorm_type,
    apply_transformation = apply_transformation,
    transformations = transformations,
    save_each_step = save_each_step,
    datasetName_full = args$dataset_name,
    do_qc = do_qc,
    qc_filtered_percentile = args$qc_filtered_percentile
  )
  
  if (verbose) {
    cat("\nCount matrix generation completed successfully!\n")
    cat("Output saved to:", args$save_dir, "\n")
  }
}
