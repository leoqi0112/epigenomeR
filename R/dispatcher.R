suppressPackageStartupMessages({
  if (!requireNamespace("argparse", quietly = TRUE))
    install.packages("argparse", repos = "https://cloud.r-project.org")
  library(argparse)
})

# Source all function modules
source("R/count_matrix_function_with_qc.R")

truthy <- function(x) { tolower(as.character(x)) %in% c("t","true","1","yes","y") }
# ============================================================================
# Main parser setup
# ============================================================================
parser <- ArgumentParser(prog = "dispatcher.R", description = "Hi-Plex CUT&Tag analysis suite")

subparsers <- parser$add_subparsers(dest = "command", help = "Available commands")

# Global arguments
parser$add_argument("--threads", type = "integer", default = NA,
                    help = "Number of threads (overrides SLURM_CPUS_PER_TASK if set).")
parser$add_argument("--helpers", type = "character", default = NULL,
                    help = "Optional helpers R file to source() (qc(), apply_transformations(), etc.).")
parser$add_argument("--verbose", default = "TRUE",
                    help = "Verbose output (TRUE/FALSE)")

# ============================================================================
# Subcommand:count_matrix_function_with_qc
# build count matrix from bam file with pre-specified genomic regions
# ============================================================================
p_count <- subparsers$add_parser("count-matrix", 
                                 help = "Build count matrix from BAMs with optional QC & transforms")
p_count$add_argument("--bam_paths", required = TRUE, 
                     help = "Comma-separated BAM file paths")
p_count$add_argument("--regions", required = TRUE, 
                     help = "Integer bin size OR comma-separated region files (.bed/.tsv/.txt/.csv)")
p_count$add_argument("--save_dir", required = TRUE, 
                     help = "Output directory")
p_count$add_argument("--libnorm_type", default = "libnorm")
p_count$add_argument("--apply_transformation", default = "TRUE")
p_count$add_argument("--transformations", default = NULL, 
                     help = "Comma-separated list")
p_count$add_argument("--save_each_step", default = "TRUE")
p_count$add_argument("--dataset_name", default = NULL)
p_count$add_argument("--do_qc", default = "FALSE")
p_count$add_argument("--qc_filtered_percentile", type = "double", default = 0.25)

# ============================================================================
# Parse and dispatch
# ============================================================================
args <- parser$parse_args()

# Set threads if specified
if (!is.na(args$threads)) {
  Sys.setenv(SLURM_CPUS_PER_TASK = args$threads)
}

# Source optional helpers file
if (!is.null(args$helpers) && file.exists(args$helpers)) {
  source(args$helpers)
}

verbose <- truthy(args$verbose)

# ============================================================================
# Dispatch to appropriate command handler
# ============================================================================

if (is.null(args$command)) {
  parser$print_help()
  quit(status = 1)
}

tryCatch({
  
  if (args$command == "count-matrix") {
    source("R/handlers/count_matrix_handler.R")
    handle_count_matrix(args, verbose)
  } else {
    stop(paste("Unknown command:", args$command))
  }

}, error = function(e) {
  cat("Error:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})