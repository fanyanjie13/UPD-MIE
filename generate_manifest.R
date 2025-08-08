#!/usr/bin/env Rscript

# Script to generate manifest.txt from batch file
# Usage example:
#   Rscript generate_manifest.R P500_trio_batch.txt manifest.txt /data/vcf
# Or from within the scripts directory:
#   Rscript generate_manifest.R ../input/P500_trio_batch.txt ../input/manifest.txt /data/vcf

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Set defaults and handle arguments
if (length(args) < 1) {
  stop("Usage: Rscript generate_manifest.R input_batch [output_dir] [vcf_base_dir]\n",
       "Example 1: Rscript generate_manifest.R P500_trio_batch.txt\n",
       "Example 2: Rscript generate_manifest.R P500_trio_batch.txt manifest.txt\n",
       "Example 3: Rscript generate_manifest.R P500_trio_batch.txt manifest.txt /storagelocal/newpipeline\n",
       "\nDefaults:\n",
       "  output_dir = current directory\n",
       "  vcf_base_dir = /storagelocal/newpipeline")
}

# Get input file
input_file <- args[1]

# Get or set output filename/path
if (length(args) >= 2) {
    # If second argument is a directory, use it as output directory
    if (dir.exists(args[2]) || grepl("/$", args[2])) {
        output_dir <- args[2]
        output_file <- file.path(output_dir, sprintf("manifest_%s.txt", tools::file_path_sans_ext(basename(input_file))))
    } else {
        # If second argument is a filename, use it directly
        output_file <- args[2]
        output_dir <- dirname(output_file)
    }
} else {
    # Default: create manifest file in current directory
    output_dir <- getwd()
    output_file <- file.path(output_dir, sprintf("manifest_%s.txt", tools::file_path_sans_ext(basename(input_file))))
}

# Get or set VCF base directory
vcf_base_dir <- if (length(args) >= 3) args[3] else "/storagelocal/newpipeline"

# Convert relative paths to absolute paths
if (!grepl("^/", output_dir)) {
    output_dir <- file.path(getwd(), output_dir)
    output_file <- file.path(getwd(), output_file)
}
if (!grepl("^/", input_file)) {
    input_file <- file.path(getwd(), input_file)
}

# Construct output filename from input filename
input_basename <- tools::file_path_sans_ext(basename(input_file))
output_file <- file.path(output_dir, sprintf("manifest_%s.txt", input_basename))

# Check if input file exists
if (!file.exists(input_file)) {
  stop("Input file does not exist: ", input_file)
}

# Create output directory if it doesn't exist
out_dir <- dirname(output_file)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Read the batch file
batch_data <- try(read.csv(input_file, header = FALSE, stringsAsFactors = FALSE))
if (inherits(batch_data, "try-error")) {
  stop("Error reading batch file. Ensure it's a comma-separated file with 7 columns.")
}
colnames(batch_data) <- c("proband", "father", "mother", "proband_proj", "father_proj", "mother_proj", "type")

# Function to construct VCF path and check existence
construct_vcf_path <- function(sample_id, project) {
  vcf_path <- file.path(vcf_base_dir, project, sample_id, paste0(sample_id, ".flt.vcf"))
  if (!file.exists(vcf_path)) {
    warning(sprintf("VCF file does not exist: %s", vcf_path))
  }
  return(vcf_path)
}

# Generate manifest data
manifest_data <- data.frame(
  probandID = batch_data$proband,
  fatherID = batch_data$father,
  motherID = batch_data$mother,
  probandvcf = sapply(1:nrow(batch_data), function(i) construct_vcf_path(batch_data$proband[i], batch_data$proband_proj[i])),
  fathervcf = sapply(1:nrow(batch_data), function(i) construct_vcf_path(batch_data$father[i], batch_data$father_proj[i])),
  mothervcf = sapply(1:nrow(batch_data), function(i) construct_vcf_path(batch_data$mother[i], batch_data$mother_proj[i]))
)

# Write manifest file with tab separation (no header)
write.table(manifest_data, file = output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Create Dell storage version of manifest data
dell_manifest_data <- manifest_data
dell_manifest_data$probandvcf <- file.path("/mnt/Storage7/yfan/UPD_trio/exchange", paste0(dell_manifest_data$probandID, ".flt.vcf"))
dell_manifest_data$fathervcf <- file.path("/mnt/Storage7/yfan/UPD_trio/exchange", paste0(dell_manifest_data$fatherID, ".flt.vcf"))
dell_manifest_data$mothervcf <- file.path("/mnt/Storage7/yfan/UPD_trio/exchange", paste0(dell_manifest_data$motherID, ".flt.vcf"))

# Write Dell storage manifest file (no header)
dell_output_file <- file.path(dirname(output_file), paste0("manifest_", tools::file_path_sans_ext(basename(input_file)), "_dell.txt"))
write.table(dell_manifest_data, file = dell_output_file, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Print summary
cat("\nManifest generation summary:\n")
cat(sprintf("Input batch file: %s\n", input_file))
cat(sprintf("Output manifest: %s\n", output_file))
cat(sprintf("Dell storage manifest: %s\n", dell_output_file))
cat(sprintf("VCF base directory: %s\n", vcf_base_dir))
cat(sprintf("Total trios processed: %d\n", nrow(manifest_data)))

# Check for any missing VCF files
missing_vcfs <- !sapply(c(as.character(manifest_data$probandvcf), 
                         as.character(manifest_data$fathervcf), 
                         as.character(manifest_data$mothervcf)), file.exists)
if (any(missing_vcfs)) {
  cat("\nWARNING: Some VCF files are missing. Please check the paths.\n")
}
