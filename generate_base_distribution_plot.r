library(dplyr)
library(ggplot2)

# Load configuration
source("config.R")

# Define the folder containing the input files
folder_path <- file.path(WORKDIR, "output")

# Create output directories if they don't exist
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STATS_DIR, recursive = TRUE, showWarnings = FALSE)

# Get a list of all files in the folder
file_list <- list.files(path = folder_path, pattern = ".*\\.homREFed\\.table$", full.names = TRUE)

# Print total number of files to process
total_files <- length(file_list)
cat(sprintf("Found %d files to process\n", total_files))

# Initialize an empty list to store the data from each file
data_list <- list()

# Loop through each file in the file_list
for (i in seq_along(file_list)) {
  file_path <- file_list[i]
  
  # Print progress
  cat(sprintf("Processing file %d of %d: %s\n", i, total_files, basename(file_path)))
  
  # Read all lines from the file
  all_lines <- tryCatch(
    {
      readLines(file_path)
    },
    error = function(e) {
      message(paste("Error reading file:", file_path))
      message(conditionMessage(e))
      NULL
    }
  )
  
  # Check if the file was read successfully
  if (!is.null(all_lines)) {
    # Find the start of the first table (line with "chr pat_hUPD pat_iUPD mat_hUPD mat_iUPD BPI informative")
    table1_start <- grep("^chr\\spat_hUPD\\spat_iUPD\\smat_hUPD\\smat_iUPD\\sBPI\\sinformative", all_lines)
    
    # Find the start of the second table (line with "chr UA_P UI_P UA_M UI_M BPI informative")
    table2_start <- grep("^chr\\sUA_P\\sUI_P\\sUA_M\\sUI_M\\sBPI\\sinformative", all_lines)
    
    # Check if both tables were found
    if (length(table1_start) > 0 && length(table2_start) > 0) {
      # Extract lines for the first table (excluding the header line of the second table)
      table1_lines <- all_lines[table1_start:(table2_start - 1)]
      
      # Read the first table
      table1 <- read.table(text = table1_lines, header = TRUE, sep = "", fill = TRUE, blank.lines.skip = TRUE, stringsAsFactors = FALSE)
      
      # Remove rows with NA in chr
      table1 <- table1[!is.na(table1$chr), ]
      
      # Remove rows with "Tally" or "chr" in the chr column
      table1 <- table1 %>%
        filter(!chr %in% c("Tally", "chr"))
      
      # Select only the required columns
      required_cols <- c("chr", "BPI", "informative", "pat_hUPD", "pat_iUPD", "mat_hUPD", "mat_iUPD")
      table1_subset <- table1[, required_cols]
      
      # Convert relevant columns to numeric
      cols_to_convert <- c("BPI", "informative", "pat_hUPD", "pat_iUPD", "mat_hUPD", "mat_iUPD")
      for (col in cols_to_convert) {
        table1_subset[[col]] <- as.numeric(table1_subset[[col]])
      }
      
      data_list[[file_path]] <- table1_subset
      cat("  Successfully processed\n")
    } else {
      message(paste("Skipping file:", file_path, "- Unable to locate both tables."))
    }
  } else {
    message(paste("Skipping file:", file_path, "- Unable to read file or missing required columns."))
  }
}

# Print summary of processed files
cat(sprintf("\nProcessing complete. Successfully processed %d out of %d files.\n", 
            length(data_list), total_files))

# --- Calculate percentages and prepare for merging ---
library(dplyr)

# Create an empty list to store the processed data frames
combined_data_list <- list()

# Iterate through the data_list, calculate percentages, and prepare for merging
for (file_name in names(data_list)) {
  df <- data_list[[file_name]]
  
  # Calculate the five percentage values
  df <- df %>%
    mutate(
      BPI_pct = ifelse(informative == 0, 0, (BPI / informative) * 100),
      pat_hUPD_pct = ifelse(informative == 0, 0, (pat_hUPD / informative) * 100),
      pat_iUPD_pct = ifelse(informative == 0, 0, (pat_iUPD / informative) * 100),
      mat_hUPD_pct = ifelse(informative == 0, 0, (mat_hUPD / informative) * 100),
      mat_iUPD_pct = ifelse(informative == 0, 0, (mat_iUPD / informative) * 100)
    )
  
  # Remove the original columns
  df <- df %>%
    select(-BPI, -informative, -pat_hUPD, -pat_iUPD, -mat_hUPD, -mat_iUPD)
  
  # Extract ID from the file name (NGS_*)
  id <- sub(".*(NGS_[0-9]+).*", "\\1", file_name)
  
  # Rename the percentage columns using the extracted ID
  colnames(df) <- c("chr", paste0(c("BPI", "pat_hUPD", "pat_iUPD", "mat_hUPD", "mat_iUPD"), "_pct_", id))
  
  # Append the processed df to combined_data_list
  combined_data_list <- append(combined_data_list, list(df))
}

# Merge all data frames by 'chr'
combined_data <- Reduce(function(x, y) merge(x, y, by = "chr", all = TRUE), combined_data_list)

# Replace NAs with 0 - should be done AFTER conversion to numeric
combined_data[is.na(combined_data)] <- 0

# --- Reshaping ---
library(reshape2)

# Function to reshape data for a given column prefix
reshape_data <- function(data, col_prefix) {
  melted_data <- data %>%
    select(chr, starts_with(col_prefix)) %>%
    melt(id.vars = "chr", variable.name = "id_file", value.name = col_prefix)
  
  # Remove prefix and "pct_" from id_file for cleaner sample IDs
  melted_data$id <- sub(paste0(col_prefix, "_"), "", melted_data$id_file)
  
  return(melted_data)
}

# Reshape data for each of the five values
melted_data_bpi <- reshape_data(combined_data, "BPI_pct")
melted_data_pathupd <- reshape_data(combined_data, "pat_hUPD_pct")
melted_data_patiupd <- reshape_data(combined_data, "pat_iUPD_pct")
melted_data_mathupd <- reshape_data(combined_data, "mat_hUPD_pct")
melted_data_matiupd <- reshape_data(combined_data, "mat_iUPD_pct")

# --- Separate data for all_chrs and other chromosomes ---
library(dplyr)

# Function to separate data for all_chrs and other chromosomes
separate_data <- function(data) {
  all_chrs_data <- data %>% filter(chr == "all_chrs")
  other_chrs_data <- data %>% filter(chr != "all_chrs")
  return(list(all_chrs_data = all_chrs_data, other_chrs_data = other_chrs_data))
}

# Separate data for each of the five values
separated_data_bpi <- separate_data(melted_data_bpi)
separated_data_pathupd <- separate_data(melted_data_pathupd)
separated_data_patiupd <- separate_data(melted_data_patiupd)
separated_data_mathupd <- separate_data(melted_data_mathupd)
separated_data_matiupd <- separate_data(melted_data_matiupd)

# --- Save "other_chrs_data" and "all_chrs_data" ---
write.csv(separated_data_bpi$other_chrs_data, file = file.path(STATS_DIR, "other_chrs_data_bpi.csv"), row.names = FALSE)
write.csv(separated_data_bpi$all_chrs_data, file = file.path(STATS_DIR, "all_chrs_data_bpi.csv"), row.names = FALSE)
write.csv(separated_data_pathupd$other_chrs_data, file = file.path(STATS_DIR, "other_chrs_data_pat_hUPD.csv"), row.names = FALSE)
write.csv(separated_data_pathupd$all_chrs_data, file = file.path(STATS_DIR, "all_chrs_data_pat_hUPD.csv"), row.names = FALSE)
write.csv(separated_data_patiupd$other_chrs_data, file = file.path(STATS_DIR, "other_chrs_data_pat_iUPD.csv"), row.names = FALSE)
write.csv(separated_data_patiupd$all_chrs_data, file = file.path(STATS_DIR, "all_chrs_data_pat_iUPD.csv"), row.names = FALSE)
write.csv(separated_data_mathupd$other_chrs_data, file = file.path(STATS_DIR, "other_chrs_data_mat_hUPD.csv"), row.names = FALSE)
write.csv(separated_data_mathupd$all_chrs_data, file = file.path(STATS_DIR, "all_chrs_data_mat_hUPD.csv"), row.names = FALSE)
write.csv(separated_data_matiupd$other_chrs_data, file = file.path(STATS_DIR, "other_chrs_data_mat_iUPD.csv"), row.names = FALSE)
write.csv(separated_data_matiupd$all_chrs_data, file = file.path(STATS_DIR, "all_chrs_data_mat_iUPD.csv"), row.names = FALSE)