library(dplyr)
library(ggplot2)

# Define the input and output directories
folder_path <- "/data/yfan/upd_trio/temp/input"
output_dir <- "/data/yfan/upd_trio_routine/output/stats"

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Get a list of all files in the folder
file_list <- list.files(path = folder_path, full.names = TRUE)

# Initialize an empty list to store the data from each file
data_list <- list()

# Loop through each file in the file_list
for (file_path in file_list) {
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
    } else {
      message(paste("Skipping file:", file_path, "- Unable to locate both tables."))
    }
  } else {
    message(paste("Skipping file:", file_path, "- Unable to read file or missing required columns."))
  }
}
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


# --- Calculate Percentiles Across All Samples ---
library(dplyr)
library(tidyr)

# Convert 'chr' to an ordered factor
combined_data$chr <- factor(combined_data$chr, levels = c(1:22, "all_chrs"))

# Function to calculate percentiles for a given column name
calculate_percentiles <- function(data, col_prefix) {
  data %>%
    select(chr, starts_with(col_prefix)) %>%
    pivot_longer(cols = starts_with(col_prefix), names_to = "sample", values_to = "value") %>%
    group_by(chr) %>%
    summarise(
      median = median(value, na.rm = TRUE),
      p99 = quantile(value, 0.99, na.rm = TRUE),
      p95 = quantile(value, 0.95, na.rm = TRUE),
      p75 = quantile(value, 0.75, na.rm = TRUE),
      p25 = quantile(value, 0.25, na.rm = TRUE),
      p05 = quantile(value, 0.05, na.rm = TRUE),
      p01 = quantile(value, 0.01, na.rm = TRUE)
    )
}

# Calculate percentiles for each of the five values
percentile_table_bpi <- calculate_percentiles(combined_data, "BPI_pct_")
percentile_table_pathupd <- calculate_percentiles(combined_data, "pat_hUPD_pct_")
percentile_table_patiupd <- calculate_percentiles(combined_data, "pat_iUPD_pct_")
percentile_table_mathupd <- calculate_percentiles(combined_data, "mat_hUPD_pct_")
percentile_table_matiupd <- calculate_percentiles(combined_data, "mat_iUPD_pct_")

# --- Save the percentile_tables ---
write.csv(percentile_table_bpi, file = file.path(output_dir, "percentile_table_bpi.csv"), row.names = FALSE)
write.csv(percentile_table_pathupd, file = file.path(output_dir, "percentile_table_pat_hUPD.csv"), row.names = FALSE)
write.csv(percentile_table_patiupd, file = file.path(output_dir, "percentile_table_pat_iUPD.csv"), row.names = FALSE)
write.csv(percentile_table_mathupd, file = file.path(output_dir, "percentile_table_mat_hUPD.csv"), row.names = FALSE)
write.csv(percentile_table_matiupd, file = file.path(output_dir, "percentile_table_mat_iUPD.csv"), row.names = FALSE)

# --- Save the other_chrs_data files ---
# First, prepare the data in the correct format
other_chrs_data <- combined_data %>%
  select(chr, starts_with("BPI_pct_")) %>%
  pivot_longer(cols = starts_with("BPI_pct_"), 
               names_to = "id", 
               values_to = "BPI_pct") %>%
  mutate(id = sub("BPI_pct_", "", id))

# Save each metric's data
write.csv(other_chrs_data, file.path(output_dir, "other_chrs_data_bpi.csv"), row.names = FALSE)

other_chrs_data <- combined_data %>%
  select(chr, starts_with("pat_hUPD_pct_")) %>%
  pivot_longer(cols = starts_with("pat_hUPD_pct_"), 
               names_to = "id", 
               values_to = "pat_hUPD_pct") %>%
  mutate(id = sub("pat_hUPD_pct_", "", id))
write.csv(other_chrs_data, file.path(output_dir, "other_chrs_data_pat_hUPD.csv"), row.names = FALSE)

other_chrs_data <- combined_data %>%
  select(chr, starts_with("pat_iUPD_pct_")) %>%
  pivot_longer(cols = starts_with("pat_iUPD_pct_"), 
               names_to = "id", 
               values_to = "pat_iUPD_pct") %>%
  mutate(id = sub("pat_iUPD_pct_", "", id))
write.csv(other_chrs_data, file.path(output_dir, "other_chrs_data_pat_iUPD.csv"), row.names = FALSE)

other_chrs_data <- combined_data %>%
  select(chr, starts_with("mat_hUPD_pct_")) %>%
  pivot_longer(cols = starts_with("mat_hUPD_pct_"), 
               names_to = "id", 
               values_to = "mat_hUPD_pct") %>%
  mutate(id = sub("mat_hUPD_pct_", "", id))
write.csv(other_chrs_data, file.path(output_dir, "other_chrs_data_mat_hUPD.csv"), row.names = FALSE)

other_chrs_data <- combined_data %>%
  select(chr, starts_with("mat_iUPD_pct_")) %>%
  pivot_longer(cols = starts_with("mat_iUPD_pct_"), 
               names_to = "id", 
               values_to = "mat_iUPD_pct") %>%
  mutate(id = sub("mat_iUPD_pct_", "", id))
write.csv(other_chrs_data, file.path(output_dir, "other_chrs_data_matiupd.csv"), row.names = FALSE)