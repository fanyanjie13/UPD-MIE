# Load required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)

# Load configuration
source("config.R")

# --- Load the percentile_tables (Base Distribution) ---
percentile_table_bpi <- read.csv(file.path(STATS_DIR, "percentile_table_bpi.csv"))
percentile_table_pathupd <- read.csv(file.path(STATS_DIR, "percentile_table_pat_hUPD.csv"))
percentile_table_patiupd <- read.csv(file.path(STATS_DIR, "percentile_table_pat_iUPD.csv"))
percentile_table_mathupd <- read.csv(file.path(STATS_DIR, "percentile_table_mat_hUPD.csv"))
percentile_table_matiupd <- read.csv(file.path(STATS_DIR, "percentile_table_mat_iUPD.csv"))

# --- Function to process homREFed.table files ---
process_homREFed_table <- function(file_path) {
  # Read all lines from the file
  all_lines <- readLines(file_path)
  
  # Find the start of the first table
  table1_start <- grep("^chr\\spat_hUPD\\spat_iUPD\\smat_hUPD\\smat_iUPD\\sBPI\\sinformative", all_lines)
  table2_start <- grep("^chr\\sUA_P\\sUI_P\\sUA_M\\sUI_M\\sBPI\\sinformative", all_lines)
  
  # Extract lines for the first table
  table1_lines <- all_lines[table1_start:(table2_start - 1)]
  
  # Read the first table
  table1 <- read.table(text = table1_lines, header = TRUE, sep = "", fill = TRUE, 
                      blank.lines.skip = TRUE, stringsAsFactors = FALSE)
  
  # Clean the data
  table1 <- table1[!is.na(table1$chr), ]
  table1 <- table1 %>% filter(!chr %in% c("Tally", "chr"))
  
  # Select required columns
  required_cols <- c("chr", "BPI", "informative", "pat_hUPD", "pat_iUPD", "mat_hUPD", "mat_iUPD")
  table1_subset <- table1[, required_cols]
  
  # Convert to numeric
  cols_to_convert <- c("BPI", "informative", "pat_hUPD", "pat_iUPD", "mat_hUPD", "mat_iUPD")
  for (col in cols_to_convert) {
    table1_subset[[col]] <- as.numeric(table1_subset[[col]])
  }
  
  # Calculate percentages
  table1_subset <- table1_subset %>%
    mutate(
      BPI_pct = ifelse(informative == 0, 0, (BPI / informative) * 100),
      pat_hUPD_pct = ifelse(informative == 0, 0, (pat_hUPD / informative) * 100),
      pat_iUPD_pct = ifelse(informative == 0, 0, (pat_iUPD / informative) * 100),
      mat_hUPD_pct = ifelse(informative == 0, 0, (mat_hUPD / informative) * 100),
      mat_iUPD_pct = ifelse(informative == 0, 0, (mat_iUPD / informative) * 100)
    )
  
  # Extract sample ID
  sample_id <- sub(".*(NGS_[0-9]+).*", "\\1", basename(file_path))
  
  # Add sample ID
  table1_subset$id <- sample_id
  
  return(table1_subset)
}

# --- Process all homREFed.table files and update other_chrs_data ---
update_reference_data <- function() {
  # Get list of homREFed.table files
  file_list <- list.files(path = input_folder, pattern = ".*\\.homREFed\\.table$", full.names = TRUE)
  
  # Print total number of files to process
  total_files <- length(file_list)
  cat(sprintf("Found %d files to process\n", total_files))
  
  # Process each file
  all_data <- do.call(rbind, lapply(seq_along(file_list), function(i) {
    file_path <- file_list[i]
    cat(sprintf("Processing file %d of %d: %s\n", i, total_files, basename(file_path)))
    result <- process_homREFed_table(file_path)
    cat("  Successfully processed\n")
    return(result)
  }))
  
  # Split data by metric and save to respective files
  write.csv(all_data[, c("id", "chr", "BPI_pct")], file.path(STATS_DIR, "other_chrs_data_bpi.csv"), row.names=FALSE)
  write.csv(all_data[, c("id", "chr", "pat_hUPD_pct")], file.path(STATS_DIR, "other_chrs_data_pat_hUPD.csv"), row.names=FALSE)
  write.csv(all_data[, c("id", "chr", "pat_iUPD_pct")], file.path(STATS_DIR, "other_chrs_data_pat_iUPD.csv"), row.names=FALSE)
  write.csv(all_data[, c("id", "chr", "mat_hUPD_pct")], file.path(STATS_DIR, "other_chrs_data_mat_hUPD.csv"), row.names=FALSE)
  write.csv(all_data[, c("id", "chr", "mat_iUPD_pct")], file.path(STATS_DIR, "other_chrs_data_mat_iUPD.csv"), row.names=FALSE)
  
  return(all_data)
}

# --- Load other_chrs_data (Aggregated Sample Data) ---
load_data <- function(value_name) {
  other_chrs_data <- read.csv(file.path(STATS_DIR, paste0("other_chrs_data_", value_name, ".csv")))
  return(other_chrs_data)
}

# --- Process input data ---
input_folder <- OUTPUT_DIR  # Using the path from config.R

# First update the reference data
cat("Processing homREFed.table files and updating reference data...\n")
update_reference_data()

# Then load the updated data
data_bpi <- load_data("bpi")
data_pathupd <- load_data("pat_hUPD")
data_patiupd <- load_data("pat_iUPD")
data_mathupd <- load_data("mat_hUPD")
data_matiupd <- load_data("mat_iUPD")

# --- Optimized Plotting ---
# --- Function to create a single plot (Modified for chr order and y-axis range) ---
create_single_plot <- function(percentile_table, other_chrs_data, sample_id, value_name) {
  
  # Filter data for the specific sample and chr 1-22
  sample_data <- other_chrs_data %>%
    filter(id == sample_id, chr != "all_chrs")
  
  # Convert 'chr' to an ordered factor to control plotting order
  sample_data$chr <- factor(sample_data$chr, levels = as.character(1:22))
  percentile_table$chr <- factor(percentile_table$chr, levels = as.character(1:22))
  
  ggplot(percentile_table %>% filter(chr != "all_chrs"), aes(x = chr)) +
    geom_boxplot(
      aes(
        ymin = p05,
        lower = p25,
        middle = median,
        upper = p75,
        ymax = p95
      ),
      stat = "identity",
      alpha = 0.7
    ) +
    geom_point(data = sample_data, aes(x = chr, y = .data[[value_name]]), color = "red", size = 3) +
    labs(
      title = paste0(sample_id, " - ", value_name),
      x = "Chromosome",
      y = paste(value_name)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
    ) +
    coord_cartesian(ylim = c(0, 100))  # Set y-axis range from 0 to 100
}

# --- Function to create a grid of plots for a sample ---
create_grid_plot <- function(sample_id) {
  p_bpi <- create_single_plot(percentile_table_bpi, data_bpi, sample_id, "BPI_pct")
  p_pathupd <- create_single_plot(percentile_table_pathupd, data_pathupd, sample_id, "pat_hUPD_pct")
  p_patiupd <- create_single_plot(percentile_table_patiupd, data_patiupd, sample_id, "pat_iUPD_pct")
  p_mathupd <- create_single_plot(percentile_table_mathupd, data_mathupd, sample_id, "mat_hUPD_pct")
  p_matiupd <- create_single_plot(percentile_table_matiupd, data_matiupd, sample_id, "mat_iUPD_pct")
  
  # Save individual plots
  ggsave(file.path(PLOTS_DIR, paste0(sample_id, "_BPI_pct.png")), p_bpi, width = 8, height = 6)
  ggsave(file.path(PLOTS_DIR, paste0(sample_id, "_pat_hUPD_pct.png")), p_pathupd, width = 8, height = 6)
  ggsave(file.path(PLOTS_DIR, paste0(sample_id, "_pat_iUPD_pct.png")), p_patiupd, width = 8, height = 6)
  ggsave(file.path(PLOTS_DIR, paste0(sample_id, "_mat_hUPD_pct.png")), p_mathupd, width = 8, height = 6)
  ggsave(file.path(PLOTS_DIR, paste0(sample_id, "_mat_iUPD_pct.png")), p_matiupd, width = 8, height = 6)
  
  # --- Create a dummy PDF device to suppress output ---  
  pdf(NULL)    
  
  # Arrange the plots in a grid
  combined_plot <- grid.arrange(
    p_bpi, 
    arrangeGrob(p_pathupd, p_patiupd, p_mathupd, p_matiupd, ncol = 2, nrow = 2), # Arrange the smaller plots in a 2x2 grid
    ncol = 2, 
    widths = c(2, 2) # Adjust the widths to make the BPI plot larger
  )
  
  # --- Close the dummy PDF device ---  
  dev.off()   
  
  # --- Save the combined plot ---
  output_filename <- file.path(PLOTS_DIR, paste0(sample_id, "_grid_plot.png"))
  ggsave(output_filename, combined_plot, width = 16, height = 6)
  
  cat("Plot saved to:", output_filename, "\n")
}

# --- Get a list of all sample IDs ---
file_list <- list.files(path = input_folder, pattern = ".*\\.homREFed\\.table$", full.names = TRUE)
sample_ids <- sub(".*(NGS_[0-9]+).*", "\\1", file_list)

# Print total number of samples to process
total_samples <- length(sample_ids)
cat(sprintf("\nGenerating plots for %d samples...\n", total_samples))

# --- Loop through all sample IDs and generate grid plots ---
for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  cat(sprintf("\nProcessing sample %d of %d: %s\n", i, total_samples, sample_id))
  create_grid_plot(sample_id)
}

cat(sprintf("\nAll plots have been generated in the %s directory\n", PLOTS_DIR))

