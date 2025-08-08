# combined_chr_bpi.R
# Summarize BPI values across all samples for each chromosome

# Load configuration
source("config.R")

# Load required libraries
library(dplyr)
library(readr)

# Read the BPI data (long format: id, chr, BPI_pct)
data <- read_csv(file.path(STATS_DIR, "other_chrs_data_bpi.csv"))

# Calculate summary statistics per chromosome (across all samples)
summary_table <- data %>%
  group_by(chr) %>%
  summarise(
    mean_ratio = mean(BPI_pct, na.rm = TRUE),
    sd_ratio = sd(BPI_pct, na.rm = TRUE),
    n = n(),
    se = sd(BPI_pct, na.rm = TRUE) / sqrt(n()),
    margin_of_error = qt(0.975, df = n() - 1) * se,
    lower_ci = mean(BPI_pct, na.rm = TRUE) - margin_of_error,
    upper_ci = mean(BPI_pct, na.rm = TRUE) + margin_of_error
  ) %>%
  ungroup()

# Save the summary table
write_csv(summary_table, file.path(STATS_DIR, "combined_chr_bpi_summary.csv"))

# Also save the per-sample, per-chromosome BPI table
bpi_each_sample <- data %>%
  select(chr, file_id = id, BPI_pct)
write_csv(bpi_each_sample, file.path(STATS_DIR, "combined_chr_bpi_eachSample.csv"))

cat("Summary table saved as", file.path(STATS_DIR, "combined_chr_bpi_summary.csv"), "\n")
cat("Per-sample BPI table saved as", file.path(STATS_DIR, "combined_chr_bpi_eachSample.csv"), "\n")