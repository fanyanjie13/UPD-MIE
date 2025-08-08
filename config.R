# Configuration for R scripts
WORKDIR <- "/data/yfan/upd_trio_routine"
SCRIPTS_DIR <- file.path(WORKDIR, "scripts")
INPUT_DIR <- file.path(WORKDIR, "input")
OUTPUT_DIR <- file.path(WORKDIR, "output")
PLOTS_DIR <- file.path(OUTPUT_DIR, "plots")
STATS_DIR <- file.path(OUTPUT_DIR, "stats")

# Create necessary directories if they don't exist
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(STATS_DIR, recursive = TRUE, showWarnings = FALSE)
