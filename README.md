# UPD Analysis Pipeline

Automated pipeline for UPD (Uniparental Disomy) analysis of trio exome data using UPDio, with statistical analysis and visualization.
## Prerequisites
- R (with packages: dplyr, ggplot2, gridExtra)
- Perl with VCFtools
- bcftools, tabix, bgzip
- UPDio toolkit

## Directory Structure

```
project_root/
├── input/    # VCF files and manifest
├── output/   # Results and plots
└── scripts/  # Analysis scripts
```
## Setup Instructions

1. Create the directory structure:
    ```bash
    mkdir -p ./project_root/{scripts,input,output/{plots,stats},temp,tools}
    ```
2. Install and configure UPDio and all required tools. Make sure all paths in your configuration files are set relative to your project directory for portability.
3. Prepare a manifest file (`input/manifest.txt`) with tab-separated fields:
    ```
    probandID    fatherID    motherID    probandvcf               fathervcf               mothervcf
    NGS_001      NGS_002     NGS_003     ./input/proband.vcf      ./input/father.vcf      ./input/mother.vcf
    ```
## Pipeline Components

### Scripts
- `updio_trioWES_new.sh` - Main UPDio analysis
- `generate_base_distribution_plot.r` - Generate base distribution plots
- `generate_percentile_rank.R` - Calculate percentile rankings
- `higlight_on_base_plot.R` - Create highlight plots comparing samples
- `combined_chr_bpi.R` - Generate combined chromosome BPI analysis

### Configuration Files
- `config.R` - Required for R analysis scripts. Contains paths for data input/output and must be in the same directory as the R scripts.
- `config.ps1` - (Optional) For PowerShell automation. Contains environment variables and paths.
- `setup_workspace.ps1` - (Optional) Script to help set up the directory structure.
## Running the Pipeline

1. Navigate to the scripts directory:
    ```bash
    cd ./scripts
    ```
2. Run UPDio analysis:
    ```bash
    bash updio_trioWES_new.sh
    ```
3. Run R analysis pipeline:
    ```bash
    Rscript generate_base_distribution_plot.r
    Rscript generate_percentile_rank.R
    Rscript higlight_on_base_plot.R
    Rscript combined_chr_bpi.R
    ```
## Output Files

### Plots Directory (`output/plots/`)
- `*_BPI_pct.png` - BPI percentage plots
- `*_pat_hUPD_pct.png` - Paternal hUPD percentage plots
- `*_pat_iUPD_pct.png` - Paternal iUPD percentage plots
- `*_mat_hUPD_pct.png` - Maternal hUPD percentage plots
- `*_mat_iUPD_pct.png` - Maternal iUPD percentage plots
- `*_grid_plot.png` - Combined grid plots

### Stats Directory (`output/stats/`)
- `other_chrs_data_*.csv` - Chromosome-specific statistics
- `all_chrs_data_*.csv` - Genome-wide statistics
- `combined_chr_bpi_*.csv` - Combined BPI analysis results
## Troubleshooting

- Ensure all required R packages are installed:
    ```R
    install.packages(c("dplyr", "ggplot2", "gridExtra", "reshape2"))
    ```
- Verify all input VCF files exist and are accessible.
- Make sure all paths in configuration files are correct and relative to the project root for portability.

## Notes

- All paths in manifest.txt should be relative or absolute as appropriate for your environment.
- Input VCF files will be automatically compressed and indexed if needed.
- Temporary files are cleaned up after processing.
- The pipeline generates both individual and combined plots for each analysis type.
UPD Analysis Pipeline
Automated pipeline for UPD (Uniparental Disomy) analysis of trio exome data using UPDio, with statistical analysis and visualization.

