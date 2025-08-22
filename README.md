# Flow Cytometry Analysis Pipeline

A comprehensive R package for automated flow cytometry analysis with interactive visualization tools, specifically designed for congenic marker analysis and engraftment studies.

## Overview

This pipeline provides a complete workflow for flow cytometry data analysis, from FlowJo workspace import to publication-ready visualizations. It features interactive data exploration, automated congenic marker detection, statistical analysis, and advanced visualization tools including UMAP dimensionality reduction.

### Key Features

- **Interactive Data Analysis**: Menu-driven selection of cell populations, markers, and analysis parameters
- **Congenic Marker Analysis**: Automated detection and analysis of CD45.1, CD45.2, CD90.1, CD90.2, and CD45.1.2 markers
- **Engraftment Analysis**: Specialized tools for donor-recipient engraftment studies
- **Advanced Visualizations**: MFI heatmaps, paired comparisons, density plots, and UMAP analysis
- **Statistical Testing**: Built-in statistical analysis with multiple test options
- **Session Management**: Persistent plot sessions with save/export functionality

## Installation

### Prerequisites

This package requires R (≥ 4.0.0) and the following dependencies:

```r
# Install required packages
install.packages(c(
  "here", "tidyverse", "ggpubr", "rstatix", "scales", 
  "glue", "purrr", "patchwork", "viridis", "umap",
  "RColorBrewer", "circlize"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("CytoML", "flowCore", "flowWorkspace", "ComplexHeatmap"))
```

### Package Installation

```r
# Source the analysis functions directly from GitHub
source("https://github.com/NathanKieswetter/R_Flow_helpers/blob/main/Congenics_flow_semi-automated.R")
```

## Project Structure and Setup

### Folder Structure

The pipeline uses the `here` package for robust file path management. Your project should be organized as follows:

```
your-project/
├── your-project.Rproj          # RStudio project file (PLACE HERE)
├── data/                       # Raw data files
│   ├── your_workspace.wsp      # FlowJo workspace file
│   └── FCS_files/             # Directory containing .fcs files
├── scripts/                   # Analysis scripts
│   └── analysis.R             # Your analysis code
└── out/                       # Output files
    ├── plots/                 # Generated plots
    ├── data/                  # Processed data
    └── GatingSet/             # Saved gating sets
```

### Setting Up the HERE Package

**IMPORTANT**: Place your RStudio project file (`.Rproj`) in the root directory of your project. The `here` package uses this file to establish the project root for all relative file paths.

The pipeline automatically creates the required folder structure:

```r
library(here)

# Automatic folder creation
if(!dir.exists(here::here("data"))) {
  dir.create(here::here("data"), recursive = TRUE)
}
if(!dir.exists(here::here("out"))) {
  dir.create(here::here("out"), recursive = TRUE)
}
if(!dir.exists(here::here("scripts"))) {
  dir.create(here::here("scripts"), recursive = TRUE)
}
if(!dir.exists(here::here("out/GatingSet"))) {
  dir.create(here::here("out/GatingSet"), recursive = TRUE)
}
```

### File Path Usage

Always use `here()` for file paths to ensure reproducibility across different systems:

```r
# Correct usage
gs <- setup_flowjo_workspace(
  xml_path = here("data", "your_workspace.wsp"),
  fcs_path = here("data", "FCS_files")
)

# Avoid absolute paths like:
# xml_path = "/Users/yourname/Desktop/project/data/workspace.wsp"
```

## Quick Start Guide

### Basic Workflow

```r
# 1. Load required libraries
library(here)
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(scales)
library(glue)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(umap)
library(viridis)
library(purrr)
library(patchwork)

# 2. Source the analysis functions
source("https://github.com/NathanKieswetter/R_Flow_helpers/blob/main/Congenics_flow_semi-automated.R")

# 3. Set up your gating set
gs <- setup_flowjo_workspace(
  xml_path = here("data", "your_workspace.wsp"),
  fcs_path = here("data", "your_fcs_files/")
)

# 4. Visualize gating hierarchy
plot(gs, fontsize = 15, bool = TRUE)

# 5. Run interactive analysis
congenics_results <- analyze_flow_data_auto(gs)
```

## Core Functions

### Data Import and Setup

#### `setup_flowjo_workspace(xml_path, fcs_path, keywords)`
Imports FlowJo workspace and creates a GatingSet object.

**Parameters:**
- `xml_path`: Path to FlowJo workspace (.wsp) file
- `fcs_path`: Path to directory containing .fcs files
- `keywords`: Sample metadata keywords to extract (default: `c("$WELLID", "GROUPNAME")`)

### Interactive Analysis

#### `analyze_flow_data_auto(gs, ...)`
Main interactive analysis function with automatic congenic detection.

**Features:**
- Interactive node/population selection
- Automatic parent node detection
- Marker/channel selection with compensation detection
- Statistical analysis options
- Automatic congenic marker identification

### Data Processing

#### `data_clean_custom(data)`
Cleans and standardizes data output.

**Features:**
- Removes special characters from column names
- Standardizes node naming
- Handles unstained sample detection and removal
- Works with both single data frames and lists

#### `assign_metadata_menu(df)`
Interactive metadata assignment for experimental design.

**Features:**
- WT/KO marker assignment
- Recipient marker identification
- Custom metadata addition (timepoint, batch, sex, etc.)
- Genotype classification

### Visualization Functions

#### `create_engraftment_plot(data, ...)`
Creates publication-ready engraftment ratio plots.

**Features:**
- Interactive statistical method selection
- Normalization options (to spleen or other control tissues)
- Paired statistical comparisons
- Customizable styling

#### `create_mfi_heatmaps_interactive_enhanced(mfi_data)`
Advanced MFI heatmap creation with statistical testing.

**Features:**
- Interactive marker and congenic selection
- Multiple scaling methods (Z-score, log, percentile)
- Statistical testing integration
- Grouped or combined tissue analysis

#### `create_paired_comparison_plots(data)`
Creates paired comparison plots for population frequencies.

**Features:**
- Interactive statistical test selection
- Flexible faceting options
- Paired t-test and Wilcoxon test support
- Automatic significance annotation

#### `analyze_flow_umap_enhanced(gs)`
Comprehensive UMAP analysis with interactive visualization.

**Features:**
- Sample exclusion options
- Event counting and assessment
- External metadata import
- Interactive visualization menu
- Session-based plot management

## Example Workflow

### Basic Setup and Data Import

```r
# 1. Load required libraries
library(here)
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(scales)
library(glue)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(umap)
library(viridis)
library(purrr)
library(patchwork)

# 2. Source the analysis functions
source("https://github.com/NathanKieswetter/R_Flow_helpers/blob/main/Congenics_flow_semi-automated.R")

# 3. Setup folder structure (automatic - creates data/, out/, scripts/ directories)
# Folders will be created automatically when needed

# 4. Import FlowJo workspace
gs <- setup_flowjo_workspace(
  xml_path = here("data", "Exp2_D7.wsp"),
  fcs_path = here("data", "FCS_D7/")
)

# 5. Save GatingSet for future use
save_gatingset(gs, "experiment_D7_gatingset.rds")

# 6. Assess gating hierarchy
plot(gs, fontsize = 15, bool = TRUE)
```

### Interactive Data Analysis

```r
# 7. Run interactive analyses for different populations
# Each analysis will guide you through node selection, parent mapping, and channel selection
congenics_results <- analyze_flow_data_auto(gs)           # For congenic marker analysis
KLRG1_CD127_results <- analyze_flow_data_auto(gs)         # For activation markers
CD69_CD103_results <- analyze_flow_data_auto(gs)          # For tissue residence markers

# 8. Save analysis results
save_analysis_results(congenics_results, "congenics_results.rds")
save_analysis_results(KLRG1_CD127_results, "KLRG1_CD127_results.rds")
save_analysis_results(CD69_CD103_results, "CD69_CD103_results.rds")

# 9. Clean up data (removes special characters, handles unstained samples)
congenics_results <- data_clean_custom(congenics_results)
KLRG1_CD127_results <- data_clean_custom(KLRG1_CD127_results)
CD69_CD103_results <- data_clean_custom(CD69_CD103_results)
```

### Metadata Assignment and Experimental Design

```r
# 10. Add experimental metadata (interactive menus for WT/KO assignment, timepoints, etc.)
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)

# The metadata assignment includes:
# - WT/KO marker identification
# - Recipient marker assignment
# - Custom metadata (timepoint, batch, sex, etc.)
```

### Visualization and Analysis

```r
# 11. Create engraftment analysis with statistical testing
plot_engraftment <- create_engraftment_plot(congenics_results_with_genotype)
print(plot_engraftment)

# Save engraftment plot (plots automatically save to out/plots/)
if(!is.null(plot_engraftment$plot)) {
  ggsave(here("out", "plots", "engraftment_analysis.png"), 
         plot_engraftment$plot, width = 10, height = 6, dpi = 300)
}

# 12. Generate interactive MFI heatmaps with statistical testing
mfi_heatmaps <- create_mfi_heatmaps_interactive_enhanced(congenics_results$mfi)
# This function provides:
# - Interactive congenic and marker selection
# - Multiple scaling methods
# - Statistical testing options
# - Automatic export to out/plots/

# 13. Create paired comparison plots
KLRG1_CD127_plots <- create_paired_comparison_plots(KLRG1_CD127_results$counts)
CD69_CD103_plots <- create_paired_comparison_plots(CD69_CD103_results$counts)

# Access individual plots:
# KLRG1_CD127_plots[["Spleen"]]        # View specific tissue
# CD69_CD103_plots[["CD103+CD69+"]]    # View specific population

# 14. Advanced UMAP analysis (interactive with session management)
umap_results <- analyze_flow_umap_enhanced(gs)
# Features include:
# - Sample exclusion (unstained, controls)
# - Interactive visualization menu
# - Plot session management
# - Multiple export options
```

### Data Export and Session Management

```r
# 15. Export processed data
# Export count/frequency data
write_csv(congenics_results$counts, here("out", "data", "congenics_counts.csv"))
write_csv(congenics_results$mfi, here("out", "data", "congenics_mfi.csv"))

# Export UMAP data
export_umap_data(umap_results, "umap_analysis.csv")

# 16. Save complete session
session_name <- save_session(gs, congenics_results, "experiment_D7_session")

# 17. Manage saved files interactively
manage_saved_files()  # Interactive menu for loading/deleting saved data
```

### Loading Previous Sessions

```r
# Load a previously saved session
loaded_session <- load_session("experiment_D7_session")
gs <- loaded_session$gs
results <- loaded_session$results

# Or load individual components
gs <- load_gatingset("experiment_D7_gatingset.rds")
congenics_results <- load_analysis_results("congenics_results.rds")

# List available files
list_available_gatingsets()
```

### Advanced Workflow Options

```r
# 18. Statistical analysis export
if(!is.null(mfi_heatmaps) && exists("stats_results")) {
  export_stats_results(stats_results, "mfi_statistical_analysis.csv")
}

# 19. Continue UMAP visualization session
restart_visualization_session(umap_results)

# 20. Export all session plots at once
export_all_session_plots()  # Saves all plots to out/plots/

# 21. Create quick plots from saved results
quick_plot <- create_quick_umap_plot(umap_results, "GROUPNAME")
ggsave(here("out", "plots", "quick_tissue_plot.png"), quick_plot)
```

## Configuration Options

### Statistical Testing
The pipeline supports multiple statistical tests:
- Paired/unpaired t-tests
- Mann-Whitney U test
- Wilcoxon signed-rank test
- ANOVA with post-hoc comparisons
- Kruskal-Wallis test

### Visualization Options
- **Color scales**: Viridis, RColorBrewer palettes
- **Plot types**: Scatter, density, heatmap, paired comparisons
- **Export formats**: PNG, PDF, SVG, JPEG
- **Statistical annotations**: P-values, significance stars, custom formatting

### Data Transformation
- **Flow cytometry**: asinh transformation (recommended)
- **Alternative**: log10, square root, or no transformation
- **Scaling**: Z-score, percentile, or raw values

## Troubleshooting

### Common Issues

1. **"Node not found" errors**: Ensure your FlowJo workspace contains the expected populations
2. **File path issues**: Always use `here()` for file paths and ensure your `.Rproj` file is in the project root
3. **Missing FCS files**: Verify that .fcs files are in the specified directory and match workspace samples
4. **Memory issues with UMAP**: Reduce `max_cells_per_sample` parameter for large datasets

### Data Requirements

- **FlowJo workspace**: Compatible .wsp files (FlowJo 10+)
- **FCS files**: Standard .fcs format with proper channel naming
- **Compensation**: Compensated parameters should contain "Comp" in channel names
- **Sample naming**: Consistent naming between workspace and FCS files

## Output Files

The pipeline generates several types of output:

### Data Files
- `*_counts.csv`: Population counts and frequencies
- `*_mfi.csv`: Mean fluorescence intensity data
- `*_umap_data.csv`: UMAP coordinates with marker expression

### Plots
- `engraftment_*.png`: Engraftment ratio plots
- `heatmap_*.png`: MFI heatmaps
- `paired_comparison_*.png`: Statistical comparison plots
- `umap_*.png`: UMAP visualizations

### Session Data
- `gating_set.rds`: Saved GatingSet object
- `analysis_results.rds`: Complete analysis results

## Contributing

To contribute to this package:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## Citation

If you use this pipeline in your research, please cite:

```
Nathan Scott Kieswetter
```

## License

[Your license information here]

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: natescottkieswetter@gmail.com

---

**Note**: This pipeline is specifically designed for congenic marker analysis in flow cytometry data. For other applications, some functions may require modification.
