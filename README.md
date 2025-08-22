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

```r
# Load libraries and setup
library(here)
# ... (load all required libraries)

# Source analysis functions
source("https://github.com/NathanKieswetter/R_Flow_helpers/blob/main/Congenics_flow_semi-automated.R")

# Setup folder structure (automatic)
# Folders will be created if they don't exist

# Import FlowJo workspace
gs <- setup_flowjo_workspace(
  xml_path = here("data", "Exp2_D7.wsp"),
  fcs_path = here("data", "FCS_D7/")
)

# Assess gating hierarchy
plot(gs, fontsize = 15, bool = TRUE)

# Run interactive analyses for different populations
congenics_results <- analyze_flow_data_auto(gs)
KLRG1_CD127_results <- analyze_flow_data_auto(gs)
CD69_CD103_results <- analyze_flow_data_auto(gs)

# Clean up data
congenics_results <- data_clean_custom(congenics_results)
KLRG1_CD127_results <- data_clean_custom(KLRG1_CD127_results)
CD69_CD103_results <- data_clean_custom(CD69_CD103_results)

# Add experimental metadata
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)

# Create engraftment analysis
plot_engraftment <- create_engraftment_plot(congenics_results_with_genotype)
print(plot_engraftment)

# Generate MFI heatmaps
mfi_heatmaps <- create_mfi_heatmaps_interactive_enhanced(congenics_results$mfi)

# Create paired comparison plots
KLRG1_CD127_plots <- create_paired_comparison_plots(KLRG1_CD127_results$counts)
CD69_CD103_plots <- create_paired_comparison_plots(CD69_CD103_results$counts)

# Advanced UMAP analysis
umap_results <- analyze_flow_umap_enhanced(gs)
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
[Your citation information here]
```

## License

[Your license information here]

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: [your contact information]

---

**Note**: This pipeline is specifically designed for congenic marker analysis in flow cytometry data. For other applications, some functions may require modification.
