# Flow Cytometry Analysis Pipeline

A comprehensive R package for automated flow cytometry analysis with interactive visualization tools, specifically designed for congenic marker analysis and engraftment studies.

## Overview

This pipeline provides a complete workflow for flow cytometry data analysis, from FlowJo workspace import to publication-ready visualizations. It features interactive data exploration, automated congenic marker detection, statistical analysis, and advanced visualization tools including UMAP dimensionality reduction.

### Key Features

- **Interactive Data Analysis**: Menu-driven selection of cell populations, markers, and analysis parameters
- **Congenic Marker Analysis**: Automated detection and analysis of CD45.1, CD45.2, CD90.1, CD90.2, and CD45.1.2 markers
- **Engraftment Analysis**: Specialized tools for donor-recipient engraftment studies with statistical testing
- **Advanced Visualizations**: MFI heatmaps with statistics, paired comparisons, density plots, and UMAP analysis
- **Statistical Testing**: Built-in statistical analysis with multiple test options and significance annotations
- **Session Management**: Persistent plot sessions with save/export functionality
- **Data Cleaning**: Automated unstained sample detection and removal
- **Enhanced Interactive Menus**: Comprehensive navigation with back options throughout analysis

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
source("https://raw.githubusercontent.com/yourusername/your-repo/main/your-script.R")
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

# 2. Source the analysis functions
source("path/to/your/analysis_script.R")

# 3. Set up your gating set
gs <- setup_flowjo_workspace(
  xml_path = here("data", "your_workspace.wsp"),
  fcs_path = here("data", "your_fcs_files/")
)

# 4. Run interactive analysis with automatic congenic detection
congenics_results <- analyze_flow_data_auto(gs)

# 5. Clean up data (removes unstained samples, cleans column names)
congenics_results <- data_clean_custom(congenics_results)

# 6. Add experimental metadata
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)
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

#### `analyze_flow_data_auto(gs, add_congenics, congenic_candidates)`
Main analysis function with automatic congenic column addition.

**Parameters:**
- `gs`: GatingSet object
- `add_congenics`: Automatically add congenic marker column (default: TRUE)
- `congenic_candidates`: List of congenic markers to detect (default: `c("CD45.1", "CD45.2", "CD90.1", "CD90.2", "CD45.1.2")`)

**Features:**
- Interactive node selection with regex and leaf name support
- Full navigation with back buttons throughout the process
- Parent node management with hierarchical and congenic-specific options
- Interactive channel selection for MFI analysis
- Automatic congenic marker detection and classification

#### `analyze_flow_data(gs, node_selection, parent_selection, channels)`
Core analysis function with full interactive navigation.

**Features:**
- **Node Selection**: Interactive, regex pattern, or leaf name selection
- **Parent Selection**: Auto-hierarchical, common parent, individual, or congenic-specific
- **Channel Selection**: Specific markers, all channels, or compensated channels only
- **Back Navigation**: Navigate backwards through any step of the analysis
- **Preview & Confirmation**: Review selections before proceeding

### Data Processing and Cleaning

#### `data_clean_custom(data)`
Enhanced data cleaning with interactive unstained sample detection.

**Features:**
- **Unstained Sample Detection**: Automatically finds samples with "unstained", "no stain", etc.
- **Interactive Removal**: User choice to remove or keep detected samples
- **Column Cleaning**: Removes "$" symbols and cleans NodeShort formatting
- **Works with Lists**: Handles both single data frames and analysis result lists

#### `assign_metadata_menu(df)`
Interactive metadata assignment for experimental design.

**Features:**
- **Genotype Assignment**: WT/KO marker selection from detected congenics
- **Recipient Marking**: Optional recipient marker identification
- **Custom Metadata**: Add timepoint, batch, sex, or custom columns
- **Menu-Driven Interface**: Easy-to-use selection menus

### Visualization Functions

#### `create_paired_comparison_plots(data)`
Generate statistical comparison plots with interactive test selection.

**Features:**
- **Statistical Tests**: Paired/unpaired t-test, Wilcoxon signed-rank test
- **Faceting Options**: By tissue (GROUPNAME), by cell population (NodeShort), or no faceting
- **Automatic Layout**: Smart plot organization based on faceting choice
- **Statistical Annotations**: P-values and significance stars

#### `create_engraftment_plot(data, ...)`
Specialized engraftment analysis with comprehensive statistical options.

**Features:**
- **Interactive Statistics**: Choose from paired t-test, unpaired t-test, ANOVA
- **Comparison Types**: All pairwise comparisons or vs. control group
- **Normalization Options**: Interactive tissue selection for normalization
- **Statistical Annotations**: Proper bracket placement and p-value formatting

#### `create_mfi_heatmaps_interactive_enhanced(mfi_data)`
Advanced MFI heatmap creation with statistical testing.

**Features:**
- **Interactive Selection**: Congenics, markers, grouping options via menus
- **Statistical Testing**: Two-group and multi-group tests with post-hoc analysis
- **Scaling Methods**: Raw, Z-score, log, sqrt, percentile scaling
- **Significance Overlay**: P-values and stars directly on heatmaps
- **Export Functions**: Save statistical results to CSV

### UMAP Analysis

#### `analyze_flow_umap_enhanced(gs)`
Comprehensive UMAP analysis with enhanced interactive features.

**Features:**
- **Sample Exclusion**: Remove unstained controls, compensation samples
- **External Metadata**: Import additional sample information from CSV
- **Interactive Visualization**: Multiple plot types with session management
- **Plot Management**: Save, export, and organize multiple plots
- **Robust Error Handling**: Handles data transformation issues automatically

**Interactive Visualization Menu:**
- Basic UMAP plots
- Color by variable (markers or metadata)
- Facet by categorical variables
- Multi-marker heatmap overlays
- Density plots
- Session management (save/load plots)

### Statistical Analysis

#### Enhanced Statistical Testing Options

**Two-Group Tests:**
- Paired t-test (recommended for engraftment data)
- Unpaired t-test
- Mann-Whitney U test (non-parametric)
- Wilcoxon signed-rank test (non-parametric, paired)

**Multi-Group Tests:**
- One-way ANOVA with post-hoc comparisons
- Kruskal-Wallis test (non-parametric)
- Repeated measures ANOVA
- Friedman test (non-parametric, paired)

**Features:**
- Interactive test selection with recommendations
- Multiple comparison correction (Bonferroni, FDR)
- Effect size calculations
- Comprehensive results export

## Example Workflow

### Complete Analysis Pipeline

```r
# 1. Setup and Data Import
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

# Source the analysis functions
source("path/to/your/script.R")

# Import FlowJo workspace
gs <- setup_flowjo_workspace(
  xml_path = here("data", "Exp2_D7.wsp"),
  fcs_path = here("data", "FCS_D7/")
)

# 2. Interactive Data Analysis
# Run main analysis with automatic congenic detection
congenics_results <- analyze_flow_data_auto(gs)

# Run additional analyses for other markers
KLRG1_CD127_results <- analyze_flow_data(gs)
CD69_CD103_results <- analyze_flow_data(gs)

# 3. Data Cleaning and Preprocessing
# Clean all datasets (removes unstained, cleans names)
congenics_results <- data_clean_custom(congenics_results)
KLRG1_CD127_results <- data_clean_custom(KLRG1_CD127_results)
CD69_CD103_results <- data_clean_custom(CD69_CD103_results)

# 4. Metadata Assignment
# Add experimental design information
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)

# 5. Statistical Analysis and Visualization
# Create engraftment analysis
engraftment_plot <- create_engraftment_plot(
  congenics_results_with_genotype,
  interactive = TRUE,
  interactive_stats = TRUE
)
print(engraftment_plot$plot)

# Generate MFI heatmaps with statistics
mfi_heatmaps <- create_mfi_heatmaps_interactive_enhanced(congenics_results$mfi)

# Create paired comparison plots
KLRG1_plots <- create_paired_comparison_plots(KLRG1_CD127_results$counts)
CD69_plots <- create_paired_comparison_plots(CD69_CD103_results$counts)

# 6. UMAP Analysis
umap_results <- analyze_flow_umap_enhanced(gs)

# 7. Data Export
# Export processed data
write_csv(congenics_results$counts, here("out", "data", "congenics_counts.csv"))
write_csv(congenics_results$mfi, here("out", "data", "congenics_mfi.csv"))

# Export statistical results
stats_results <- export_stats_results(mfi_heatmaps$statistics, "mfi_stats.csv")

# Export UMAP data
export_umap_data(umap_results, "umap_analysis.csv")
export_all_session_plots("umap_plots", format = "png")
```

### Advanced Features

#### Interactive Node Selection
The pipeline provides multiple ways to select cell populations:

```r
# During analysis, you'll see options like:
# 1. Select from list (interactive)
# 2. Use regex pattern  
# 3. Use leaf names
# 4. Exit/Cancel

# Example regex patterns:
# - "CD8.*pos" to find CD8+ populations
# - "CD45\\.1" to find CD45.1 populations
# - "Live.*CD3" to find Live CD3+ populations
```

#### Parent Node Management
Sophisticated parent node selection with multiple strategies:

```r
# Available parent selection methods:
# 1. Auto (hierarchical parent) - automatically detects hierarchy
# 2. Same parent for all nodes - apply one parent to all selected nodes
# 3. Individual parent for each node - customize each node's parent
# 4. Congenic-specific - automatic detection for CD45.1, CD45.2, etc.
```

#### Channel Selection for MFI
Flexible marker selection for MFI analysis:

```r
# Options include:
# 1. Select specific channels (interactive list)
# 2. Use all channels
# 3. Use all compensated channels (recommended)
# 4. Skip MFI analysis
```

## Advanced Configuration

### Congenic Marker Detection
The pipeline automatically detects congenic markers in your data:

```r
# Default congenic candidates
congenic_candidates <- c("CD45.1", "CD45.2", "CD90.1", "CD90.2", "CD45.1.2")

# Custom congenic detection
results <- analyze_flow_data_auto(
  gs, 
  congenic_candidates = c("CD45.1", "CD45.2", "Thy1.1", "Thy1.2")
)
```

### Statistical Testing Configuration

```r
# Configure statistical tests for heatmaps
stats_config <- list(
  test_type = "t_test",           # or "wilcox_test", "anova", etc.
  alpha = 0.05,                   # significance level
  sig_display = "both",           # "p_values", "stars", "both", "none"
  multiple_correction = "bonferroni"
)
```

### UMAP Parameters

```r
# Customize UMAP parameters during analysis
# - Number of neighbors (default: 15)
# - Minimum distance (default: 0.1)
# - Data transformation: asinh (recommended), log10, sqrt, none
# - Maximum cells per sample (default: 5000)
```

## Output Files

The pipeline generates comprehensive output:

### Data Files
- `*_counts.csv`: Population counts and frequencies with congenic classification
- `*_mfi.csv`: Mean fluorescence intensity data
- `*_umap_data.csv`: UMAP coordinates with marker expression
- `*_stats_results.csv`: Statistical analysis results

### Plots
- `engraftment_*.png`: Statistical engraftment ratio plots
- `heatmap_*.png`: MFI heatmaps with statistical annotations
- `paired_comparison_*.png`: Statistical comparison plots
- `umap_*.png`: UMAP visualizations
- `multi_marker_*.png`: Multi-marker heatmap overlays

## Troubleshooting

### Common Issues

1. **"Node not found" errors**: 
   - Use the interactive node selection to see available populations
   - Check that your FlowJo workspace contains the expected gates

2. **File path issues**: 
   - Always use `here()` for file paths
   - Ensure your `.Rproj` file is in the project root directory

3. **Missing FCS files**: 
   - Verify that .fcs files are in the specified directory
   - Check that sample names match between workspace and files

4. **Memory issues with UMAP**: 
   - Reduce `max_cells_per_sample` parameter
   - Use sample exclusion to remove unnecessary samples

5. **Unstained sample detection**: 
   - The pipeline automatically detects samples with "unstained", "no stain", etc.
   - You can choose to remove or keep these samples during cleaning

### Data Requirements

- **FlowJo workspace**: Compatible .wsp files (FlowJo 10+)
- **FCS files**: Standard .fcs format
- **Compensation**: Compensated parameters should contain "Comp" in channel names
- **Sample naming**: Consistent naming between workspace and FCS files

### Performance Tips

- **Large datasets**: Use sample exclusion and cell subsampling for UMAP
- **Memory management**: Clean up large objects after analysis
- **Plot sessions**: Use UMAP session management to avoid recomputing plots

## Contributing

To contribute to this package:

1. Fork the repository
2. Create a feature branch
3. Make your changes with comprehensive testing
4. Update documentation
5. Submit a pull request

## Citation

If you use this pipeline in your research, please cite:

```
Flow Cytometry Analysis Pipeline with Interactive Visualization Tools. 2025.
GitHub: [repository URL]
```

## License

[Your license information here]

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: [your contact information]

---

**Note**: This pipeline provides a comprehensive, interactive workflow for flow cytometry analysis with particular strength in congenic marker analysis, statistical testing, and advanced visualization. The system includes robust error handling, comprehensive navigation options, and publication-ready output generation.
