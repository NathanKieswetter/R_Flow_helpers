# Flow Cytometry Analysis Pipeline for Adoptive Transfer Experiments

A comprehensive R package for automated flow cytometry analysis with interactive visualization tools, specifically designed for congenic marker analysis and engraftment studies.

**NB**: _Full disclosure, these tools are VERY much still in development and have MANY bugs that I'm working through!_

## Overview

This pipeline provides a comprehensive workflow for flow cytometry data analysis, from importing data into the FlowJo workspace to creating publication-ready visualizations. It features interactive data exploration, automated detection of congenic markers, statistical analysis, and advanced visualization tools, including UMAP dimensionality reduction and clustering.

**NB**: This is a data analysis toolkit that requires **fully compensated and gated data**. 

### Key Features

- **Enhanced Interactive Setup**: Automated package installation, intelligent keyword detection, and project structure management
- **Congenic Marker Analysis**: Automated detection and analysis of CD45.1, CD45.2, CD90.1, CD90.2, and CD45.1.2 markers
- **Engraftment Analysis**: Specialized tools for donor-recipient engraftment studies with comprehensive statistical testing
- **Advanced Visualizations**: MFI heatmaps with statistics, paired comparisons, density plots, and UMAP analysis
- **Statistical Testing**: Built-in statistical analysis with multiple test options and significance annotations
- **Enhanced Session Management**: Persistent plot sessions with comprehensive save/export functionality
- **Data Cleaning**: Automated unstained sample detection and removal with interactive confirmation
- **UMAP Integration**: Advanced dimensionality reduction with clustering capabilities
- **Clustering Analysis**: FlowSOM, ConsensusClusterPlus, K-means, hierarchical, and DBSCAN clustering methods

## Installation

### Prerequisites

This package automatically installs and loads all required dependencies:

**CRAN Packages:**
```r
# Automatically installed via the pipeline
c("here", "tidyverse", "ggpubr", "rstatix", "scales", 
  "glue", "circlize", "RColorBrewer", "umap", "cluster", 
  "factoextra", "patchwork")
```

**Bioconductor Packages:**
```r
# Automatically installed via the pipeline
c("CytoML", "flowCore", "flowWorkspace", "CATALYST", 
  "FlowSOM", "ConsensusClusterPlus", "ComplexHeatmap")
```

### Package Installation

```r
# Source the analysis functions directly from GitHub
source("https://raw.githubusercontent.com/NathanKieswetter/R_Flow_helpers/main/Congenics_flow_semi-automated.R")
```

## Project Structure and Setup

### Enhanced Folder Structure

The pipeline automatically creates a comprehensive project structure:

```
your-project/
├── your-project.Rproj          # RStudio project file (PLACE HERE)
├── data/                       # Raw and processed data
│   ├── raw/                    # Original data files
│   ├── processed/              # Cleaned datasets
│   ├── metadata/               # Sample metadata files
│   ├── statistics/             # Statistical results
│   ├── exported_dataframes/    # Analysis outputs
│   └── umap_data/             # UMAP analysis results
├── scripts/                   # Analysis scripts
│   └── analysis.R             # Your analysis code
└── out/                       # Output files
    ├── plots/                 # Generated plots
    │   ├── paired_comparisons/ # Statistical comparison plots
    │   ├── engraftment/       # Engraftment analysis plots
    │   ├── heatmaps/         # MFI and other heatmaps
    │   │   └── mfi/          # MFI-specific heatmaps
    │   ├── umap/             # UMAP visualizations
    │   └── exploratory/      # Exploratory plots
    ├── tables/               # Summary tables
    ├── reports/              # Analysis reports
    ├── sessions/             # Session data
    └── GatingSet/            # Saved gating sets
```
### Setting Up the HERE Package

**IMPORTANT**: Place your RStudio project file (`.Rproj`) in the root directory of your project. The `here` package uses this file to establish the project root for all relative file paths.Always use `here()` for file paths to ensure reproducibility across different systems.

### Automatic Setup

```r
# Source the pipeline (includes automatic package installation)
source("https://raw.githubusercontent.com/NathanKieswetter/R_Flow_helpers/main/Congenics_flow_semi-automated.R")

# Create enhanced project structure
setup_enhanced_project_structure()
```

## Quick Start Guide

### Basic Workflow

```r
# 1. Source the analysis functions (auto-installs packages)
source("https://raw.githubusercontent.com/NathanKieswetter/R_Flow_helpers/main/Congenics_flow_semi-automated.R")

# 2. Create enhanced project structure
setup_enhanced_project_structure()

# 3. Set up your gating set with interactive keyword selection
gs <- setup_flowjo_workspace_interactive(
  xml_path = here("data", "your_workspace.wsp"),
  fcs_path = here("data", "your_fcs_files/")
)

# 4. Run enhanced interactive analysis with automatic congenic detection
congenics_results <- analyze_flow_data_auto_enhanced(gs)

# 5. Clean up data (removes unstained samples, cleans column names)
congenics_results <- data_clean_custom(congenics_results)

# 6. Add experimental metadata and define analysis factors
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)
```

## Core Functions

### Enhanced Data Import and Setup

#### `setup_flowjo_workspace_interactive(xml_path, fcs_path, keywords, sample_limit)`
Enhanced FlowJo workspace import with intelligent keyword detection and selection.

**Features:**
- **Interactive Keyword Selection**: Smart detection of biologically relevant keywords
- **Keyword Preview**: Shows sample distributions and metadata examples
- **Pattern-Based Detection**: Automatic identification of tissue types, treatments, timepoints
- **Validation**: Comprehensive data quality checks

#### `setup_enhanced_project_structure()`
Creates a comprehensive folder structure for an organized analysis workflow.

**Features:**
- **Automatic Folder Creation**: Creates all necessary directories
- **Organized Output Structure**: Separate folders for different analysis types
- **Session Management**: Dedicated folders for plot sessions and reports

### Enhanced Interactive Analysis

#### `analyze_flow_data_auto_enhanced(gs, add_congenics, define_factors)`
Main analysis function with enhanced workflow and factor definition.

**Features:**
- **Automatic Congenic Detection**: Intelligent identification of congenic markers
- **Factor Definition**: Interactive setup of tissue_factor and pairing_factor
- **Enhanced Navigation**: Full back navigation throughout the process
- **Data Validation**: Comprehensive error checking and user guidance
- **Automatic Saving**: Optional interactive saving of all analysis outputs

#### `define_analysis_factors(df)`
Interactive definition of standardized analysis factors for generalizable workflows.

**Features:**
- **Tissue Factor Selection**: Choose grouping variables (tissue, treatment, condition)
- **Pairing Factor Setup**: Define paired analysis variables (mouse ID, well position)
- **Biological Relevance Detection**: Smart suggestions based on data content
- **Pattern Extraction**: Extract factors from sample names using regex patterns

### Enhanced Data Processing and Cleaning

#### `data_clean_custom(data, auto_save)`
Advanced data cleaning with interactive sample exclusion and quality control.

**Features:**
- **Intelligent Unstained Detection**: Finds samples with "unstained", "no stain", etc.
- **Interactive Removal**: User choice with preview of samples to be removed
- **Batch Processing**: Applies cleaning decisions to all datasets consistently
- **Column Standardization**: Removes problematic characters and standardizes naming

#### `assign_metadata_menu_enhanced(df, include_factor_definition)`
Comprehensive metadata assignment with factor definition integration.

**Features:**
- **Factor Definition Integration**: Seamless integration with analysis factor setup
- **Enhanced Genotype Assignment**: Smart congenic marker detection and assignment
- **Custom Metadata Addition**: Timepoint, batch, sex, and custom column support
- **Data Validation**: Comprehensive checks for data consistency and completeness

### Advanced Visualization Functions

#### `create_mfi_heatmaps_interactive_enhanced(mfi_data, auto_save)`
State-of-the-art MFI heatmap creation with statistical testing and multiple display options.

**Features:**
- **Interactive Selection**: Comprehensive menus for congenics, markers, and grouping
- **Statistical Integration**: Built-in statistical testing with significance overlays
- **Multiple Scaling Methods**: Raw, Z-score, log, sqrt, percentile scaling options
- **Advanced Statistics**: Two-group and multi-group tests with post-hoc analysis
- **Export Functions**: Automated saving of statistical results and visualizations

#### `create_engraftment_plot(data, interactive, interactive_stats)`
Enhanced engraftment analysis with comprehensive statistical options and normalization.

**Features:**
- **Interactive Statistics**: Choose from paired t-test, unpaired t-test, ANOVA
- **Flexible Normalization**: Interactive tissue selection for ratio calculations
- **Comparison Types**: All pairwise or vs. control group comparisons
- **Enhanced Visualization**: Proper error bars, significance annotations, and formatting

#### `create_paired_comparison_plots_enhanced(data, tissue_col, pairing_col)`
Advanced paired comparison plots using standardized factor definitions.

**Features:**
- **Factor Integration**: Uses tissue_factor and pairing_factor for consistent analysis
- **Enhanced Statistical Testing**: Multiple test options with proper pairing detection
- **Flexible Faceting**: By tissue, by cell population, or no faceting
- **Automatic Layout**: Smart plot organization based on data structure

### UMAP Analysis and Clustering

#### `analyze_flow_umap_enhanced(gs, keywords)`
Comprehensive UMAP analysis with enhanced sample management and visualization.

**Features:**
- **Sample Exclusion**: Interactive removal of unstained controls and compensation samples
- **External Metadata Import**: CSV import with robust merging and validation
- **Enhanced Visualization**: Multiple plot types with session management
- **Data Transformation**: asinh, log10, sqrt, or no transformation options
- **Persistent Sessions**: Save and manage multiple plots with export functionality

#### `perform_clustering_analysis(umap_results, method)`
Advanced clustering analysis with multiple algorithms and comprehensive evaluation.

**Available Methods:**
- **FlowSOM**: Self-organizing maps optimized for flow cytometry
- **FlowSOM + ConsensusClusterPlus**: Most robust option with consensus validation
- **K-means**: Classical clustering with elbow method optimization
- **Hierarchical**: Dendrogram-based clustering with silhouette analysis
- **DBSCAN**: Density-based clustering for complex shapes

**Features:**
- **Interactive Method Selection**: Guided selection with method recommendations
- **Automatic Parameter Optimization**: Data-driven parameter selection
- **Cluster Characterization**: MFI analysis and signature marker identification
- **Differential Expression**: Statistical comparison between clusters
- **Comprehensive Export**: Results, plots, and summary reports

## Example Workflow

### Complete Analysis Pipeline

```r
# 1. Setup and Data Import
source("https://raw.githubusercontent.com/NathanKieswetter/R_Flow_helpers/main/Congenics_flow_semi-automated.R")
setup_enhanced_project_structure()

# 2. Import FlowJo workspace with interactive keyword selection
gs <- setup_flowjo_workspace_interactive(
  xml_path = here("data", "Exp2_D7.wsp"),
  fcs_path = here("data", "FCS_D7/")
)

# 3. Interactive Data Analysis with Enhanced Features
congenics_results <- analyze_flow_data_auto_enhanced(gs)
KLRG1_CD127_results <- analyze_flow_data_auto_enhanced(gs)
CD69_CD103_results <- analyze_flow_data_auto_enhanced(gs)

# 4. Enhanced Data Cleaning and Preprocessing
congenics_results <- data_clean_custom(congenics_results)
KLRG1_CD127_results <- data_clean_custom(KLRG1_CD127_results)
CD69_CD103_results <- data_clean_custom(CD69_CD103_results)

# 5. Metadata Assignment with Factor Definition
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)

# 6. Advanced Statistical Analysis and Visualization

# Enhanced engraftment analysis
plot_engraftment <- create_engraftment_plot(congenics_results_with_genotype)
print(plot_engraftment)

# Interactive MFI heatmaps with statistics
mfi_heatmaps <- create_mfi_heatmaps_interactive_enhanced(congenics_results$mfi)

# Enhanced paired comparison plots
KLRG1_plots <- create_paired_comparison_plots_enhanced(KLRG1_CD127_results$counts)
CD69_plots <- create_paired_comparison_plots_enhanced(CD69_CD103_results$counts)

# 7. UMAP Analysis with Clustering
umap_results <- analyze_flow_umap_enhanced(gs)

# 8. Clustering Analysis (Multiple Methods Available)
clustering_results <- perform_clustering_analysis(umap_results)

# 9. Advanced Data Export and Session Management
# All analysis functions include interactive saving options
# Session data is automatically managed for continued analysis
```

### Multiple Experiment Workflow

```r
# Setup multiple GatingSets
gs1 <- setup_flowjo_workspace_with_metadata(
  xml_path = here("data", "Exp1_D7.wsp"),
  fcs_path = here("data", "FCS_D7/"),
  add_metadata = TRUE
)

gs2 <- setup_flowjo_workspace_with_metadata(
  xml_path = here("data", "Exp2_D7.wsp"),
  fcs_path = here("data", "FCS_D7/"),
  add_metadata = TRUE
)

# Check compatibility and merge
compatibility <- check_gating_set_compatibility(gs1, gs2)
if(compatibility$compatible) {
  gs_merged <- merge_gating_sets_with_tracking(gs1, gs2, "Experiment1", "Experiment2")
}

# Continue with standard workflow
congenics_results <- analyze_flow_data_auto_enhanced(gs_merged)
```

## Advanced Configuration

### Enhanced Statistical Testing Configuration

```r
# Configure statistical tests for heatmaps
stats_config <- list(
  test_type = "t_test",           # or "wilcox_test", "anova", etc.
  alpha = 0.05,                   # significance level
  sig_display = "both",           # "p_values", "stars", "both", "none"
  multiple_correction = "bonferroni",
  pairing = "paired"              # "paired" or "unpaired"
)
```

### UMAP and Clustering Parameters

```r
# Enhanced UMAP configuration during analysis
# - Sample exclusion patterns: "unstained", "compensation", "controls"
# - External metadata import: CSV files with sample matching
# - Transformation options: asinh (recommended), log10, sqrt, none
# - Clustering integration: Automatic transition to clustering analysis
```

### Factor Definition System

```r
# The pipeline uses standardized factors for consistent analysis:
# - tissue_factor: Groups samples by tissue/treatment/condition
# - pairing_factor: Identifies paired samples for statistical analysis

# These factors are automatically created and used throughout all analyses
# ensuring consistent and reproducible workflows across different datasets
```

## Output Files

The enhanced pipeline generates comprehensive output with organized folder structure:

### Data Files
- `*_counts.csv`: Population counts and frequencies with congenic classification
- `*_mfi.csv`: Mean fluorescence intensity data with marker information
- `*_umap_data.csv`: UMAP coordinates with marker expression and metadata
- `*_clustering_results.csv`: Cluster assignments and characterization
- `*_stats_results.csv`: Comprehensive statistical analysis results

### Visualizations
- `engraftment_*.png`: Statistical engraftment ratio plots with error bars
- `heatmap_mfi_*.png`: MFI heatmaps with statistical annotations
- `paired_comparison_*.png`: Enhanced statistical comparison plots
- `umap_*.png`: UMAP visualizations with various color schemes
- `clustering_*.png`: Cluster analysis plots and characterization
- `multi_marker_*.png`: Multi-marker heatmap overlays

### Reports
- `analysis_report.txt`: Comprehensive analysis summary
- `cluster_analysis_report.txt`: Detailed clustering results
- `statistical_summary.csv`: All statistical test results

## Troubleshooting

### Common Issues and Solutions

1. **Package Installation Issues**:
   - The pipeline automatically handles package installation
   - If issues persist, manually install BiocManager: `install.packages("BiocManager")`

2. **"Node not found" errors**:
   - Use interactive node selection to see available populations
   - Check FlowJo workspace gating consistency

3. **File path issues**:
   - Always use `here()` for file paths
   - Ensure `.Rproj` file is in project root directory

4. **Memory issues with UMAP/Clustering**:
   - Use sample exclusion to remove unnecessary samples
   - Reduce `max_cells_per_sample` parameter in UMAP analysis

5. **Factor Definition Issues**:
   - Ensure sample names contain identifiable patterns
   - Use interactive factor definition for custom setups

### Performance Optimization

- **Large Datasets**: Use sample exclusion and cell subsampling
- **Memory Management**: The pipeline includes automatic memory optimization
- **Session Management**: Use persistent sessions to avoid recomputing analyses
- **Export Organization**: Automated file organization prevents data loss

## Advanced Features

### Session Management
- **Plot Sessions**: Persistent storage of all visualizations
- **Analysis Continuity**: Resume analysis from any point
- **Export Management**: Organized export of all results

### Data Integration
- **Multiple Experiments**: Merge compatible GatingSets
- **External Metadata**: Robust CSV import with validation
- **Factor Standardization**: Consistent analysis variables across datasets

### Statistical Robustness
- **Multiple Test Correction**: Bonferroni, FDR, and other methods
- **Effect Size Calculations**: Cohen's d and other effect size measures
- **Confidence Intervals**: Proper uncertainty quantification

## Contributing

To contribute to this package:

1. Fork the repository
2. Create a feature branch with comprehensive testing
3. Update documentation and examples
4. Submit a pull request with detailed description

## Citation

If you use this pipeline in your research, please cite:

```
Enhanced Flow Cytometry Analysis Pipeline with Interactive Visualization 
and Clustering Tools. 2024. GitHub: https://github.com/NathanKieswetter/R_Flow_helpers
```

## License

This project is licensed under the MIT License - see the repository for details.

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub: https://github.com/NathanKieswetter/R_Flow_helpers
- Check the troubleshooting section above
- Review the example workflow for proper usage patterns

---

**Note**: This enhanced pipeline provides a comprehensive, production-ready workflow for flow cytometry analysis with particular strength in congenic marker analysis, statistical testing, UMAP analysis, and clustering. The system includes robust error handling, comprehensive navigation options, session management, and publication-ready output generation suitable for both exploratory analysis and final publication.
