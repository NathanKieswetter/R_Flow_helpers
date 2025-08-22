# Flow Cytometry Analysis Pipeline

A comprehensive R package for automated flow cytometry analysis with interactive visualization tools, specifically designed for congenic marker analysis and engraftment studies.

## Overview

This pipeline provides a complete workflow for flow cytometry data analysis, from FlowJo workspace import to publication-ready visualizations. It features interactive data exploration, automated congenic marker detection, statistical analysis, and advanced visualization tools including UMAP dimensionality reduction.

### Key Features

- **Interactive Data Analysis**: Menu-driven selection of cell populations, markers, and analysis parameters
- **Enhanced Metadata Configuration**: Flexible metadata detection and mapping system
- **Interactive Channel Annotation**: Comprehensive marker naming and annotation tools
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
source("https://raw.githubusercontent.com/NathanKieswetter/R_Flow_helpers/main/Congenics_flow_semi-automated.R")
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

## Enhanced Features (NEW)

### Metadata Configuration System

The pipeline now includes a flexible metadata detection and configuration system that automatically adapts to different experimental setups.

#### Automatic Metadata Detection

```r
# Detect metadata structure automatically
metadata_config <- detect_metadata_structure(gs, interactive = TRUE)

# Preview available metadata
preview_metadata_structure(gs)

# Validate configuration
validate_metadata_config(gs, metadata_config)

# Save configuration for reuse
save_metadata_config(metadata_config, "my_experiment_config.rds")

# Load saved configuration
loaded_config <- load_metadata_config("my_experiment_config.rds")
```

#### Flexible Analysis with Custom Metadata

```r
# Run analysis with custom metadata configuration
congenics_results <- analyze_flow_data_flexible(
  gs, 
  metadata_config = metadata_config,
  interactive_annotation = TRUE  # Enable channel annotation
)

# Or use the enhanced auto analysis
congenics_results <- analyze_flow_data_auto_enhanced(
  gs, 
  metadata_config = metadata_config,
  add_congenics = TRUE
)
```

### Interactive Channel Annotation System

For datasets where marker names weren't properly defined during acquisition, the pipeline provides comprehensive annotation tools.

#### Channel Annotation Workflow

```r
# Check current channel-marker mapping
current_lookup <- get_marker_lookup_enhanced(gs, use_saved_annotations = TRUE)

# Run interactive annotation (comprehensive menu system)
annotated_lookup <- annotate_channels_interactive(gs, save_annotation = TRUE)

# Load previously saved annotations
saved_lookup <- load_channel_annotations(gs)

# Apply annotations during analysis
congenics_results <- analyze_flow_data_flexible(
  gs,
  interactive_annotation = TRUE  # Triggers annotation workflow if needed
)
```

#### Annotation Features

- **Pattern-based annotation**: Bulk annotate channels using regex patterns
- **Common panel loading**: Pre-defined marker panels (T-cell, PBMC, tissue analysis)
- **Data preview**: View expression statistics to help identify markers
- **Search functionality**: Find markers by keyword
- **Save/load system**: Reuse annotations across sessions

### Flexible Plotting Functions

Enhanced plotting functions that work with any metadata structure:

```r
# Flexible paired comparison plots
plots <- create_flexible_paired_plots(data, metadata_config = metadata_config)

# Flexible engraftment plots
engraftment_plot <- create_flexible_engraftment_plot(
  data, 
  metadata_config = metadata_config,
  interactive = TRUE,
  interactive_stats = TRUE
)

# Enhanced UMAP with flexible metadata
umap_results <- analyze_flow_umap_flexible(gs, metadata_config = metadata_config)
```

## Core Functions

### Data Import and Setup

#### `setup_flowjo_workspace(xml_path, fcs_path, keywords)`
Imports FlowJo workspace and creates a GatingSet object.

**Parameters:**
- `xml_path`: Path to FlowJo workspace (.wsp) file
- `fcs_path`: Path to directory containing .fcs files
- `keywords`: Sample metadata keywords to extract (default: `c("$WELLID", "GROUPNAME")`)

### Enhanced Metadata Functions (NEW)

#### `detect_metadata_structure(gs, custom_keywords, interactive)`
Automatically detects and configures metadata mappings.

**Parameters:**
- `gs`: GatingSet object
- `custom_keywords`: Additional custom keyword mappings
- `interactive`: Enable interactive refinement (default: TRUE)

**Features:**
- Auto-detection of common metadata patterns
- Interactive refinement and validation
- Support for custom experimental designs
- Comprehensive error handling

#### `annotate_channels_interactive(gs, save_annotation)`
Interactive channel annotation system for undefined markers.

**Parameters:**
- `gs`: GatingSet object
- `save_annotation`: Save annotations for future use (default: TRUE)

**Features:**
- Menu-driven annotation workflow
- Pattern-based bulk annotation
- Common marker panel loading
- Data preview for marker identification
- Search and filtering capabilities

### Interactive Analysis

#### `analyze_flow_data_auto_enhanced(gs, ...)`
Enhanced main analysis function with flexible metadata support.

**New Features:**
- Automatic metadata detection
- Optional channel annotation
- Congenic marker detection
- Statistical analysis integration
- Customizable metadata mappings

#### `analyze_flow_data_flexible(gs, metadata_config, ...)`
Completely flexible analysis function that adapts to any metadata structure.

**Parameters:**
- `gs`: GatingSet object
- `metadata_config`: Custom metadata configuration
- `node_selection`: Node selection method (default: "interactive")
- `interactive_annotation`: Enable channel annotation workflow

### Data Processing

#### `data_clean_custom(data)`
Enhanced data cleaning with unstained sample detection.

**Features:**
- Interactive unstained sample removal
- Column name standardization
- Comprehensive data validation
- Works with both single data frames and lists

#### `assign_metadata_menu(df)`
Interactive metadata assignment for experimental design.

**Features:**
- WT/KO marker assignment
- Recipient marker identification
- Custom metadata addition (timepoint, batch, sex, etc.)
- Genotype classification

### Enhanced Visualization Functions

#### `create_mfi_heatmaps_interactive_enhanced(mfi_data)`
Advanced MFI heatmap creation with statistical testing.

**New Features:**
- Interactive congenic and marker selection
- Multiple scaling methods (Z-score, log, percentile)
- Comprehensive statistical testing suite
- Significance annotation on heatmaps
- Export functionality

#### `create_flexible_paired_plots(data, metadata_config)`
Flexible paired comparison plots that work with any metadata structure.

**Features:**
- Automatic metadata adaptation
- Multiple statistical test options
- Flexible grouping variables
- Publication-ready output

#### `create_flexible_engraftment_plot(data, metadata_config, ...)`
Enhanced engraftment analysis with flexible metadata support.

**Features:**
- Interactive statistical method selection
- Flexible normalization options
- Multiple comparison types
- Comprehensive error handling

#### `analyze_flow_umap_enhanced(gs)`
Comprehensive UMAP analysis with enhanced features.

**Features:**
- Interactive sample exclusion
- External metadata import
- Session-based plot management
- Enhanced visualization menu
- Robust error handling

## Example Workflow

### Enhanced Setup and Data Import

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

# 3. Setup folder structure (automatic)
ensure_output_dirs()

# 4. Import FlowJo workspace
gs <- setup_flowjo_workspace(
  xml_path = here("data", "Exp2_D7.wsp"),
  fcs_path = here("data", "FCS_D7/")
)

# 5. Configure metadata (NEW)
metadata_config <- detect_metadata_structure(gs, interactive = TRUE)
save_metadata_config(metadata_config, "experiment_D7_metadata.rds")

# 6. Annotate channels if needed (NEW)
preview_metadata_structure(gs)
annotated_lookup <- annotate_channels_interactive(gs, save_annotation = TRUE)

# 7. Save GatingSet for future use
save_gatingset(gs, "experiment_D7_gatingset.rds")

# 8. Assess gating hierarchy
plot(gs, fontsize = 15, bool = TRUE)
```

### Enhanced Interactive Data Analysis

```r
# 9. Run enhanced flexible analyses
congenics_results <- analyze_flow_data_flexible(
  gs, 
  metadata_config = metadata_config,
  interactive_annotation = TRUE,
  add_congenics = TRUE
)

# 10. Additional marker analyses with same configuration
KLRG1_CD127_results <- analyze_flow_data_flexible(gs, metadata_config = metadata_config)
CD69_CD103_results <- analyze_flow_data_flexible(gs, metadata_config = metadata_config)

# 11. Save analysis results
save_analysis_results(congenics_results, "congenics_results.rds")
save_analysis_results(KLRG1_CD127_results, "KLRG1_CD127_results.rds")
save_analysis_results(CD69_CD103_results, "CD69_CD103_results.rds")

# 12. Clean up data
congenics_results <- data_clean_custom(congenics_results)
KLRG1_CD127_results <- data_clean_custom(KLRG1_CD127_results)
CD69_CD103_results <- data_clean_custom(CD69_CD103_results)
```

### Enhanced Metadata Assignment and Visualization

```r
# 13. Add experimental metadata (enhanced interface)
congenics_results_with_genotype <- assign_metadata_menu(congenics_results$counts)

# 14. Create enhanced engraftment analysis
plot_engraftment <- create_flexible_engraftment_plot(
  congenics_results_with_genotype,
  metadata_config = metadata_config,
  interactive = TRUE,
  interactive_stats = TRUE
)
print(plot_engraftment)

# 15. Generate enhanced MFI heatmaps with statistics
mfi_heatmaps <- create_mfi_heatmaps_interactive_enhanced(congenics_results$mfi)
# Export statistical results
stats_results <- export_stats_results(mfi_heatmaps$statistics, "mfi_statistical_analysis.csv")

# 16. Create flexible paired comparison plots
KLRG1_CD127_plots <- create_flexible_paired_plots(
  KLRG1_CD127_results$counts, 
  metadata_config = metadata_config
)
CD69_CD103_plots <- create_flexible_paired_plots(
  CD69_CD103_results$counts,
  metadata_config = metadata_config
)

# 17. Enhanced UMAP analysis with flexible metadata
umap_results <- analyze_flow_umap_flexible(gs, metadata_config = metadata_config)
```

### Enhanced Data Export and Session Management

```r
# 18. Export processed data with metadata configuration
write_csv(congenics_results$counts, here("out", "data", "congenics_counts.csv"))
write_csv(congenics_results$mfi, here("out", "data", "congenics_mfi.csv"))

# Save metadata configuration with results
saveRDS(metadata_config, here("out", "data", "metadata_config.rds"))

# 19. Export UMAP data and session
export_umap_data(umap_results, "umap_analysis.csv")
export_all_session_plots()

# 20. Save complete enhanced session
session_name <- save_session(gs, congenics_results, "enhanced_experiment_D7_session")

# 21. Advanced session management
manage_saved_files()  # Interactive menu for loading/deleting saved data
```

### Loading and Continuing Enhanced Sessions

```r
# Load previous session
loaded_session <- load_session("enhanced_experiment_D7_session")
gs <- loaded_session$gs
results <- loaded_session$results

# Load saved metadata configuration
metadata_config <- load_metadata_config("experiment_D7_metadata.rds")

# Load saved channel annotations
annotated_lookup <- load_channel_annotations(gs, "channel_annotations.rds")

# Continue UMAP visualization with enhanced features
restart_visualization_session(umap_results)

# Create quick plots with flexible metadata
quick_plot <- create_quick_umap_plot(umap_results, metadata_config$tissue)
ggsave(here("out", "plots", "quick_tissue_plot.png"), quick_plot)
```

## Enhanced Configuration Options

### Metadata Configuration

The pipeline supports flexible metadata mapping:

```r
# Custom metadata mappings
custom_config <- list(
  tissue = "Location",           # Map tissue to "Location" column
  sample_id = "SampleName",      # Map sample ID to "SampleName" column
  time_point = "Day",            # Map timepoint to "Day" column
  treatment = "Condition",       # Map treatment to "Condition" column
  batch = "Experiment",          # Map batch to "Experiment" column
  subject = "MouseID"            # Map subject to "MouseID" column
)

# Apply custom configuration
results <- analyze_flow_data_flexible(gs, metadata_config = custom_config)
```

### Channel Annotation

**Common Marker Panels:**
- Mouse T-cell Basic (CD3, CD4, CD8, CD45, Live/Dead)
- Human PBMC Basic (CD3, CD4, CD8, CD19, CD14, CD45, Live/Dead)
- Mouse Tissue Analysis (CD45, CD3, CD11b, F4/80, Ly6G, Live/Dead)

**Annotation Methods:**
- **Interactive selection**: Point-and-click marker assignment
- **Pattern matching**: Bulk annotation using regex patterns
- **Panel loading**: Apply pre-defined marker panels
- **Data preview**: Statistical analysis to help identify markers

### Statistical Testing

Enhanced statistical testing with comprehensive options:
- **Two-group tests**: Paired/unpaired t-tests, Mann-Whitney U, Wilcoxon signed-rank
- **Multi-group tests**: ANOVA, Kruskal-Wallis, repeated measures ANOVA, Friedman
- **Multiple comparisons**: Bonferroni, FDR, Tukey HSD
- **Effect size calculations**: Cohen's d, eta-squared
- **Power analysis**: Sample size recommendations

### Visualization Options

**Enhanced Features:**
- **Interactive menus**: Streamlined selection workflows
- **Session management**: Persistent plot storage and retrieval
- **Export options**: Multiple formats (PNG, PDF, SVG, JPEG)
- **Statistical annotations**: P-values, significance stars, effect sizes
- **Color schemes**: Viridis, RColorBrewer palettes, custom gradients

## Troubleshooting

### Common Issues

1. **"Node not found" errors**: Ensure your FlowJo workspace contains the expected populations
2. **File path issues**: Always use `here()` for file paths and ensure your `.Rproj` file is in the project root
3. **Missing FCS files**: Verify that .fcs files are in the specified directory and match workspace samples
4. **Memory issues with UMAP**: Reduce `max_cells_per_sample` parameter for large datasets
5. **Metadata mapping errors**: Use `preview_metadata_structure()` to check available columns
6. **Channel annotation issues**: Use `get_marker_lookup_enhanced()` to verify current mappings

### Enhanced Data Requirements

- **FlowJo workspace**: Compatible .wsp files (FlowJo 10+)
- **FCS files**: Standard .fcs format with proper channel naming
- **Compensation**: Compensated parameters should contain "Comp" in channel names
- **Sample naming**: Consistent naming between workspace and FCS files
- **Metadata structure**: Any column-based metadata format supported

### Performance Optimization

- **Large datasets**: Use sample exclusion and cell subsampling
- **Memory management**: Process data in batches for very large experiments
- **Visualization**: Use session management to avoid recomputing plots
- **Statistical testing**: Consider multiple testing correction for large marker panels

## Enhanced Output Files

The pipeline generates comprehensive output with enhanced organization:

### Data Files
- `*_counts.csv`: Population counts and frequencies
- `*_mfi.csv`: Mean fluorescence intensity data
- `*_umap_data.csv`: UMAP coordinates with marker expression
- `*_stats_results.csv`: Statistical analysis results
- `metadata_config.rds`: Saved metadata configuration
- `channel_annotations.csv/rds`: Saved channel annotations

### Plots
- `engraftment_*.png`: Enhanced engraftment ratio plots
- `heatmap_*.png`: MFI heatmaps with statistical annotations
- `paired_comparison_*.png`: Statistical comparison plots
- `umap_*.png`: UMAP visualizations with session management
- `multi_marker_*.png`: Multi-marker heatmap overlays

### Session Data
- `*_gatingset.rds`: Saved GatingSet object
- `*_results.rds`: Complete analysis results
- `*_session.rds`: Complete analysis session with metadata

### Configuration Files
- `metadata_config.rds`: Reusable metadata mappings
- `channel_annotations.rds`: Reusable channel annotations
- `analysis_config.rds`: Complete analysis configuration

## Contributing

To contribute to this package:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Update documentation
6. Submit a pull request

## Citation

If you use this pipeline in your research, please cite:

```
Nathan Scott Kieswetter. Enhanced Interactive Flow Cytometry Analysis Pipeline. 2025.
```

## License

[Your license information here]

## Support

For questions, issues, or feature requests:
- Open an issue on GitHub
- Contact: natescottkieswetter@gmail.com

---

**Note**: This enhanced pipeline now supports flexible experimental designs and provides comprehensive tools for metadata management, channel annotation, and statistical analysis. The system automatically adapts to different data structures while maintaining backward compatibility with existing workflows.
