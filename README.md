# Flow Cytometry Analysis Pipeline (Semi-specific to Aurora for now)

A comprehensive R package for interactive flow cytometry data analysis with FlowJo workspace integration, featuring automated population analysis, statistical testing, visualization, and UMAP dimensionality reduction.

## Overview

This pipeline provides a complete workflow for flow cytometry analysis including:
- **FlowJo workspace integration** with automated gating set creation
- **Interactive population selection** with regex and hierarchical support
- **Automated frequency and MFI extraction** with custom parent mapping
- **Statistical analysis** with multiple test options and corrections
- **Advanced visualizations** including heatmaps, paired comparisons, and engraftment plots
- **UMAP analysis** for single-cell dimensionality reduction
- **Congenic marker analysis** for transplantation studies

## Installation

### Dependencies

```r
# CRAN packages
install.packages(c(
  "tidyverse", "here", "ggpubr", "rstatix", "scales", "glue",
  "RColorBrewer", "umap", "ComplexHeatmap", "circlize"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("CytoML", "flowCore", "flowWorkspace"))
```

### Loading the Package

```r
source("path/to/flowcyt_pipeline.R")
```

## Quick Start

### Basic Workflow

```r
# 1. Setup FlowJo workspace
gs <- setup_flowjo_workspace(
  xml_path = here("data/your_workspace.wsp"),
  fcs_path = here("data/fcs_files/"),
  keywords = c("$WELLID", "GROUPNAME")
)

# 2. Interactive analysis with automatic congenic detection
results <- analyze_flow_data_auto(gs)

# 3. Clean up data (removes unstained samples, cleans names)
results <- data_clean_custom(results)

# 4. Add metadata annotations
results_annotated <- assign_metadata_menu(results$counts)

# 5. Create visualizations
plots <- create_paired_comparison_plots(results_annotated)
```

## Core Functions

### Data Extraction

#### `setup_flowjo_workspace(xml_path, fcs_path, keywords)`
Initialize FlowJo workspace and create gating set.

**Parameters:**
- `xml_path`: Path to FlowJo workspace file (.wsp)
- `fcs_path`: Directory containing FCS files  
- `keywords`: FCS keywords to extract (default: `c("$WELLID", "GROUPNAME")`)

#### `analyze_flow_data_auto(..., add_congenics = TRUE)`
Main interactive analysis function with automatic congenic detection.

**Features:**
- Interactive node selection with back navigation
- Hierarchical parent mapping options
- MFI channel selection with compensation detection
- Automatic congenic marker identification
- Returns counts/frequencies and MFI data

### Data Processing

#### `data_clean_custom(data)`
Comprehensive data cleaning function.

**Operations:**
- Removes `$` symbols from column names
- Cleans NodeShort formatting
- Interactive unstained sample detection and removal
- Handles both single data frames and lists

#### `assign_metadata_menu(df)`
Interactive metadata assignment with genotype classification.

**Features:**
- Wild-type (WT) marker selection
- Knock-out (KO) marker assignment  
- Recipient marker identification
- Custom metadata addition (timepoint, batch, sex, etc.)

### Statistical Analysis

#### `create_paired_comparison_plots(data)`
Generate paired comparison plots with statistical testing.

**Options:**
- Paired t-test or Wilcoxon signed-rank test
- Flexible faceting (by tissue or population)
- Automatic p-value adjustment
- Publication-ready formatting

### Advanced Visualizations

#### `create_mfi_heatmaps_interactive_enhanced(mfi_data)`
Comprehensive MFI heatmap creation with statistical overlays.

**Features:**
- Interactive marker and congenic selection
- Multiple scaling methods (Z-score, log, percentile)
- Statistical testing integration (t-tests, ANOVA, non-parametric)
- Significance annotations on heatmaps
- Export functionality

#### `create_engraftment_plot(data, ...)`
Specialized engraftment analysis for transplantation studies.

**Parameters:**
- `ko_marker`, `wt_marker`: Congenic markers for ratio calculation
- `spleen_group`: Normalization tissue (default: "Spleen")
- `interactive`: Enable interactive parameter selection
- `add_stats`: Include statistical comparisons

### UMAP Analysis

#### `analyze_flow_umap_enhanced(gs)`
Interactive single-cell UMAP analysis workflow.

**Features:**
- Sample exclusion (unstained controls, etc.)
- Population event counting
- External metadata import
- Interactive marker selection
- Persistent visualization session
- Multiple plot types and export options

## Example Workflows

### Congenic Transplantation Analysis

```r
# Setup and extract data
gs <- setup_flowjo_workspace(xml_path, fcs_path)
results <- analyze_flow_data_auto(gs, add_congenics = TRUE)
results <- data_clean_custom(results)

# Add genotype information
counts_with_genotype <- assign_metadata_menu(results$counts)

# Create paired comparison plots
plots <- create_paired_comparison_plots(counts_with_genotype)

# Engraftment analysis
engraftment_plot <- create_engraftment_plot(
  data = counts_with_genotype,
  ko_marker = "CD45.1.2",
  wt_marker = "CD45.1",
  interactive = TRUE
)
```

### MFI Heatmap Analysis

```r
# Interactive heatmap with statistics
heatmaps <- create_mfi_heatmaps_interactive_enhanced(results$mfi)

# Export statistical results
stats_results <- heatmaps$statistics
export_df <- export_stats_results(stats_results, "mfi_statistics.csv")
```

### Single-Cell UMAP Analysis

```r
# Comprehensive UMAP workflow
umap_results <- analyze_flow_umap_enhanced(gs)

# Continue visualization session
interactive_visualization_menu()

# Export results
export_umap_data(umap_results, "umap_results.csv")
export_all_session_plots("umap_plots/")
```

## Key Features

### Interactive Navigation
- **Back navigation** throughout all menus
- **Pattern matching** for population selection (regex support)
- **Preview options** before final selections
- **Error handling** with retry mechanisms

### Flexible Parent Mapping
- **Automatic hierarchical** parent detection
- **Congenic-specific** mapping for transplantation studies
- **Individual parent** assignment per population
- **Custom parent** selection with validation

### Statistical Testing
- **Multiple test options**: t-tests, ANOVA, non-parametric alternatives
- **Paired and unpaired** designs
- **Multiple comparison corrections**
- **Effect size calculations**
- **Publication-ready output**

### Data Quality Control
- **Event counting** and population assessment
- **Sample exclusion** workflows
- **Missing data handling**
- **Validation checks** throughout pipeline

### Export and Reproducibility
- **Session management** for continued analysis
- **Plot storage** and batch export
- **Data export** in multiple formats
- **Parameter logging** for reproducibility

## Congenic Marker Support

The pipeline automatically detects and processes these congenic markers:
- `CD45.1`, `CD45.2`, `CD45.1.2`
- `CD90.1`, `CD90.2`
- Custom marker support via configuration

### Engraftment Calculations
```r
# Automatic ratio calculation: log2(KO:WT)
# Optional normalization to control tissue (e.g., spleen)
# Statistical comparisons across tissues
```

## Output Structure

### Main Results Object
```r
results <- list(
  counts = tibble(),      # Population frequencies and counts
  mfi = tibble(),         # MFI data by marker and population  
  nodes = character(),    # Selected population paths
  parent_mapping = list(), # Custom parent assignments
  channels = character()  # Selected MFI channels
)
```

### Enhanced Results (with congenics)
```r
# Automatic congenic column added to both counts and mfi
# Genotype classification (donor_WT, donor_KO, Recipient)
# Metadata integration (timepoint, batch, sex, custom fields)
```

## Troubleshooting

### Common Issues

1. **No populations found**: Check FlowJo workspace file path and FCS directory
2. **Missing channels**: Verify compensation has been applied in FlowJo
3. **Sample naming**: Ensure FCS filenames match workspace sample names
4. **Memory issues**: Use sample/cell limits for large datasets

### Performance Optimization

```r
# For large datasets, use sampling
extract_single_cell_data(
  gs, node, markers,
  sample_limit = 50,      # Limit number of samples
  max_cells_per_sample = 5000  # Subsample cells per sample
)
```

## Citation

If you use this pipeline in your research, please cite:

```
FlowCyt Analysis Pipeline
Author: [Your Name]
Year: 2024
GitHub: [Repository URL]
```

## License

[Specify your license here]

## Support

For questions, issues, or feature requests:
- GitHub Issues: [Repository URL]/issues
- Email: [Your email]

## Changelog

### Version 1.0.0
- Initial release with full FlowJo integration
- Interactive analysis workflows
- Statistical testing and visualization
- UMAP single-cell analysis
- Congenic marker support
