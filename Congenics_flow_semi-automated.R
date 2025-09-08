#=============================================================================
# Package Installation and Loading
#=============================================================================

# Function to check and install CRAN packages
install_cran_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing CRAN packages:", paste(new_packages, collapse = ", "), "\n")
    install.packages(new_packages, dependencies = TRUE)
  }
}

# Function to check and install Bioconductor packages
install_bioc_packages <- function(packages) {
  # First ensure BiocManager is installed
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) {
    cat("Installing Bioconductor packages:", paste(new_packages, collapse = ", "), "\n")
    BiocManager::install(new_packages, ask = FALSE, update = FALSE)
  }
}

# Define package lists
cran_packages <- c(
  "here",
  "tidyverse", 
  "ggpubr",
  "rstatix",
  "scales",
  "glue",
  "circlize",
  "RColorBrewer",
  "umap",
  "cluster",
  "factoextra",
  "patchwork"
)

bioc_packages <- c(
  "CytoML",
  "flowCore", 
  "flowWorkspace",
  "CATALYST",
  "FlowSOM",
  "ConsensusClusterPlus",
  "ComplexHeatmap"
)

# Check and install packages
cat("Checking and installing required packages...\n")
install_cran_packages(cran_packages)
install_bioc_packages(bioc_packages)

# Load all packages
cat("Loading packages...\n")

# Load CRAN packages
library(here)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(scales)
library(glue)
library(circlize)
library(RColorBrewer)
library(umap)
library(cluster)
library(factoextra)
library(patchwork)

# Load Bioconductor packages
library(CytoML)
library(flowCore)
library(flowWorkspace)
library(CATALYST)
library(FlowSOM)
library(ConsensusClusterPlus)
library(ComplexHeatmap)

cat("All packages loaded successfully!\n")

#=============================================================================
# Folder structure helper functions
#=============================================================================

setup_enhanced_project_structure <- function() {
  
  # Your original folders
  if(!dir.exists(here::here("data"))) {
    cat(sprintf("Creating folder %s\n", here::here("data")))
    dir.create(here::here("data"), recursive = TRUE)
  }
  if(!dir.exists(here::here("out"))) {
    cat(sprintf("Creating folder %s\n", here::here("out")))
    dir.create(here::here("out"), recursive = TRUE)
  }
  if(!dir.exists(here::here("scripts"))) {
    cat(sprintf("Creating folder %s\n", here::here("scripts")))
    dir.create(here::here("scripts"), recursive = TRUE)
  }
  if(!dir.exists(here::here("out/GatingSet"))) {
    cat(sprintf("Creating folder %s\n", here::here("out/GatingSet")))
    dir.create(here::here("out/GatingSet"), recursive = TRUE)
  }
  
  # Additional enhanced folders
  enhanced_folders <- c(
    # Data subdirectories
    "data/raw",
    "data/processed", 
    "data/metadata",
    "data/statistics",
    "data/exported_dataframes",
    "data/umap_data",
    
    # Output subdirectories
    "out/plots",
    "out/plots/paired_comparisons",
    "out/plots/engraftment", 
    "out/plots/heatmaps",
    "out/plots/heatmaps/mfi",
    "out/plots/umap",
    "out/plots/exploratory",
    "out/tables",
    "out/reports",
    "out/sessions"
  )
  
  for(folder in enhanced_folders) {
    full_path <- here::here(folder)
    if(!dir.exists(full_path)) {
      cat(sprintf("Creating folder %s\n", full_path))
      dir.create(full_path, recursive = TRUE)
    }
  }
  
  cat("Enhanced project structure created successfully!\n")
}

# ============================================================================
# SAVE UTILITY FUNCTIONS
# ============================================================================

interactive_save_plot <- function(plot_object, suggested_name = "plot", plot_type = "general", 
                                  auto_save = FALSE, save_folder = NULL) {
  if(auto_save) return(NULL)
  
  # Determine save folder
  if(is.null(save_folder)) {
    save_folder <- switch(plot_type,
                          "paired_comparison" = here::here("out/plots/paired_comparisons"),
                          "engraftment" = here::here("out/plots/engraftment"),
                          "mfi_heatmap" = here::here("out/plots/heatmaps/mfi"),
                          "umap" = here::here("out/plots/umap"),
                          "exploratory" = here::here("out/plots/exploratory"),
                          here::here("out/plots")
    )
  }
  
  cat("\n=== Save Plot ===\n")
  save_choice <- menu(c("Save plot", "Don't save"), title = "Would you like to save this plot?")
  
  if(save_choice == 1) {
    default_name <- paste0(suggested_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    filename <- readline(paste0("Enter filename (default: ", default_name, "): "))
    if(filename == "") filename <- default_name
    
    formats <- c("png", "pdf", "svg", "jpg")
    format_choice <- menu(formats, title = "Choose file format:")
    if(format_choice == 0) format_choice <- 1
    selected_format <- formats[format_choice]
    
    if(!str_detect(filename, paste0("\\.", selected_format, "$"))) {
      filename <- paste0(filename, ".", selected_format)
    }
    
    width_input <- readline("Width in inches (default 10): ")
    height_input <- readline("Height in inches (default 8): ")
    dpi_input <- readline("DPI (default 300): ")
    
    width <- if(width_input == "") 10 else as.numeric(width_input)
    height <- if(height_input == "") 8 else as.numeric(height_input)
    dpi <- if(dpi_input == "") 300 else as.numeric(dpi_input)
    
    if(is.na(width)) width <- 10
    if(is.na(height)) height <- 8
    if(is.na(dpi)) dpi <- 300
    
    filepath <- file.path(save_folder, filename)
    
    tryCatch({
      ggsave(filepath, plot = plot_object, width = width, height = height, dpi = dpi)
      cat("Plot saved to:", filepath, "\n")
      return(filepath)
    }, error = function(e) {
      cat("Error saving plot:", e$message, "\n")
      return(NULL)
    })
  }
  return(NULL)
}

interactive_save_dataframe <- function(df, suggested_name = "data", data_type = "general", 
                                       auto_save = FALSE, save_folder = NULL) {
  if(auto_save) return(NULL)
  
  if(is.null(save_folder)) {
    save_folder <- switch(data_type,
                          "statistics" = here::here("data/statistics"),
                          "processed" = here::here("data/processed"),
                          "metadata" = here::here("data/metadata"),
                          "umap" = here::here("data/umap_data"),
                          here::here("data/exported_dataframes")
    )
  }
  
  cat("\n=== Save Data ===\n")
  save_choice <- menu(c("Save data", "Don't save"), title = "Would you like to save this data?")
  
  if(save_choice == 1) {
    default_name <- paste0(suggested_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    filename <- readline(paste0("Enter filename (default: ", default_name, "): "))
    if(filename == "") filename <- default_name
    
    formats <- c("csv", "xlsx", "rds", "tsv")
    format_choice <- menu(formats, title = "Choose file format:")
    if(format_choice == 0) format_choice <- 1
    selected_format <- formats[format_choice]
    
    if(!str_detect(filename, paste0("\\.", selected_format, "$"))) {
      filename <- paste0(filename, ".", selected_format)
    }
    
    filepath <- file.path(save_folder, filename)
    
    tryCatch({
      switch(selected_format,
             "csv" = write_csv(df, filepath),
             "xlsx" = {
               if(require(openxlsx, quietly = TRUE)) {
                 write.xlsx(df, filepath)
               } else {
                 filepath <- str_replace(filepath, "\\.xlsx$", ".csv")
                 write_csv(df, filepath)
               }
             },
             "rds" = saveRDS(df, filepath),
             "tsv" = write_tsv(df, filepath)
      )
      cat("Data saved to:", filepath, "\n")
      return(filepath)
    }, error = function(e) {
      cat("Error saving data:", e$message, "\n")
      return(NULL)
    })
  }
  return(NULL)
}

# ============================================================================
# SECTION 1: CORE SETUP AND DATA LOADING
# ============================================================================

# Function to extract and preview keywords from FCS files
preview_fcs_keywords <- function(fcs_path, sample_limit = 5) {
  cat("Scanning FCS files for available keywords...\n")
  
  # Get list of FCS files
  fcs_files <- list.files(fcs_path, pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)
  
  if(length(fcs_files) == 0) {
    stop("No FCS files found in the specified path: ", fcs_path)
  }
  
  cat("Found", length(fcs_files), "FCS files\n")
  
  # Limit number of files to sample for keyword extraction (for performance)
  files_to_check <- if(length(fcs_files) > sample_limit) {
    sample(fcs_files, sample_limit)
  } else {
    fcs_files
  }
  
  cat("Sampling", length(files_to_check), "files for keyword analysis...\n")
  
  # Extract keywords from sampled files
  all_keywords <- list()
  
  for(i in seq_along(files_to_check)) {
    file_path <- files_to_check[i]
    file_name <- basename(file_path)
    
    cat("Reading keywords from:", file_name, "\n")
    
    tryCatch({
      # Read FCS file
      fcs_data <- flowCore::read.FCS(file_path, transformation = FALSE)
      
      # Extract keywords
      keywords <- flowCore::keyword(fcs_data)
      
      # Store keywords with file info
      all_keywords[[file_name]] <- keywords
      
    }, error = function(e) {
      cat("  Warning: Could not read", file_name, ":", e$message, "\n")
    })
  }
  
  if(length(all_keywords) == 0) {
    stop("Could not extract keywords from any FCS files")
  }
  
  # Find common keywords across files
  keyword_names <- unique(unlist(lapply(all_keywords, names)))
  
  # Create summary of keywords with examples
  keyword_summary <- map_dfr(keyword_names, function(kw) {
    # Get values for this keyword across files
    values <- map_chr(all_keywords, function(file_keywords) {
      if(kw %in% names(file_keywords)) {
        val <- file_keywords[[kw]]
        
        # Handle NULL values first
        if(is.null(val)) {
          return("(empty)")
        }
        
        # Handle different value types safely
        tryCatch({
          # Check if it's a vector with multiple elements
          if(length(val) > 1) {
            # For vectors, show first few elements
            if(length(val) > 3) {
              collapsed_val <- paste(c(as.character(head(val, 3)), "..."), collapse = ", ")
            } else {
              collapsed_val <- paste(as.character(val), collapse = ", ")
            }
            return(paste0("[vector:", length(val), "] ", collapsed_val))
          }
          
          # Handle single values
          if(length(val) == 1) {
            # Convert to character safely
            char_val <- as.character(val)
            
            # Check for empty/NA after conversion
            if(is.na(char_val) || char_val == "" || char_val == " " || char_val == "NA") {
              return("(empty)")
            } else {
              return(char_val)
            }
          }
          
          # Handle zero-length values
          if(length(val) == 0) {
            return("(empty)")
          }
          
          # Fallback
          return("(unknown)")
          
        }, error = function(e) {
          # If any error occurs in processing, return a safe value
          return(paste0("(error: ", class(val)[1], ")"))
        })
        
      } else {
        return("(missing)")
      }
    })
    
    # Calculate presence statistics
    non_empty_values <- values[values != "(empty)" & values != "(missing)"]
    unique_values <- unique(non_empty_values)
    
    # Create proper example values string
    example_values_str <- if(length(unique_values) > 0) {
      if(length(unique_values) <= 3) {
        paste(unique_values, collapse = ", ")
      } else {
        paste(c(head(unique_values, 3), "..."), collapse = ", ")
      }
    } else {
      "(no examples)"
    }
    
    tibble(
      keyword = kw,
      present_in_files = sum(values != "(missing)"),
      total_files = length(values),
      non_empty_values = length(non_empty_values),
      unique_values = length(unique_values),
      example_values = example_values_str,
      all_values = list(values),
      all_unique_values = list(unique_values)
    )
  }) %>%
    arrange(desc(present_in_files), desc(non_empty_values), keyword)
  
  return(list(
    keyword_summary = keyword_summary,
    file_keywords = all_keywords,
    files_sampled = files_to_check
  ))
}

# Function to display keyword information in a user-friendly way
display_keyword_info <- function(keyword_data, show_all = FALSE) {
  
  summary_df <- keyword_data$keyword_summary
  
  if(nrow(summary_df) == 0) {
    cat("No keywords found that are present in ALL FCS files\n")
    cat("This suggests inconsistent FCS file formats\n")
    return()
  }
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n")
  cat("KEYWORDS PRESENT IN ALL SAMPLED FILES\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Files analyzed:", length(keyword_data$files_sampled), "\n")
  cat("Keywords present in ALL files:", nrow(summary_df), "\n\n")
  
  # Categorize keywords for better display
  useful_keywords <- summary_df %>% 
    dplyr::filter(non_empty_values == total_files, unique_values > 1) %>%
    arrange(desc(unique_values))
  
  constant_keywords <- summary_df %>% 
    dplyr::filter(non_empty_values == total_files, unique_values == 1)
  
  mostly_empty <- summary_df %>% 
    dplyr::filter(non_empty_values < total_files)
  
  # Display useful keywords (varying values)
  if(nrow(useful_keywords) > 0) {
    cat("🎯 USEFUL KEYWORDS (varying values across files):\n")
    cat("These are most likely to be useful for analysis\n\n")
    
    display_df <- useful_keywords %>%
      select(keyword, unique_values, example_values) %>%
      slice_head(n = if(show_all) nrow(useful_keywords) else 15)
    
    for(i in 1:nrow(display_df)) {
      row <- display_df[i, ]
      cat(sprintf("%-30s | %2d unique values | Examples: %s\n", 
                  row$keyword, row$unique_values, row$example_values))
    }
    
    if(!show_all && nrow(useful_keywords) > 15) {
      cat(sprintf("... and %d more useful keywords (use option 2 to see all)\n", 
                  nrow(useful_keywords) - 15))
    }
    cat("\n")
  }
  
  # Display constant keywords
  if(nrow(constant_keywords) > 0) {
    cat("📊 CONSTANT KEYWORDS (same value in all files):\n")
    cat("These might be useful for experiment metadata\n\n")
    
    display_df <- constant_keywords %>%
      select(keyword, example_values) %>%
      slice_head(n = if(show_all) nrow(constant_keywords) else 10)
    
    for(i in 1:nrow(display_df)) {
      row <- display_df[i, ]
      cat(sprintf("%-30s | Constant value: %s\n", 
                  row$keyword, row$example_values))
    }
    
    if(!show_all && nrow(constant_keywords) > 10) {
      cat(sprintf("... and %d more constant keywords\n", nrow(constant_keywords) - 10))
    }
    cat("\n")
  }
  
  # Display mostly empty keywords
  if(nrow(mostly_empty) > 0) {
    cat("⚠️  MOSTLY EMPTY KEYWORDS (some empty values):\n")
    cat("Use with caution\n\n")
    
    display_df <- mostly_empty %>%
      select(keyword, non_empty_values, total_files, example_values) %>%
      slice_head(n = if(show_all) nrow(mostly_empty) else 5)
    
    for(i in 1:nrow(display_df)) {
      row <- display_df[i, ]
      cat(sprintf("%-30s | %2d/%2d non-empty | Examples: %s\n", 
                  row$keyword, row$non_empty_values, row$total_files, row$example_values))
    }
    
    if(!show_all && nrow(mostly_empty) > 5) {
      cat(sprintf("... and %d more mostly empty keywords\n", nrow(mostly_empty) - 5))
    }
    cat("\n")
  }
  
  cat("💡 TIP: Use option 6 to see detailed examples for any keyword\n")
}

# Function to show detailed examples for a specific keyword
show_keyword_examples <- function(keyword_data, keyword_name) {
  
  if(!keyword_name %in% keyword_data$keyword_summary$keyword) {
    cat("Keyword", keyword_name, "not found in available keywords\n")
    return()
  }
  
  keyword_info <- keyword_data$keyword_summary %>% 
    dplyr::filter(keyword == keyword_name)
  
  values <- keyword_info$all_values[[1]]
  unique_vals <- keyword_info$all_unique_values[[1]]
  file_names <- names(keyword_data$file_keywords)
  
  cat("\n", paste(rep("-", 80), collapse = ""), "\n")
  cat("DETAILED EXAMPLES FOR:", keyword_name, "\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  
  cat("Present in ALL", keyword_info$total_files, "files\n")
  cat("Non-empty values:", keyword_info$non_empty_values, "\n")
  cat("Unique values:", keyword_info$unique_values, "\n\n")
  
  # Show all unique values if reasonable number
  if(length(unique_vals) <= 20 && length(unique_vals) > 1) {
    cat("🎯 ALL UNIQUE VALUES:\n")
    iwalk(unique_vals, ~cat(sprintf("  %2d. %s\n", .y, .x)))
    cat("\n")
  } else if(length(unique_vals) > 20) {
    cat("🎯 SAMPLE OF UNIQUE VALUES (", length(unique_vals), "total):\n")
    sample_vals <- if(length(unique_vals) > 20) head(unique_vals, 20) else unique_vals
    iwalk(sample_vals, ~cat(sprintf("  %2d. %s\n", .y, .x)))
    cat("  ... and", length(unique_vals) - 20, "more unique values\n\n")
  } else if(length(unique_vals) == 1) {
    cat("🎯 CONSTANT VALUE: ", unique_vals[1], "\n\n")
  }
  
  cat("📁 VALUES FROM EACH FILE:\n")
  for(i in seq_along(values)) {
    # Truncate very long values for display
    display_value <- values[i]
    if(nchar(display_value) > 60) {
      display_value <- paste0(str_trunc(display_value, 57), "...")
    }
    cat(sprintf("%-35s: %s\n", str_trunc(file_names[i], 33), display_value))
  }
  cat("\n")
}

# Improved biological relevance detection
detect_biological_relevance <- function(keyword_data) {
  
  # Define more specific biological patterns with word boundaries
  bio_patterns <- list(
    tissues = list(
      pattern = "\\b(liver|spleen|kidney|lymph_node|LN|medLN|cLN|gut|fat|tumor|
      skin|blood|PBMC|SPL|cIEL|cLN|lymph node|cervical|intestine|colon|SI|SG|
      salivary gland|Spleen|SPL|Sp|Peripheralblood|PB|PBMCs|Bonemarrow|BM|Thymus|
      Thy|Inguinallymphnode|iLN|Axillarylymphnode|aLN|Brachiallymphnode|
      bLN|Cervicallymphnode|cLN|Mandibularlymphnode|mandLN|Superficialcervicallymphnode|sCLN
      |Popliteallymphnode|pLN|Mesentericallymphnode|mLN|Mediastinallymphnode|medLN
      |Para-aortic/lumbarlymphnode|paLN|Draininglymphnode|dLN|Skin-draininglymphnodes
      |Peyer’spatches|PP|Smallintestine|SI|Smallintestinalintraepitheliallymphocytes
      |siIEL|Intraepitheliallymphocytes|IEL|Largeintestinalintraepitheliallymphocytes
      |liIEL|Intestinallaminapropria|siLP|cLP|Laminaproprialymphocytes|LPL
      |Colon/largeintestine|LI|Nasal-associatedlymphoidtissue|NALT|Bronchus/bronchus-associatedlymphoidtissue
      |BALT|Lung|Bronchoalveolarlavage|BAL|Liver|LIV|Peritonealcavity/peritonealexudatecells
      |PerC|PEC|Skin|Earskin|Dorsalskin|Visceraladiposetissue|VAT|Epididymalwhiteadiposetissue
      |eWAT|Subcutaneousadiposetissue|sWAT|Tumor|Tumor-infiltratinglymphocytes|TILs|Brain
      |Centralnervoussystem|CNS|Spinalcord|SC|Meninges|Kidney|Heart|Salivarygland
      |SG|Eye/oculartissues|Eye)\\b",
      description = "tissue types"
    ),
    treatments = list(
      pattern = "\\b(stim|stimulated|control|unstained|KO|WT|wildtype|congenic
      |Tamoxifen|tamox|LPS|PMA|SEB|gp33|listeria|LCMV)\\b",
      description = "treatments"
    ),
    timepoints = list(
      pattern = "\\b(D[0-9]+|DPI|day|days)\\b",
      description = "timepoints"
    ),
    well_positions = list(
      pattern = "\\b[A-H][0-9]{1,2}\\b",
      description = "well positions"
    ),
    mouse_info = list(
      pattern = "\\b(M[0-9]+|mouse|mice)\\b",
      description = "mouse identifiers"
    ),
    sample_info = list(
      pattern = "\\b(tissue_type|sample_type|group|condition|treatment)\\b",
      description = "sample metadata"
    )
  )
  
  # Check each keyword for biological relevance
  bio_relevant_indices <- c()
  relevance_info <- list()
  
  for(i in 1:nrow(keyword_data$keyword_summary)) {
    row <- keyword_data$keyword_summary[i, ]
    found_patterns <- character(0)
    
    # Check keyword name and values against each pattern
    for(pattern_name in names(bio_patterns)) {
      pattern <- bio_patterns[[pattern_name]]$pattern
      description <- bio_patterns[[pattern_name]]$description
      
      # Check keyword name
      keyword_match <- str_detect(tolower(row$keyword), regex(pattern, ignore_case = TRUE))
      
      # Check example values
      values_match <- FALSE
      if(!is.na(row$example_values) && row$example_values != "" && row$example_values != "(no examples)") {
        values_match <- str_detect(tolower(row$example_values), regex(pattern, ignore_case = TRUE))
      }
      
      # Check all unique values for more thorough matching
      if(!values_match && "all_unique_values" %in% names(row)) {
        all_vals <- row$all_unique_values[[1]]
        if(length(all_vals) > 0) {
          values_match <- any(str_detect(tolower(all_vals), regex(pattern, ignore_case = TRUE)))
        }
      }
      
      if(keyword_match || values_match) {
        found_patterns <- c(found_patterns, description)
      }
    }
    
    if(length(found_patterns) > 0) {
      bio_relevant_indices <- c(bio_relevant_indices, i)
      relevance_info[[row$keyword]] <- found_patterns
    }
  }
  
  return(list(
    indices = bio_relevant_indices,
    relevance_info = relevance_info
  ))
}

# Main interactive keyword selection function
select_keywords_interactive <- function(fcs_path, sample_limit = 5) {
  
  cat("=== Interactive Keyword Selection for FlowJo Workspace ===\n")
  cat("This will help you choose the best keywords for your analysis\n\n")
  
  # Extract keyword information
  cat("Step 1: Analyzing FCS files for available keywords...\n")
  keyword_data <- preview_fcs_keywords(fcs_path, sample_limit)
  
  selected_keywords <- character(0)
  
  while(TRUE) {
    cat("\n", paste(rep("=", 60), collapse = ""), "\n")
    cat("KEYWORD SELECTION MENU\n")
    cat(paste(rep("=", 60), collapse = ""), "\n")
    
    if(length(selected_keywords) > 0) {
      cat("Currently selected keywords:", paste(selected_keywords, collapse = ", "), "\n\n")
    } else {
      cat("No keywords selected yet\n\n")
    }
    
    cat("Options:\n")
    cat("1. Show keyword overview\n")
    cat("2. Show all keywords (detailed)\n")
    cat("3. Select keywords from complete list\n")
    cat("4. Select from recommended keywords (Best for 1st pass)\n")
    cat("5. Search keywords by pattern\n")
    cat("6. View detailed examples for a keyword\n")
    cat("7. Remove selected keywords\n")
    cat("8. Finish selection\n")
    
    choice <- readline("Choose option (1-8): ")
    
    if(choice == "1") {
      # Show overview
      display_keyword_info(keyword_data, show_all = FALSE)
      readline("Press Enter to continue...")
      
    } else if(choice == "2") {
      # Show all keywords
      display_keyword_info(keyword_data, show_all = TRUE)
      readline("Press Enter to continue...")
      
    } else if(choice == "3") {
      # Select from all keywords
      all_keywords <- keyword_data$keyword_summary$keyword
      cat("\nAll available keywords:\n")
      iwalk(all_keywords, ~cat(sprintf("%d. %s\n", .y, .x)))
      
      cat("\nEnter keyword numbers (space or comma-separated): ")
      selection <- readline()
      
      if(selection != "") {
        # Parse selection
        tryCatch({
          if(grepl("\\s", selection)) {
            selected_indices <- as.numeric(str_trim(str_split(selection, "\\s+")[[1]]))
          } else {
            selected_indices <- as.numeric(str_trim(str_split(selection, ",")[[1]]))
          }
          
          selected_indices <- selected_indices[!is.na(selected_indices)]
          valid_indices <- selected_indices[selected_indices >= 1 & selected_indices <= length(all_keywords)]
          
          if(length(valid_indices) > 0) {
            new_keywords <- all_keywords[valid_indices]
            selected_keywords <- unique(c(selected_keywords, new_keywords))
            cat("Added keywords:", paste(new_keywords, collapse = ", "), "\n")
          } else {
            cat("No valid selections made\n")
          }
        }, error = function(e) {
          cat("Invalid input format\n")
        })
      }
      
    } else if(choice == "4") {
      # Improved recommended keywords with better biological relevance detection
      
      # Detect biologically relevant keywords
      bio_detection <- detect_biological_relevance(keyword_data)
      bio_relevant_keywords <- keyword_data$keyword_summary[bio_detection$indices, ] %>%
        arrange(desc(unique_values), desc(non_empty_values))
      
      # Also get other useful keywords (varying values) as backup
      other_useful <- keyword_data$keyword_summary %>%
        dplyr::filter(
          non_empty_values == total_files, 
          unique_values > 1,
          # Exclude those already in bio_relevant
          !keyword %in% bio_relevant_keywords$keyword
        ) %>%
        arrange(desc(unique_values)) %>%
        head(10)
      
      # Combine recommendations with bio-relevant first
      if(nrow(bio_relevant_keywords) > 0) {
        cat("\n🧬 BIOLOGICALLY RELEVANT KEYWORDS:\n")
        cat("These keywords contain tissue, treatment, or experimental information\n\n")
        
        for(i in 1:nrow(bio_relevant_keywords)) {
          row <- bio_relevant_keywords[i, ]
          
          # Get relevance information
          relevance_info <- bio_detection$relevance_info[[row$keyword]]
          highlight_info <- if(length(relevance_info) > 0) {
            paste0(" [detected: ", paste(relevance_info, collapse = ", "), "]")
          } else {
            ""
          }
          
          cat(sprintf("%2d. %-25s | %2d values | %s%s\n", 
                      i, row$keyword, row$unique_values, 
                      row$example_values, highlight_info))
        }
        
        recommended <- bio_relevant_keywords$keyword
        
      } else {
        cat("\n🎯 NO BIOLOGICALLY RELEVANT KEYWORDS FOUND\n")
        cat("Showing other useful keywords instead:\n\n")
        recommended <- character(0)
      }
      
      # Add other useful keywords if we have space or no bio-relevant ones found
      if(nrow(other_useful) > 0) {
        start_num <- length(recommended) + 1
        
        if(length(recommended) > 0) {
          cat(sprintf("\n🔬 OTHER USEFUL KEYWORDS (%d more):\n", nrow(other_useful)))
        } else {
          cat("\n🔬 USEFUL KEYWORDS (varying values):\n")
        }
        
        for(i in 1:nrow(other_useful)) {
          row <- other_useful[i, ]
          cat(sprintf("%2d. %-25s | %2d values | %s\n", 
                      start_num + i - 1, row$keyword, row$unique_values, row$example_values))
        }
        
        recommended <- c(recommended, other_useful$keyword)
      }
      
      if(length(recommended) == 0) {
        cat("No suitable keywords found for recommendation\n")
        next
      }
      
      cat(sprintf("\nTotal recommendations: %d keywords\n", length(recommended)))
      cat("Enter keyword numbers (space or comma-separated): ")
      selection <- readline()
      
      if(selection != "") {
        # Parse selection (same logic as above)
        tryCatch({
          if(grepl("\\s", selection)) {
            selected_indices <- as.numeric(str_trim(str_split(selection, "\\s+")[[1]]))
          } else {
            selected_indices <- as.numeric(str_trim(str_split(selection, ",")[[1]]))
          }
          
          selected_indices <- selected_indices[!is.na(selected_indices)]
          valid_indices <- selected_indices[selected_indices >= 1 & selected_indices <= length(recommended)]
          
          if(length(valid_indices) > 0) {
            new_keywords <- recommended[valid_indices]
            selected_keywords <- unique(c(selected_keywords, new_keywords))
            cat("Added keywords:", paste(new_keywords, collapse = ", "), "\n")
            
            # Show what biological relevance was detected
            for(kw in new_keywords) {
              if(kw %in% names(bio_detection$relevance_info)) {
                relevance_info <- bio_detection$relevance_info[[kw]]
                if(length(relevance_info) > 0) {
                  cat("  └─", kw, "detected as:", paste(relevance_info, collapse = ", "), "\n")
                }
              }
            }
          } else {
            cat("No valid selections made\n")
          }
        }, error = function(e) {
          cat("Invalid input format\n")
        })
      }
      
    } else if(choice == "5") {
      # Search by pattern
      pattern <- readline("Enter search pattern (case-insensitive): ")
      if(pattern != "") {
        matches <- grep(pattern, keyword_data$keyword_summary$keyword, 
                        ignore.case = TRUE, value = TRUE)
        
        if(length(matches) > 0) {
          cat("Matching keywords:\n")
          iwalk(matches, ~cat(sprintf("%d. %s\n", .y, .x)))
          
          cat("\nEnter keyword numbers to add (space or comma-separated): ")
          selection <- readline()
          
          if(selection != "") {
            # Parse selection (same logic as above)
            tryCatch({
              if(grepl("\\s", selection)) {
                selected_indices <- as.numeric(str_trim(str_split(selection, "\\s+")[[1]]))
              } else {
                selected_indices <- as.numeric(str_trim(str_split(selection, ",")[[1]]))
              }
              
              selected_indices <- selected_indices[!is.na(selected_indices)]
              valid_indices <- selected_indices[selected_indices >= 1 & selected_indices <= length(matches)]
              
              if(length(valid_indices) > 0) {
                new_keywords <- matches[valid_indices]
                selected_keywords <- unique(c(selected_keywords, new_keywords))
                cat("Added keywords:", paste(new_keywords, collapse = ", "), "\n")
              } else {
                cat("No valid selections made\n")
              }
            }, error = function(e) {
              cat("Invalid input format\n")
            })
          }
        } else {
          cat("No keywords match pattern:", pattern, "\n")
        }
      }
      
    } else if(choice == "6") {
      # View detailed examples
      keyword_name <- readline("Enter keyword name: ")
      if(keyword_name != "") {
        show_keyword_examples(keyword_data, keyword_name)
        readline("Press Enter to continue...")
      }
      
    } else if(choice == "7") {
      # Remove keywords
      if(length(selected_keywords) == 0) {
        cat("No keywords selected to remove\n")
      } else {
        cat("Currently selected keywords:\n")
        iwalk(selected_keywords, ~cat(sprintf("%d. %s\n", .y, .x)))
        
        cat("\nEnter numbers to remove (space or comma-separated): ")
        selection <- readline()
        
        if(selection != "") {
          tryCatch({
            if(grepl("\\s", selection)) {
              remove_indices <- as.numeric(str_trim(str_split(selection, "\\s+")[[1]]))
            } else {
              remove_indices <- as.numeric(str_trim(str_split(selection, ",")[[1]]))
            }
            
            remove_indices <- remove_indices[!is.na(remove_indices)]
            valid_indices <- remove_indices[remove_indices >= 1 & remove_indices <= length(selected_keywords)]
            
            if(length(valid_indices) > 0) {
              removed_keywords <- selected_keywords[valid_indices]
              selected_keywords <- selected_keywords[-valid_indices]
              cat("Removed keywords:", paste(removed_keywords, collapse = ", "), "\n")
            } else {
              cat("No valid selections made\n")
            }
          }, error = function(e) {
            cat("Invalid input format\n")
          })
        }
      }
      
    } else if(choice == "8") {
      # Finish selection
      break
      
    } else {
      cat("Invalid choice. Please select 1-8.\n")
    }
  }
  
  # Final summary
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("FINAL KEYWORD SELECTION\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  if(length(selected_keywords) == 0) {
    cat("No keywords selected. Using default keywords: pairing_factor, tissue_factor\n")
    selected_keywords <- c("pairing_factor", "tissue_factor")
  } else {
    cat("Selected", length(selected_keywords), "keywords:\n")
    iwalk(selected_keywords, ~cat(sprintf("  %d. %s\n", .y, .x)))
  }
  
  return(selected_keywords)
}

# ============================================================================
# FUNCTIONS TO ADD METADATA TO EXISTING GATING SET
# ============================================================================

# Function to add metadata to a GatingSet interactively
add_metadata_to_gatingset <- function(gs) {
  
  cat("=== Add Metadata to GatingSet ===\n")
  
  # Get current phenoData
  current_pd <- pData(gs)
  
  cat("Current sample information:\n")
  cat("Samples:", nrow(current_pd), "\n")
  cat("Current columns:", paste(names(current_pd), collapse = ", "), "\n\n")
  
  # Show first few samples
  cat("Current data preview:\n")
  print(head(current_pd, 3))
  
  # Get sample names for reference
  sample_names <- rownames(current_pd)
  
  while(TRUE) {
    cat("\n=== Metadata Addition Options ===\n")
    cat("1. Add timepoint\n")
    cat("2. Add batch information\n") 
    cat("3. Add treatment/condition\n")
    cat("4. Add custom column with single value (all samples)\n")
    cat("5. Add custom column with individual values per sample\n")
    cat("6. Import metadata from CSV file\n")
    cat("7. Preview current metadata\n")
    cat("8. Finish and return GatingSet\n")
    
    choice <- readline("Choose option (1-8): ")
    
    if(choice == "8") {
      cat("Metadata addition complete!\n")
      break
    }
    
    if(choice == "1") {
      # Add timepoint
      timepoint_value <- readline("Enter timepoint (e.g., D0, Day7, 2weeks): ")
      if(timepoint_value != "") {
        current_pd$timepoint <- timepoint_value
        cat("Added timepoint:", timepoint_value, "to all samples\n")
      }
      
    } else if(choice == "2") {
      # Add batch
      batch_value <- readline("Enter batch identifier (e.g., Batch1, 20231201): ")
      if(batch_value != "") {
        current_pd$batch <- batch_value
        cat("Added batch:", batch_value, "to all samples\n")
      }
      
    } else if(choice == "3") {
      # Add treatment/condition
      treatment_value <- readline("Enter treatment/condition (e.g., Control, LPS, Vehicle): ")
      if(treatment_value != "") {
        current_pd$treatment <- treatment_value
        cat("Added treatment:", treatment_value, "to all samples\n")
      }
      
    } else if(choice == "4") {
      # Add custom column with single value
      col_name <- readline("Enter column name: ")
      if(col_name != "" && !col_name %in% names(current_pd)) {
        col_value <- readline(paste("Enter value for", col_name, ": "))
        if(col_value != "") {
          current_pd[[col_name]] <- col_value
          cat("Added", col_name, ":", col_value, "to all samples\n")
        }
      } else if(col_name %in% names(current_pd)) {
        cat("Column already exists. Choose a different name.\n")
      }
      
    } else if(choice == "5") {
      # Add custom column with individual values
      col_name <- readline("Enter column name: ")
      if(col_name != "" && !col_name %in% names(current_pd)) {
        cat("Enter values for each sample:\n")
        values <- character(length(sample_names))
        
        for(i in seq_along(sample_names)) {
          sample <- sample_names[i]
          value <- readline(paste0("Value for ", sample, " (", i, "/", length(sample_names), "): "))
          values[i] <- if(value == "") NA_character_ else value
        }
        
        current_pd[[col_name]] <- values
        cat("Added individual values for", col_name, "\n")
      } else if(col_name %in% names(current_pd)) {
        cat("Column already exists. Choose a different name.\n")
      }
      
    } else if(choice == "6") {
      # Import from CSV
      csv_path <- readline("Enter path to CSV file (must have 'Sample' column): ")
      
      if(file.exists(csv_path)) {
        tryCatch({
          metadata_df <- read_csv(csv_path, show_col_types = FALSE)
          
          if(!"Sample" %in% names(metadata_df)) {
            cat("Error: CSV must contain 'Sample' column matching sample names\n")
            next
          }
          
          # Show preview
          cat("CSV preview:\n")
          print(head(metadata_df))
          
          # Check which samples match
          matching_samples <- intersect(metadata_df$Sample, sample_names)
          missing_samples <- setdiff(sample_names, metadata_df$Sample)
          
          cat("Matching samples:", length(matching_samples), "\n")
          if(length(missing_samples) > 0) {
            cat("Missing samples:", paste(head(missing_samples, 5), collapse = ", "))
            if(length(missing_samples) > 5) cat("... and", length(missing_samples) - 5, "more")
            cat("\n")
          }
          
          confirm <- readline("Proceed with CSV import? (y/n): ")
          if(tolower(confirm) == "y") {
            # Prepare metadata for joining
            metadata_df <- metadata_df %>%
              column_to_rownames("Sample")
            
            # Only keep columns not already in current_pd
            new_cols <- setdiff(names(metadata_df), names(current_pd))
            if(length(new_cols) > 0) {
              # Add new columns to current_pd
              for(col in new_cols) {
                current_pd[[col]] <- metadata_df[rownames(current_pd), col]
              }
              cat("Added columns from CSV:", paste(new_cols, collapse = ", "), "\n")
            } else {
              cat("No new columns to add (all columns already exist)\n")
            }
          }
          
        }, error = function(e) {
          cat("Error reading CSV:", e$message, "\n")
        })
      } else {
        cat("File not found:", csv_path, "\n")
      }
      
    } else if(choice == "7") {
      # Preview current metadata
      cat("\nCurrent metadata:\n")
      cat("Columns:", paste(names(current_pd), collapse = ", "), "\n")
      print(head(current_pd))
      readline("Press Enter to continue...")
    }
    
    # Update the GatingSet with modified phenoData
    pData(gs) <- current_pd
  }
  
  # Show final summary
  final_pd <- pData(gs)
  cat("\n=== Final Metadata Summary ===\n")
  cat("Total samples:", nrow(final_pd), "\n")
  cat("Total columns:", ncol(final_pd), "\n")
  cat("Column names:", paste(names(final_pd), collapse = ", "), "\n")
  
  return(gs)
}

# Function to add metadata from pattern matching (e.g., extract from sample names)
add_metadata_from_patterns <- function(gs) {
  
  cat("=== Extract Metadata from Sample Names ===\n")
  
  current_pd <- pData(gs)
  sample_names <- rownames(current_pd)
  
  cat("Current sample names:\n")
  print(head(sample_names, 10))
  if(length(sample_names) > 10) {
    cat("... and", length(sample_names) - 10, "more samples\n")
  }
  
  cat("\nCommon patterns to extract:\n")
  cat("1. Extract timepoint (e.g., 'D0', 'Day7', 'Week2' from sample names)\n")
  cat("2. Extract treatment (e.g., 'Ctrl', 'LPS', 'Vehicle' from sample names)\n")
  cat("3. Extract mouse/subject ID (e.g., 'M1', 'Mouse01' from sample names)\n")
  cat("4. Custom regex pattern\n")
  cat("5. Back to main menu\n")
  
  pattern_choice <- readline("Choose option (1-5): ")
  
  if(pattern_choice == "5") return(gs)
  
  if(pattern_choice == "1") {
    # Extract timepoint
    cat("Common timepoint patterns:\n")
    cat("a. D + number (D0, D7, D14)\n")
    cat("b. Day + number (Day0, Day7)\n") 
    cat("c. Week + number (Week1, Week2)\n")
    cat("d. Custom pattern\n")
    
    tp_choice <- readline("Choose timepoint pattern (a-d): ")
    
    pattern <- switch(tp_choice,
                      "a" = "D\\d+",
                      "b" = "Day\\d+", 
                      "c" = "Week\\d+",
                      "d" = {
                        custom_pattern <- readline("Enter regex pattern for timepoint: ")
                        if(custom_pattern == "") NULL else custom_pattern
                      },
                      NULL)
    
    if(!is.null(pattern)) {
      extracted_tp <- str_extract(sample_names, pattern)
      
      # Show preview
      preview_df <- data.frame(
        Sample = sample_names,
        Extracted_Timepoint = extracted_tp
      )
      cat("Preview of extracted timepoints:\n")
      print(head(preview_df, 10))
      
      # Show summary
      tp_counts <- table(extracted_tp, useNA = "ifany")
      cat("\nTimepoint counts:\n")
      print(tp_counts)
      
      confirm <- readline("Add extracted timepoints? (y/n): ")
      if(tolower(confirm) == "y") {
        current_pd$timepoint <- extracted_tp
        pData(gs) <- current_pd
        cat("Timepoint column added successfully!\n")
      }
    }
    
  } else if(pattern_choice == "4") {
    # Custom regex pattern
    cat("Enter custom regex pattern and column name\n")
    col_name <- readline("Column name for extracted data: ")
    pattern <- readline("Regex pattern to extract: ")
    
    if(col_name != "" && pattern != "") {
      extracted_data <- str_extract(sample_names, pattern)
      
      # Show preview
      preview_df <- data.frame(
        Sample = sample_names,
        Extracted = extracted_data
      )
      names(preview_df)[2] <- col_name
      
      cat("Preview of extracted data:\n")
      print(head(preview_df, 10))
      
      confirm <- readline("Add extracted data? (y/n): ")
      if(tolower(confirm) == "y") {
        current_pd[[col_name]] <- extracted_data
        pData(gs) <- current_pd
        cat("Column", col_name, "added successfully!\n")
      }
    }
  }
  
  return(gs)
}

# Enhanced setup function that includes metadata addition
setup_flowjo_workspace_with_metadata <- function(xml_path, fcs_path, keywords = NULL, 
                                                 sample_limit = 5, add_metadata = TRUE) {
  
  # First, set up the basic workspace
  gs <- setup_flowjo_workspace_interactive(xml_path, fcs_path, keywords, sample_limit)
  
  # Then add metadata if requested
  if(add_metadata) {
    cat("\n=== Would you like to add metadata to the GatingSet? ===\n")
    add_meta <- menu(c("Yes, add metadata interactively", 
                       "Yes, extract from sample name patterns",
                       "No, skip metadata addition"), 
                     title = "Add metadata?")
    
    if(add_meta == 1) {
      gs <- add_metadata_to_gatingset(gs)
    } else if(add_meta == 2) {
      gs <- add_metadata_from_patterns(gs)
    }
  }
  
  return(gs)
}


# Enhanced setup function with interactive keyword selection

setup_flowjo_workspace_interactive <- function(xml_path, fcs_path, keywords = NULL, sample_limit = 5) {
  
  cat("=== Enhanced FlowJo Workspace Setup ===\n")
  
  # Interactive keyword selection if not provided
  if(is.null(keywords)) {
    cat("No keywords specified. Starting interactive keyword selection...\n")
    keywords <- select_keywords_interactive(fcs_path, sample_limit)
  } else {
    cat("Using provided keywords:", paste(keywords, collapse = ", "), "\n")
  }
  
  cat("\nSetting up FlowJo workspace...\n")
  
  # Setup workspace with selected keywords
  ws <- open_flowjo_xml(xml_path)
  gs <- flowjo_to_gatingset(
    ws, 
    keywords = keywords,
    keywords.source = "FCS", 
    path = fcs_path, 
    extend_val = -10000
  )
  
  cat("Workspace setup complete!\n")
  cat("Keywords imported:", paste(keywords, collapse = ", "), "\n")
  cat("Samples loaded:", length(sampleNames(gs)), "\n")
  
  return(gs)
}

# Convenience function for just previewing keywords without selection
preview_keywords_only <- function(fcs_path, sample_limit = 5) {
  cat("=== FCS Keyword Preview ===\n")
  keyword_data <- preview_fcs_keywords(fcs_path, sample_limit)
  display_keyword_info(keyword_data, show_all = FALSE)
  
  cat("\nFor detailed exploration, use: select_keywords_interactive(fcs_path)\n")
  
  return(keyword_data)
}

#Extract data/counts
# Completely bulletproof extract_counts_freqs_base function
extract_counts_freqs_base <- function(gs, nodes, parent_mapping = NULL, keywords = NULL) {
  
  # Helper to safely get counts for a node in a sample
  safe_get_count <- function(gh, node) {
    tryCatch({
      if(node %in% gh_get_pop_paths(gh)) {
        list(
          count = gh_pop_get_count(gh, node),
          parent_count = gh_pop_get_count(gh, gh_pop_get_parent(gh, node)),
          success = TRUE
        )
      } else {
        list(success = FALSE)
      }
    }, error = function(e) list(success = FALSE))
  }
  
  # Extract counts for all node-sample combinations
  results <- map_dfr(nodes, function(node) {
    map_dfr(seq_along(gs), function(i) {
      sample_name <- sampleNames(gs)[i]
      gh <- gs[[i]]
      
      result <- safe_get_count(gh, node)
      if(!result$success) return(tibble())
      
      tibble(
        Sample = sample_name,
        Node = node,
        Subpop = result$count,  # Name it Subpop directly to avoid any conflicts
        ParentCount = result$parent_count
      )
    })
  })
  
  if(nrow(results) == 0) {
    stop("No counts extracted. Check that nodes exist in the gating set.")
  }
  
  # Get ALL available sample metadata
  pd <- pData(gs) %>% 
    rownames_to_column("Sample")
  
  # Debug: Check for column name conflicts BEFORE join
  cat("Results columns before join:", paste(names(results), collapse = ", "), "\n")
  cat("Metadata columns:", paste(names(pd), collapse = ", "), "\n")
  
  # Check for conflicts and rename them proactively
  conflicting_cols <- intersect(names(results), names(pd))
  conflicting_cols <- setdiff(conflicting_cols, "Sample")  # Sample is expected to match
  
  if(length(conflicting_cols) > 0) {
    cat("⚠️  Column name conflicts detected:", paste(conflicting_cols, collapse = ", "), "\n")
    cat("Renaming metadata columns to avoid conflicts\n")
    
    # Rename conflicting columns in metadata with .meta suffix
    for(col in conflicting_cols) {
      old_name <- col
      new_name <- paste0(col, ".meta")
      names(pd)[names(pd) == old_name] <- new_name
      cat("Renamed metadata column:", old_name, "->", new_name, "\n")
    }
  }
  
  # Verify results structure before join
  cat("Verifying 'Subpop' column exists:", "Subpop" %in% names(results), "\n")
  
  # Perform the join
  cat("Performing left_join...\n")
  results <- results %>%
    left_join(pd, by = "Sample")
  
  cat("Results columns after join:", paste(names(results), collapse = ", "), "\n")
  
  # Verify Subpop column still exists after join
  if(!"Subpop" %in% names(results)) {
    stop("ERROR: Subpop column disappeared after join! This indicates a serious column conflict.")
  }
  
  # Calculate frequencies with custom parents if provided
  if(!is.null(parent_mapping)) {
    parent_results <- map_dfr(names(parent_mapping), function(node) {
      parent_node <- parent_mapping[node]
      
      map_dfr(seq_along(gs), function(i) {
        sample_name <- sampleNames(gs)[i]
        gh <- gs[[i]]
        
        result <- safe_get_count(gh, parent_node)
        if(!result$success) return(tibble())
        
        tibble(
          Sample = sample_name,
          Node = node,
          CustomParentCount = result$count
        )
      })
    })
    
    results <- results %>%
      left_join(parent_results, by = c("Sample", "Node")) %>%
      mutate(
        Freq = case_when(
          !is.na(CustomParentCount) & CustomParentCount > 0 ~ (Subpop / CustomParentCount) * 100,
          ParentCount > 0 ~ (Subpop / ParentCount) * 100,
          TRUE ~ 0
        )
      )
  } else {
    results <- results %>%
      mutate(Freq = if_else(ParentCount > 0, (Subpop / ParentCount) * 100, 0))
  }
  
  results %>% mutate(NodeShort = basename(Node))
}

# Similarly, fix MFI extraction to not assume standardized columns
extract_mfi_base <- function(gs, nodes, channels = NULL, summary_fun = median, keywords = NULL) {
  if(is.null(channels)) {
    channels <- select_channels_interactive(gs)
    if(length(channels) == 0) return(tibble())
  } else if(length(channels) == 1 && tolower(channels) == "all") {
    channels <- get_marker_lookup(gs)$colname
  } else if(length(channels) == 1 && tolower(channels) == "compensated") {
    lookup <- get_marker_lookup(gs)
    channels <- lookup %>%
      dplyr::filter(str_detect(colname, "(?i)comp")) %>%
      pull(colname)
    
    if(length(channels) == 0) {
      warning("No compensated channels found. Skipping MFI analysis.")
      return(tibble())
    }
  }
  
  # Validate channels
  lookup <- get_marker_lookup(gs)
  missing_channels <- setdiff(channels, lookup$colname)
  if(length(missing_channels) > 0) {
    stop("Channels not found: ", paste(missing_channels, collapse = ", "))
  }
  
  # Extract MFI for all node-sample-channel combinations
  results <- map_dfr(nodes, function(node) {
    map_dfr(seq_along(gs), function(i) {
      sample_name <- sampleNames(gs)[i]
      gh <- gs[[i]]
      
      if(!node %in% gh_get_pop_paths(gh)) return(tibble())
      
      ff <- tryCatch({
        gh_pop_get_data(gh, node)
      }, error = function(e) return(NULL))
      
      if(is.null(ff)) return(tibble())
      
      expr_data <- if(inherits(ff, "cytoframe")) exprs(ff) else exprs(ff)
      df <- as.data.frame(expr_data)
      
      available_channels <- intersect(channels, colnames(df))
      if(length(available_channels) == 0) return(tibble())
      
      map_dfr(available_channels, function(ch) {
        tibble(
          Sample = sample_name,
          Node = node,
          colname = ch,
          MFI = summary_fun(df[[ch]], na.rm = TRUE)
        )
      })
    })
  })
  
  if(nrow(results) == 0) return(tibble())
  
  results <- results %>%
    left_join(lookup, by = "colname") %>%
    select(Sample, Node, colname, marker, MFI)
  
  # Join with ALL sample metadata (don't assume specific columns exist)
  pd <- pData(gs) %>% 
    rownames_to_column("Sample")
  
  results <- results %>%
    left_join(pd, by = "Sample") %>%
    mutate(NodeShort = basename(Node))
  
  return(results)
}

# Now create the corrected main analysis function
analyze_flow_data_corrected <- function(gs, 
                                        node_selection = "interactive", 
                                        parent_selection = "interactive",
                                        channels = NULL,
                                        summary_fun = median) {
  
  cat("=== Flow Cytometry Analysis Pipeline ===\n")
  
  # Get nodes and parents first
  nodes <- if(node_selection == "interactive") {
    select_nodes_interactive(gs)
  } else {
    node_selection
  }
  
  if(is.null(nodes) || length(nodes) == 0) {
    cat("No nodes selected. Exiting.\n")
    return(NULL)
  }
  
  parent_mapping <- if(parent_selection == "interactive") {
    parent_result <- select_parents_interactive(nodes, gs)
    if(is.character(parent_result) && length(parent_result) == 1 && parent_result == "BACK") {
      cat("Going back to node selection...\n")
      return(NULL)
    }
    parent_result
  } else {
    parent_selection
  }
  
  # Extract data using base functions (no standardized columns assumed)
  cat("Extracting counts and frequencies...\n")
  counts <- extract_counts_freqs_base(gs, nodes, parent_mapping)
  
  cat("Extracting MFI data...\n")
  mfi <- extract_mfi_base(gs, nodes, channels, summary_fun)
  
  return(list(
    counts = counts, 
    mfi = mfi, 
    nodes = nodes, 
    parent_mapping = parent_mapping,
    channels = channels
  ))
}

# Corrected enhanced analysis function with proper workflow order
analyze_flow_data_auto_enhanced <- function(gs, 
                                            add_congenics = TRUE, 
                                            congenic_candidates = c("CD45.1", "CD45.2", "CD90.1", "CD90.2", "CD45.1.2"),
                                            define_factors = TRUE,
                                            auto_save = FALSE,
                                            ...) {
  
  # Step 1: Extract raw data (no standardized columns assumed)
  results <- analyze_flow_data_corrected(gs, ...)
  
  if(is.null(results)) return(NULL)
  
  # Step 2: Add congenics if requested
  if(add_congenics) {
    results <- add_congenics_column(results, candidates = congenic_candidates)
  }
  
  # Step 3: Define analysis factors AFTER data extraction
  if(define_factors) {
    cat("\n=== DEFINING ANALYSIS FACTORS ===\n")
    cat("Setting up standardized factor columns for downstream analyses.\n\n")
    
    # Define factors using the counts data
    if(!is.null(results$counts) && nrow(results$counts) > 0) {
      cat("Defining factors using counts data...\n")
      counts_factors <- define_analysis_factors(results$counts)
      results$counts <- counts_factors$data
      
      # Store the factor mapping info
      factor_info <- list(
        tissue_factor = "tissue_factor",
        pairing_factor = "pairing_factor", 
        original_tissue_col = counts_factors$original_tissue_col,
        original_pairing_col = counts_factors$original_pairing_col,
        pairing_enabled = counts_factors$pairing_enabled
      )
    }
    
    # Apply the SAME factor definitions to MFI data - FIXED VERSION
    if(!is.null(results$mfi) && nrow(results$mfi) > 0 && exists("factor_info")) {
      cat("Applying same factor definitions to MFI data...\n")
      
      # For MFI data, we need to handle the fact that it has multiple rows per sample
      # First, create a sample-level metadata table with the factor definitions
      sample_factors <- results$counts %>%
        select(Sample, tissue_factor, pairing_factor) %>%
        distinct()
      
      # Remove any existing tissue_factor and pairing_factor columns from MFI data to avoid conflicts
      results$mfi <- results$mfi %>%
        select(-any_of(c("tissue_factor", "pairing_factor")))
      
      # Join the factor definitions to MFI data
      results$mfi <- results$mfi %>%
        left_join(sample_factors, by = "Sample")
      
      cat("✅ Applied factor definitions to MFI data\n")
    }
    
    # Store factor information in results
    if(exists("factor_info")) {
      results$factor_info <- factor_info
    }
    
    cat("\n✅ Factor definition complete!\n")
    cat("All subsequent analyses will use:\n")
    cat("• tissue_factor for grouping\n")
    cat("• pairing_factor for paired analyses\n")
  }
  
  # Interactive save option
  if(!auto_save) {
    cat("\n=== Save Analysis Results ===\n")
    
    # Save counts data
    if(!is.null(results$counts)) {
      save_counts <- menu(c("Save counts data", "Don't save"), title = "Save cell counts data?")
      if(save_counts == 1) {
        interactive_save_dataframe(results$counts, "flow_counts", "processed")
      }
    }
    
    # Save MFI data
    if(!is.null(results$mfi)) {
      save_mfi <- menu(c("Save MFI data", "Don't save"), title = "Save MFI data?")
      if(save_mfi == 1) {
        interactive_save_dataframe(results$mfi, "flow_mfi", "processed")
      }
    }
  }
  
  return(results)
}

# ============================================================================
# FUNCTIONS TO MERGE TWO GATING SETS
# ============================================================================

# Function to merge two GatingSets
merge_gating_sets <- function(gs1, gs2, conflict_resolution = "interactive") {
  
  cat("=== Merging GatingSets ===\n")
  cat("GatingSet 1 - Samples:", length(sampleNames(gs1)), "\n")
  cat("GatingSet 2 - Samples:", length(sampleNames(gs2)), "\n")
  
  # Check for sample name conflicts
  samples1 <- sampleNames(gs1)
  samples2 <- sampleNames(gs2)
  overlapping_samples <- intersect(samples1, samples2)
  
  if(length(overlapping_samples) > 0) {
    cat("\nWarning: Overlapping sample names found:\n")
    cat(paste(head(overlapping_samples, 10), collapse = ", "))
    if(length(overlapping_samples) > 10) {
      cat("... and", length(overlapping_samples) - 10, "more")
    }
    cat("\n")
    
    if(conflict_resolution == "interactive") {
      cat("\nConflict resolution options:\n")
      cat("1. Rename samples in GS2 (add suffix)\n")
      cat("2. Skip overlapping samples from GS2\n")
      cat("3. Cancel merge\n")
      
      choice <- readline("Choose option (1-3): ")
      
      if(choice == "3") {
        stop("Merge cancelled by user")
      } else if(choice == "1") {
        # Rename overlapping samples in gs2
        suffix <- readline("Enter suffix for GS2 samples (default: '_GS2'): ")
        if(suffix == "") suffix <- "_GS2"
        
        new_names2 <- samples2
        new_names2[samples2 %in% overlapping_samples] <- paste0(samples2[samples2 %in% overlapping_samples], suffix)
        sampleNames(gs2) <- new_names2
        
        cat("Renamed", length(overlapping_samples), "samples in GS2\n")
        
      } else if(choice == "2") {
        # Remove overlapping samples from gs2
        keep_indices <- !samples2 %in% overlapping_samples
        gs2 <- gs2[keep_indices]
        
        cat("Removed", length(overlapping_samples), "overlapping samples from GS2\n")
        cat("GS2 now has", length(sampleNames(gs2)), "samples\n")
      }
    }
  }
  
  # Check gating tree compatibility
  cat("\nChecking gating tree compatibility...\n")
  
  # Get population paths from both GatingSets
  pops1 <- gh_get_pop_paths(gs1[[1]])
  pops2 <- gh_get_pop_paths(gs2[[1]])
  
  common_pops <- intersect(pops1, pops2)
  unique_to_gs1 <- setdiff(pops1, pops2)
  unique_to_gs2 <- setdiff(pops2, pops1)
  
  cat("Common populations:", length(common_pops), "\n")
  cat("Unique to GS1:", length(unique_to_gs1), "\n")
  cat("Unique to GS2:", length(unique_to_gs2), "\n")
  
  if(length(unique_to_gs1) > 0 || length(unique_to_gs2) > 0) {
    cat("\nWarning: Gating trees are not identical\n")
    
    if(length(unique_to_gs1) > 0) {
      cat("Populations only in GS1:\n")
      cat(paste(head(unique_to_gs1, 5), collapse = "\n"))
      if(length(unique_to_gs1) > 5) cat("\n... and", length(unique_to_gs1) - 5, "more\n")
    }
    
    if(length(unique_to_gs2) > 0) {
      cat("Populations only in GS2:\n")
      cat(paste(head(unique_to_gs2, 5), collapse = "\n"))
      if(length(unique_to_gs2) > 5) cat("\n... and", length(unique_to_gs2) - 5, "more\n")
    }
    
    cat("\nNote: Merged GatingSet will contain all populations, but some samples may not have all gates\n")
  }
  
  # Check phenoData compatibility
  cat("\nChecking metadata compatibility...\n")
  pd1 <- pData(gs1)
  pd2 <- pData(gs2)
  
  cols1 <- names(pd1)
  cols2 <- names(pd2)
  
  common_cols <- intersect(cols1, cols2)
  unique_to_pd1 <- setdiff(cols1, cols2)
  unique_to_pd2 <- setdiff(cols2, cols1)
  
  cat("Common metadata columns:", length(common_cols), "\n")
  cat("Unique to GS1 metadata:", length(unique_to_pd1), "\n")
  cat("Unique to GS2 metadata:", length(unique_to_pd2), "\n")
  
  # Perform the merge using rbind2
  cat("\nPerforming merge...\n")
  
  tryCatch({
    merged_gs <- rbind2(gs1, gs2)
    
    cat("Merge successful!\n")
    cat("Merged GatingSet contains:\n")
    cat("- Samples:", length(sampleNames(merged_gs)), "\n")
    cat("- Populations:", length(gh_get_pop_paths(merged_gs[[1]])), "\n")
    cat("- Metadata columns:", ncol(pData(merged_gs)), "\n")
    
    return(merged_gs)
    
  }, error = function(e) {
    cat("Error during merge:", e$message, "\n")
    cat("This may be due to incompatible gating structures or channel differences\n")
    return(NULL)
  })
}

# Function to merge multiple GatingSets
merge_multiple_gating_sets <- function(gs_list, conflict_resolution = "interactive") {
  
  if(length(gs_list) < 2) {
    stop("Need at least 2 GatingSets to merge")
  }
  
  cat("=== Merging Multiple GatingSets ===\n")
  cat("Number of GatingSets:", length(gs_list), "\n")
  
  # Show summary of each GatingSet
  for(i in seq_along(gs_list)) {
    cat(sprintf("GS%d - Samples: %d\n", i, length(sampleNames(gs_list[[i]]))))
  }
  
  # Start with first GatingSet and progressively merge others
  merged_gs <- gs_list[[1]]
  
  for(i in 2:length(gs_list)) {
    cat(sprintf("\nMerging GS%d...\n", i))
    merged_gs <- merge_gating_sets(merged_gs, gs_list[[i]], conflict_resolution)
    
    if(is.null(merged_gs)) {
      stop(sprintf("Failed to merge GS%d", i))
    }
  }
  
  cat("\nFinal merged GatingSet:\n")
  cat("- Total samples:", length(sampleNames(merged_gs)), "\n")
  cat("- Populations:", length(gh_get_pop_paths(merged_gs[[1]])), "\n")
  
  return(merged_gs)
}

# Function to add distinguishing metadata before merging
add_source_metadata <- function(gs, source_name) {
  
  current_pd <- pData(gs)
  current_pd$source_gs <- source_name
  pData(gs) <- current_pd
  
  cat("Added source metadata:", source_name, "to", nrow(current_pd), "samples\n")
  
  return(gs)
}

# Enhanced merge with source tracking
merge_gating_sets_with_tracking <- function(gs1, gs2, gs1_name = "GS1", gs2_name = "GS2", 
                                            conflict_resolution = "interactive") {
  
  # Add source tracking metadata
  gs1_tracked <- add_source_metadata(gs1, gs1_name)
  gs2_tracked <- add_source_metadata(gs2, gs2_name) 
  
  # Perform merge
  merged_gs <- merge_gating_sets(gs1_tracked, gs2_tracked, conflict_resolution)
  
  if(!is.null(merged_gs)) {
    # Show source distribution
    source_counts <- table(pData(merged_gs)$source_gs)
    cat("\nSource distribution in merged GatingSet:\n")
    print(source_counts)
  }
  
  return(merged_gs)
}

# Function to check compatibility before attempting merge
check_gating_set_compatibility <- function(gs1, gs2) {
  
  cat("=== GatingSet Compatibility Check ===\n")
  
  # Check sample overlap
  samples1 <- sampleNames(gs1)
  samples2 <- sampleNames(gs2)
  overlapping <- intersect(samples1, samples2)
  
  cat("Sample name overlap:", length(overlapping), "samples\n")
  
  # Check channel compatibility
  channels1 <- colnames(gh_pop_get_data(gs1[[1]]))
  channels2 <- colnames(gh_pop_get_data(gs2[[1]]))
  
  common_channels <- intersect(channels1, channels2)
  missing_in_gs2 <- setdiff(channels1, channels2)
  missing_in_gs1 <- setdiff(channels2, channels1)
  
  cat("Common channels:", length(common_channels), "\n")
  cat("Channels only in GS1:", length(missing_in_gs2), "\n")
  cat("Channels only in GS2:", length(missing_in_gs1), "\n")
  
  # Check gating tree similarity
  pops1 <- gh_get_pop_paths(gs1[[1]])
  pops2 <- gh_get_pop_paths(gs2[[1]])
  
  common_pops <- intersect(pops1, pops2)
  similarity <- length(common_pops) / max(length(pops1), length(pops2)) * 100
  
  cat("Gating tree similarity:", round(similarity, 1), "%\n")
  
  # Overall compatibility assessment
  compatible <- TRUE
  issues <- character(0)
  
  if(length(overlapping) > 0) {
    issues <- c(issues, "Sample name conflicts")
    compatible <- FALSE
  }
  
  if(length(missing_in_gs1) > 0 || length(missing_in_gs2) > 0) {
    issues <- c(issues, "Channel differences")
  }
  
  if(similarity < 80) {
    issues <- c(issues, "Major gating tree differences")
  }
  
  if(compatible) {
    cat("Status: Compatible for merging\n")
  } else {
    cat("Status: Issues detected\n")
    cat("Issues:", paste(issues, collapse = ", "), "\n")
  }
  
  return(list(
    compatible = compatible,
    issues = issues,
    sample_overlap = length(overlapping),
    channel_similarity = length(common_channels) / max(length(channels1), length(channels2)) * 100,
    gating_similarity = similarity
  ))
}

# Usage examples:

# Example 1: Basic merge
# merged_gs <- merge_gating_sets(gs1, gs2)

# Example 2: Merge with source tracking
# merged_gs <- merge_gating_sets_with_tracking(gs1, gs2, "Experiment1", "Experiment2")

# Example 3: Check compatibility first
# compatibility <- check_gating_set_compatibility(gs1, gs2)
# if(compatibility$compatible) {
#   merged_gs <- merge_gating_sets(gs1, gs2)
# }

# Example 4: Merge multiple GatingSets
# gs_list <- list(gs1, gs2, gs3)
# merged_gs <- merge_multiple_gating_sets(gs_list)

# ============================================================================
# NODE SELECTION AND RESOLUTION
# ============================================================================

# Get all unique population paths
get_all_populations <- function(gs) {
  unique(unlist(map(gs, gh_get_pop_paths)))
}

# Interactive node selection with regex support and back navigation
select_nodes_interactive <- function(gs) {
  pops <- get_all_populations(gs)
  
  while(TRUE) {
    cat("\n=== Node Selection ===\n")
    cat("1. Select from list (interactive)\n")
    cat("2. Use regex pattern\n") 
    cat("3. Use leaf names\n")
    cat("4. Exit/Cancel\n")
    
    choice <- readline("Choose method (1-4): ")
    
    if(choice == "4") {
      stop("Node selection cancelled by user.")
    }
    
    result <- switch(choice,
                     "1" = {
                       cat("\nAvailable populations:\n")
                       sel <- select.list(pops, multiple = TRUE, title = "Select populations (Cancel to go back)")
                       if(length(sel) == 0) {
                         cat("No populations selected. Going back...\n")
                         next
                       }
                       sel
                     },
                     
                     "2" = {
                       while(TRUE) {
                         cat("\nRegex Pattern Selection\n")
                         pattern <- readline("Enter regex pattern (e.g., 'CD8a_pos', 'CD45\\.1') or 'back': ")
                         
                         if(tolower(pattern) == "back") break
                         if(pattern == "") {
                           cat("Empty pattern. Please try again.\n")
                           next
                         }
                         
                         matches <- grep(pattern, pops, value = TRUE, ignore.case = TRUE)
                         if(length(matches) == 0) {
                           cat("No matches found. Try a different pattern.\n")
                           retry <- readline("Try again? (y/n): ")
                           if(tolower(retry) != "y") break
                           next
                         }
                         
                         cat(sprintf("Found %d matches:\n", length(matches)))
                         iwalk(matches, ~cat(sprintf("%d. %s\n", .y, .x)))
                         
                         while(TRUE) {
                           indices <- readline("Enter numbers (comma-separated), 'all', or 'back': ")
                           if(tolower(indices) == "back") break
                           if(indices == "all") return(matches)
                           
                           tryCatch({
                             selected_idx <- as.numeric(str_trim(str_split(indices, ",")[[1]]))
                             if(any(is.na(selected_idx)) || any(selected_idx < 1) || any(selected_idx > length(matches))) {
                               cat("Invalid selection. Please try again.\n")
                               next
                             }
                             return(matches[selected_idx])
                           }, error = function(e) {
                             cat("Invalid input format. Please try again.\n")
                           })
                         }
                       }
                       NULL
                     },
                     
                     "3" = {
                       while(TRUE) {
                         cat("\nLeaf Name Selection\n")
                         input <- readline("Enter leaf names (comma-separated) or 'back': ")
                         if(tolower(input) == "back") break
                         if(input == "") {
                           cat("No leaf names provided. Please try again.\n")
                           next
                         }
                         
                         leaf_names <- str_trim(str_split(input, ",")[[1]])
                         matches <- resolve_nodes_by_leaf(gs, leaf_names)
                         
                         if(length(matches) == 0) {
                           cat("No matches found for those leaf names.\n")
                           retry <- readline("Try again? (y/n): ")
                           if(tolower(retry) != "y") break
                           next
                         }
                         
                         cat(sprintf("Found %d matches:\n", length(matches)))
                         iwalk(matches, ~cat(sprintf("  %s\n", .x)))
                         
                         confirm <- readline("Use these matches? (y/n): ")
                         if(tolower(confirm) == "y") return(matches)
                       }
                       NULL
                     },
                     
                     {
                       cat("Invalid choice. Please select 1-4.\n")
                       NULL
                     }
    )
    
    if(!is.null(result)) return(result)
  }
}

# Resolve nodes by leaf name
resolve_nodes_by_leaf <- function(gs, leaf_names, exact = TRUE) {
  pops <- get_all_populations(gs)
  
  get_leaf <- function(path) {
    parts <- str_split(path, "/")[[1]]
    tail(parts[parts != ""], 1)
  }
  
  leaves <- map_chr(pops, get_leaf)
  
  if(exact) {
    pops[leaves %in% leaf_names]
  } else {
    matches <- map(leaf_names, ~pops[str_detect(leaves, .x, ignore_case = TRUE)]) %>%
      set_names(leaf_names)
    unlist(matches, use.names = FALSE)
  }
}

# Get hierarchical parent of a node
get_hierarchical_parent <- function(node_path) {
  # Split the path by "/" and remove empty elements
  parts <- str_split(node_path, "/")[[1]]
  parts <- parts[parts != ""]
  
  # If there's only one part or it's the root, return "/"
  if(length(parts) <= 1) {
    return("/")
  }
  
  # Remove the last part to get the parent
  parent_parts <- parts[-length(parts)]
  
  # Reconstruct the parent path
  parent_path <- paste(parent_parts, collapse = "/")
  
  # Add leading "/" if it doesn't start with one
  if(!str_starts(parent_path, "/")) {
    parent_path <- paste0("/", parent_path)
  }
  
  return(parent_path)
}
# ============================================================================
# PARENT NODE MANAGEMENT
# ============================================================================

# Interactive parent selection with back navigation
select_parents_interactive <- function(nodes, gs) {
  # local helper: extract leaf name from a path
  get_leaf <- function(path) {
    parts <- str_split(path, "/")[[1]]
    tail(parts[parts != ""], 1)
  }
  
  # local helper: prefix before underscore
  leaf_prefix <- function(leaf) {
    sub("^([^_]+).*", "\1", leaf)
  }
  
  while (TRUE) {
    cat("\n=== Parent Node Selection ===\n")
    cat("Selected nodes:\n")
    iwalk(nodes, ~cat(sprintf("  %d. %s\n", .y, basename(.x))))
    
    cat("\nAvailable methods:\n")
    cat("1. Auto (hierarchical parent)\n")
    cat("2. Same parent for all nodes\n")
    cat("3. Individual parent for each node\n")
    cat("4. Congenic-specific (CD45.1, CD45.2, CD90.1, CD90.2, CD45.1.2)\n")
    cat("5. Back to node selection\n")
    
    method <- readline("Choose method (1-5): ")
    
    if (method == "5") {
      return("BACK")
    }
    
    result <- switch(method,
                     "1" = {
                       # Auto hierarchical
                       cat("\nUsing automatic hierarchical parent detection...\n")
                       mapping <- set_names(map_chr(nodes, get_hierarchical_parent), nodes)
                       list(mapping = mapping, method = "auto")
                     },
                     
                     "2" = {
                       # Same parent for all
                       while (TRUE) {
                         cat("\nSelect common parent for all nodes:\n")
                         all_pops <- get_all_populations(gs)
                         parent <- select.list(all_pops, multiple = FALSE, title = "Select common parent (Cancel to go back)")
                         
                         if (length(parent) == 0) {
                           cat("No parent selected.\n")
                           retry <- readline("Try again? (y/n/back): ")
                           if (tolower(retry) == "back") break
                           if (tolower(retry) != "y") break
                           next
                         }
                         
                         mapping <- set_names(rep(parent, length(nodes)), nodes)
                         return(list(mapping = mapping, method = "common"))
                       }
                       NULL
                     },
                     
                     "3" = {
                       # Individual parents
                       all_pops <- get_all_populations(gs)
                       mapping <- character(length(nodes))
                       names(mapping) <- nodes
                       
                       for (i in seq_along(nodes)) {
                         node <- nodes[i]
                         
                         while (TRUE) {
                           cat(sprintf("\n--- Node %d/%d: %s ---\n", i, length(nodes), basename(node)))
                           
                           # Filter to likely parents (hierarchically related)
                           likely_parents <- all_pops[
                             map_lgl(all_pops, ~ str_detect(node, fixed(.x))) &
                               nchar(all_pops) < nchar(node) &
                               all_pops != node
                           ]
                           
                           if (length(likely_parents) > 0) {
                             cat("Suggested parents:\n")
                             iwalk(likely_parents, ~cat(sprintf("%d. %s\n", .y, .x)))
                             cat(sprintf("%d. Choose from all populations\n", length(likely_parents) + 1))
                             cat(sprintf("%d. Back to previous node\n", length(likely_parents) + 2))
                             
                             choice <- readline("Select parent: ")
                             
                             if (choice == as.character(length(likely_parents) + 2)) {
                               if (i == 1) {
                                 cat("Already at first node. Going back to method selection.\n")
                                 return("BACK_TO_METHOD")
                               }
                               # Go back to previous node
                               i <- i - 1
                               break
                             } else if (choice == as.character(length(likely_parents) + 1)) {
                               # Choose from all populations
                               parent <- select.list(all_pops, title = paste("Parent for", basename(node)))
                               if (length(parent) == 0) next
                               mapping[node] <- parent
                               break
                             } else {
                               choice_num <- as.numeric(choice)
                               if (!is.na(choice_num) && choice_num >= 1 && choice_num <= length(likely_parents)) {
                                 mapping[node] <- likely_parents[choice_num]
                                 break
                               } else {
                                 cat("Invalid choice. Please try again.\n")
                               }
                             }
                           } else {
                             cat("No hierarchically related parents found.\n")
                             cat("1. Choose from all populations\n")
                             cat("2. Back to previous node\n")
                             
                             choice <- readline("Select option: ")
                             if (choice == "2") {
                               if (i == 1) {
                                 cat("Already at first node. Going back to method selection.\n")
                                 return("BACK_TO_METHOD")
                               }
                               i <- i - 1
                               break
                             } else if (choice == "1") {
                               parent <- select.list(all_pops, title = paste("Parent for", basename(node)))
                               if (length(parent) == 0) next
                               mapping[node] <- parent
                               break
                             }
                           }
                         }
                       }
                       
                       list(mapping = mapping, method = "individual")
                     },
                     
                     "4" = {
                       # Congenic-specific (updated)
                       cat("\nUsing congenic-specific mapping...\n")
                       candidates <- c("CD45.1", "CD45.2", "CD90.1", "CD90.2", "CD45.1.2")
                       all_pops <- get_all_populations(gs)
                       
                       # Precompute leaves and prefixes
                       pop_leaves <- map_chr(all_pops, get_leaf)
                       pop_prefixes <- map_chr(pop_leaves, leaf_prefix)
                       
                       mapping <- set_names(character(length(nodes)), nodes)
                       
                       for (node in nodes) {
                         node_leaf <- get_leaf(node)
                         node_pref <- leaf_prefix(node_leaf)
                         
                         matches <- character(0)
                         
                         # Case 1: node leaf prefix directly equals one of the congenic candidates
                         if (node_pref %in% candidates) {
                           matches <- all_pops[pop_prefixes == node_pref]
                         }
                         
                         # Case 2: node path contains a congenic string anywhere (e.g. '/.../CD45.1_KO/...')
                         if (length(matches) == 0) {
                           found_cands <- candidates[map_lgl(candidates, ~ str_detect(node, fixed(.x)) || str_detect(node, fixed(paste0(.x, "_"))))]
                           if (length(found_cands) > 0) {
                             matches <- all_pops[pop_prefixes %in% found_cands]
                           }
                         }
                         
                         # If multiple matches prefer the most specific (longest path)
                         if (length(matches) > 0) {
                           mapping[node] <- matches[which.max(nchar(matches))]
                         } else {
                           # Fallback to hierarchical parent detection
                           mapping[node] <- get_hierarchical_parent(node)
                         }
                       }
                       
                       list(mapping = mapping, method = "congenic")
                     },
                     
                     {
                       cat("Invalid choice. Please select 1-5.\n")
                       NULL
                     }
    )
    
    # Handle results with preview
    if (!is.null(result)) {
      if (is.character(result) && length(result) == 1 && result == "BACK_TO_METHOD") {
        next
      }
      
      if (is.list(result) && !is.null(result$mapping)) {
        # Show preview and get user confirmation
        preview_result <- preview_parent_mapping(nodes, result$mapping, gs)
        
        if (preview_result == "ACCEPT") {
          cat("\nParent mapping accepted!\n")
          return(result$mapping)
        } else if (preview_result == "BACK") {
          cat("\nReturning to parent selection method...\n")
          next
        } else if (preview_result == "CANCEL") {
          cat("\nParent selection cancelled.\n")
          return("BACK")
        }
      }
      
      if (is.character(result) && length(result) > 1) {
        # For direct mapping returns, also show preview
        preview_result <- preview_parent_mapping(nodes, result, gs)
        
        if (preview_result == "ACCEPT") {
          return(result)
        } else if (preview_result == "BACK") {
          next
        } else if (preview_result == "CANCEL") {
          return("BACK")
        }
      }
    }
  }
}

# Preview parent mapping with back navigation
preview_parent_mapping <- function(nodes, parent_mapping, gs) {
  all_pops <- get_all_populations(gs)
  
  while(TRUE) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("PARENT MAPPING PREVIEW\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    
    iwalk(parent_mapping, function(parent, node) {
      cat(sprintf("\nNode: %s\n", basename(node)))
      cat(sprintf("  Full path: %s\n", node))
      cat(sprintf("  Parent: %s\n", parent))
      cat(sprintf("  Node exists: %s\n", ifelse(node %in% all_pops, "✓", "✗")))
      cat(sprintf("  Parent exists: %s\n", ifelse(parent %in% all_pops, "✓", "✗")))
      
      if(parent %in% all_pops && node %in% all_pops) {
        is_hierarchical <- str_detect(node, fixed(parent)) && parent != node
        cat(sprintf("  Hierarchical: %s\n", ifelse(is_hierarchical, "✓", "✗")))
      }
    })
    
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat("Options:\n")
    cat("1. Accept mapping and continue\n")
    cat("2. Go back to modify parent selection\n")
    cat("3. Cancel analysis\n")
    
    choice <- readline("Choose (1-3): ")
    
    switch(choice,
           "1" = return("ACCEPT"),
           "2" = return("BACK"), 
           "3" = return("CANCEL"),
           {
             cat("Invalid choice. Please select 1-3.\n")
           }
    )
  }
}

# ============================================================================
# DATA EXTRACTION
# ============================================================================

# Get marker lookup table
get_marker_lookup <- function(gs) {
  fh <- gh_pop_get_data(gs[[1]], "/")
  params <- pData(parameters(fh))
  
  tibble(
    colname = params$name,
    marker = if_else(is.na(params$desc) | params$desc == "", params$name, params$desc)
  )
}

# Select channels for MFI analysis with back navigation
select_channels_interactive <- function(gs) {
  lookup <- get_marker_lookup(gs)
  
  # Identify compensated channels
  comp_channels <- lookup %>%
    dplyr::filter(str_detect(colname, "(?i)comp")) %>%
    pull(colname)
  
  while(TRUE) {
    cat("\n=== Channel Selection for MFI Analysis ===\n")
    cat("Available channels:", nrow(lookup), "\n")
    if(length(comp_channels) > 0) {
      cat("Compensated channels found:", length(comp_channels), "\n")
    }
    
    cat("\nOptions:\n")
    cat("1. Select specific channels\n")
    cat("2. Use all channels\n")
    if(length(comp_channels) > 0) {
      cat("3. Use all compensated channels (containing 'Comp')\n")
      cat("4. Skip MFI analysis\n")
      cat("5. Back\n")
    } else {
      cat("3. Skip MFI analysis\n")
      cat("4. Back\n")
    }
    
    choice <- readline(paste0("Choose (1-", if(length(comp_channels) > 0) "5" else "4", "): "))
    
    if(length(comp_channels) > 0) {
      # Menu with compensated channels option
      switch(choice,
             "1" = {
               choices <- paste0(lookup$colname, " :: ", lookup$marker)
               sel <- select.list(choices, multiple = TRUE, title = "Select channels for MFI (Cancel to go back)")
               
               if(length(sel) == 0) {
                 cat("No channels selected.\n")
                 retry <- readline("Try again? (y/n): ")
                 if(tolower(retry) != "y") next
                 next
               }
               
               selected_channels <- map_chr(sel, ~str_split(.x, " :: ")[[1]][1])
               
               cat(sprintf("Selected %d channels:\n", length(selected_channels)))
               iwalk(selected_channels, ~cat(sprintf("  %d. %s\n", .y, .x)))
               
               confirm <- readline("Confirm selection? (y/n): ")
               if(tolower(confirm) == "y") return(selected_channels)
             },
             
             "2" = {
               cat("Using all available channels for MFI analysis.\n")
               return("all")
             },
             
             "3" = {
               cat(sprintf("Using %d compensated channels:\n", length(comp_channels)))
               comp_lookup <- lookup %>% dplyr::filter(colname %in% comp_channels)
               iwalk(comp_lookup$colname, ~cat(sprintf("  %d. %s :: %s\n", .y, .x, 
                                                       comp_lookup$marker[comp_lookup$colname == .x])))
               
               confirm <- readline("Confirm compensated channels? (y/n): ")
               if(tolower(confirm) == "y") return(comp_channels)
             },
             
             "4" = {
               cat("Skipping MFI analysis.\n")
               return(character(0))
             },
             
             "5" = {
               return("BACK")
             },
             
             {
               cat("Invalid choice. Please select 1-5.\n")
             }
      )
    } else {
      # Menu without compensated channels option
      switch(choice,
             "1" = {
               choices <- paste0(lookup$colname, " :: ", lookup$marker)
               sel <- select.list(choices, multiple = TRUE, title = "Select channels for MFI (Cancel to go back)")
               
               if(length(sel) == 0) {
                 cat("No channels selected.\n")
                 retry <- readline("Try again? (y/n): ")
                 if(tolower(retry) != "y") next
                 next
               }
               
               selected_channels <- map_chr(sel, ~str_split(.x, " :: ")[[1]][1])
               
               cat(sprintf("Selected %d channels:\n", length(selected_channels)))
               iwalk(selected_channels, ~cat(sprintf("  %d. %s\n", .y, .x)))
               
               confirm <- readline("Confirm selection? (y/n): ")
               if(tolower(confirm) == "y") return(selected_channels)
             },
             
             "2" = {
               cat("Using all available channels for MFI analysis.\n")
               return("all")
             },
             
             "3" = {
               cat("Skipping MFI analysis.\n")
               return(character(0))
             },
             
             "4" = {
               return("BACK")
             },
             
             {
               cat(sprintf("Invalid choice. Please select 1-%s.\n", if(length(comp_channels) > 0) "5" else "4"))
             }
      )
    }
  }
}

# Extract counts and frequencies
extract_counts_freqs <- function(gs, nodes, parent_mapping = NULL, keywords = c("pairing_factor", "tissue_factor")) {
  
  # Helper to safely get counts for a node in a sample
  safe_get_count <- function(gh, node) {
    tryCatch({
      if(node %in% gh_get_pop_paths(gh)) {
        list(
          count = gh_pop_get_count(gh, node),
          parent_count = gh_pop_get_count(gh, gh_pop_get_parent(gh, node)),
          success = TRUE
        )
      } else {
        list(success = FALSE)
      }
    }, error = function(e) list(success = FALSE))
  }
  
  # Extract counts for all node-sample combinations
  results <- map_dfr(nodes, function(node) {
    map_dfr(seq_along(gs), function(i) {
      sample_name <- sampleNames(gs)[i]
      gh <- gs[[i]]
      
      result <- safe_get_count(gh, node)
      if(!result$success) return(tibble())
      
      tibble(
        Sample = sample_name,
        Node = node,
        Count = result$count,
        ParentCount = result$parent_count
      )
    })
  })
  
  if(nrow(results) == 0) {
    stop("No counts extracted. Check that nodes exist in the gating set.")
  }
  
  # Join with sample metadata
  pd <- pData(gs) %>% 
    rownames_to_column("Sample") %>%
    select(Sample, all_of(keywords))
  
  results <- results %>%
    left_join(pd, by = "Sample") %>%
    rename(Subpop = Count)
  
  # Calculate frequencies with custom parents if provided
  if(!is.null(parent_mapping)) {
    # Get custom parent counts
    parent_results <- map_dfr(names(parent_mapping), function(node) {
      parent_node <- parent_mapping[node]
      
      map_dfr(seq_along(gs), function(i) {
        sample_name <- sampleNames(gs)[i]
        gh <- gs[[i]]
        
        result <- safe_get_count(gh, parent_node)
        if(!result$success) return(tibble())
        
        tibble(
          Sample = sample_name,
          Node = node,
          CustomParentCount = result$count
        )
      })
    })
    
    results <- results %>%
      left_join(parent_results, by = c("Sample", "Node")) %>%
      mutate(
        Freq = case_when(
          !is.na(CustomParentCount) & CustomParentCount > 0 ~ (Subpop / CustomParentCount) * 100,
          ParentCount > 0 ~ (Subpop / ParentCount) * 100,
          TRUE ~ 0
        )
      )
  } else {
    results <- results %>%
      mutate(Freq = if_else(ParentCount > 0, (Subpop / ParentCount) * 100, 0))
  }
  
  results %>% mutate(NodeShort = basename(Node))
}

# Extract MFI data
extract_mfi <- function(gs, nodes, channels = NULL, summary_fun = median, keywords = c("pairing_factor", "tissue_factor")) {
  if(is.null(channels)) {
    channels <- select_channels_interactive(gs)
    if(length(channels) == 0) return(tibble())
  } else if(length(channels) == 1 && tolower(channels) == "all") {
    channels <- get_marker_lookup(gs)$colname
  } else if(length(channels) == 1 && tolower(channels) == "compensated") {
    # Handle compensated channels selection
    lookup <- get_marker_lookup(gs)
    channels <- lookup %>%
      dplyr::filter(str_detect(colname, "(?i)comp")) %>%
      pull(colname)
    
    if(length(channels) == 0) {
      warning("No compensated channels found. Skipping MFI analysis.")
      return(tibble())
    }
  }
  
  # Validate channels
  lookup <- get_marker_lookup(gs)
  missing_channels <- setdiff(channels, lookup$colname)
  if(length(missing_channels) > 0) {
    stop("Channels not found: ", paste(missing_channels, collapse = ", "))
  }
  
  # Extract MFI for all node-sample-channel combinations
  results <- map_dfr(nodes, function(node) {
    map_dfr(seq_along(gs), function(i) {
      sample_name <- sampleNames(gs)[i]
      gh <- gs[[i]]
      
      # Check if node exists in this sample
      if(!node %in% gh_get_pop_paths(gh)) return(tibble())
      
      # Get data
      ff <- tryCatch({
        gh_pop_get_data(gh, node)
      }, error = function(e) return(NULL))
      
      if(is.null(ff)) return(tibble())
      
      # Extract expression data
      expr_data <- if(inherits(ff, "cytoframe")) exprs(ff) else exprs(ff)
      df <- as.data.frame(expr_data)
      
      # Calculate MFI for available channels
      available_channels <- intersect(channels, colnames(df))
      if(length(available_channels) == 0) return(tibble())
      
      map_dfr(available_channels, function(ch) {
        tibble(
          Sample = sample_name,
          Node = node,
          colname = ch,
          MFI = summary_fun(df[[ch]], na.rm = TRUE)
        )
      })
    })
  })
  
  if(nrow(results) == 0) return(tibble())
  
  # Add marker information
  results <- results %>%
    left_join(lookup, by = "colname") %>%
    select(Sample, Node, colname, marker, MFI)
  
  # Join with sample metadata (same as in extract_counts_freqs)
  pd <- pData(gs) %>% 
    rownames_to_column("Sample") %>%
    select(Sample, all_of(keywords))
  
  results <- results %>%
    left_join(pd, by = "Sample") %>%
    mutate(NodeShort = basename(Node))
  
  return(results)
}

# ============================================================================
# MAIN ANALYSIS FUNCTION WITH FULL BACK NAVIGATION
# ============================================================================



# ---------------------------------------------------------------------------
# CONGENIC EXTRACTION HELPERS
# ---------------------------------------------------------------------------
# Returns the best-matching congenic candidate found in any path segment of `node`.
# Matching checks the substring before the first '_' in each path segment and
# compares it to the candidate list. If multiple candidates match, the longest
# (most specific) candidate is returned. If none match, NA is returned.
get_candidate_from_node <- function(node, candidates) {
  parts <- str_split(node, "/")[[1]]
  parts <- parts[parts != ""]
  if(length(parts) == 0) return(NA_character_)
  
  prefixes <- map_chr(parts, ~ sub("^([^_]+).*", "\\1", .x))
  
  # find intersection while preserving specificity (prefer longest match)
  hits <- intersect(prefixes, candidates)
  if(length(hits) == 0) return(NA_character_)
  hits[which.max(nchar(hits))]
}

# Adds a `congenics` column to results$counts and results$mfi (if present)
# based on Node path. Default candidates include the congenics used elsewhere.
add_congenics_column <- function(results, candidates = c("CD45.1", "CD45.2", "CD90.1", "CD90.2", "CD45.1.2")) {
  if(!is.null(results$counts) && nrow(results$counts) > 0) {
    results$counts <- results$counts %>%
      mutate(congenics = map_chr(Node, ~ get_candidate_from_node(.x, candidates)))
  }
  
  if(!is.null(results$mfi) && nrow(results$mfi) > 0) {
    results$mfi <- results$mfi %>%
      mutate(congenics = map_chr(Node, ~ get_candidate_from_node(.x, candidates)))
  }
  
  results
}

# Example: call this after running analyze_flow_data()
# results <- analyze_flow_data(gs)
# results <- add_congenics_column(results)

#==============================================================================
# Additional data metadata annotation for downstream analysis
#==============================================================================
assign_metadata_menu <- function(df) {
  
  # Check what congenic markers are present
  available_markers <- unique(df$congenics)
  available_markers <- available_markers[!is.na(available_markers)]
  
  if(length(available_markers) < 1) {
    stop("No congenic markers found in the data")
  }
  
  # WT selection menu
  cat("=== Select Wild Type (WT) Marker ===\n")
  wt_choice <- menu(available_markers, title = "Choose WT marker:")
  if(wt_choice == 0) stop("WT marker selection cancelled")
  wt_marker <- available_markers[wt_choice]
  
  # KO selection menu  
  remaining_markers <- available_markers[-wt_choice]
  ko_marker <- NULL
  
  if(length(remaining_markers) > 0) {
    cat("\n=== Select Knock Out (KO) Marker (Optional) ===\n")
    ko_options <- c(remaining_markers, "Skip KO assignment")
    ko_choice <- menu(ko_options, title = "Choose KO marker:")
    
    if(ko_choice > 0 && ko_choice <= length(remaining_markers)) {
      ko_marker <- remaining_markers[ko_choice]
      remaining_markers <- remaining_markers[-ko_choice]
    }
  }
  
  # Recipient selection menu
  recipient_marker <- NULL
  
  if(length(remaining_markers) > 0) {
    cat("\n=== Select Recipient Marker (Optional) ===\n")
    recipient_options <- c(remaining_markers, "Skip Recipient assignment")
    recipient_choice <- menu(recipient_options, title = "Choose Recipient marker:")
    
    if(recipient_choice > 0 && recipient_choice <= length(remaining_markers)) {
      recipient_marker <- remaining_markers[recipient_choice]
    }
  }
  
  # Create genotype column
  df_updated <- df %>%
    mutate(genotype = case_when(
      congenics == wt_marker ~ "donor_WT",
      TRUE ~ "Other"
    ))
  
  # Add KO assignments if ko_marker was selected
  if(!is.null(ko_marker)) {
    df_updated <- df_updated %>%
      mutate(genotype = case_when(
        congenics == ko_marker ~ "donor_KO",
        TRUE ~ genotype
      ))
  }
  
  # Add Recipient assignments if recipient_marker was selected
  if(!is.null(recipient_marker)) {
    df_updated <- df_updated %>%
      mutate(genotype = case_when(
        congenics == recipient_marker ~ "Recipient",
        TRUE ~ genotype
      ))
  }
  
  # Display genotype results
  cat("\n=== Genotype Assignment Complete ===\n")
  cat(paste("WT:", wt_marker, "\n"))
  if(!is.null(ko_marker)) cat(paste("KO:", ko_marker, "\n"))
  if(!is.null(recipient_marker)) cat(paste("Recipient:", recipient_marker, "\n"))
  cat("\nGenotype counts:\n")
  print(table(df_updated$genotype, useNA = "ifany"))
  
  # Metadata addition section
  cat("\n=== Add Custom Metadata ===\n")
  
  metadata_complete <- FALSE
  
  while(!metadata_complete) {
    metadata_options <- c(
      "Add timepoint",
      "Add batch",
      "Add sex",
      "Add custom metadata column",
      "Finish and return data"
    )
    
    cat("\n--- Metadata Menu ---\n")
    metadata_choice <- menu(metadata_options, title = "Choose metadata option:")
    
    if(metadata_choice == 0) {
      cat("Metadata addition cancelled. Returning data with genotype assignments only.\n")
      break
    }
    
    if(metadata_choice == 1) {
      # Add timepoint
      cat("\n=== Add Timepoint ===\n")
      timepoint_value <- readline(prompt = "Enter timepoint value: ")
      if(timepoint_value != "") {
        df_updated <- df_updated %>%
          mutate(timepoint = timepoint_value)
        cat(paste("Added timepoint:", timepoint_value, "\n"))
      } else {
        cat("Timepoint addition skipped (empty input).\n")
      }
      
    } else if(metadata_choice == 2) {
      # Add batch
      cat("\n=== Add Batch ===\n")
      batch_value <- readline(prompt = "Enter batch value: ")
      if(batch_value != "") {
        df_updated <- df_updated %>%
          mutate(batch = batch_value)
        cat(paste("Added batch:", batch_value, "\n"))
      } else {
        cat("Batch addition skipped (empty input).\n")
      }
      
    } else if(metadata_choice == 3) {
      # Add sex
      cat("\n=== Add Sex ===\n")
      sex_options <- c("m", "f", "Enter custom value")
      sex_choice <- menu(sex_options, title = "Choose sex:")
      
      if(sex_choice == 0) {
        cat("Sex addition cancelled.\n")
      } else if(sex_choice == 1) {
        df_updated <- df_updated %>%
          mutate(sex = "m")
        cat("Added sex: m\n")
      } else if(sex_choice == 2) {
        df_updated <- df_updated %>%
          mutate(sex = "f")
        cat("Added sex: f\n")
      } else if(sex_choice == 3) {
        sex_value <- readline(prompt = "Enter custom sex value: ")
        if(sex_value != "") {
          df_updated <- df_updated %>%
            mutate(sex = sex_value)
          cat(paste("Added sex:", sex_value, "\n"))
        } else {
          cat("Sex addition skipped (empty input).\n")
        }
      }
      
    } else if(metadata_choice == 4) {
      # Add custom metadata column
      cat("\n=== Add Custom Metadata Column ===\n")
      column_name <- readline(prompt = "Enter column name: ")
      
      if(column_name != "" && !column_name %in% names(df_updated)) {
        column_value <- readline(prompt = paste("Enter value for", column_name, ": "))
        if(column_value != "") {
          df_updated <- df_updated %>%
            mutate(!!sym(column_name) := column_value)
          cat(paste("Added", column_name, ":", column_value, "\n"))
        } else {
          cat("Custom column addition skipped (empty value).\n")
        }
      } else if(column_name == "") {
        cat("Custom column addition skipped (empty column name).\n")
      } else {
        cat("Column name already exists. Please choose a different name.\n")
      }
      
    } else if(metadata_choice == 5) {
      # Finish and return data
      metadata_complete <- TRUE
    }
  }
  
  # Display final results
  cat("\n=== Final Results ===\n")
  cat("Columns in dataset:\n")
  cat(paste(names(df_updated), collapse = ", "), "\n")
  cat("\nDataset summary:\n")
  cat(paste("Rows:", nrow(df_updated), "\n"))
  cat(paste("Columns:", ncol(df_updated), "\n"))
  
  return(df_updated)
}

# ============================================================================
# INTERACTIVE FACTOR DEFINITION FOR GENERALIZABLE ANALYSIS
# ============================================================================

#' Interactive definition of tissue and pairing factors
#' This function allows users to specify which columns to use for tissue grouping
#' and paired analysis, making the code more generalizable
define_analysis_factors <- function(df) {
  
  cat("=== Analysis Factor Definition ===\n")
  cat("This step will define the key variables for your analysis:\n")
  cat("• Tissue Factor: Groups samples by tissue/treatment/condition\n")
  cat("• Pairing Factor: Identifies paired samples for statistical analysis\n\n")
  
  # Show available columns
  available_cols <- names(df)
  cat("Available columns in your data:\n")
  iwalk(available_cols, ~cat(sprintf("  %d. %s\n", .y, .x)))
  
  # Preview data
  cat("\nData preview (first 5 rows):\n")
  print(head(df %>% select(1:min(6, ncol(df))), 5))
  
  # ==========================================================================
  # TISSUE FACTOR DEFINITION
  # ==========================================================================
  
  tissue_factor_col <- NULL
  
  while(is.null(tissue_factor_col)) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("TISSUE FACTOR DEFINITION\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    cat("The tissue factor groups your samples for analysis.\n")
    cat("Common examples:\n")
    cat("• tissue_factor (tissue types: Spleen, Liver, etc.)\n")
    cat("• Treatment (Control, Drug_A, Drug_B, etc.)\n")
    cat("• Timepoint (Day_0, Day_7, Day_14, etc.)\n")
    cat("• Genotype (WT, KO, etc.)\n\n")
    
    # Show unique values for potential tissue factors with safer processing
    cat("Column preview with unique values:\n")
    for(i in seq_along(available_cols)) {
      col_name <- available_cols[i]
      
      # Safely get unique values
      tryCatch({
        col_data <- df[[col_name]]
        
        # Handle different data types safely
        if(is.list(col_data)) {
          cat(sprintf("%2d. %-20s | List column (not suitable for grouping)\n", i, col_name))
          next
        }
        
        # Convert to character to handle factors and other types
        col_data_char <- as.character(col_data)
        unique_vals <- unique(col_data_char[!is.na(col_data_char)])
        
        # Only show if reasonable number of unique values (likely categorical)
        if(length(unique_vals) <= 20 && length(unique_vals) > 1) {
          cat(sprintf("%2d. %-20s | %d groups: %s\n", 
                      i, col_name, length(unique_vals),
                      paste(head(unique_vals, 4), collapse = ", ")))
          if(length(unique_vals) > 4) cat(sprintf("%-24s   ... and %d more\n", "", length(unique_vals) - 4))
        } else if(length(unique_vals) == 1) {
          cat(sprintf("%2d. %-20s | 1 group: %s\n", i, col_name, unique_vals[1]))
        } else {
          cat(sprintf("%2d. %-20s | %d unique values (likely continuous)\n", i, col_name, length(unique_vals)))
        }
        
      }, error = function(e) {
        cat(sprintf("%2d. %-20s | Error processing column\n", i, col_name))
      })
    }
    
    tissue_choice <- readline("\nEnter column number or name for TISSUE FACTOR: ")
    
    # Parse choice
    selected_tissue_factor <- NULL
    if(grepl("^\\d+$", tissue_choice)) {
      choice_num <- as.numeric(tissue_choice)
      if(!is.na(choice_num) && choice_num >= 1 && choice_num <= length(available_cols)) {
        selected_tissue_factor <- available_cols[choice_num]
      }
    } else if(tissue_choice %in% available_cols) {
      selected_tissue_factor <- tissue_choice
    }
    
    if(is.null(selected_tissue_factor)) {
      cat("❌ Invalid selection. Please try again.\n")
      next
    }
    
    # Validate the selected column
    selected_col_data <- df[[selected_tissue_factor]]
    
    if(is.list(selected_col_data)) {
      cat("❌ Selected column contains list data and cannot be used for grouping. Please choose another column.\n")
      next
    }
    
    # Show tissue factor summary with safer processing
    tryCatch({
      # Get the selected column data and check its structure more thoroughly
      selected_data <- df[[selected_tissue_factor]]
      
      # More robust validation
      if(is.null(selected_data)) {
        stop("Selected column is NULL")
      }
      
      # Handle different data types more safely
      if(is.data.frame(selected_data) || is.matrix(selected_data) || 
         (is.list(selected_data) && !is.factor(selected_data))) {
        stop("Selected column contains complex data structure (not a simple vector)")
      }
      
      # Convert to character more safely
      if(is.factor(selected_data)) {
        char_data <- as.character(selected_data)
      } else if(is.numeric(selected_data) || is.logical(selected_data)) {
        char_data <- as.character(selected_data)
      } else if(is.character(selected_data)) {
        char_data <- selected_data
      } else {
        # Try to coerce to character, but wrap in additional error handling
        char_data <- tryCatch({
          as.character(selected_data)
        }, error = function(e) {
          stop("Cannot convert selected column to character: ", e$message)
        })
      }
      
      # Filter out NA and empty values
      valid_data <- char_data[!is.na(char_data) & char_data != "" & char_data != "NA"]
      
      if(length(valid_data) == 0) {
        stop("No valid (non-empty, non-NA) values found in selected column")
      }
      
      # Create summary using base R to avoid dplyr issues
      unique_vals <- unique(valid_data)
      counts <- table(valid_data)
      
      # Create a simple summary
      temp_df <- data.frame(
        group = names(counts),
        n_samples = as.numeric(counts),
        stringsAsFactors = FALSE
      )
      temp_df <- temp_df[order(temp_df$n_samples, decreasing = TRUE), ]
      names(temp_df)[1] <- selected_tissue_factor
      
      cat(sprintf("\n✅ Selected tissue factor: %s\n", selected_tissue_factor))
      cat("Groups found:\n")
      print(temp_df)
      
      confirm <- readline("Confirm this tissue factor? (y/n): ")
      if(tolower(confirm) == "y") {
        tissue_factor_col <- selected_tissue_factor
      } else {
        cat("Selection cancelled. Choose again.\n")
      }
      
    }, error = function(e) {
      cat("❌ Error processing selected column:", e$message, "\n")
      cat("Please choose a different column.\n")
    })
  }
  
  # ==========================================================================
  # PAIRING FACTOR DEFINITION
  # ==========================================================================
  
  pairing_factor_col <- NULL
  
  while(is.null(pairing_factor_col)) {
    cat("\n", paste(rep("=", 70), collapse = ""), "\n")
    cat("PAIRING FACTOR DEFINITION\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    
    cat("The pairing factor identifies which samples come from the same subject/mouse.\n")
    cat("This is crucial for paired statistical analyses.\n\n")
    
    cat("Pairing options:\n")
    cat("1. Use existing column (e.g., Mouse_ID, Subject_ID)\n")
    cat("2. Extract from WELL ID (row letter = same mouse)\n")
    cat("3. Create custom pairing based on pattern\n")
    cat("4. No pairing (unpaired analysis only)\n")
    
    pairing_option <- readline("Choose pairing option (1-4): ")
    
    if(pairing_option == "1") {
      # Use existing column - with safer processing
      cat("\nColumns that might contain pairing information:\n")
      
      pairing_candidates <- character(0)
      counter <- 1
      
      for(col_name in available_cols) {
        tryCatch({
          col_data <- df[[col_name]]
          
          # Skip list columns
          if(is.list(col_data)) {
            next
          }
          
          # Convert to character safely
          col_data_char <- as.character(col_data)
          unique_vals <- unique(col_data_char[!is.na(col_data_char) & col_data_char != ""])
          n_unique <- length(unique_vals)
          n_total <- nrow(df)
          
          # Good pairing candidates: fewer unique values than total rows, but more than 1
          if(n_unique > 1 && n_unique < n_total && n_unique >= (n_total * 0.1)) {
            avg_samples_per_group <- round(n_total / n_unique, 1)
            cat(sprintf("%2d. %-20s | %d groups, avg %.1f samples/group\n", 
                        counter, col_name, n_unique, avg_samples_per_group))
            pairing_candidates <- c(pairing_candidates, col_name)
            counter <- counter + 1
          }
        }, error = function(e) {
          # Skip problematic columns silently
        })
      }
      
      if(length(pairing_candidates) == 0) {
        cat("No obvious pairing candidates found. Showing all columns:\n")
        iwalk(available_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
        pairing_candidates <- available_cols
      }
      
      pair_choice <- readline("Enter column number or name for PAIRING: ")
      
      # Parse pairing choice
      if(grepl("^\\d+$", pair_choice)) {
        choice_num <- as.numeric(pair_choice)
        if(!is.na(choice_num) && choice_num >= 1 && choice_num <= length(pairing_candidates)) {
          candidate_col <- pairing_candidates[choice_num]
          
          # Validate it's not a list column
          if(!is.list(df[[candidate_col]])) {
            pairing_factor_col <- candidate_col
          } else {
            cat("❌ Selected column contains list data. Please choose another.\n")
          }
        }
      } else if(pair_choice %in% available_cols && !is.list(df[[pair_choice]])) {
        pairing_factor_col <- pair_choice
      }
      
    } else if(pairing_option == "2") {
      # Extract from WELL ID - similar logic as before but with safer processing
      cat("\nWELL ID Pairing Options:\n")
      cat("a. Row letter pairing (A1,A2,A3 = Mouse 1; B1,B2,B3 = Mouse 2)\n")
      cat("b. Custom pattern extraction\n")
      
      # Look for well ID columns
      wellid_cols <- available_cols[grepl("well|WellID|Well_ID|\\$WELLID|\\$LOCATIONID|LOCATIONID", available_cols, ignore.case = TRUE)]
      
      if(length(wellid_cols) > 0) {
        cat("\nFound potential well ID columns:\n")
        iwalk(wellid_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
        
        # Show examples safely
        for(col in wellid_cols) {
          tryCatch({
            if(!is.list(df[[col]])) {
              sample_vals <- head(unique(as.character(df[[col]])), 5)
              cat(sprintf("  %s examples: %s\n", col, paste(sample_vals, collapse = ", ")))
            }
          }, error = function(e) {
            cat(sprintf("  %s: Error reading examples\n", col))
          })
        }
        
        wellid_choice <- readline("Enter well ID column number or name: ")
        
        # Parse well ID choice
        selected_wellid <- NULL
        if(grepl("^\\d+$", wellid_choice)) {
          choice_num <- as.numeric(wellid_choice)
          if(!is.na(choice_num) && choice_num >= 1 && choice_num <= length(wellid_cols)) {
            selected_wellid <- wellid_cols[choice_num]
          }
        } else if(wellid_choice %in% wellid_cols) {
          selected_wellid <- wellid_choice
        }
        
        if(!is.null(selected_wellid) && !is.list(df[[selected_wellid]])) {
          cat("\nPairing method:\n")
          cat("a. Row letter (A1,A2→A; B1,B2→B)\n")
          cat("b. Custom pattern\n")
          
          method_choice <- readline("Choose method (a/b): ")
          
          if(tolower(method_choice) == "a") {
            # Create row letter pairing
            pairing_factor_col <- "pairing_factor_created"
            df <- df %>%
              mutate(pairing_factor_created = str_extract(as.character(.data[[selected_wellid]]), "^[A-H]"))
            
            # Show pairing preview using base R to avoid dplyr issues
            tryCatch({
              pairing_data <- df$pairing_factor_created
              valid_pairing <- pairing_data[!is.na(pairing_data) & pairing_data != ""]
              
              if(length(valid_pairing) > 0) {
                pairing_counts <- table(valid_pairing)
                pairing_preview <- data.frame(
                  pairing_factor_created = names(pairing_counts),
                  n_samples = as.numeric(pairing_counts),
                  stringsAsFactors = FALSE
                )
                pairing_preview <- pairing_preview[order(pairing_preview$pairing_factor_created), ]
                
                cat("\nPairing groups created from row letters:\n")
                print(pairing_preview)
              } else {
                cat("\nWarning: No valid pairing groups created from row letters.\n")
              }
            }, error = function(e) {
              cat("\nWarning: Could not display pairing preview, but pairing factor was created.\n")
            })
          }
        }
        
      } else {
        cat("No well ID columns found. Please choose option 1 or 3.\n")
        next
      }
      
    } else if(pairing_option == "4") {
      # No pairing
      pairing_factor_col <- "none"
      cat("✅ No pairing selected. Only unpaired analyses will be available.\n")
    } else {
      cat("Option 3 (custom pairing) not fully implemented. Please choose 1, 2, or 4.\n")
    }
    
    if(is.null(pairing_factor_col)) {
      cat("❌ Invalid selection. Please try again.\n")
      next
    }
    
    # Confirm pairing factor (unless "none")
    if(pairing_factor_col != "none") {
      # Show pairing summary safely
      if(pairing_factor_col %in% names(df)) {
        tryCatch({
          # Get the pairing column data and validate it
          pairing_data <- df[[pairing_factor_col]]
          
          # More robust validation for pairing data
          if(is.null(pairing_data)) {
            stop("Pairing column is NULL")
          }
          
          if(is.data.frame(pairing_data) || is.matrix(pairing_data) || 
             (is.list(pairing_data) && !is.factor(pairing_data))) {
            stop("Pairing column contains complex data structure")
          }
          
          # Convert to character safely
          if(is.factor(pairing_data)) {
            char_pairing_data <- as.character(pairing_data)
          } else if(is.numeric(pairing_data) || is.logical(pairing_data)) {
            char_pairing_data <- as.character(pairing_data)
          } else if(is.character(pairing_data)) {
            char_pairing_data <- pairing_data
          } else {
            char_pairing_data <- tryCatch({
              as.character(pairing_data)
            }, error = function(e) {
              stop("Cannot convert pairing column to character: ", e$message)
            })
          }
          
          # Filter out invalid values
          valid_pairing_data <- char_pairing_data[!is.na(char_pairing_data) & 
                                                    char_pairing_data != "" & 
                                                    char_pairing_data != "NA"]
          
          if(length(valid_pairing_data) == 0) {
            stop("No valid values found in pairing column")
          }
          
          # Create summary using base R
          pairing_counts <- table(valid_pairing_data)
          temp_pairing_df <- data.frame(
            group = names(pairing_counts),
            n_samples = as.numeric(pairing_counts),
            stringsAsFactors = FALSE
          )
          temp_pairing_df <- temp_pairing_df[order(temp_pairing_df$n_samples, decreasing = TRUE), ]
          names(temp_pairing_df)[1] <- pairing_factor_col
          
          cat(sprintf("\n✅ Selected pairing factor: %s\n", pairing_factor_col))
          cat("Pairing groups:\n")
          print(head(temp_pairing_df, 10))
          if(nrow(temp_pairing_df) > 10) {
            cat("... and", nrow(temp_pairing_df) - 10, "more groups\n")
          }
          
          confirm <- readline("Confirm this pairing factor? (y/n): ")
          if(tolower(confirm) != "y") {
            pairing_factor_col <- NULL
            cat("Selection cancelled. Choose again.\n")
          }
        }, error = function(e) {
          cat("❌ Error processing pairing factor:", e$message, "\n")
          pairing_factor_col <- NULL
        })
      }
    }
  }
  
  # ==========================================================================
  # CREATE FINAL FACTOR COLUMNS
  # ==========================================================================
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CREATING ANALYSIS FACTORS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  # Create tissue_factor column (duplicate and rename for clarity)
  df <- df %>%
    mutate(tissue_factor = as.character(.data[[tissue_factor_col]]))
  
  cat("✅ Created 'tissue_factor' column based on:", tissue_factor_col, "\n")
  
  # Create pairing_factor column if not "none"
  if(pairing_factor_col != "none") {
    df <- df %>%
      mutate(pairing_factor = as.character(.data[[pairing_factor_col]]))
    cat("✅ Created 'pairing_factor' column based on:", pairing_factor_col, "\n")
  } else {
    df <- df %>%
      mutate(pairing_factor = "no_pairing")
    cat("✅ Created 'pairing_factor' column with value: no_pairing\n")
  }
  
  # Final summary with safer processing
  cat("\n=== FACTOR DEFINITION COMPLETE ===\n")
  cat("Analysis factors created:\n")
  
  tryCatch({
    tissue_final <- df %>% 
      dplyr::filter(!is.na(tissue_factor), tissue_factor != "") %>%
      count(tissue_factor, name = "n") %>% 
      arrange(desc(n))
    cat("\nTissue Factor Groups:\n")
    print(tissue_final)
    
    if(pairing_factor_col != "none") {
      pairing_final <- df %>% 
        dplyr::filter(!is.na(pairing_factor), pairing_factor != "") %>%
        count(pairing_factor, name = "n") %>% 
        arrange(desc(n))
      cat("\nPairing Factor Groups:\n")
      print(head(pairing_final, 15))
      if(nrow(pairing_final) > 15) {
        cat("... and", nrow(pairing_final) - 15, "more pairing groups\n")
      }
    } else {
      cat("\nPairing Factor: No pairing (unpaired analyses only)\n")
    }
  }, error = function(e) {
    cat("Warning: Could not display factor summaries, but factors were created.\n")
  })
  
  return(list(
    data = df,
    tissue_factor = "tissue_factor",
    pairing_factor = "pairing_factor",
    original_tissue_col = tissue_factor_col,
    original_pairing_col = if(pairing_factor_col != "none") pairing_factor_col else NULL,
    pairing_enabled = pairing_factor_col != "none"
  ))
}


# Enhanced metadata assignment with factor definition
assign_metadata_menu_enhanced <- function(df, include_factor_definition = TRUE) {
  
  # First, define analysis factors if requested
  if(include_factor_definition) {
    factor_results <- define_analysis_factors(df)
    df <- factor_results$data
    cat("\nFactor definition complete. Continuing with metadata assignment...\n")
  }
  
  # Check what congenic markers are present
  available_markers <- unique(df$congenics)
  available_markers <- available_markers[!is.na(available_markers)]
  
  if(length(available_markers) < 1) {
    stop("No congenic markers found in the data")
  }
  
  # WT selection menu
  cat("\n=== Select Wild Type (WT) Marker ===\n")
  wt_choice <- menu(available_markers, title = "Choose WT marker:")
  if(wt_choice == 0) stop("WT marker selection cancelled")
  wt_marker <- available_markers[wt_choice]
  
  # KO selection menu  
  remaining_markers <- available_markers[-wt_choice]
  ko_marker <- NULL
  
  if(length(remaining_markers) > 0) {
    cat("\n=== Select Knock Out (KO) Marker (Optional) ===\n")
    ko_options <- c(remaining_markers, "Skip KO assignment")
    ko_choice <- menu(ko_options, title = "Choose KO marker:")
    
    if(ko_choice > 0 && ko_choice <= length(remaining_markers)) {
      ko_marker <- remaining_markers[ko_choice]
      remaining_markers <- remaining_markers[-ko_choice]
    }
  }
  
  # Recipient selection menu
  recipient_marker <- NULL
  
  if(length(remaining_markers) > 0) {
    cat("\n=== Select Recipient Marker (Optional) ===\n")
    recipient_options <- c(remaining_markers, "Skip Recipient assignment")
    recipient_choice <- menu(recipient_options, title = "Choose Recipient marker:")
    
    if(recipient_choice > 0 && recipient_choice <= length(remaining_markers)) {
      recipient_marker <- remaining_markers[recipient_choice]
    }
  }
  
  # Create genotype column
  df_updated <- df %>%
    mutate(genotype = case_when(
      congenics == wt_marker ~ "donor_WT",
      TRUE ~ "Other"
    ))
  
  # Add KO assignments if ko_marker was selected
  if(!is.null(ko_marker)) {
    df_updated <- df_updated %>%
      mutate(genotype = case_when(
        congenics == ko_marker ~ "donor_KO",
        TRUE ~ genotype
      ))
  }
  
  # Add Recipient assignments if recipient_marker was selected
  if(!is.null(recipient_marker)) {
    df_updated <- df_updated %>%
      mutate(genotype = case_when(
        congenics == recipient_marker ~ "Recipient",
        TRUE ~ genotype
      ))
  }
  
  # Display genotype results
  cat("\n=== Genotype Assignment Complete ===\n")
  cat(paste("WT:", wt_marker, "\n"))
  if(!is.null(ko_marker)) cat(paste("KO:", ko_marker, "\n"))
  if(!is.null(recipient_marker)) cat(paste("Recipient:", recipient_marker, "\n"))
  cat("\nGenotype counts:\n")
  print(table(df_updated$genotype, useNA = "ifany"))
  
  # Continue with rest of original metadata assignment function...
  # [Rest of the original assign_metadata_menu function continues here]
  
  # Metadata addition section
  cat("\n=== Add Custom Metadata ===\n")
  
  metadata_complete <- FALSE
  
  while(!metadata_complete) {
    metadata_options <- c(
      "Add timepoint",
      "Add batch",
      "Add sex",
      "Add custom metadata column",
      "Finish and return data"
    )
    
    cat("\n--- Metadata Menu ---\n")
    metadata_choice <- menu(metadata_options, title = "Choose metadata option:")
    
    if(metadata_choice == 0) {
      cat("Metadata addition cancelled. Returning data with genotype assignments only.\n")
      break
    }
    
    if(metadata_choice == 1) {
      # Add timepoint
      cat("\n=== Add Timepoint ===\n")
      timepoint_value <- readline(prompt = "Enter timepoint value: ")
      if(timepoint_value != "") {
        df_updated <- df_updated %>%
          mutate(timepoint = timepoint_value)
        cat(paste("Added timepoint:", timepoint_value, "\n"))
      } else {
        cat("Timepoint addition skipped (empty input).\n")
      }
      
    } else if(metadata_choice == 2) {
      # Add batch
      cat("\n=== Add Batch ===\n")
      batch_value <- readline(prompt = "Enter batch value: ")
      if(batch_value != "") {
        df_updated <- df_updated %>%
          mutate(batch = batch_value)
        cat(paste("Added batch:", batch_value, "\n"))
      } else {
        cat("Batch addition skipped (empty input).\n")
      }
      
    } else if(metadata_choice == 3) {
      # Add sex
      cat("\n=== Add Sex ===\n")
      sex_options <- c("m", "f", "Enter custom value")
      sex_choice <- menu(sex_options, title = "Choose sex:")
      
      if(sex_choice == 0) {
        cat("Sex addition cancelled.\n")
      } else if(sex_choice == 1) {
        df_updated <- df_updated %>%
          mutate(sex = "m")
        cat("Added sex: m\n")
      } else if(sex_choice == 2) {
        df_updated <- df_updated %>%
          mutate(sex = "f")
        cat("Added sex: f\n")
      } else if(sex_choice == 3) {
        sex_value <- readline(prompt = "Enter custom sex value: ")
        if(sex_value != "") {
          df_updated <- df_updated %>%
            mutate(sex = sex_value)
          cat(paste("Added sex:", sex_value, "\n"))
        } else {
          cat("Sex addition skipped (empty input).\n")
        }
      }
      
    } else if(metadata_choice == 4) {
      # Add custom metadata column
      cat("\n=== Add Custom Metadata Column ===\n")
      column_name <- readline(prompt = "Enter column name: ")
      
      if(column_name != "" && !column_name %in% names(df_updated)) {
        column_value <- readline(prompt = paste("Enter value for", column_name, ": "))
        if(column_value != "") {
          df_updated <- df_updated %>%
            mutate(!!sym(column_name) := column_value)
          cat(paste("Added", column_name, ":", column_value, "\n"))
        } else {
          cat("Custom column addition skipped (empty value).\n")
        }
      } else if(column_name == "") {
        cat("Custom column addition skipped (empty column name).\n")
      } else {
        cat("Column name already exists. Please choose a different name.\n")
      }
      
    } else if(metadata_choice == 5) {
      # Finish and return data
      metadata_complete <- TRUE
    }
  }
  
  # Display final results
  cat("\n=== Final Results ===\n")
  cat("Columns in dataset:\n")
  cat(paste(names(df_updated), collapse = ", "), "\n")
  cat("\nDataset summary:\n")
  cat(paste("Rows:", nrow(df_updated), "\n"))
  cat(paste("Columns:", ncol(df_updated), "\n"))
  
  # Summary of analysis factors if they exist
  if("tissue_factor" %in% names(df_updated)) {
    cat("\nTissue Factor Summary:\n")
    tissue_summary <- df_updated %>% count(tissue_factor, name = "n") %>% arrange(desc(n))
    print(tissue_summary)
  }
  
  if("tissue_factor" %in% names(df_updated)) {
    pairing_summary <- df_updated %>% count(tissue_factor, name = "n") %>% arrange(desc(n))
    cat("\nPairing Factor Summary:\n")
    print(head(pairing_summary, 10))
    if(nrow(pairing_summary) > 10) {
      cat("... and", nrow(pairing_summary) - 10, "more pairing groups\n")
    }
  }
  
  return(df_updated)
}

#============================================================================
# Data Cleanup Functions
#============================================================================
data_clean_custom <- function(data, auto_save = FALSE) {
  
  # Function to handle unstained sample detection and removal
  handle_unstained_samples <- function(df, data_type = "data") {
    if(!"Sample" %in% names(df)) {
      return(df)
    }
    
    unstained_patterns <- c("unstained", "no stain", "no_stain")
    pattern <- paste(unstained_patterns, collapse = "|")
    
    unstained_rows <- df %>%
      mutate(row_id = row_number()) %>%
      dplyr::filter(str_detect(tolower(Sample), pattern)) %>%
      pull(row_id)
    
    if(length(unstained_rows) == 0) {
      return(df)
    }
    
    # Only show interactive prompt for the first dataset
    if(!exists(".unstained_decision", envir = .GlobalEnv)) {
      cat("\n=== Unstained Samples Detected in", data_type, "===\n")
      cat(paste("Found", length(unstained_rows), "unstained sample(s)\n"))
      
      unstained_data <- df %>% slice(unstained_rows)
      if("NodeShort" %in% names(df)) {
        print(unstained_data %>% select(Sample, NodeShort) %>% head(10))
      } else {
        print(unstained_data %>% select(Sample) %>% head(10))
      }
      
      action_options <- c(
        "Remove unstained samples from ALL datasets",
        "Keep unstained samples in ALL datasets"
      )
      
      user_choice <- menu(action_options, title = "What would you like to do with unstained samples?")
      
      # Store decision globally
      .unstained_decision <<- if(user_choice == 1) "remove" else "keep"
      
      if(.unstained_decision == "remove") {
        cat("Will remove unstained samples from all datasets\n")
      } else {
        cat("Will keep unstained samples in all datasets\n")
      }
    }
    
    # Apply the stored decision
    if(.unstained_decision == "remove") {
      df_cleaned <- df %>% slice(-unstained_rows)
      cat("Removed", length(unstained_rows), "unstained samples from", data_type, "\n")
      return(df_cleaned)
    } else {
      return(df)
    }
  }
  
  # Function to clean a single data frame
  clean_single_df <- function(df, data_type = "data") {
    names(df) <- str_remove_all(names(df), "\\$")
    
    if("NodeShort" %in% names(df)) {
      df <- df %>%
        mutate(
          NodeShort = NodeShort %>%
            str_remove(".*:") %>%
            str_remove_all(" ") %>%
            str_remove_all(",") %>%
            str_remove_all(":")
        )
    }
    
    df <- handle_unstained_samples(df, data_type)
    return(df)
  }
  
  # Clear any existing decision at the start
  if(exists(".unstained_decision", envir = .GlobalEnv)) {
    rm(.unstained_decision, envir = .GlobalEnv)
  }
  
  # Check if input is a list
  if(is.list(data) && !is.data.frame(data)) {
    # Apply cleaning function to each element in the list
    cleaned_data <- imap(data, ~ {
      if(is.data.frame(.x)) {
        data_type <- if(.y == "counts") "counts data" else if(.y == "mfi") "MFI data" else .y
        clean_single_df(.x, data_type)
      } else {
        .x
      }
    })
    
    # Clean up the global decision variable
    if(exists(".unstained_decision", envir = .GlobalEnv)) {
      rm(.unstained_decision, envir = .GlobalEnv)
    }
    
    # Interactive save option
    if(!auto_save) {
      cat("\n=== Save Cleaned Data ===\n")
      save_choice <- menu(c("Save cleaned data", "Don't save"), title = "Save cleaned datasets?")
      
      if(save_choice == 1) {
        if("counts" %in% names(cleaned_data)) {
          interactive_save_dataframe(cleaned_data$counts, "cleaned_counts", "processed")
        }
        if("mfi" %in% names(cleaned_data)) {
          interactive_save_dataframe(cleaned_data$mfi, "cleaned_mfi", "processed")
        }
      }
    }
    
    return(cleaned_data)
    
  } else if(is.data.frame(data)) {
    result <- clean_single_df(data, "single dataset")
    
    # Clean up the global decision variable
    if(exists(".unstained_decision", envir = .GlobalEnv)) {
      rm(.unstained_decision, envir = .GlobalEnv)
    }
    
    # Interactive save option
    if(!auto_save) {
      cat("\n=== Save Cleaned Data ===\n")
      save_choice <- menu(c("Save cleaned data", "Don't save"), title = "Save cleaned dataset?")
      
      if(save_choice == 1) {
        interactive_save_dataframe(result, "cleaned_data", "processed")
      }
    }
    
    return(result)
  } else {
    warning("Input is neither a list nor a data frame. Returning unchanged.")
    return(data)
  }
}

#=============================================================================
# Plotting helper functions
#=============================================================================
# Helper function to create individual subgroup plots
create_subgroup_plot <- function(df, test_type = "t_test", facet_var = "tissue_factor") {
  # Check if we have enough data for statistical testing
  if (length(unique(df$congenics)) < 2) {
    warning(paste("Subgroup has less than 2 'congenics' levels for comparison. Skipping statistical test."))
    stat_test_results <- NULL
  } else {
    # Perform statistical test based on test_type and facet_var
    if (facet_var == "NodeShort") {
      # When faceting by NodeShort, group by NodeShort for stats
      stat_test_results <- df %>%
        group_by(NodeShort) %>%
        dplyr::filter(n_distinct(congenics) == 2) %>%
        {
          if (test_type == "t_test") {
            t_test(., Freq ~ congenics, paired = TRUE)
          } else if (test_type == "wilcox_test") {
            wilcox_test(., Freq ~ congenics, paired = TRUE)
          }
        } %>%
        add_xy_position(x = "congenics") %>%
        mutate(p.format = scales::pvalue(p, accuracy = 0.001, add_p = TRUE))
    } else {
      # When faceting by tissue_factor or no faceting, group by tissue_factor
      stat_test_results <- df %>%
        group_by(tissue_factor) %>%
        dplyr::filter(n_distinct(congenics) == 2) %>%
        {
          if (test_type == "t_test") {
            t_test(., Freq ~ congenics, paired = TRUE)
          } else if (test_type == "wilcox_test") {
            wilcox_test(., Freq ~ congenics, paired = TRUE)
          }
        } %>%
        add_xy_position(x = "congenics") %>%
        mutate(p.format = scales::pvalue(p, accuracy = 0.001, add_p = TRUE))
    }
  }
  
  # Calculate y-axis limits
  max_freq <- max(df$Freq, na.rm = TRUE)
  max_y_position <- if (!is.null(stat_test_results) && nrow(stat_test_results) > 0) {
    max(max_freq, stat_test_results$y.position, na.rm = TRUE)
  } else {
    max_freq
  }
  upper_y_limit <- max_y_position * 1.2
  
  # Create the base plot
  p <- ggplot(df, aes(x = congenics, y = Freq, group = pairing_factor)) +
    geom_line(color = "black", linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
    geom_point(shape = 21, fill = "white", color = "black", size = 3, show.legend = FALSE)
  
  # Add faceting if requested
  if (facet_var == "tissue_factor") {
    p <- p + facet_wrap(~ tissue_factor, scales = "free_y", strip.position = "top")
  } else if (facet_var == "NodeShort") {
    p <- p + facet_wrap(~ NodeShort, scales = "free_y", strip.position = "top")
  }
  
  # Dynamic title and y-axis label based on faceting
  if (facet_var == "tissue_factor") {
    # When faceting by tissue, title should show the NodeShort
    plot_title <- paste("Cell Population:", unique(df$NodeShort))
    y_label <- paste0(unique(df$NodeShort), " (%)")
  } else if (facet_var == "NodeShort") {
    # When faceting by NodeShort, title should show the tissue
    plot_title <- paste("Tissue:", unique(df$tissue_factor))
    y_label <- "Cell Frequency (%)"
  } else {
    # No faceting
    plot_title <- paste("Subgroup:", unique(df$NodeShort))
    y_label <- paste0(unique(df$NodeShort), " (%)")
  }
  
  # Continue with plot formatting
  p <- p +
    scale_y_continuous(
      limits = c(0, upper_y_limit),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      title = plot_title,
      x = "Congenic Marker",
      y = y_label
    ) +
    theme_minimal(base_size = 14) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 0.4),
      axis.line.y = element_line(color = "black", linewidth = 0.4),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      strip.placement = "outside"
    )
  
  # Add statistical annotations if available
  if (!is.null(stat_test_results) && nrow(stat_test_results) > 0) {
    p <- p + stat_pvalue_manual(
      stat_test_results,
      label = "p.format",
      tip.length = 0,
      bracket.size = 0.5,
      hide.ns = FALSE
    )
  }
  
  return(p)
}
#===============================================================================
# Main Paired compariosn plots interactive function
#===============================================================================
create_paired_comparison_plots <- function(data, auto_save = FALSE) {
  
  # Load required libraries
  require(dplyr)
  require(purrr)
  require(ggplot2)
  require(rstatix)
  require(ggpubr)
  require(scales)
  
  # Validate inputs
  if(missing(data) || !is.data.frame(data)) {
    stop("Please provide a valid data frame")
  }
  
  required_cols <- c("NodeShort", "tissue_factor", "pairing_factor", "congenics", "Freq")
  missing_cols <- setdiff(required_cols, names(data))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
  }
  
  # Interactive test selection
  test_choice <- menu(c("Paired t-test", "Wilcoxon signed-rank test"),
                      title = "Choose statistical test:")
  if(test_choice == 0) test_choice <- 1
  test_type <- c("t_test", "wilcox_test")[test_choice]
  
  # Interactive faceting selection
  facet_choice <- menu(c("tissue_factor (tissue)", "NodeShort (subpopulation of interest)", "No faceting"),
                       title = "Choose how to facet the plots:")
  if(facet_choice == 0) facet_choice <- 1
  facet_var <- c("tissue_factor", "NodeShort", "none")[facet_choice]
  
  # Print information about the analysis
  test_description <- switch(test_type, 
                             't_test' = 'paired t-test', 
                             'wilcox_test' = 'Wilcoxon signed-rank test')
  facet_description <- switch(facet_var,
                              'tissue_factor' = 'tissue',
                              'NodeShort' = 'subpopulation of interest',
                              'none' = 'no faceting (separate plots)')
  
  # Choose split variable based on what you want separate plots for
  if(facet_var == "tissue_factor") {
    split_var <- "NodeShort"
    plot_description <- "cell population(s)"
  } else if(facet_var == "NodeShort") {
    split_var <- "tissue_factor"  
    plot_description <- "tissue(s)"
  } else {
    split_var <- c("NodeShort", "tissue_factor")
    plot_description <- "NodeShort-tissue combination(s)"
  }
  
  cat("Creating paired comparison plots using", test_description, "\n")
  cat("Faceting by", facet_description, "\n")
  
  # Calculate number of plots that will be created
  if(facet_var == "none") {
    n_plots <- data %>% 
      distinct(NodeShort, tissue_factor) %>% 
      nrow()
    cat("Processing", n_plots, plot_description, "\n\n")
  } else {
    n_plots <- n_distinct(data[[split_var]])
    cat("Processing", n_plots, plot_description, "\n\n")
  }
  
  # Helper function to create individual subgroup plots
  create_subgroup_plot <- function(df, test_type = "t_test", facet_var = "tissue_factor") {
    # Check if we have enough data for statistical testing
    if(length(unique(df$congenics)) < 2) {
      warning(paste("Subgroup has less than 2 'congenics' levels for comparison. Skipping statistical test."))
      stat_test_results <- NULL
    } else {
      # Perform statistical test based on test_type and facet_var
      if(facet_var == "NodeShort") {
        # When faceting by NodeShort, group by NodeShort for stats
        stat_test_results <- df %>%
          group_by(NodeShort) %>%
          dplyr::filter(n_distinct(congenics) == 2) %>%
          {
            if(test_type == "t_test") {
              t_test(., Freq ~ congenics, paired = TRUE)
            } else if(test_type == "wilcox_test") {
              wilcox_test(., Freq ~ congenics, paired = TRUE)
            }
          } %>%
          add_xy_position(x = "congenics") %>%
          mutate(p.format = scales::pvalue(p, accuracy = 0.001, add_p = TRUE))
      } else {
        # When faceting by tissue_factor or no faceting, group by tissue_factor
        stat_test_results <- df %>%
          group_by(tissue_factor) %>%
          dplyr::filter(n_distinct(congenics) == 2) %>%
          {
            if(test_type == "t_test") {
              t_test(., Freq ~ congenics, paired = TRUE)
            } else if(test_type == "wilcox_test") {
              wilcox_test(., Freq ~ congenics, paired = TRUE)
            }
          } %>%
          add_xy_position(x = "congenics") %>%
          mutate(p.format = scales::pvalue(p, accuracy = 0.001, add_p = TRUE))
      }
    }
    
    # Calculate y-axis limits
    max_freq <- max(df$Freq, na.rm = TRUE)
    max_y_position <- if(!is.null(stat_test_results) && nrow(stat_test_results) > 0) {
      max(max_freq, stat_test_results$y.position, na.rm = TRUE)
    } else {
      max_freq
    }
    upper_y_limit <- max_y_position * 1.2
    
    # Create the base plot
    p <- ggplot(df, aes(x = congenics, y = Freq, group = pairing_factor)) +
      geom_line(color = "black", linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
      geom_point(shape = 21, fill = "white", color = "black", size = 3, show.legend = FALSE)
    
    # Add faceting if requested
    if(facet_var == "tissue_factor") {
      p <- p + facet_wrap(~ tissue_factor, scales = "free_y", strip.position = "top")
    } else if(facet_var == "NodeShort") {
      p <- p + facet_wrap(~ NodeShort, scales = "free_y", strip.position = "top")
    }
    
    # Dynamic title and y-axis label based on faceting
    if(facet_var == "tissue_factor") {
      plot_title <- paste("Cell Population:", unique(df$NodeShort))
      y_label <- paste0(unique(df$NodeShort), " (%)")
    } else if(facet_var == "NodeShort") {
      plot_title <- paste("Tissue:", unique(df$tissue_factor))
      y_label <- "Cell Frequency (%)"
    } else {
      plot_title <- paste("Subgroup:", unique(df$NodeShort))
      y_label <- paste0(unique(df$NodeShort), " (%)")
    }
    
    # Continue with plot formatting
    p <- p +
      scale_y_continuous(
        limits = c(0, upper_y_limit),
        expand = expansion(mult = c(0.02, 0.05))
      ) +
      labs(
        title = plot_title,
        x = "Congenic Marker",
        y = y_label
      ) +
      theme_minimal(base_size = 14) +
      theme(
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.placement = "outside"
      )
    
    # Add statistical annotations if available
    if(!is.null(stat_test_results) && nrow(stat_test_results) > 0) {
      p <- p + stat_pvalue_manual(
        stat_test_results,
        label = "p.format",
        tip.length = 0,
        bracket.size = 0.5,
        hide.ns = FALSE
      )
    }
    
    return(p)
  }
  
  # Create plots based on the correct splitting logic
  if(facet_var == "tissue_factor") {
    plots <- data %>%
      group_split(NodeShort, .keep = TRUE) %>%
      map(~ create_subgroup_plot(.x, test_type = test_type, facet_var = "tissue_factor")) %>%
      set_names(map_chr(data %>% distinct(NodeShort) %>% pull(NodeShort), as.character))
  } else if(facet_var == "NodeShort") {
    plots <- data %>%
      group_split(tissue_factor, .keep = TRUE) %>%
      map(~ create_subgroup_plot(.x, test_type = test_type, facet_var = "NodeShort")) %>%
      set_names(map_chr(data %>% distinct(tissue_factor) %>% pull(tissue_factor), as.character))
  } else {
    plots <- data %>%
      group_split(NodeShort, tissue_factor, .keep = TRUE) %>%
      map(~ create_subgroup_plot(.x, test_type = test_type, facet_var = "none")) %>%
      set_names(map_chr(1:length(.), ~ {
        df <- data %>% group_split(NodeShort, tissue_factor, .keep = TRUE) %>% .[[.x]]
        paste0(unique(df$NodeShort), "_", unique(df$tissue_factor))
      }))
  }
  
  # Print summary
  plot_count <- length(plots)
  cat("Successfully created", plot_count, "plot(s)", "\n")
  cat("Available plots:", paste(names(plots), collapse = ", "), "\n")
  
  # Interactive save option
  if(!auto_save && length(plots) > 0) {
    cat("\n=== Save Paired Comparison Plots ===\n")
    save_choice <- menu(c("Save all plots", "Save selected plots", "Don't save"), 
                        title = "Save options for paired comparison plots:")
    
    if(save_choice == 1) {
      # Save all plots
      saved_count <- 0
      iwalk(plots, function(plot, plot_name) {
        suggested_name <- paste0("paired_comparison_", make.names(plot_name))
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        filename <- paste0(suggested_name, "_", timestamp, ".png")
        filepath <- file.path(here::here("out/plots/paired_comparisons"), filename)
        
        tryCatch({
          ggsave(filepath, plot = plot, width = 10, height = 8, dpi = 300)
          cat("Saved:", filename, "\n")
          saved_count <<- saved_count + 1
        }, error = function(e) {
          cat("Error saving", plot_name, ":", e$message, "\n")
        })
      })
      cat("Saved", saved_count, "out of", length(plots), "plots\n")
      
    } else if(save_choice == 2) {
      # Save selected plots
      cat("\nAvailable plots:\n")
      iwalk(names(plots), ~cat(sprintf("%d. %s\n", .y, .x)))
      
      selection <- readline("Enter plot numbers to save (space or comma-separated): ")
      if(selection != "") {
        tryCatch({
          if(grepl("\\s", selection)) {
            selected_indices <- as.numeric(str_trim(str_split(selection, "\\s+")[[1]]))
          } else {
            selected_indices <- as.numeric(str_trim(str_split(selection, ",")[[1]]))
          }
          
          selected_indices <- selected_indices[!is.na(selected_indices)]
          valid_indices <- selected_indices[selected_indices >= 1 & selected_indices <= length(plots)]
          
          if(length(valid_indices) > 0) {
            selected_plots <- plots[valid_indices]
            selected_names <- names(plots)[valid_indices]
            
            saved_count <- 0
            iwalk(selected_plots, function(plot, idx) {
              plot_name <- selected_names[idx]
              suggested_name <- paste0("paired_comparison_", make.names(plot_name))
              filepath <- interactive_save_plot(plot, suggested_name, "paired_comparison")
              if(!is.null(filepath)) saved_count <<- saved_count + 1
            })
            
            cat("Saved", saved_count, "selected plots\n")
          }
          
        }, error = function(e) {
          cat("Invalid input format\n")
        })
      }
    }
  }
  
  return(plots)
}

# Example usage:
# plots <- create_paired_comparison_plots(your_data)  # User will be prompted to choose test type
# plots[["CD103+CD69+"]]  # View specific cell population plot (when faceting by tissue_factor)
# plots[["Spleen"]]       # View specific tissue plot (when faceting by NodeShort)


#===============================================================================
# Enhanced MFI heatmaps with additional interactive options
#===============================================================================
# 
# Enhanced MFI Heatmap Code with Statistical Testing - Cleaned Version
# Optimized for create_mfi_heatmaps_interactive_enhanced() workflow

# Required libraries
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(grid)

# ===== INTERACTIVE SELECTION FUNCTIONS =====

select_congenics_interactive <- function(mfi_data) {
  available_congenics <- mfi_data %>%
    dplyr::filter(!is.na(congenics)) %>%
    pull(congenics) %>%
    unique() %>%
    sort()
  
  if (length(available_congenics) == 0) {
    stop("No valid congenic values found in data")
  }
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("CONGENIC SELECTION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  cat("Available congenics in your data:\n")
  iwalk(available_congenics, function(congenic, idx) {
    sample_count <- mfi_data %>%
      dplyr::filter(congenics == congenic) %>%
      pull(Sample) %>%
      n_distinct()
    
    cat(sprintf("%2d. %s (%d samples)\n", idx, congenic, sample_count))
  })
  
  cat(sprintf("%2d. Select all congenics\n", length(available_congenics) + 1))
  cat(sprintf("%2d. Show congenic distribution by tissue\n", length(available_congenics) + 2))
  
  while (TRUE) {
    choice <- readline("\nEnter congenic numbers (space or comma-separated, e.g., '1 3 5' or '1,3,5') or option: ")
    
    if (choice == as.character(length(available_congenics) + 1)) {
      selected_congenics <- available_congenics
      cat(sprintf("✓ Selected all %d congenics\n", length(selected_congenics)))
      return(selected_congenics)
    }
    
    if (choice == as.character(length(available_congenics) + 2)) {
      cat("\n--- Congenic Distribution by Tissue ---\n")
      distribution <- mfi_data %>%
        dplyr::filter(!is.na(congenics)) %>%
        count(tissue_factor, congenics) %>%
        pivot_wider(names_from = congenics, values_from = n, values_fill = 0)
      
      print(distribution)
      cat("\nPress Enter to continue...")
      readline()
      next
    }
    
    tryCatch({
      if (grepl("\\s", choice)) {
        selected_indices <- as.numeric(str_trim(str_split(choice, "\\s+")[[1]]))
      } else {
        selected_indices <- as.numeric(str_trim(str_split(choice, ",")[[1]]))
      }
      
      selected_indices <- selected_indices[!is.na(selected_indices)]
      
      if (length(selected_indices) > 0 && all(selected_indices >= 1 & selected_indices <= length(available_congenics))) {
        selected_congenics <- available_congenics[selected_indices]
        cat(sprintf("✓ Selected %d congenics: %s\n", 
                    length(selected_congenics), 
                    paste(selected_congenics, collapse = ", ")))
        return(selected_congenics)
      } else {
        cat("Invalid selection. Please enter numbers within the valid range.\n")
      }
    }, error = function(e) {
      cat("Invalid format. Please enter space or comma-separated numbers (e.g., '1 3 5' or '1,3,5').\n")
    })
  }
}

select_grouping_option <- function(mfi_data) {
  available_groups <- unique(mfi_data$tissue_factor)
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("GROUPING OPTIONS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  cat("1. Separate heatmaps for each tissue (current behavior)\n")
  cat("2. Combined heatmap showing all tissues together\n")
  cat("3. Show tissue distribution first\n")
  
  while (TRUE) {
    choice <- readline("\nChoose grouping option (1-3): ")
    
    if (choice == "1") {
      return(list(type = "separate", groups = available_groups))
    } else if (choice == "2") {
      return(list(type = "combined", groups = available_groups))
    } else if (choice == "3") {
      cat("\n--- Tissue Distribution ---\n")
      distribution <- mfi_data %>%
        dplyr::filter(!is.na(congenics)) %>%
        count(tissue_factor, congenics) %>%
        arrange(tissue_factor, congenics)
      
      cat("Available tissues:", paste(available_groups, collapse = ", "), "\n")
      print(distribution)
      
      marker_dist <- mfi_data %>%
        dplyr::filter(!is.na(congenics)) %>%
        count(tissue_factor, marker) %>%
        group_by(tissue_factor) %>%
        summarise(
          n_markers = n_distinct(marker),
          total_measurements = sum(n),
          .groups = "drop"
        )
      
      cat("\n--- Markers per Tissue ---\n")
      print(marker_dist)
      
      cat("\nPress Enter to continue...")
      readline()
      next
    } else {
      cat("Invalid choice. Please enter 1, 2, or 3.\n")
    }
  }
}

select_markers_interactive <- function(mfi_data) {
  available_markers <- mfi_data %>%
    dplyr::filter(!is.na(congenics)) %>%
    pull(marker) %>%
    unique() %>%
    sort()
  
  if (length(available_markers) == 0) {
    stop("No valid markers found in data")
  }
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("MARKER SELECTION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  cat("Available markers in your data:\n")
  iwalk(available_markers, function(marker, idx) {
    measurement_count <- mfi_data %>%
      dplyr::filter(marker == .env$marker, !is.na(congenics)) %>%
      nrow()
    
    cat(sprintf("%2d. %s (%d measurements)\n", idx, marker, measurement_count))
  })
  
  cat(sprintf("%2d. Select all markers\n", length(available_markers) + 1))
  cat(sprintf("%2d. Show marker distribution by tissue\n", length(available_markers) + 2))
  cat(sprintf("%2d. Search markers by keyword\n", length(available_markers) + 3))
  
  while (TRUE) {
    choice <- readline("\nEnter marker numbers (space or comma-separated, e.g., '1 3 5' or '1,3,5') or option: ")
    
    if (choice == as.character(length(available_markers) + 1)) {
      selected_markers <- available_markers
      cat(sprintf("✓ Selected all %d markers\n", length(selected_markers)))
      return(selected_markers)
    }
    
    if (choice == as.character(length(available_markers) + 2)) {
      cat("\n--- Marker Distribution by Tissue ---\n")
      distribution <- mfi_data %>%
        dplyr::filter(!is.na(congenics)) %>%
        count(tissue_factor, marker) %>%
        pivot_wider(names_from = tissue_factor, values_from = n, values_fill = 0)
      
      print(distribution)
      cat("\nPress Enter to continue...")
      readline()
      next
    }
    
    if (choice == as.character(length(available_markers) + 3)) {
      keyword <- readline("Enter search keyword (case-insensitive): ")
      if (keyword != "") {
        matching_markers <- available_markers[grepl(keyword, available_markers, ignore.case = TRUE)]
        if (length(matching_markers) > 0) {
          cat(sprintf("Found %d matching markers:\n", length(matching_markers)))
          iwalk(matching_markers, function(marker, idx) {
            original_idx <- which(available_markers == marker)
            cat(sprintf("%2d. %s (original #%d)\n", original_idx, marker, original_idx))
          })
          cat("\nYou can now select these by their original numbers.\n")
        } else {
          cat("No markers found matching your search term.\n")
        }
      }
      cat("\nPress Enter to continue...")
      readline()
      next
    }
    
    tryCatch({
      if (grepl("\\s", choice)) {
        selected_indices <- as.numeric(str_trim(str_split(choice, "\\s+")[[1]]))
      } else {
        selected_indices <- as.numeric(str_trim(str_split(choice, ",")[[1]]))
      }
      
      selected_indices <- selected_indices[!is.na(selected_indices)]
      
      if (length(selected_indices) > 0 && all(selected_indices >= 1 & selected_indices <= length(available_markers))) {
        selected_markers <- available_markers[selected_indices]
        cat(sprintf("✓ Selected %d markers: %s\n", 
                    length(selected_markers), 
                    paste(selected_markers, collapse = ", ")))
        return(selected_markers)
      } else {
        cat("Invalid selection. Please enter numbers within the valid range.\n")
      }
    }, error = function(e) {
      cat("Invalid format. Please enter space or comma-separated numbers (e.g., '1 3 5' or '1,3,5').\n")
    })
  }
}

select_scaling_method_interactive <- function() {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("SCALING METHOD SELECTION\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  cat("Available scaling methods:\n")
  cat("1. No scaling (raw MFI values)\n")
  cat("   → Use when: Markers have similar dynamic ranges, absolute values matter\n\n")
  
  cat("2. Row scaling (Z-score per marker)\n")
  cat("   → Use when: Markers have very different expression ranges\n\n")
  
  cat("3. Global scaling (Z-score across all data)\n")
  cat("   → Use when: Want to see overall patterns while preserving relative differences\n\n")
  
  cat("4. Log transformation\n")
  cat("   → Use when: Data spans several orders of magnitude, has multiplicative effects\n\n")
  
  cat("5. Square root transformation\n")
  cat("   → Use when: Data has moderate skewness, want gentler transformation than log\n\n")
  
  cat("6. Percentile scaling\n")
  cat("   → Use when: Want to focus on relative rankings, ignore extreme outliers\n\n")
  
  cat("7. Show data diagnostic\n")
  cat("   → Use when: Unsure which scaling to choose, want to examine data first\n\n")
  
  while (TRUE) {
    choice <- readline("\nChoose scaling method (1-7): ")
    
    if (choice == "1") {
      return(list(scale_method = "none", aggregation_method = "mean"))
    } else if (choice == "2") {
      return(list(scale_method = "row", aggregation_method = "mean"))
    } else if (choice == "3") {
      return(list(scale_method = "global", aggregation_method = "mean"))
    } else if (choice == "4") {
      log_choice <- readline("Enter log base (default 2): ")
      log_base <- ifelse(log_choice == "", 2, as.numeric(log_choice))
      if (is.na(log_base) || log_base <= 0) {
        cat("Invalid log base. Using default base 2.\n")
        log_base <- 2
      }
      return(list(scale_method = "log", aggregation_method = "mean", log_base = log_base))
    } else if (choice == "5") {
      return(list(scale_method = "sqrt", aggregation_method = "mean"))
    } else if (choice == "6") {
      perc_choice <- readline("Enter percentile range (space or comma-separated, default '0.05 0.95'): ")
      if (perc_choice == "") {
        percentile_range <- c(0.05, 0.95)
      } else {
        tryCatch({
          if (grepl("\\s", perc_choice)) {
            percentile_range <- as.numeric(str_trim(str_split(perc_choice, "\\s+")[[1]]))
          } else {
            percentile_range <- as.numeric(str_trim(str_split(perc_choice, ",")[[1]]))
          }
          
          if (length(percentile_range) != 2 || any(percentile_range < 0) || any(percentile_range > 1)) {
            stop("Invalid range")
          }
          if (percentile_range[1] >= percentile_range[2]) {
            stop("Lower percentile must be less than upper percentile")
          }
        }, error = function(e) {
          cat("Invalid percentile range. Using default 0.05, 0.95.\n")
          percentile_range <<- c(0.05, 0.95)
        })
      }
      return(list(scale_method = "percentile", aggregation_method = "mean", 
                  percentile_range = percentile_range))
    } else if (choice == "7") {
      return("DIAGNOSTIC")
    } else {
      cat("Invalid choice. Please enter 1-7.\n")
    }
  }
}

# ===== STATISTICAL TESTING FUNCTIONS =====

select_statistical_test_interactive <- function(mfi_data, selected_congenics, selected_markers) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("STATISTICAL TESTING OPTIONS\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  n_groups <- length(selected_congenics)
  
  cat("Number of congenic groups selected:", n_groups, "\n")
  cat("Groups:", paste(selected_congenics, collapse = ", "), "\n\n")
  
  if (n_groups < 2) {
    cat("⚠️  Statistical testing requires at least 2 groups.\n")
    cat("Please select more congenics in the previous step.\n")
    return(list(perform_stats = FALSE))
  }
  
  cat("Available statistical tests:\n")
  cat("1. No statistical testing (heatmap only)\n")
  cat("2. Show data distribution and recommendations\n")
  
  if (n_groups == 2) {
    cat("3. Unpaired t-test (independent samples)\n")
    cat("4. Paired t-test (matched samples)\n") 
    cat("5. Mann-Whitney U test (non-parametric, unpaired)\n")
    cat("6. Wilcoxon signed-rank test (non-parametric, paired)\n")
  } else {
    cat("3. One-way ANOVA (parametric)\n")
    cat("4. Kruskal-Wallis test (non-parametric)\n")
    cat("5. Repeated measures ANOVA (paired design)\n")
    cat("6. Friedman test (non-parametric, paired)\n")
  }
  
  while (TRUE) {
    choice <- readline(paste("\nChoose statistical test (1-6): "))
    
    if (choice == "1") {
      return(list(perform_stats = FALSE))
    } else if (choice == "2") {
      show_data_distribution(mfi_data, selected_congenics, selected_markers)
      cat("\nPress Enter to continue...")
      readline()
      next
    } else if (choice == "3") {
      if (n_groups == 2) {
        return(configure_two_group_test("t.test", "unpaired", selected_congenics, mfi_data))
      } else {
        return(configure_multi_group_test("anova", selected_congenics, mfi_data))
      }
    } else if (choice == "4") {
      if (n_groups == 2) {
        return(configure_two_group_test("t.test", "paired", selected_congenics, mfi_data))
      } else {
        return(configure_multi_group_test("kruskal", selected_congenics, mfi_data))
      }
    } else if (choice == "5") {
      if (n_groups == 2) {
        return(configure_two_group_test("mannwhitney", "unpaired", selected_congenics, mfi_data))
      } else {
        return(configure_multi_group_test("rm_anova", selected_congenics, mfi_data))
      }
    } else if (choice == "6") {
      if (n_groups == 2) {
        return(configure_two_group_test("wilcoxon", "paired", selected_congenics, mfi_data))
      } else {
        return(configure_multi_group_test("friedman", selected_congenics, mfi_data))
      }
    } else {
      cat("Invalid choice. Please enter a number between 1 and 6.\n")
    }
  }
}

configure_two_group_test <- function(test_type, pairing, selected_congenics, mfi_data) {
  cat(sprintf("\n--- Configuring %s (%s) ---\n", 
              switch(test_type,
                     "t.test" = "T-test",
                     "mannwhitney" = "Mann-Whitney U test",
                     "wilcoxon" = "Wilcoxon signed-rank test"),
              pairing))
  
  tissue_col <- "tissue_factor"  # Default tissue column
  sig_display <- select_significance_display()
  alpha <- get_alpha_level()
  
  pairing_var <- NULL
  if (pairing == "paired") {
    cat("\nFor paired analysis, samples need to be matched.\n")
    cat("Default: [pairing_factor]; Other Common pairing variables: Sample, NodeShort, Subject_ID, etc.\n")
    pairing_var <- readline("Enter column name for pairing (default: 'Sample'): ")
    if (pairing_var == "") pairing_var <- "Sample"
  }
  
  return(list(
    perform_stats = TRUE,
    test_type = test_type,
    pairing = pairing,
    pairing_var = pairing_var,
    alpha = alpha,
    sig_display = sig_display,
    groups = selected_congenics,
    post_hoc = FALSE,
    tissue_col = tissue_col
  ))
}

configure_multi_group_test <- function(test_type, selected_congenics, mfi_data) {
  cat(sprintf("\n--- Configuring %s ---\n", 
              switch(test_type,
                     "anova" = "One-way ANOVA",
                     "kruskal" = "Kruskal-Wallis test", 
                     "rm_anova" = "Repeated measures ANOVA",
                     "friedman" = "Friedman test")))
  
  tissue_col <- "tissue_factor"  # Default tissue column
  sig_display <- select_significance_display()
  alpha <- get_alpha_level()
  
  post_hoc <- FALSE
  if (test_type %in% c("anova", "kruskal")) {
    cat("\nPost-hoc testing options:\n")
    cat("1. No post-hoc tests\n")
    cat("2. Perform pairwise comparisons\n")
    
    ph_choice <- readline("Choose post-hoc option (1-2): ")
    if (ph_choice == "2") {
      post_hoc <- TRUE
      cat("Will perform pairwise comparisons with multiple comparison correction.\n")
    }
  }
  
  pairing_var <- NULL
  if (test_type %in% c("rm_anova", "friedman")) {
    cat("\nFor repeated measures analysis, samples need to be matched.\n")
    cat("Common pairing variables: Sample, NodeShort, Subject_ID, etc.\n")
    pairing_var <- readline("Enter column name for pairing (default: 'Sample'): ")
    if (pairing_var == "") pairing_var <- "Sample"
  }
  
  return(list(
    perform_stats = TRUE,
    test_type = test_type,
    pairing_var = pairing_var,
    alpha = alpha,
    sig_display = sig_display,
    groups = selected_congenics,
    post_hoc = post_hoc,
    tissue_col = tissue_col
  ))
}

select_significance_display <- function() {
  cat("\nSignificance display options:\n")
  cat("1. P-values (exact numbers, e.g., 0.032)\n")
  cat("2. Significance stars (*, **, ***)\n")
  cat("3. Both p-values and stars\n")
  cat("4. No significance overlay (stats in separate output)\n")
  
  while (TRUE) {
    choice <- readline("Choose display option (1-4): ")
    
    if (choice == "1") {
      return("p_values")
    } else if (choice == "2") {
      return("stars")
    } else if (choice == "3") {
      return("both")
    } else if (choice == "4") {
      return("none")
    } else {
      cat("Invalid choice. Please enter 1-4.\n")
    }
  }
}

get_alpha_level <- function() {
  cat("\nSignificance level options:\n")
  cat("1. α = 0.05 (standard)\n")
  cat("2. α = 0.01 (strict)\n")
  cat("3. α = 0.001 (very strict)\n")
  cat("4. Custom alpha level\n")
  
  while (TRUE) {
    choice <- readline("Choose significance level (1-4): ")
    
    if (choice == "1") {
      return(0.05)
    } else if (choice == "2") {
      return(0.01)
    } else if (choice == "3") {
      return(0.001)
    } else if (choice == "4") {
      custom_alpha <- readline("Enter custom alpha (e.g., 0.025): ")
      alpha_val <- as.numeric(custom_alpha)
      if (!is.na(alpha_val) && alpha_val > 0 && alpha_val < 1) {
        return(alpha_val)
      } else {
        cat("Invalid alpha value. Using default 0.05.\n")
        return(0.05)
      }
    } else {
      cat("Invalid choice. Please enter 1-4.\n")
    }
  }
}

# ===== DATA DIAGNOSTIC FUNCTION =====

diagnose_mfi_data <- function(mfi_data) {
  cat("\n=== MFI DATA DIAGNOSTIC ===\n")
  
  cat("Dataset dimensions:", nrow(mfi_data), "rows x", ncol(mfi_data), "columns\n")
  cat("Column names:", paste(names(mfi_data), collapse = ", "), "\n\n")
  
  required_cols <- c("marker", "MFI", "tissue_factor", "congenics", "NodeShort")
  missing_cols <- setdiff(required_cols, names(mfi_data))
  if (length(missing_cols) > 0) {
    cat("⚠️  Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  } else {
    cat("✓ All required columns present\n")
  }
  
  cat("\n--- Data Completeness ---\n")
  na_summary <- mfi_data %>%
    summarise(across(everything(), ~sum(is.na(.))))
  print(na_summary)
  
  if ("MFI" %in% names(mfi_data)) {
    cat("\n--- MFI Distribution ---\n")
    mfi_stats <- summary(mfi_data$MFI)
    print(mfi_stats)
    
    cat("MFI range:", min(mfi_data$MFI, na.rm = TRUE), "to", 
        max(mfi_data$MFI, na.rm = TRUE), "\n")
    
    zero_neg <- sum(mfi_data$MFI <= 0, na.rm = TRUE)
    if (zero_neg > 0) {
      cat("⚠️ ", zero_neg, "zero or negative MFI values found\n")
    }
  }
  
  cat("\n--- Unique Values ---\n")
  if ("congenics" %in% names(mfi_data)) {
    congenics_unique <- unique(mfi_data$congenics[!is.na(mfi_data$congenics)])
    cat("Congenics (", length(congenics_unique), "):", paste(congenics_unique, collapse = ", "), "\n")
  }
  
  if ("tissue_factor" %in% names(mfi_data)) {
    groups_unique <- unique(mfi_data$tissue_factor[!is.na(mfi_data$tissue_factor)])
    cat("Tissues (", length(groups_unique), "):", paste(groups_unique, collapse = ", "), "\n")
  }
  
  if ("marker" %in% names(mfi_data)) {
    markers_unique <- unique(mfi_data$marker[!is.na(mfi_data$marker)])
    cat("Markers (", length(markers_unique), "):", paste(head(markers_unique, 10), collapse = ", "))
    if (length(markers_unique) > 10) cat(", ... and", length(markers_unique) - 10, "more")
    cat("\n")
  }
  
  if (all(c("tissue_factor", "congenics") %in% names(mfi_data))) {
    cat("\n--- Sample Distribution ---\n")
    sample_dist <- mfi_data %>%
      dplyr::filter(!is.na(congenics), !is.na(tissue_factor)) %>%
      count(tissue_factor, congenics) %>%
      pivot_wider(names_from = congenics, values_from = n, values_fill = 0)
    
    print(sample_dist)
  }
}

show_data_distribution <- function(mfi_data, selected_congenics, selected_markers) {
  cat("\n=== DATA DISTRIBUTION ANALYSIS ===\n")
  
  filtered_data <- mfi_data %>%
    dplyr::filter(congenics %in% selected_congenics,
           marker %in% selected_markers,
           !is.na(MFI))
  
  if (nrow(filtered_data) == 0) {
    cat("No data available for analysis.\n")
    return()
  }
  
  cat("--- Sample Sizes ---\n")
  sample_sizes <- filtered_data %>%
    group_by(congenics) %>%
    summarise(
      n_measurements = n(),
      n_unique_samples = n_distinct(Sample, na.rm = TRUE),
      .groups = "drop"
    )
  print(sample_sizes)
  
  cat("\n--- Normality Testing ---\n")
  normality_results <- filtered_data %>%
    group_by(congenics, marker) %>%
    summarise(
      n = n(),
      mean_mfi = mean(MFI, na.rm = TRUE),
      sd_mfi = sd(MFI, na.rm = TRUE),
      shapiro_p = ifelse(n >= 3 & n <= 5000, 
                         tryCatch(shapiro.test(MFI)$p.value, error = function(e) NA),
                         NA),
      .groups = "drop"
    ) %>%
    mutate(
      normal_likely = case_when(
        is.na(shapiro_p) ~ "Unknown",
        shapiro_p > 0.05 ~ "Yes",
        TRUE ~ "No"
      )
    )
  
  norm_summary <- normality_results %>%
    count(normal_likely, name = "n_marker_group_combinations") %>%
    mutate(percentage = round(n_marker_group_combinations / sum(n_marker_group_combinations) * 100, 1))
  
  cat("Normality assessment (Shapiro-Wilk test, p > 0.05 suggests normal):\n")
  print(norm_summary)
  
  cat("\n--- RECOMMENDATIONS ---\n")
  n_groups <- length(selected_congenics)
  normal_pct <- norm_summary %>% dplyr::filter(normal_likely == "Yes") %>% pull(percentage)
  if (length(normal_pct) == 0) normal_pct <- 0
  
  if (n_groups == 2) {
    cat("For 2 groups:\n")
    if (normal_pct >= 70) {
      cat("✓ Data appears mostly normal → Consider t-test\n")
    } else {
      cat("⚠️  Data may not be normal → Consider Mann-Whitney U test\n")
    }
    cat("• Use paired tests if samples are matched (same subjects/locations)\n")
    cat("• Use unpaired tests if samples are independent\n")
  } else {
    cat("For", n_groups, "groups:\n")
    if (normal_pct >= 70) {
      cat("✓ Data appears mostly normal → Consider ANOVA\n")
    } else {
      cat("⚠️  Data may not be normal → Consider Kruskal-Wallis test\n")
    }
    cat("• Use repeated measures if samples are matched\n")
    cat("• Consider post-hoc tests if significant differences found\n")
  }
  
  cat("\n--- Example Distributions (first 3 markers) ---\n")
  example_markers <- head(selected_markers, 3)
  for (marker in example_markers) {
    marker_data <- filtered_data %>% dplyr::filter(marker == !!marker)
    cat(sprintf("\n%s:\n", marker))
    marker_summary <- marker_data %>%
      group_by(congenics) %>%
      summarise(
        n = n(),
        mean = round(mean(MFI, na.rm = TRUE), 2),
        median = round(median(MFI, na.rm = TRUE), 2),
        sd = round(sd(MFI, na.rm = TRUE), 2),
        .groups = "drop"
      )
    print(marker_summary)
  }
}

# ===== STATISTICAL ANALYSIS FUNCTIONS =====

perform_statistical_analysis <- function(mfi_data, selected_congenics, selected_markers, stats_config) {
  if (!stats_config$perform_stats) {
    return(NULL)
  }
  
  cat("\nPerforming statistical analysis...\n")
  cat("Test type:", stats_config$test_type, "\n")
  
  filtered_data <- mfi_data %>%
    dplyr::filter(congenics %in% selected_congenics,
           marker %in% selected_markers,
           !is.na(MFI))
  
  if (nrow(filtered_data) == 0) {
    warning("No data available for statistical analysis")
    return(NULL)
  }
  
  unique_tissues <- unique(filtered_data$tissue_factor)
  cat("Performing tests within each tissue separately:\n")
  cat("Tissues:", paste(unique_tissues, collapse = ", "), "\n")
  
  results_list <- list()
  
  for (current_tissue in unique_tissues) {
    cat(sprintf("\nProcessing tissue: %s\n", current_tissue))
    
    tissue_data <- filtered_data %>% dplyr::filter(tissue_factor == current_tissue)
    
    if (nrow(tissue_data) == 0) {
      warning(paste("No data for tissue:", current_tissue))
      next
    }
    
    for (current_marker in selected_markers) {
      marker_tissue_data <- tissue_data %>% dplyr::filter(marker == current_marker)
      
      if (nrow(marker_tissue_data) < 2) {
        warning(paste("Insufficient data for marker:", current_marker, "in tissue:", current_tissue))
        next
      }
      
      combo_id <- paste(current_tissue, current_marker, sep = "_")
      
      test_result <- switch(stats_config$test_type,
                            "t.test" = perform_t_test(marker_tissue_data, stats_config),
                            "mannwhitney" = perform_mannwhitney_test(marker_tissue_data, stats_config),
                            "wilcoxon" = perform_wilcoxon_test(marker_tissue_data, stats_config),
                            "anova" = perform_anova_test(marker_tissue_data, stats_config),
                            "kruskal" = perform_kruskal_test(marker_tissue_data, stats_config),
                            "rm_anova" = perform_rm_anova_test(marker_tissue_data, stats_config),
                            "friedman" = perform_friedman_test(marker_tissue_data, stats_config))
      
      if (!is.null(test_result)) {
        test_result$tissue <- current_tissue
        test_result$marker <- current_marker
        results_list[[combo_id]] <- test_result
        
        cat(sprintf("  - %s: p = %.4f\n", current_marker, test_result$p_value))
      }
    }
  }
  
  stats_summary <- compile_stats_summary(results_list, stats_config)
  
  return(list(
    individual_results = results_list,
    summary = stats_summary,
    config = stats_config
  ))
}

# Individual test functions
perform_t_test <- function(marker_data, stats_config) {
  groups <- stats_config$groups
  group1_data <- marker_data %>% dplyr::filter(congenics == groups[1]) %>% pull(MFI)
  group2_data <- marker_data %>% dplyr::filter(congenics == groups[2]) %>% pull(MFI)
  
  if (length(group1_data) == 0 || length(group2_data) == 0) {
    return(NULL)
  }
  
  if (stats_config$pairing == "paired") {
    if (is.null(stats_config$pairing_var)) {
      warning("Pairing variable not specified for paired t-test")
      return(NULL)
    }
    
    paired_data <- marker_data %>%
      select(all_of(c(stats_config$pairing_var, "congenics", "MFI"))) %>%
      pivot_wider(names_from = congenics, values_from = MFI) %>%
      dplyr::filter(!is.na(.data[[groups[1]]]), !is.na(.data[[groups[2]]]))
    
    if (nrow(paired_data) < 2) {
      warning("Insufficient paired data for t-test")
      return(NULL)
    }
    
    test_result <- t.test(paired_data[[groups[1]]], paired_data[[groups[2]]], 
                          paired = TRUE, conf.level = 1 - stats_config$alpha)
    
    return(list(
      test_name = "Paired t-test",
      p_value = test_result$p.value,
      statistic = test_result$statistic,
      df = test_result$parameter,
      method = test_result$method,
      n_pairs = nrow(paired_data),
      mean_diff = test_result$estimate
    ))
    
  } else {
    test_result <- t.test(group1_data, group2_data, 
                          var.equal = FALSE, conf.level = 1 - stats_config$alpha)
    
    return(list(
      test_name = "Unpaired t-test",
      p_value = test_result$p.value,
      statistic = test_result$statistic,
      df = test_result$parameter,
      method = test_result$method,
      n_group1 = length(group1_data),
      n_group2 = length(group2_data),
      mean_diff = test_result$estimate[1] - test_result$estimate[2]
    ))
  }
}

perform_mannwhitney_test <- function(marker_data, stats_config) {
  groups <- stats_config$groups
  group1_data <- marker_data %>% dplyr::filter(congenics == groups[1]) %>% pull(MFI)
  group2_data <- marker_data %>% dplyr::filter(congenics == groups[2]) %>% pull(MFI)
  
  if (length(group1_data) == 0 || length(group2_data) == 0) {
    return(NULL)
  }
  
  test_result <- wilcox.test(group1_data, group2_data, 
                             paired = FALSE, conf.level = 1 - stats_config$alpha)
  
  return(list(
    test_name = "Mann-Whitney U test",
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    method = test_result$method,
    n_group1 = length(group1_data),
    n_group2 = length(group2_data)
  ))
}

perform_wilcoxon_test <- function(marker_data, stats_config) {
  groups <- stats_config$groups
  
  if (is.null(stats_config$pairing_var)) {
    warning("Pairing variable not specified for Wilcoxon signed-rank test")
    return(NULL)
  }
  
  paired_data <- marker_data %>%
    select(all_of(c(stats_config$pairing_var, "congenics", "MFI"))) %>%
    pivot_wider(names_from = congenics, values_from = MFI) %>%
    dplyr::filter(!is.na(.data[[groups[1]]]), !is.na(.data[[groups[2]]]))
  
  if (nrow(paired_data) < 2) {
    warning("Insufficient paired data for Wilcoxon test")
    return(NULL)
  }
  
  test_result <- wilcox.test(paired_data[[groups[1]]], paired_data[[groups[2]]], 
                             paired = TRUE, conf.level = 1 - stats_config$alpha)
  
  return(list(
    test_name = "Wilcoxon signed-rank test",
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    method = test_result$method,
    n_pairs = nrow(paired_data)
  ))
}

perform_anova_test <- function(marker_data, stats_config) {
  anova_data <- marker_data %>%
    select(congenics, MFI) %>%
    dplyr::filter(!is.na(MFI))
  
  if (nrow(anova_data) < 3) {
    warning("Insufficient data for ANOVA")
    return(NULL)
  }
  
  aov_result <- aov(MFI ~ congenics, data = anova_data)
  aov_summary <- summary(aov_result)
  
  result <- list(
    test_name = "One-way ANOVA",
    p_value = aov_summary[[1]][["Pr(>F)"]][1],
    f_statistic = aov_summary[[1]][["F value"]][1],
    df_between = aov_summary[[1]][["Df"]][1],
    df_within = aov_summary[[1]][["Df"]][2],
    method = "One-way ANOVA"
  )
  
  group_ns <- anova_data %>% count(congenics, name = "n")
  for (i in 1:nrow(group_ns)) {
    result[[paste0("n_", group_ns$congenics[i])]] <- group_ns$n[i]
  }
  
  if (stats_config$post_hoc && result$p_value < stats_config$alpha) {
    posthoc_result <- TukeyHSD(aov_result)
    result$posthoc <- posthoc_result
  }
  
  return(result)
}

perform_kruskal_test <- function(marker_data, stats_config) {
  anova_data <- marker_data %>%
    select(congenics, MFI) %>%
    dplyr::filter(!is.na(MFI))
  
  if (nrow(anova_data) < 3) {
    warning("Insufficient data for Kruskal-Wallis test")
    return(NULL)
  }
  
  test_result <- kruskal.test(MFI ~ congenics, data = anova_data)
  
  result <- list(
    test_name = "Kruskal-Wallis test",
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    df = test_result$parameter,
    method = test_result$method
  )
  
  group_ns <- anova_data %>% count(congenics, name = "n")
  for (i in 1:nrow(group_ns)) {
    result[[paste0("n_", group_ns$congenics[i])]] <- group_ns$n[i]
  }
  
  if (stats_config$post_hoc && result$p_value < stats_config$alpha) {
    tryCatch({
      posthoc_result <- pairwise.wilcox.test(anova_data$MFI, anova_data$congenics, 
                                             p.adjust.method = "bonferroni")
      result$posthoc <- posthoc_result
    }, error = function(e) {
      warning("Could not perform post-hoc tests")
    })
  }
  
  return(result)
}

perform_rm_anova_test <- function(marker_data, stats_config) {
  if (is.null(stats_config$pairing_var)) {
    warning("Pairing variable not specified for repeated measures ANOVA")
    return(NULL)
  }
  
  rm_data <- marker_data %>%
    select(all_of(c(stats_config$pairing_var, "congenics", "MFI"))) %>%
    dplyr::filter(!is.na(MFI)) %>%
    group_by(.data[[stats_config$pairing_var]]) %>%
    dplyr::filter(n_distinct(congenics) == length(stats_config$groups)) %>%
    ungroup()
  
  if (nrow(rm_data) < 3) {
    warning("Insufficient paired data for repeated measures ANOVA")
    return(NULL)
  }
  
  aov_formula <- as.formula(paste("MFI ~ congenics + Error(", stats_config$pairing_var, ")"))
  aov_result <- aov(aov_formula, data = rm_data)
  aov_summary <- summary(aov_result)
  
  if (length(aov_summary) >= 2 && "Pr(>F)" %in% names(aov_summary[[2]][[1]])) {
    p_value <- aov_summary[[2]][[1]][["Pr(>F)"]][1]
    f_stat <- aov_summary[[2]][[1]][["F value"]][1]
    df1 <- aov_summary[[2]][[1]][["Df"]][1]
    df2 <- aov_summary[[2]][[1]][["Df"]][2]
  } else {
    warning("Could not extract statistics from repeated measures ANOVA")
    return(NULL)
  }
  
  return(list(
    test_name = "Repeated measures ANOVA",
    p_value = p_value,
    f_statistic = f_stat,
    df_between = df1,
    df_within = df2,
    method = "Repeated measures ANOVA",
    n_subjects = n_distinct(rm_data[[stats_config$pairing_var]])
  ))
}

perform_friedman_test <- function(marker_data, stats_config) {
  if (is.null(stats_config$pairing_var)) {
    warning("Pairing variable not specified for Friedman test")
    return(NULL)
  }
  
  friedman_data <- marker_data %>%
    select(all_of(c(stats_config$pairing_var, "congenics", "MFI"))) %>%
    dplyr::filter(!is.na(MFI)) %>%
    pivot_wider(names_from = congenics, values_from = MFI) %>%
    dplyr::filter(if_all(all_of(stats_config$groups), ~ !is.na(.)))
  
  if (nrow(friedman_data) < 3) {
    warning("Insufficient paired data for Friedman test")
    return(NULL)
  }
  
  test_matrix <- as.matrix(friedman_data[, stats_config$groups])
  test_result <- friedman.test(test_matrix)
  
  return(list(
    test_name = "Friedman test",
    p_value = test_result$p.value,
    statistic = test_result$statistic,
    df = test_result$parameter,
    method = test_result$method,
    n_subjects = nrow(friedman_data)
  ))
}

compile_stats_summary <- function(results_list, stats_config) {
  if (length(results_list) == 0) {
    return(NULL)
  }
  
  summary_df <- map_dfr(results_list, function(result) {
    data.frame(
      tissue = as.character(result$tissue),
      marker = as.character(result$marker),
      test_name = result$test_name,
      p_value = result$p_value,
      significant = result$p_value < stats_config$alpha,
      stars = case_when(
        result$p_value < 0.001 ~ "***",
        result$p_value < 0.01 ~ "**", 
        result$p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      stringsAsFactors = FALSE
    )
  }, .id = "tissue_marker_combo")
  
  summary_df <- summary_df %>%
    mutate(
      p_formatted = case_when(
        p_value < 0.001 ~ "< 0.001",
        p_value < 0.01 ~ sprintf("%.3f", p_value),
        TRUE ~ sprintf("%.3f", p_value)
      )
    )
  
  return(summary_df)
}

display_stats_summary <- function(stats_results) {
  if (is.null(stats_results) || is.null(stats_results$summary)) {
    return()
  }
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("STATISTICAL TESTING RESULTS (BY TISSUE)\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")
  
  summary_df <- stats_results$summary
  config <- stats_results$config
  
  cat("Test performed:", summary_df$test_name[1], "\n")
  cat("Significance level: α =", config$alpha, "\n")
  cat("Total tissue-marker combinations tested:", nrow(summary_df), "\n")
  
  n_significant <- sum(summary_df$significant)
  cat("Significant results:", n_significant, "/", nrow(summary_df), 
      sprintf("(%.1f%%)\n", n_significant/nrow(summary_df)*100))
  
  cat("\n--- Summary by Tissue ---\n")
  tissue_summary <- summary_df %>%
    group_by(tissue) %>%
    summarise(
      total_tests = n(),
      significant_tests = sum(significant),
      pct_significant = round(mean(significant) * 100, 1),
      .groups = "drop"
    )
  
  for (i in 1:nrow(tissue_summary)) {
    row <- tissue_summary[i, ]
    cat(sprintf("%-20s: %d/%d significant (%.1f%%)\n", 
                row$tissue, row$significant_tests, row$total_tests, row$pct_significant))
  }
  
  cat("\n--- Detailed Results by Tissue ---\n")
  
  summary_display <- summary_df %>%
    arrange(tissue, p_value) %>%
    mutate(
      result_summary = case_when(
        significant ~ paste0("✓ ", p_formatted, " ", stars),
        TRUE ~ paste0("  ", p_formatted, " ", stars)
      )
    )
  
  current_tissue <- ""
  for (i in 1:nrow(summary_display)) {
    row <- summary_display[i, ]
    
    if (row$tissue != current_tissue) {
      cat(sprintf("\n%s:\n", row$tissue))
      current_tissue <- row$tissue
    }
    
    cat(sprintf("  %-23s %s\n", row$marker, row$result_summary))
  }
  
  if (nrow(summary_df) > 1) {
    cat("\n⚠️  Note: Multiple comparisons performed across tissues and markers.\n")
    cat("   Consider adjusting for multiple testing if performing many comparisons.\n")
  }
  
  cat("\nLegend: *** p < 0.001, ** p < 0.01, * p < 0.05\n")
  cat("Each test was performed independently within each tissue.\n")
}

# ===== HEATMAP CREATION FUNCTIONS =====

create_mfi_heatmaps_with_stats <- function(mfi_data, 
                                           selected_congenics = NULL,
                                           selected_markers = NULL,
                                           grouping_option = list(type = "separate"),
                                           scale_method = "none",
                                           aggregation_method = "mean",
                                           cluster_rows = TRUE, 
                                           cluster_cols = TRUE,
                                           show_row_names = TRUE,
                                           show_column_names = TRUE,
                                           show_cell_values = FALSE,
                                           log_base = 2,
                                           percentile_range = c(0.05, 0.95),
                                           stats_config = NULL) {
  
  required_cols <- c("marker", "MFI", "tissue_factor", "congenics", "NodeShort")
  missing_cols <- setdiff(required_cols, names(mfi_data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
  }
  
  mfi_clean <- mfi_data %>%
    dplyr::filter(!is.na(congenics))
  
  if (!is.null(selected_congenics)) {
    mfi_clean <- mfi_clean %>%
      dplyr::filter(congenics %in% selected_congenics)
    cat("Filtered to selected congenics:", paste(selected_congenics, collapse = ", "), "\n")
  }
  
  if (!is.null(selected_markers)) {
    mfi_clean <- mfi_clean %>%
      dplyr::filter(marker %in% selected_markers)
    cat("Filtered to selected markers:", paste(selected_markers, collapse = ", "), "\n")
  }
  
  if (nrow(mfi_clean) == 0) {
    stop("No data remaining after filtering")
  }
  
  stats_results <- NULL
  if (!is.null(stats_config) && stats_config$perform_stats) {
    stats_results <- perform_statistical_analysis(mfi_clean, selected_congenics, 
                                                  selected_markers, stats_config)
    
    if (!is.null(stats_results)) {
      display_stats_summary(stats_results)
    }
  }
  
  agg_fun <- if (aggregation_method == "mean") mean else median
  legend_suffix <- if (aggregation_method == "mean") "Mean" else "Median"
  
  if (grouping_option$type == "separate") {
    return(create_separate_heatmaps_with_stats(mfi_clean, scale_method, aggregation_method, 
                                               agg_fun, legend_suffix, cluster_rows, cluster_cols,
                                               show_row_names, show_column_names, show_cell_values,
                                               log_base, percentile_range, stats_results, stats_config))
  } else {
    return(create_combined_heatmap_with_stats(mfi_clean, scale_method, aggregation_method,
                                              agg_fun, legend_suffix, cluster_rows, cluster_cols,
                                              show_row_names, show_column_names, show_cell_values,
                                              log_base, percentile_range, stats_results, stats_config))
  }
}

create_separate_heatmaps_with_stats <- function(mfi_clean, scale_method, aggregation_method,
                                                agg_fun, legend_suffix, cluster_rows, cluster_cols,
                                                show_row_names, show_column_names, show_cell_values,
                                                log_base, percentile_range, stats_results, stats_config) {
  
  group_names <- unique(mfi_clean$tissue_factor)
  cat("Creating separate heatmaps for", length(group_names), "tissue(s):", paste(group_names, collapse = ", "), "\n")
  
  heatmap_list <- list()
  
  for (group in group_names) {
    cat("Processing tissue:", group, "\n")
    
    group_data <- mfi_clean %>%
      dplyr::filter(tissue_factor == group)
    
    if (nrow(group_data) == 0) {
      warning(paste("No data found for tissue:", group))
      next
    }
    
    heatmap_matrix <- group_data %>%
      group_by(marker, congenics) %>%
      summarise(agg_MFI = agg_fun(MFI, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = congenics, values_from = agg_MFI) %>%
      column_to_rownames("marker") %>%
      as.matrix()
    
    if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
      warning(paste("Empty matrix for tissue:", group))
      next
    }
    
    result <- create_single_heatmap_with_stats(heatmap_matrix, as.character(group), scale_method, legend_suffix,
                                               cluster_rows, cluster_cols, show_row_names, 
                                               show_column_names, show_cell_values, log_base, 
                                               percentile_range, stats_results, stats_config)
    
    heatmap_list[[group]] <- result$heatmap
    
    cat(sprintf("  - Created heatmap with %d markers and %d congenics\n", 
                nrow(heatmap_matrix), ncol(heatmap_matrix)))
  }
  
  return(heatmap_list)
}

create_combined_heatmap_with_stats <- function(mfi_clean, scale_method, aggregation_method,
                                               agg_fun, legend_suffix, cluster_rows, cluster_cols,
                                               show_row_names, show_column_names, show_cell_values,
                                               log_base, percentile_range, stats_results, stats_config) {
  
  cat("Creating combined heatmap for all tissues\n")
  
  unique_tissues <- sort(unique(mfi_clean$tissue_factor))
  unique_congenics <- sort(unique(mfi_clean$congenics))
  
  cat("Tissues found:", paste(unique_tissues, collapse = ", "), "\n")
  cat("Congenics found:", paste(unique_congenics, collapse = ", "), "\n")
  
  ordered_columns <- c()
  tissue_groups <- list()
  
  for (tissue in unique_tissues) {
    tissue_cols <- c()
    for (congenic in unique_congenics) {
      col_name <- paste(tissue, congenic, sep = "_")
      ordered_columns <- c(ordered_columns, col_name)
      tissue_cols <- c(tissue_cols, col_name)
    }
    tissue_groups[[tissue]] <- tissue_cols
  }
  
  heatmap_data <- mfi_clean %>%
    group_by(marker, tissue_factor, congenics) %>%
    summarise(agg_MFI = agg_fun(MFI, na.rm = TRUE), .groups = "drop") %>%
    mutate(tissue_congenic = paste(tissue_factor, congenics, sep = "_"))
  
  heatmap_matrix <- heatmap_data %>%
    select(marker, tissue_congenic, agg_MFI) %>%
    pivot_wider(names_from = tissue_congenic, values_from = agg_MFI) %>%
    column_to_rownames("marker") %>%
    as.matrix()
  
  missing_cols <- setdiff(ordered_columns, colnames(heatmap_matrix))
  if (length(missing_cols) > 0) {
    missing_matrix <- matrix(NA, nrow = nrow(heatmap_matrix), ncol = length(missing_cols))
    colnames(missing_matrix) <- missing_cols
    rownames(missing_matrix) <- rownames(heatmap_matrix)
    heatmap_matrix <- cbind(heatmap_matrix, missing_matrix)
  }
  
  heatmap_matrix <- heatmap_matrix[, ordered_columns, drop = FALSE]
  
  if (nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
    stop("Empty combined matrix")
  }
  
  cat("Matrix dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix), "\n")
  
  col_data <- data.frame(
    tissue_congenic = colnames(heatmap_matrix),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      tissue = str_extract(tissue_congenic, "^[^_]+"),
      congenic = str_extract(tissue_congenic, "[^_]+$")
    ) %>%
    column_to_rownames("tissue_congenic")
  
  if (length(unique_tissues) <= 11) {
    tissue_colors <- RColorBrewer::brewer.pal(max(3, length(unique_tissues)), "Set3")[1:length(unique_tissues)]
  } else {
    tissue_colors <- rainbow(length(unique_tissues))
  }
  names(tissue_colors) <- unique_tissues
  
  congenic_colors <- c("CD45.1" = "lightblue", "CD45.2" = "lightcoral", "CD45.1.2" = "lightgreen")
  congenic_colors <- congenic_colors[names(congenic_colors) %in% unique_congenics]
  if (length(unique_congenics) > length(congenic_colors)) {
    additional_colors <- rainbow(length(unique_congenics) - length(congenic_colors))
    names(additional_colors) <- setdiff(unique_congenics, names(congenic_colors))
    congenic_colors <- c(congenic_colors, additional_colors)
  }
  
  tissue_splits <- rep(unique_tissues, each = length(unique_congenics))
  names(tissue_splits) <- colnames(heatmap_matrix)
  
  # Create cell-level significance annotations if available
  sig_matrix <- NULL
  if (!is.null(stats_results) && !is.null(stats_config) && 
      stats_config$sig_display != "none") {
    
    sig_matrix <- create_combined_significance_matrix(stats_results, rownames(heatmap_matrix), 
                                                      unique_tissues, stats_config, colnames(heatmap_matrix))
  }
  
  ha_col <- HeatmapAnnotation(
    Tissue = col_data$tissue,
    Congenic = col_data$congenic,
    col = list(
      Tissue = tissue_colors,
      Congenic = congenic_colors
    ),
    gap = unit(1, "mm"),
    annotation_name_gp = gpar(fontsize = 10),
    which = "column"
  )
  
  scaled_result <- apply_scaling_method(heatmap_matrix, scale_method, legend_suffix, 
                                        log_base, percentile_range)
  
  # Create cell function that shows both values and significance
  cell_function <- NULL
  if (show_cell_values || !is.null(sig_matrix)) {
    cell_function <- function(j, i, x, y, width, height, fill) {
      # Show cell values if requested
      if (show_cell_values && !is.na(scaled_result$matrix[i, j])) {
        grid.text(sprintf("%.2f", scaled_result$matrix[i, j]), 
                  x, y - unit(2, "mm"), gp = gpar(fontsize = 8))
      }
      
      # Show significance annotations
      if (!is.null(sig_matrix) && sig_matrix[i, j] != "") {
        y_pos <- if (show_cell_values) y + unit(2, "mm") else y
        grid.text(sig_matrix[i, j], 
                  x, y_pos, gp = gpar(fontsize = 7, fontface = "bold", col = "black"))
      }
    }
  }
  
  ht <- Heatmap(
    scaled_result$matrix,
    name = scaled_result$legend_title,
    col = scaled_result$color_function,
    cluster_rows = cluster_rows,
    cluster_columns = FALSE,
    column_split = tissue_splits,
    column_gap = unit(1, "mm"),
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    column_names_gp = gpar(fontsize = 9, angle = 45),
    column_title = "All Tissues Combined",
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title = "Markers",
    row_title_gp = gpar(fontsize = 12),
    top_annotation = ha_col,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    ),
    cell_fun = cell_function
  )
  
  cat(sprintf("Created combined heatmap with %d markers and %d tissue-congenic combinations\n", 
              nrow(heatmap_matrix), ncol(heatmap_matrix)))
  
  return(list("Combined" = ht))
}

create_single_heatmap_with_stats <- function(heatmap_matrix, title, scale_method, legend_suffix,
                                             cluster_rows, cluster_cols, show_row_names, 
                                             show_column_names, show_cell_values, log_base, 
                                             percentile_range, stats_results, stats_config) {
  
  scaled_result <- apply_scaling_method(heatmap_matrix, scale_method, legend_suffix, 
                                        log_base, percentile_range)
  
  # Create cell-level significance annotations if available
  sig_matrix <- NULL
  if (!is.null(stats_results) && !is.null(stats_config) && 
      stats_config$sig_display != "none") {
    
    sig_annotations <- create_tissue_significance_matrix(stats_results, rownames(heatmap_matrix), 
                                                         title, stats_config)
    
    if (!is.null(sig_annotations) && length(sig_annotations) > 0) {
      # Create a matrix for cell annotations (for separate heatmaps, we show significance on all cells for significant markers)
      sig_matrix <- matrix("", nrow = nrow(heatmap_matrix), ncol = ncol(heatmap_matrix))
      rownames(sig_matrix) <- rownames(heatmap_matrix)
      colnames(sig_matrix) <- colnames(heatmap_matrix)
      
      # Fill in significance for markers that have significant results
      for (marker in names(sig_annotations)) {
        if (marker %in% rownames(sig_matrix) && sig_annotations[marker] != "") {
          sig_matrix[marker, ] <- sig_annotations[marker]
        }
      }
    }
  }
  
  # Create cell function that shows both values and significance
  cell_function <- NULL
  if (show_cell_values || !is.null(sig_matrix)) {
    cell_function <- function(j, i, x, y, width, height, fill) {
      # Show cell values if requested
      if (show_cell_values && !is.na(scaled_result$matrix[i, j])) {
        grid.text(sprintf("%.2f", scaled_result$matrix[i, j]), 
                  x, y - unit(2, "mm"), gp = gpar(fontsize = 8))
      }
      
      # Show significance annotations
      if (!is.null(sig_matrix) && sig_matrix[i, j] != "") {
        y_pos <- if (show_cell_values) y + unit(2, "mm") else y
        grid.text(sig_matrix[i, j], 
                  x, y_pos, gp = gpar(fontsize = 8, fontface = "bold", col = "black"))
      }
    }
  }
  
  ht <- Heatmap(
    scaled_result$matrix,
    name = scaled_result$legend_title,
    col = scaled_result$color_function,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_cols,
    show_row_names = show_row_names,
    show_column_names = show_column_names,
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title = "Markers",
    row_title_gp = gpar(fontsize = 12),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 10)
    ),
    cell_fun = cell_function
  )
  
  return(list(heatmap = ht, matrix = scaled_result$matrix))
}

# ===== SIGNIFICANCE ANNOTATION FUNCTIONS =====

create_combined_significance_matrix <- function(stats_results, selected_markers, tissues, stats_config, column_names) {
  if (is.null(stats_results) || is.null(stats_results$summary)) {
    cat("No stats results available for combined heatmap\n")
    return(NULL)
  }
  
  summary_df <- stats_results$summary
  summary_df$tissue <- as.character(summary_df$tissue)
  
  cat("Creating cell-specific significance annotations\n")
  cat("Available tissues:", paste(unique(summary_df$tissue), collapse = ", "), "\n")
  
  # Create a matrix to store significance for each cell
  sig_matrix <- matrix("", nrow = length(selected_markers), ncol = length(column_names))
  rownames(sig_matrix) <- selected_markers
  colnames(sig_matrix) <- column_names
  
  # For each statistical test result, determine which cells should show significance
  for (i in 1:nrow(summary_df)) {
    test_result <- summary_df[i, ]
    marker <- test_result$marker
    tissue <- test_result$tissue
    
    if (!marker %in% selected_markers) next
    
    # Find columns that correspond to this tissue
    tissue_columns <- column_names[grepl(paste0("^", tissue, "_"), column_names)]
    
    # If significant, mark the tissue columns with the significance indicator
    if (test_result$significant) {
      sig_text <- switch(stats_config$sig_display,
                         "p_values" = test_result$p_formatted,
                         "stars" = test_result$stars,
                         "both" = paste(test_result$p_formatted, test_result$stars),
                         "")
      
      # Add significance annotation to all columns for this tissue
      for (col in tissue_columns) {
        if (col %in% colnames(sig_matrix)) {
          sig_matrix[marker, col] <- sig_text
        }
      }
    }
  }
  
  cat(sprintf("Created cell-specific annotations for %d markers across %d columns\n", 
              nrow(sig_matrix), ncol(sig_matrix)))
  
  return(sig_matrix)
}

create_tissue_significance_matrix <- function(stats_results, selected_markers, current_tissue, stats_config) {
  if (is.null(stats_results) || is.null(stats_results$summary)) {
    return(NULL)
  }
  
  summary_df <- stats_results$summary
  summary_df$tissue <- as.character(summary_df$tissue)
  current_tissue <- as.character(current_tissue)
  
  tissue_summary <- summary_df %>%
    dplyr::filter(tissue == current_tissue)
  
  if (nrow(tissue_summary) == 0) {
    return(NULL)
  }
  
  if (stats_config$sig_display == "p_values") {
    annotation_values <- setNames(tissue_summary$p_formatted, tissue_summary$marker)
  } else if (stats_config$sig_display == "stars") {
    annotation_values <- setNames(tissue_summary$stars, tissue_summary$marker)
  } else if (stats_config$sig_display == "both") {
    annotation_values <- setNames(
      paste0(tissue_summary$p_formatted, " ", tissue_summary$stars), 
      tissue_summary$marker
    )
  } else {
    return(NULL)
  }
  
  return(annotation_values[names(annotation_values) %in% selected_markers])
}

# ===== SCALING FUNCTION =====

apply_scaling_method <- function(heatmap_matrix, scale_method, legend_suffix, 
                                 log_base, percentile_range) {
  
  scaled_matrix <- heatmap_matrix
  legend_title <- paste(legend_suffix, "MFI")
  
  if (scale_method == "row") {
    scaled_matrix <- t(scale(t(heatmap_matrix)))
    legend_title <- "Z-score"
  } else if (scale_method == "global") {
    global_mean <- mean(heatmap_matrix, na.rm = TRUE)
    global_sd <- sd(as.vector(heatmap_matrix), na.rm = TRUE)
    scaled_matrix <- (heatmap_matrix - global_mean) / global_sd
    legend_title <- "Global Z-score"
  } else if (scale_method == "log") {
    min_val <- min(heatmap_matrix[heatmap_matrix > 0], na.rm = TRUE)
    constant <- min_val / 10
    scaled_matrix <- log(heatmap_matrix + constant, base = log_base)
    legend_title <- paste0("Log", log_base, " MFI")
  } else if (scale_method == "sqrt") {
    scaled_matrix <- sqrt(pmax(heatmap_matrix, 0))
    legend_title <- "√MFI"
  } else if (scale_method == "percentile") {
    lower_perc <- quantile(heatmap_matrix, percentile_range[1], na.rm = TRUE)
    upper_perc <- quantile(heatmap_matrix, percentile_range[2], na.rm = TRUE)
    scaled_matrix <- pmax(pmin(heatmap_matrix, upper_perc), lower_perc)
    scaled_matrix <- (scaled_matrix - lower_perc) / (upper_perc - lower_perc)
    legend_title <- paste0("Percentile (", percentile_range[1]*100, "-", percentile_range[2]*100, "%)")
  }
  
  if (scale_method %in% c("row", "global")) {
    max_val <- max(abs(scaled_matrix), na.rm = TRUE)
    if (max_val == 0) max_val <- 1
    col_fun <- colorRamp2(c(-max_val, 0, max_val), c("blue", "white", "red"))
  } else if (scale_method == "percentile") {
    col_fun <- colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
  } else {
    min_val <- min(scaled_matrix, na.rm = TRUE)
    max_val <- max(scaled_matrix, na.rm = TRUE)
    
    if (min_val == max_val) {
      col_fun <- colorRamp2(c(min_val - 0.1, min_val, min_val + 0.1), 
                            c("lightblue", "white", "lightcoral"))
    } else {
      if (scale_method %in% c("log", "sqrt")) {
        col_fun <- colorRamp2(
          c(min_val, 
            min_val + 0.3 * (max_val - min_val),
            min_val + 0.6 * (max_val - min_val),
            max_val),
          c("darkblue", "lightblue", "yellow", "red")
        )
      } else {
        col_fun <- colorRamp2(
          seq(min_val, max_val, length.out = 9),
          rev(RColorBrewer::brewer.pal(9, "Spectral"))
        )
      }
    }
  }
  
  return(list(
    matrix = scaled_matrix,
    legend_title = legend_title,
    color_function = col_fun
  ))
}

# ===== MAIN INTERACTIVE FUNCTION =====

create_mfi_heatmaps_interactive_enhanced <- function(mfi_data, auto_save = FALSE, ...) {
  
  cat("=== Enhanced Interactive MFI Heatmap Creation with Statistics ===\n")
  
  while(TRUE) {
    # Step 1: Select congenics
    selected_congenics <- select_congenics_interactive(mfi_data)
    if(is.null(selected_congenics)) {
      cat("Congenic selection cancelled.\n")
      return(NULL)
    }
    
    # Step 2: Select grouping option
    grouping_option <- select_grouping_option(mfi_data)
    
    # Step 3: Select markers
    filtered_data <- mfi_data %>%
      dplyr::filter(!is.na(congenics), congenics %in% selected_congenics)
    
    selected_markers <- select_markers_interactive(filtered_data)
    if(is.null(selected_markers)) {
      cat("Marker selection cancelled.\n")
      return(NULL)
    }
    
    # Step 4: Select statistical testing
    stats_config <- select_statistical_test_interactive(mfi_data, selected_congenics, selected_markers)
    
    # Step 5: Get scaling method
    scaling_params <- select_scaling_method_interactive()
    
    if(is.null(scaling_params)) {
      cat("Scaling method selection cancelled.\n")
      return(NULL)
    }
    
    if(is.character(scaling_params) && scaling_params == "DIAGNOSTIC") {
      cat("\nRunning data diagnostic...\n")
      diagnose_mfi_data(mfi_data)
      cat("\nPress Enter to continue...")
      readline()
      next
    }
    
    # Step 6: Create heatmaps
    cat("\nCreating heatmaps with your selections...\n")
    cat("- Congenics:", paste(selected_congenics, collapse = ", "), "\n")
    cat("- Markers:", length(selected_markers), "selected\n")
    cat("- Grouping:", grouping_option$type, "\n")
    cat("- Scaling:", scaling_params$scale_method, "\n")
    if(stats_config$perform_stats) {
      cat("- Statistics:", stats_config$test_type, "\n")
    } else {
      cat("- Statistics: None\n")
    }
    
    # Combine all parameters
    all_params <- c(
      list(
        mfi_data = mfi_data,
        selected_congenics = selected_congenics,
        selected_markers = selected_markers,
        grouping_option = grouping_option,
        stats_config = stats_config
      ),
      scaling_params,
      list(...)
    )
    
    heatmaps <- do.call(create_mfi_heatmaps_with_stats, all_params)
    
    if(length(heatmaps) > 0) {
      cat(sprintf("\nDisplaying heatmap: %s\n", names(heatmaps)[1]))
      draw(heatmaps[[1]])
      
      # Interactive save option
      if(!auto_save) {
        cat("\n=== Save MFI Heatmaps ===\n")
        save_choice <- menu(c("Save all heatmaps", "Save selected heatmaps", "Don't save"), 
                            title = "Save options for MFI heatmaps:")
        
        if(save_choice == 1) {
          # Save all heatmaps
          saved_count <- 0
          iwalk(heatmaps, function(heatmap, heatmap_name) {
            tryCatch({
              timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
              suggested_name <- paste0("mfi_heatmap_", make.names(heatmap_name))
              filename <- paste0(suggested_name, "_", timestamp, ".png")
              filepath <- file.path(here::here("out/plots/heatmaps/mfi"), filename)
              
              png(filepath, width = 12, height = 10, units = "in", res = 300)
              draw(heatmap)
              dev.off()
              
              cat("Heatmap saved to:", filepath, "\n")
              saved_count <<- saved_count + 1
              
            }, error = function(e) {
              cat("Error saving heatmap", heatmap_name, ":", e$message, "\n")
            })
          })
          cat("Saved", saved_count, "out of", length(heatmaps), "heatmaps\n")
          
        } else if(save_choice == 2) {
          # Save selected heatmaps
          cat("\nAvailable heatmaps:\n")
          iwalk(names(heatmaps), ~cat(sprintf("%d. %s\n", .y, .x)))
          
          selection <- readline("Enter heatmap numbers to save (space or comma-separated): ")
          if(selection != "") {
            tryCatch({
              if(grepl("\\s", selection)) {
                selected_indices <- as.numeric(str_trim(str_split(selection, "\\s+")[[1]]))
              } else {
                selected_indices <- as.numeric(str_trim(str_split(selection, ",")[[1]]))
              }
              
              selected_indices <- selected_indices[!is.na(selected_indices)]
              valid_indices <- selected_indices[selected_indices >= 1 & selected_indices <= length(heatmaps)]
              
              if(length(valid_indices) > 0) {
                selected_heatmaps <- heatmaps[valid_indices]
                selected_names <- names(heatmaps)[valid_indices]
                
                saved_count <- 0
                iwalk(selected_heatmaps, function(heatmap, idx) {
                  heatmap_name <- selected_names[idx]
                  
                  # Interactive save for each selected heatmap
                  cat(sprintf("\n--- Saving heatmap: %s ---\n", heatmap_name))
                  default_name <- paste0("mfi_heatmap_", make.names(heatmap_name), "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
                  filename <- readline(paste0("Enter filename (default: ", default_name, "): "))
                  if(filename == "") filename <- default_name
                  
                  if(!str_detect(filename, "\\.(png|pdf)$")) {
                    format_choice <- menu(c("PNG", "PDF"), title = "Choose format:")
                    extension <- if(format_choice == 2) ".pdf" else ".png"
                    filename <- paste0(filename, extension)
                  }
                  
                  filepath <- file.path(here::here("out/plots/heatmaps/mfi"), filename)
                  
                  tryCatch({
                    if(str_detect(filename, "\\.pdf$")) {
                      pdf(filepath, width = 12, height = 10)
                    } else {
                      png(filepath, width = 12, height = 10, units = "in", res = 300)
                    }
                    draw(heatmap)
                    dev.off()
                    
                    cat("Heatmap saved to:", filepath, "\n")
                    saved_count <<- saved_count + 1
                    
                  }, error = function(e) {
                    cat("Error saving heatmap:", e$message, "\n")
                  })
                })
                
                cat("Saved", saved_count, "selected heatmaps\n")
              }
              
            }, error = function(e) {
              cat("Invalid input format\n")
            })
          }
        }
        
        # Save statistics if available
        if(stats_config$perform_stats && exists("stats_results") && !is.null(stats_results)) {
          stats_save_choice <- menu(c("Save statistical results", "Don't save statistics"), 
                                    title = "Save heatmap statistics?")
          if(stats_save_choice == 1) {
            interactive_save_dataframe(export_stats_results(stats_results), "mfi_heatmap_statistics", "statistics")
          }
        }
      }
      
      return(heatmaps)
    } else {
      cat("No heatmaps were created. Please check your selections.\n")
      return(NULL)
    }
  }
}

#===============================================================================
# Engraftment plot of congenic markers with Save (USE)
#===============================================================================

create_engraftment_plot <- function(data, 
                                    wellid_col = "pairing_factor",
                                    group_col = "tissue_factor", 
                                    marker_col = "congenics",
                                    freq_col = "Freq",
                                    ko_marker = "CD45.1.2",
                                    wt_marker = "CD45.1",
                                    spleen_group = "Spleen",
                                    fill_colors = c("Spleen" = "white", "SG" = "black", "IEL" = "grey60"),
                                    plot_title = "Engraftment Ratio by Tissue",
                                    x_label = "Tissue",
                                    y_label = "Log2(ratio KO:WT in tissues/spleen)",
                                    caption_text = "*Normalized to paired spleen controls",
                                    order_by_mean = TRUE,
                                    interactive = TRUE,
                                    add_stats = TRUE,
                                    stat_method = "paired_t_test",
                                    interactive_stats = TRUE,
                                    auto_save = FALSE) {
  
  # Load required libraries with explicit namespacing
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("dplyr package required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 package required")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("tidyr package required")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("rlang package required")
  
  # Check if required columns exist
  required_cols <- c(wellid_col, group_col, marker_col, freq_col)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if(length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", "), 
               "\nAvailable columns:", paste(names(data), collapse = ", ")))
  }
  
  # Interactive statistical method selection
  comparison_type <- "pairwise"
  control_group <- NULL
  
  # Add variables to store the user's choices
  show_all_p_values <- FALSE # Default to showing only significant ones
  show_numeric_p_values <- FALSE # Default to showing asterisks
  
  if(interactive_stats && add_stats) {
    cat("\n=== Statistical Analysis Options ===\n")
    cat("1: Paired t-test (recommended for engraftment data)\n")
    cat("2: Unpaired t-test\n") 
    cat("3: ANOVA with post-hoc comparisons\n")
    cat("4: No statistical testing\n")
    
    stat_choice <- as.integer(readline(prompt = "Select statistical method (enter number): "))
    
    if(is.na(stat_choice) || stat_choice < 1 || stat_choice > 4) {
      stop("Invalid selection. Please run the function again and choose a valid number.")
    }
    
    if(stat_choice == 1) {
      stat_method <- "paired_t_test"
      cat("Selected: Paired t-test\n")
    } else if(stat_choice == 2) {
      stat_method <- "t_test"
      cat("Selected: Unpaired t-test\n")
    } else if(stat_choice == 3) {
      stat_method <- "anova"
      cat("Selected: ANOVA with post-hoc\n")
    } else {
      add_stats <- FALSE
      cat("Selected: No statistical testing\n")
    }
    
    # If doing comparisons, ask about control group
    if(add_stats && stat_choice %in% c(1, 2, 3)) {
      cat("\n=== Comparison Options ===\n")
      cat("1: All pairwise comparisons\n")
      cat("2: Compare all groups to a control group\n")
      
      comp_choice <- as.integer(readline(prompt = "Select comparison type (enter number): "))
      
      if(is.na(comp_choice) || comp_choice < 1 || comp_choice > 2) {
        comparison_type <- "pairwise"
        cat("Invalid selection - defaulting to pairwise comparisons\n")
      } else if(comp_choice == 1) {
        comparison_type <- "pairwise"
        cat("Selected: All pairwise comparisons\n")
      } else {
        comparison_type <- "vs_control"
        
        # Select control group
        available_groups <- unique(data[[group_col]])
        cat("\nAvailable groups for control:\n")
        for(i in seq_along(available_groups)) {
          cat(paste(i, ":", available_groups[i], "\n"))
        }
        
        control_choice <- as.integer(readline(prompt = "Select control group (enter number): "))
        
        if(is.na(control_choice) || control_choice < 1 || control_choice > length(available_groups)) {
          control_group <- available_groups[1]
          cat(paste("Invalid selection - defaulting to:", control_group, "\n"))
        } else {
          control_group <- available_groups[control_choice]
          cat(paste("Selected control group:", control_group, "\n"))
        }
      }
      
      # === RESTORED PROMPT: Ask whether to show all p-values or just significant ones ===
      cat("\n=== P-Value Display Options ===\n")
      cat("1: Show all p-values\n")
      cat("2: Show only significant p-values (p < 0.05)\n")
      
      p_val_choice <- as.integer(readline(prompt = "Select display option (enter number): "))
      
      if(is.na(p_val_choice) || p_val_choice < 1 || p_val_choice > 2) {
        cat("Invalid selection - defaulting to showing significant p-values only.\n")
      } else if (p_val_choice == 1) {
        show_all_p_values <- TRUE
        cat("Selected: Show all p-values\n")
      } else {
        show_all_p_values <- FALSE
        cat("Selected: Show only significant p-values\n")
      }
      
      # === NEW PROMPT: Ask whether to show numerical or symbolic p-values ===
      cat("\n=== P-Value Format Options ===\n")
      cat("1: Show numeric p-values (e.g., p = 0.045)\n")
      cat("2: Show significance asterisks (e.g., * for p < 0.05)\n")
      
      p_format_choice <- as.integer(readline(prompt = "Select format option (enter number): "))
      
      if(is.na(p_format_choice) || p_format_choice < 1 || p_format_choice > 2) {
        cat("Invalid selection - defaulting to significance asterisks.\n")
      } else if (p_format_choice == 1) {
        show_numeric_p_values <- TRUE
        cat("Selected: Show numeric p-values\n")
      } else {
        show_numeric_p_values <- FALSE
        cat("Selected: Show significance asterisks\n")
      }
      # === END NEW PROMPT ===
    }
  }
  
  # Interactive selection of normalization tissue
  normalization_tissue <- NULL
  use_normalization <- FALSE
  
  if(interactive) {
    # Show available tissue types
    available_tissues <- unique(data[[group_col]])
    cat("\nAvailable tissues for normalization:\n")
    for(i in seq_along(available_tissues)) {
      cat(paste(i, ":", available_tissues[i], "\n"))
    }
    cat(paste(length(available_tissues) + 1, ": No normalization (use raw frequencies)\n"))
    
    # Get user choice
    choice <- as.integer(readline(prompt = "Select normalization tissue (enter number): "))
    
    if(is.na(choice) || choice < 1 || choice > length(available_tissues) + 1) {
      stop("Invalid selection. Please run the function again and choose a valid number.")
    }
    
    if(choice <= length(available_tissues)) {
      normalization_tissue <- available_tissues[choice]
      use_normalization <- TRUE
      cat(paste("Selected:", normalization_tissue, "for normalization\n"))
    } else {
      use_normalization <- FALSE
      normalization_tissue <- NULL
      cat("No normalization selected - using raw frequencies\n")
    }
  } else {
    # Non-interactive mode - use provided spleen_group or no normalization if NULL
    if(!is.null(spleen_group) && spleen_group %in% unique(data[[group_col]])) {
      use_normalization <- TRUE
      normalization_tissue <- spleen_group
    } else {
      use_normalization <- FALSE
      normalization_tissue <- NULL
    }
  }
  
  # Process the data to calculate engraftment ratios
  congenic_engraftment <- data %>%
    dplyr::select(dplyr::all_of(c(wellid_col, group_col, marker_col, freq_col))) %>%
    # Extract the first letter from pairing_factor to create matching groups
    dplyr::mutate(
      well_letter = substr(.data[[wellid_col]], 1, 1)
    )
  
  # Apply normalization if selected
  if(use_normalization) {
    # First, validate that normalization tissue exists in the data
    available_tissues <- unique(congenic_engraftment[[group_col]])
    if(!normalization_tissue %in% available_tissues) {
      stop(paste("Normalization tissue '", normalization_tissue, "' not found in data.\n",
                 "Available tissues:", paste(available_tissues, collapse = ", ")))
    }
    
    # Create a lookup table for normalization frequencies
    norm_lookup <- congenic_engraftment %>%
      dplyr::filter(.data[[group_col]] == normalization_tissue) %>%
      dplyr::select(dplyr::all_of(c(marker_col, "well_letter", freq_col)))
    
    # Rename the frequency column to norm_freq
    names(norm_lookup)[names(norm_lookup) == freq_col] <- "norm_freq"
    
    # Join normalization frequencies with the main data
    congenic_engraftment <- congenic_engraftment %>%
      dplyr::left_join(norm_lookup, by = c(marker_col, "well_letter")) %>%
      # Check for missing normalization values
      dplyr::filter(!is.na(norm_freq)) %>%
      dplyr::mutate(
        normalized_freq = ifelse(norm_freq == 0, .data[[freq_col]], .data[[freq_col]] / norm_freq)
      ) %>%
      dplyr::select(dplyr::all_of(c(wellid_col, group_col, marker_col)), normalized_freq) %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(marker_col), values_from = normalized_freq)
    
    # Update plot labels for normalized data
    y_label <- paste0("Log2(ratio KO:WT normalized to ", normalization_tissue, ")")
    caption_text <- paste0("*Normalized to paired ", normalization_tissue, " controls")
    
  } else {
    congenic_engraftment <- congenic_engraftment %>%
      dplyr::ungroup() %>%
      dplyr::select(dplyr::all_of(c(wellid_col, group_col, marker_col, freq_col))) %>%
      tidyr::pivot_wider(names_from = dplyr::all_of(marker_col), values_from = dplyr::all_of(freq_col))
    
    # Update plot labels for non-normalized data
    y_label <- "Log2(ratio KO:WT - raw frequencies)"
    caption_text <- "*Using raw frequencies (no normalization)"
  }
  
  # Calculate engraftment ratio
  congenic_engraftment <- congenic_engraftment %>%
    dplyr::mutate(
      engraftment_ratio = log2(.data[[ko_marker]] / .data[[wt_marker]])
    ) %>%
    # Filter out infinite values
    dplyr::filter(is.finite(engraftment_ratio)) %>%
    dplyr::select(dplyr::all_of(c(wellid_col, group_col)), engraftment_ratio)
  
  # Check if we have data after processing
  if(nrow(congenic_engraftment) == 0) {
    stop("No valid engraftment ratios calculated. Check your marker names and data.")
  }
  
  # Calculate summary statistics
  summary_stats <- congenic_engraftment %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_col))) %>%
    dplyr::summarise(
      mean_val = mean(engraftment_ratio, na.rm = TRUE),
      sd_val = stats::sd(engraftment_ratio, na.rm = TRUE),
      n = dplyr::n(),
      sem_val = sd_val / sqrt(n),
      .groups = 'drop'
    )
  
  # Rename the grouping column for consistency
  names(summary_stats)[1] <- group_col
  
  # Create factor levels
  all_groups <- unique(summary_stats[[group_col]])
  
  if(use_normalization && normalization_tissue %in% all_groups) {
    # Put normalization tissue first
    other_groups <- setdiff(all_groups, normalization_tissue)
    
    if(order_by_mean) {
      # Order other groups by mean (excluding normalization tissue)
      other_stats <- summary_stats %>% 
        dplyr::filter(.data[[group_col]] != normalization_tissue)
      other_groups_ordered <- other_stats[[group_col]][order(-other_stats$mean_val)]
      factor_levels <- c(normalization_tissue, other_groups_ordered)
    } else {
      factor_levels <- c(normalization_tissue, other_groups)
    }
  } else {
    # No normalization, just order by mean if requested
    if(order_by_mean) {
      factor_levels <- summary_stats[[group_col]][order(-summary_stats$mean_val)]
    } else {
      factor_levels <- all_groups
    }
  }
  
  # Apply factor levels to both summary stats and original data
  summary_stats <- summary_stats %>%
    dplyr::mutate(!!rlang::sym(group_col) := factor(.data[[group_col]], levels = factor_levels))
  
  congenic_engraftment <- congenic_engraftment %>%
    dplyr::mutate(!!rlang::sym(group_col) := factor(.data[[group_col]], levels = factor_levels))
  
  # Calculate error bar positions
  summary_stats <- summary_stats %>%
    dplyr::mutate(
      error_ymin = ifelse(mean_val >= 0, mean_val, mean_val - sem_val),
      error_ymax = ifelse(mean_val >= 0, mean_val + sem_val, mean_val)
    )
  
  # Perform statistical tests if requested
  stat_results <- NULL
  sig_data <- NULL
  
  if(add_stats && nrow(summary_stats) > 1) {
    
    cat("\n=== DEBUG: Starting statistical testing ===\n")
    cat("add_stats:", add_stats, "\n")
    cat("stat_method:", stat_method, "\n")
    cat("Number of summary groups:", nrow(summary_stats), "\n")
    
    # Prepare data for statistical testing
    stat_data <- congenic_engraftment %>%
      dplyr::filter(!is.na(engraftment_ratio) & is.finite(engraftment_ratio)) %>%
      dplyr::select(dplyr::all_of(c(wellid_col, group_col)), engraftment_ratio) %>%
      dplyr::rename(
        well_id = dplyr::all_of(wellid_col),
        group = dplyr::all_of(group_col),
        value = engraftment_ratio
      ) %>%
      # Add pairing information for paired tests
      dplyr::mutate(
        pair_id = substr(well_id, 1, 1)  # Extract letter (A, B, C, etc.)
      )
    
    cat("Stat data dimensions:", nrow(stat_data), "x", ncol(stat_data), "\n")
    cat("Unique groups:", paste(unique(stat_data$group), collapse = ", "), "\n")
    cat("Unique pairs:", paste(unique(stat_data$pair_id), collapse = ", "), "\n")
    
    # Check if we have enough data for statistics
    group_counts <- stat_data %>% 
      dplyr::count(group) %>%
      dplyr::filter(n >= 2)
    
    cat("Groups with >=2 observations:\n")
    print(group_counts)
    
    if(nrow(group_counts) >= 2) {
      
      # Filter to only include groups with enough data
      stat_data <- stat_data %>%
        dplyr::filter(group %in% group_counts$group)
      
      cat("Entering statistical testing block...\n")
      
      tryCatch({
        if(stat_method == "paired_t_test") {
          
          cat("Performing paired t-tests...\n")
          
          if(comparison_type == "pairwise") {
            groups <- unique(stat_data$group)
            cat("Groups for pairwise comparison:", paste(groups, collapse = ", "), "\n")
            stat_results <- data.frame()
            
            for(i in 1:(length(groups)-1)) {
              for(j in (i+1):length(groups)) {
                group1 <- groups[i]
                group2 <- groups[j]
                
                cat("Comparing", group1, "vs", group2, "\n")
                
                # Get data for both groups
                data1 <- stat_data %>% dplyr::filter(group == group1)
                data2 <- stat_data %>% dplyr::filter(group == group2)
                
                cat("  Group1 n =", nrow(data1), ", Group2 n =", nrow(data2), "\n")
                
                # Check if we have matching pairs
                common_pairs <- intersect(data1$pair_id, data2$pair_id)
                cat("  Common pairs:", paste(common_pairs, collapse = ", "), "\n")
                
                if(length(common_pairs) >= 3) { # Need at least 3 pairs
                  # Match pairs and perform paired t-test
                  matched_data1 <- data1 %>% 
                    dplyr::filter(pair_id %in% common_pairs) %>%
                    dplyr::arrange(pair_id)
                  matched_data2 <- data2 %>% 
                    dplyr::filter(pair_id %in% common_pairs) %>%
                    dplyr::arrange(pair_id)
                  
                  cat("  Matched pairs - Group1:", nrow(matched_data1), ", Group2:", nrow(matched_data2), "\n")
                  
                  test_result <- stats::t.test(matched_data1$value, matched_data2$value, 
                                               paired = TRUE)
                  
                  cat("  T-test p-value:", test_result$p.value, "\n")
                  
                  stat_results <- rbind(stat_results, data.frame(
                    group1 = group1,
                    group2 = group2,
                    p = test_result$p.value,
                    method = "Paired t-test",
                    stringsAsFactors = FALSE
                  ))
                } else {
                  cat("  Not enough common pairs (need >=3, have", length(common_pairs), ")\n")
                }
              }
            }
            
          } else if(comparison_type == "vs_control" && !is.null(control_group)) {
            
            cat("Performing vs control comparisons with control group:", control_group, "\n")
            
            # Compare all groups to control
            control_data <- stat_data %>% dplyr::filter(group == control_group)
            other_groups <- setdiff(unique(stat_data$group), control_group)
            
            cat("Control group data n =", nrow(control_data), "\n")
            cat("Other groups:", paste(other_groups, collapse = ", "), "\n")
            
            stat_results <- data.frame()
            
            for(test_group in other_groups) {
              cat("Comparing", test_group, "vs", control_group, "\n")
              
              test_data <- stat_data %>% dplyr::filter(group == test_group)
              common_pairs <- intersect(control_data$pair_id, test_data$pair_id)
              
              cat("  Test group n =", nrow(test_data), "\n")
              cat("  Common pairs:", paste(common_pairs, collapse = ", "), "\n")
              
              if(length(common_pairs) >= 3) {
                matched_control <- control_data %>% 
                  dplyr::filter(pair_id %in% common_pairs) %>%
                  dplyr::arrange(pair_id)
                matched_test <- test_data %>% 
                  dplyr::filter(pair_id %in% common_pairs) %>%
                  dplyr::arrange(pair_id)
                
                cat("  Matched pairs - Control:", nrow(matched_control), ", Test:", nrow(matched_test), "\n")
                
                test_result <- stats::t.test(matched_test$value, matched_control$value, 
                                             paired = TRUE)
                
                cat("  T-test p-value:", test_result$p.value, "\n")
                
                stat_results <- rbind(stat_results, data.frame(
                  group1 = test_group,
                  group2 = control_group,
                  p = test_result$p.value,
                  method = "Paired t-test vs control",
                  stringsAsFactors = FALSE
                ))
              } else {
                cat("  Not enough common pairs (need >=3, have", length(common_pairs), ")\n")
              }
            }
          }
          
          cat("Statistical results before p-adjustment:\n")
          print(stat_results)
          
          # Apply multiple testing correction
          if(nrow(stat_results) > 0) {
            stat_results$p.adj <- stats::p.adjust(stat_results$p, method = "bonferroni")
            cat("Statistical results after p-adjustment:\n")
            print(stat_results)
          }
          
        } else if(stat_method == "t_test") {
          cat("Performing unpaired t-tests...\n")
          # Unpaired t-test code here...
          
        } else if(stat_method == "anova") {
          cat("Performing ANOVA...\n")
          # ANOVA code here...
        }
        
      }, error = function(e) {
        cat("ERROR in statistical testing:", e$message, "\n")
        stat_results <- NULL
      })
    } else {
      cat("Not enough groups for statistical testing\n")
    }
    
    # === P-VALUE DISPLAY LOGIC ===
    if(!is.null(stat_results) && is.data.frame(stat_results) && nrow(stat_results) > 0) {
      
      # Filter based on 'show_all_p_values' choice
      if(!show_all_p_values) {
        stat_results <- dplyr::filter(stat_results, p.adj < 0.05)
      }
      
      # Format based on 'show_numeric_p_values' choice
      if(show_numeric_p_values) {
        sig_data <- stat_results %>%
          dplyr::mutate(
            p_display = paste("p =", formatC(p.adj, format = "g", digits = 2)),
            group1 = as.character(group1),
            group2 = as.character(group2)
          )
      } else {
        sig_data <- stat_results %>%
          dplyr::mutate(
            p_display = dplyr::case_when(
              p.adj < 0.001 ~ "***",
              p.adj < 0.01 ~ "**", 
              p.adj < 0.05 ~ "*",
              TRUE ~ ""
            ),
            group1 = as.character(group1),
            group2 = as.character(group2)
          )
      }
    }
    # === END P-VALUE DISPLAY LOGIC ===
  }
  
  # Create base plot
  p <- ggplot2::ggplot(summary_stats, ggplot2::aes(x = .data[[group_col]], y = mean_val)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data[[group_col]]), alpha = 0.8, 
                      color = "black", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = error_ymin, ymax = error_ymax), 
                           width = 0.3, color = "black", linewidth = 0.3) +
    ggplot2::geom_jitter(data = congenic_engraftment, 
                         ggplot2::aes(x = .data[[group_col]], y = engraftment_ratio), 
                         size = 2, alpha = 0.7, width = 0.2, height = 0,
                         inherit.aes = FALSE) +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = y_label,
      caption = caption_text
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.3),
      panel.border = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5)
  
  # Add fill colors if provided
  if(!is.null(fill_colors)) {
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  }
  
  # Add statistical annotations if available
  if(!is.null(sig_data) && nrow(sig_data) > 0) {
    
    cat("\n=== DEBUG: Adding", nrow(sig_data), "annotations to plot ===\n")
    # Calculate y positions for significance bars
    max_y <- max(summary_stats$error_ymax, na.rm = TRUE)
    y_step <- (max_y - min(summary_stats$error_ymin, na.rm = TRUE)) * 0.1
    
    # Add significance annotations manually
    for(i in 1:nrow(sig_data)) {
      y_pos <- max_y + (i * y_step)
      
      group1_pos <- which(levels(summary_stats[[group_col]]) == sig_data$group1[i])
      group2_pos <- which(levels(summary_stats[[group_col]]) == sig_data$group2[i])
      
      p <- p + 
        ggplot2::annotate("segment", 
                          x = group1_pos, xend = group2_pos,
                          y = y_pos, yend = y_pos,
                          color = "black", linewidth = 0.3) +
        ggplot2::annotate("text", 
                          x = (group1_pos + group2_pos) / 2,
                          y = y_pos + y_step * 0.3,
                          label = sig_data$p_display[i], # This now works because p_display is created
                          size = 3)
    }
  }
  
  # Return results
  result <- list(plot = p)
  
  if(add_stats && !is.null(stat_results)) {
    result$statistics <- stat_results
    result$summary_stats <- summary_stats
  }
  
  # If only plot requested, return just the plot
  if(!add_stats) {
    return(if(is.list(result)) result$plot else result)
  } else {
    return(result)
  }
}


# ============================================================================
# UPDATED ANALYSIS FUNCTIONS USING STANDARDIZED FACTORS
# ============================================================================

#' Updated paired comparison plots using standardized factors
create_paired_comparison_plots_enhanced <- function(data, 
                                                    tissue_col = "tissue_factor",
                                                    pairing_col = "tissue_factor",
                                                    test_type = "paired_t_test", 
                                                    facet_var = "tissue_factor") {
  
  # Interactive test selection if not specified
  if(missing(test_type)) {
    test_choice <- menu(c("Paired t-test", "Wilcoxon signed-rank test"),
                        title = "Choose statistical test:")
    test_type <- c("paired_t_test", "wilcox_test")[test_choice]
  }
  
  # Interactive faceting selection if not specified
  if(missing(facet_var)) {
    facet_choice <- menu(c("tissue_factor (tissue)", "NodeShort (subpopulation of interest)", "No faceting"),
                         title = "Choose how to facet the plots:")
    facet_var <- c("tissue_factor", "NodeShort", "none")[facet_choice]
  }
  
  # Validate required columns
  required_cols <- c("NodeShort", tissue_col, pairing_col, "congenics", "Freq")
  missing_cols <- setdiff(required_cols, names(data))
  
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
  }
  
  # Check if pairing is enabled
  pairing_enabled <- !all(data[[pairing_col]] %in% c("no_pairing", "none", NA))
  
  if(!pairing_enabled && test_type == "paired_t_test") {
    cat("⚠️  Pairing not available. Switching to unpaired t-test.\n")
    test_type <- "t_test"
  }
  
  # Print information about the analysis
  test_description <- switch(test_type, 
                             'paired_t_test' = 'paired t-test', 
                             't_test' = 'unpaired t-test',
                             'wilcox_test' = 'Wilcoxon signed-rank test')
  
  facet_description <- switch(facet_var,
                              'tissue_factor' = 'tissue',
                              'NodeShort' = 'subpopulation of interest',
                              'none' = 'no faceting (separate plots)')
  
  # Determine split variable based on faceting choice
  if(facet_var == "tissue_factor") {
    split_var <- "NodeShort"
    plot_description <- "cell population(s)"
  } else if(facet_var == "NodeShort") {
    split_var <- tissue_col
    plot_description <- "tissue(s)"
  } else {
    # No faceting - create separate plots for each NodeShort-tissue combination
    split_var <- c("NodeShort", tissue_col)
    plot_description <- "NodeShort-tissue combination(s)"
  }
  
  cat("Creating paired comparison plots using", test_description, "\n")
  cat("Faceting by", facet_description, "\n")
  
  # Calculate number of plots that will be created
  if(facet_var == "none") {
    n_plots <- data %>% 
      distinct(NodeShort, .data[[tissue_col]]) %>% 
      nrow()
    cat("Processing", n_plots, plot_description, "\n\n")
  } else {
    n_plots <- n_distinct(data[[split_var]])
    cat("Processing", n_plots, plot_description, "\n\n")
  }
  
  # Helper function to create individual subgroup plots with enhanced pairing support
  create_subgroup_plot_enhanced <- function(df, test_type = "paired_t_test", facet_var = "tissue_factor") {
    
    # Check if we have enough data for statistical testing
    if(length(unique(df$congenics)) < 2) {
      warning(paste("Subgroup has less than 2 'congenics' levels for comparison. Skipping statistical test."))
      stat_test_results <- NULL
    } else {
      
      # Perform statistical test based on test_type and facet_var
      if(facet_var == "NodeShort") {
        # When faceting by NodeShort, group by NodeShort for stats
        if(test_type == "paired_t_test" && pairing_enabled) {
          stat_test_results <- df %>%
            group_by(NodeShort) %>%
            dplyr::filter(n_distinct(congenics) == 2) %>%
            pairwise_t_test(as.formula(paste("Freq ~", "congenics")), 
                            paired = TRUE, 
                            p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else if(test_type == "wilcox_test" && pairing_enabled) {
          stat_test_results <- df %>%
            group_by(NodeShort) %>%
            dplyr::filter(n_distinct(congenics) == 2) %>%
            pairwise_wilcox_test(as.formula(paste("Freq ~", "congenics")), 
                                 paired = TRUE,
                                 p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else {
          # Unpaired tests
          stat_test_results <- df %>%
            group_by(NodeShort) %>%
            dplyr::filter(n_distinct(congenics) == 2) %>%
            pairwise_t_test(as.formula(paste("Freq ~", "congenics")), 
                            paired = FALSE,
                            p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        }
      } else {
        # When faceting by tissue_factor or no faceting, group by tissue_factor
        if(test_type == "paired_t_test" && pairing_enabled) {
          stat_test_results <- df %>%
            group_by(.data[[tissue_col]]) %>%
            dplyr::filter(n_distinct(congenics) == 2) %>%
            pairwise_t_test(as.formula(paste("Freq ~", "congenics")), 
                            paired = TRUE,
                            p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else if(test_type == "wilcox_test" && pairing_enabled) {
          stat_test_results <- df %>%
            group_by(.data[[tissue_col]]) %>%
            dplyr::filter(n_distinct(congenics) == 2) %>%
            pairwise_wilcox_test(as.formula(paste("Freq ~", "congenics")), 
                                 paired = TRUE,
                                 p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else {
          # Unpaired tests
          stat_test_results <- df %>%
            group_by(.data[[tissue_col]]) %>%
            dplyr::filter(n_distinct(congenics) == 2) %>%
            pairwise_t_test(as.formula(paste("Freq ~", "congenics")), 
                            paired = FALSE,
                            p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        }
      }
    }
    
    # Calculate y-axis limits
    max_freq <- max(df$Freq, na.rm = TRUE)
    max_y_position <- if(!is.null(stat_test_results) && nrow(stat_test_results) > 0) {
      max(max_freq, stat_test_results$y.position, na.rm = TRUE)
    } else {
      max_freq
    }
    upper_y_limit <- max_y_position * 1.2
    
    # Create the base plot using pairing_col for connecting lines
    p <- ggplot(df, aes(x = congenics, y = Freq, group = pairing_factor))+
      geom_line(color = "black", linewidth = 0.4, alpha = 0.7, show.legend = FALSE) +
      geom_point(shape = 21, fill = "white", color = "black", size = 3, show.legend = FALSE)
    
    # Add faceting if requested
    if(facet_var == tissue_col) {
      p <- p + facet_wrap(as.formula(paste("~", tissue_col)), scales = "free_y", strip.position = "top")
    } else if(facet_var == "NodeShort") {
      p <- p + facet_wrap(~ NodeShort, scales = "free_y", strip.position = "top")
    }
    
    # Dynamic title and y-axis label based on faceting
    if(facet_var == tissue_col) {
      plot_title <- paste("Cell Population:", unique(df$NodeShort))
      y_label <- paste0(unique(df$NodeShort), " (%)")
    } else if(facet_var == "NodeShort") {
      plot_title <- paste("Tissue:", unique(df[[tissue_col]]))
      y_label <- "Cell Frequency (%)"
    } else {
      plot_title <- paste("Subgroup:", unique(df$NodeShort), "in", unique(df[[tissue_col]]))
      y_label <- paste0(unique(df$NodeShort), " (%)")
    }
    
    # Continue with plot formatting
    p <- p +
      scale_y_continuous(
        limits = c(0, upper_y_limit),
        expand = expansion(mult = c(0.02, 0.05))
      ) +
      labs(
        title = plot_title,
        x = "Congenic Marker",
        y = y_label
      ) +
      theme_minimal(base_size = 14) +
      theme(
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 0.4),
        axis.line.y = element_line(color = "black", linewidth = 0.4),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        strip.placement = "outside"
      )
    
    # Add statistical annotations if available
    if(!is.null(stat_test_results) && nrow(stat_test_results) > 0) {
      p <- p + stat_pvalue_manual(
        stat_test_results,
        label = "p.format",
        tip.length = 0,
        bracket.size = 0.5,
        hide.ns = FALSE
      )
    }
    
    return(p)
  }
  
  # Create plots based on the correct splitting logic
  if(facet_var == tissue_col) {
    # Split by NodeShort to get separate plots for each cell population
    # Each plot will be faceted by tissue_factor
    plots <- data %>%
      group_split(NodeShort, .keep = TRUE) %>%
      map(~ create_subgroup_plot_enhanced(.x, test_type = test_type, facet_var = tissue_col)) %>%
      set_names(map_chr(data %>% distinct(NodeShort) %>% pull(NodeShort), as.character))
  } else if(facet_var == "NodeShort") {
    # Split by tissue_factor to get separate plots for each tissue
    # Each plot will be faceted by cell population (NodeShort)  
    plots <- data %>%
      group_split(.data[[tissue_col]], .keep = TRUE) %>%
      map(~ create_subgroup_plot_enhanced(.x, test_type = test_type, facet_var = "NodeShort")) %>%
      set_names(map_chr(data %>% distinct(.data[[tissue_col]]) %>% pull(.data[[tissue_col]]), as.character))
  } else {
    # No faceting - create separate plots for each NodeShort-tissue combination
    plots <- data %>%
      group_split(NodeShort, .data[[tissue_col]], .keep = TRUE) %>%
      map(~ create_subgroup_plot_enhanced(.x, test_type = test_type, facet_var = "none")) %>%
      set_names(map_chr(1:length(.), ~ {
        df <- data %>% group_split(NodeShort, .data[[tissue_col]], .keep = TRUE) %>% .[[.x]]
        paste0(unique(df$NodeShort), "_", unique(df[[tissue_col]]))
      }))
  }
  
  # Print summary
  plot_count <- length(plots)
  cat("Successfully created", plot_count, "plot(s)", "\n")
  cat("Available plots:", paste(names(plots), collapse = ", "), "\n")
  
  return(plots)
}


# ============================================================================
# ENHANCED INTERACTIVE UMAP ANALYSIS FOR FLOW CYTOMETRY DATA
# Fixed version with sample exclusion and streamlined workflow
# ============================================================================

get_all_populations <- function(gs) {
  if(requireNamespace("flowWorkspace", quietly = TRUE)) {
    return(flowWorkspace::gs_get_pop_paths(gs))
  } else {
    stop("flowWorkspace package required. Please load flowWorkspace.")
  }
}

get_marker_lookup <- function(gs) {
  if(requireNamespace("flowWorkspace", quietly = TRUE)) {
    # Get channel info from first sample
    first_sample <- flowWorkspace::gs_cyto_data(gs)[[1]]
    markers <- flowCore::markernames(first_sample)
    channels <- names(markers)
    
    return(tibble(
      colname = channels,
      marker = ifelse(is.na(markers), channels, markers)
    ))
  } else {
    stop("flowWorkspace package required. Please load flowWorkspace.")
  }
}

resolve_nodes_by_leaf <- function(gs, leaf_name) {
  all_pops <- get_all_populations(gs)
  matches <- all_pops[grepl(paste0("\\b", leaf_name, "\\b"), basename(all_pops), ignore.case = TRUE)]
  return(matches)
}

gh_get_pop_paths <- function(gh) {
  if(requireNamespace("flowWorkspace", quietly = TRUE)) {
    return(flowWorkspace::gs_get_pop_paths(gh))
  } else {
    stop("flowWorkspace package required")
  }
}

gh_pop_get_data <- function(gh, node) {
  if(requireNamespace("flowWorkspace", quietly = TRUE)) {
    return(flowWorkspace::gh_pop_get_data(gh, node))
  } else {
    stop("flowWorkspace package required")
  }
}

# ============================================================================
# ENHANCED SAMPLE EXCLUSION FUNCTIONALITY
# ============================================================================

interactive_sample_exclusion <- function(gs) {
  all_samples <- sampleNames(gs)
  cat("\n=== Sample Exclusion Menu ===\n")
  cat("Total samples available:", length(all_samples), "\n")
  
  # Get sample metadata for informed decisions
  if(exists("pData") && is.function(pData)) {
    tryCatch({
      pd <- pData(gs) %>% 
        rownames_to_column("Sample")
      
      cat("\nSample preview with metadata:\n")
      print(head(pd %>% select(1:min(4, ncol(pd))), 10))
    }, error = function(e) {
      cat("Could not load sample metadata\n")
      pd <- tibble(Sample = all_samples)
    })
  } else {
    pd <- tibble(Sample = all_samples)
  }
  
  while(TRUE) {
    cat("\n--- Exclusion Options ---\n")
    cat("1. Exclude samples by pattern (recommended for unstained, controls)\n")
    cat("2. Exclude specific samples (interactive selection)\n")
    cat("3. Exclude samples by metadata criteria\n")
    cat("4. Preview current sample list\n")
    cat("5. Continue with all samples\n")
    cat("6. Continue with current exclusions\n")
    
    choice <- readline("Choose option (1-6): ")
    
    if(choice == "1") {
      return(exclude_samples_by_pattern(all_samples, gs))
    } else if(choice == "2") {
      return(exclude_samples_interactive(all_samples))
    } else if(choice == "3") {
      return(exclude_samples_by_metadata(gs, pd))
    } else if(choice == "4") {
      cat("\nCurrent samples:\n")
      iwalk(all_samples, ~cat(sprintf("%d. %s\n", .y, .x)))
      cat("\nPress Enter to continue...")
      readline()
    } else if(choice == "5") {
      cat("Proceeding with all", length(all_samples), "samples\n")
      return(all_samples)
    } else if(choice == "6") {
      cat("Proceeding with current sample list\n")
      return(all_samples)
    } else {
      cat("Invalid choice. Please select 1-6.\n")
    }
  }
}

exclude_samples_by_pattern <- function(all_samples, gs) {
  cat("\n=== Pattern-Based Sample Exclusion ===\n")
  cat("This is useful for removing unstained controls, compensation controls, etc.\n")
  
  excluded_samples <- character(0)
  
  while(TRUE) {
    cat("\nCurrent samples:", length(all_samples), "\n")
    cat("Currently excluded:", length(excluded_samples), "\n")
    
    if(length(excluded_samples) > 0) {
      cat("Excluded samples:", paste(head(excluded_samples, 5), collapse = ", "))
      if(length(excluded_samples) > 5) cat("...")
      cat("\n")
    }
    
    cat("\nPattern exclusion options:\n")
    cat("1. Add exclusion pattern\n")
    cat("2. Remove exclusion pattern\n")
    cat("3. Preview samples matching pattern\n")
    cat("4. Common patterns (unstained, compensation, etc.)\n")
    cat("5. Clear all exclusions\n")
    cat("6. Finish exclusions\n")
    
    choice <- readline("Choose (1-6): ")
    
    if(choice == "1") {
      pattern <- readline("Enter pattern to exclude (case-insensitive regex): ")
      if(pattern == "") {
        cat("Empty pattern, skipping\n")
        next
      }
      
      matches <- grep(pattern, all_samples, ignore.case = TRUE, value = TRUE)
      if(length(matches) == 0) {
        cat("No samples match pattern:", pattern, "\n")
      } else {
        cat("Pattern '", pattern, "' matches", length(matches), "samples:\n")
        iwalk(matches, ~cat(sprintf("  %d. %s\n", .y, .x)))
        
        confirm <- readline("Exclude these samples? (y/n): ")
        if(tolower(confirm) == "y") {
          excluded_samples <- unique(c(excluded_samples, matches))
          all_samples <- setdiff(all_samples, matches)
          cat("Excluded", length(matches), "samples\n")
        }
      }
      
    } else if(choice == "4") {
      # Common patterns
      common_patterns <- list(
        "Unstained" = "unstained|blank|negative",
        "Compensation" = "comp|compensation",
        "Controls" = "control|ctrl",
        "Single stain" = "single|ss[_-]",
        "FMO" = "fmo|minus.one"
      )
      
      cat("\nCommon exclusion patterns:\n")
      iwalk(names(common_patterns), ~cat(sprintf("%d. %s (%s)\n", .y, .x, common_patterns[[.x]])))
      
      pattern_choice <- readline("Select pattern number or enter custom: ")
      
      selected_pattern <- NULL
      
      if(grepl("^\\d+$", pattern_choice)) {
        pattern_num <- as.numeric(pattern_choice)
        if(!is.na(pattern_num) && pattern_num >= 1 && pattern_num <= length(common_patterns)) {
          selected_pattern <- common_patterns[[pattern_num]]
        }
      } else if(pattern_choice != "") {
        selected_pattern <- pattern_choice
      }
      
      if(!is.null(selected_pattern)) {
        matches <- grep(selected_pattern, all_samples, ignore.case = TRUE, value = TRUE)
        if(length(matches) > 0) {
          cat("Pattern matches", length(matches), "samples:\n")
          iwalk(matches, ~cat(sprintf("  %d. %s\n", .y, .x)))
          
          confirm <- readline("Exclude these samples? (y/n): ")
          if(tolower(confirm) == "y") {
            excluded_samples <- unique(c(excluded_samples, matches))
            all_samples <- setdiff(all_samples, matches)
            cat("Excluded", length(matches), "samples\n")
          }
        } else {
          cat("No samples match this pattern\n")
        }
      } else {
        cat("Invalid selection\n")
      }
      
    } else if(choice == "6") {
      break
    } else {
      cat("Option not yet implemented or invalid choice\n")
    }
  }
  
  remaining_samples <- setdiff(sampleNames(gs), excluded_samples)
  cat("\nFinal sample count:", length(remaining_samples), "(excluded", length(excluded_samples), ")\n")
  
  return(remaining_samples)
}

exclude_samples_interactive <- function(all_samples) {
  cat("\n=== Interactive Sample Selection for Exclusion ===\n")
  cat("Select samples to EXCLUDE from analysis\n")
  
  selected_for_exclusion <- select.list(
    all_samples, 
    multiple = TRUE, 
    title = "Select samples to EXCLUDE (Cancel to skip exclusions)"
  )
  
  if(length(selected_for_exclusion) == 0) {
    cat("No samples selected for exclusion\n")
    return(all_samples)
  }
  
  cat("\nSamples selected for exclusion:\n")
  iwalk(selected_for_exclusion, ~cat(sprintf("%d. %s\n", .y, .x)))
  
  confirm <- readline(paste0("Exclude these ", length(selected_for_exclusion), " samples? (y/n): "))
  
  if(tolower(confirm) == "y") {
    remaining_samples <- setdiff(all_samples, selected_for_exclusion)
    cat("Excluded", length(selected_for_exclusion), "samples\n")
    cat("Remaining samples:", length(remaining_samples), "\n")
    return(remaining_samples)
  } else {
    cat("Exclusion cancelled, keeping all samples\n")
    return(all_samples)
  }
}

exclude_samples_by_metadata <- function(gs, pd) {
  cat("\n=== Metadata-Based Sample Exclusion ===\n")
  
  if(ncol(pd) <= 1) {
    cat("No metadata columns available for filtering\n")
    return(sampleNames(gs))
  }
  
  metadata_cols <- setdiff(names(pd), "Sample")
  cat("Available metadata columns:\n")
  iwalk(metadata_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
  
  col_choice <- readline("Enter column number or name: ")
  
  selected_col <- NULL
  if(grepl("^\\d+$", col_choice)) {
    col_num <- as.numeric(col_choice)
    if(!is.na(col_num) && col_num >= 1 && col_num <= length(metadata_cols)) {
      selected_col <- metadata_cols[col_num]
    }
  } else if(col_choice %in% metadata_cols) {
    selected_col <- col_choice
  }
  
  if(is.null(selected_col)) {
    cat("Invalid selection\n")
    return(sampleNames(gs))
  }
  
  unique_values <- unique(pd[[selected_col]])
  cat("\nUnique values in", selected_col, ":\n")
  iwalk(unique_values, ~cat(sprintf("%d. %s\n", .y, .x)))
  
  values_to_exclude <- select.list(
    as.character(unique_values), 
    multiple = TRUE,
    title = paste("Select", selected_col, "values to EXCLUDE")
  )
  
  if(length(values_to_exclude) == 0) {
    cat("No values selected for exclusion\n")
    return(sampleNames(gs))
  }
  
  excluded_samples <- pd %>%
    dplyr::filter(.data[[selected_col]] %in% values_to_exclude) %>%
    pull(Sample)
  
  cat("Samples to exclude based on", selected_col, ":\n")
  iwalk(excluded_samples, ~cat(sprintf("%d. %s\n", .y, .x)))
  
  confirm <- readline(paste0("Exclude these ", length(excluded_samples), " samples? (y/n): "))
  
  if(tolower(confirm) == "y") {
    remaining_samples <- setdiff(sampleNames(gs), excluded_samples)
    cat("Excluded", length(excluded_samples), "samples\n")
    cat("Remaining samples:", length(remaining_samples), "\n")
    return(remaining_samples)
  } else {
    return(sampleNames(gs))
  }
}

# ============================================================================
# INTERACTIVE NODE SELECTION (STREAMLINED)
# ============================================================================

select_single_node_for_umap <- function(gs) {
  pops <- get_all_populations(gs)
  
  while(TRUE) {
    cat("\n=== Node Selection for UMAP Analysis ===\n")
    cat("1. Select from list (interactive)\n")
    cat("2. Use regex pattern\n") 
    cat("3. Use leaf names\n")
    cat("4. Exit/Cancel\n")
    
    choice <- readline("Choose method (1-4): ")
    
    if(choice == "4") {
      stop("Node selection cancelled by user.")
    }
    
    result <- switch(choice,
                     "1" = {
                       cat("\nAvailable populations:\n")
                       sel <- select.list(pops, multiple = FALSE, title = "Select ONE population for UMAP")
                       if(length(sel) == 0) {
                         cat("No population selected. Going back...\n")
                         next
                       }
                       sel
                     },
                     
                     "2" = {
                       while(TRUE) {
                         pattern <- readline("Enter regex pattern or 'back': ")
                         
                         if(tolower(pattern) == "back") break
                         if(pattern == "") {
                           cat("Empty pattern. Please try again.\n")
                           next
                         }
                         
                         matches <- grep(pattern, pops, value = TRUE, ignore.case = TRUE)
                         if(length(matches) == 0) {
                           cat("No matches found. Try a different pattern.\n")
                           retry <- readline("Try again? (y/n): ")
                           if(tolower(retry) != "y") break
                           next
                         }
                         
                         cat(sprintf("Found %d matches:\n", length(matches)))
                         iwalk(matches, ~cat(sprintf("%d. %s\n", .y, .x)))
                         
                         while(TRUE) {
                           index <- readline("Enter number of population to select, or 'back': ")
                           if(tolower(index) == "back") break
                           
                           if(grepl("^\\d+$", index)) {
                             selected_idx <- as.numeric(index)
                             if(!is.na(selected_idx) && selected_idx >= 1 && selected_idx <= length(matches)) {
                               return(matches[selected_idx])
                             }
                           }
                           cat("Invalid selection. Please try again.\n")
                         }
                       }
                       NULL
                     },
                     
                     "3" = {
                       while(TRUE) {
                         input <- readline("Enter leaf name or 'back': ")
                         if(tolower(input) == "back") break
                         if(input == "") {
                           cat("No leaf name provided. Please try again.\n")
                           next
                         }
                         
                         matches <- resolve_nodes_by_leaf(gs, input)
                         
                         if(length(matches) == 0) {
                           cat("No matches found for that leaf name.\n")
                           retry <- readline("Try again? (y/n): ")
                           if(tolower(retry) != "y") break
                           next
                         }
                         
                         if(length(matches) == 1) {
                           return(matches[1])
                         } else {
                           cat(sprintf("Found %d matches:\n", length(matches)))
                           iwalk(matches, ~cat(sprintf("%d. %s\n", .y, .x)))
                           
                           while(TRUE) {
                             index <- readline("Enter number to select, or 'back': ")
                             if(tolower(index) == "back") break
                             
                             if(grepl("^\\d+$", index)) {
                               selected_idx <- as.numeric(index)
                               if(!is.na(selected_idx) && selected_idx >= 1 && selected_idx <= length(matches)) {
                                 return(matches[selected_idx])
                               }
                             }
                             cat("Invalid selection. Please try again.\n")
                           }
                         }
                       }
                       NULL
                     },
                     
                     {
                       cat("Invalid choice. Please select 1-4.\n")
                       NULL
                     }
    )
    
    if(!is.null(result)) return(result)
  }
}

# ============================================================================
# INTERACTIVE MARKER/CHANNEL SELECTION (STREAMLINED)
# ============================================================================

select_markers_for_umap <- function(gs) {
  lookup <- get_marker_lookup(gs)
  
  exclude_patterns <- c("FSC", "SSC", "Time", "Event", "Original", "Width", "Height")
  exclude_regex <- paste0("(?i)", paste(exclude_patterns, collapse = "|"))
  
  potential_markers <- lookup %>%
    dplyr::filter(!str_detect(colname, exclude_regex)) %>%
    dplyr::filter(!str_detect(marker, exclude_regex))
  
  comp_channels <- potential_markers %>%
    dplyr::filter(str_detect(colname, "(?i)comp")) %>%
    pull(colname)
  
  while(TRUE) {
    cat("\n=== Marker Selection for UMAP Analysis ===\n")
    cat("Available markers:", nrow(potential_markers), "\n")
    if(length(comp_channels) > 0) {
      cat("Compensated channels found:", length(comp_channels), "\n")
    }
    
    cat("\nRecommendation: Use compensated channels for best results\n")
    
    cat("\nOptions:\n")
    cat("1. Select specific markers/channels\n")
    cat("2. Use all compensated channels (recommended)\n")
    cat("3. Use all available channels\n")
    cat("4. Preview marker distributions\n")
    cat("5. Back\n")
    
    choice <- readline("Choose (1-5): ")
    
    switch(choice,
           "1" = {
             choices <- paste0(potential_markers$colname, " :: ", potential_markers$marker)
             sel <- select.list(choices, multiple = TRUE, 
                                title = "Select markers for UMAP (Cancel to go back)")
             
             if(length(sel) == 0) {
               cat("No markers selected.\n")
               retry <- readline("Try again? (y/n): ")
               if(tolower(retry) != "y") next
               next
             }
             
             if(length(sel) < 2) {
               cat("UMAP requires at least 2 dimensions. Please select more markers.\n")
               retry <- readline("Try again? (y/n): ")
               if(tolower(retry) != "y") next
               next
             }
             
             selected_channels <- map_chr(sel, ~str_split(.x, " :: ")[[1]][1])
             selected_markers <- map_chr(sel, ~str_split(.x, " :: ")[[1]][2])
             
             cat(sprintf("Selected %d markers:\n", length(selected_channels)))
             iwalk(selected_channels, ~cat(sprintf("  %d. %s (%s)\n", .y, 
                                                   selected_markers[.y], .x)))
             
             confirm <- readline("Confirm selection? (y/n): ")
             if(tolower(confirm) == "y") {
               return(list(
                 channels = selected_channels,
                 markers = selected_markers,
                 lookup = potential_markers %>% dplyr::filter(colname %in% selected_channels)
               ))
             }
           },
           
           "2" = {
             if(length(comp_channels) == 0) {
               cat("No compensated channels found. Try option 1 or 3.\n")
               next
             }
             
             selected_lookup <- potential_markers %>% dplyr::filter(colname %in% comp_channels)
             
             cat(sprintf("Using %d compensated channels:\n", nrow(selected_lookup)))
             iwalk(selected_lookup$colname, ~cat(sprintf("  %d. %s (%s)\n", .y, 
                                                         selected_lookup$marker[.y], .x)))
             
             confirm <- readline("Confirm compensated channels? (y/n): ")
             if(tolower(confirm) == "y") {
               return(list(
                 channels = selected_lookup$colname,
                 markers = selected_lookup$marker,
                 lookup = selected_lookup
               ))
             }
           },
           
           "3" = {
             cat(sprintf("Using all %d available channels:\n", nrow(potential_markers)))
             cat("First 10 markers:\n")
             iwalk(head(potential_markers$marker, 10), 
                   ~cat(sprintf("  %d. %s\n", .y, .x)))
             if(nrow(potential_markers) > 10) {
               cat(sprintf("  ... and %d more\n", nrow(potential_markers) - 10))
             }
             
             confirm <- readline("Confirm all channels? (y/n): ")
             if(tolower(confirm) == "y") {
               return(list(
                 channels = potential_markers$colname,
                 markers = potential_markers$marker,
                 lookup = potential_markers
               ))
             }
           },
           
           "4" = {
             cat("\n--- Marker Distribution ---\n")
             cat("Total potential markers:", nrow(potential_markers), "\n")
             cat("Compensated channels:", length(comp_channels), "\n")
             
             cat("\nSample of available markers:\n")
             sample_markers <- potential_markers %>% 
               slice_sample(n = min(20, nrow(potential_markers)))
             iwalk(sample_markers$marker, 
                   ~cat(sprintf("  %s (%s)\n", .x, sample_markers$colname[.y])))
             
             cat("\nPress Enter to continue...")
             readline()
             next
           },
           
           "5" = {
             return("BACK")
           },
           
           {
             cat("Invalid choice. Please select 1-5.\n")
           }
    )
  }
}

# ============================================================================
# ENHANCED SINGLE-CELL DATA EXTRACTION WITH SAMPLE FILTERING
# ============================================================================

get_population_event_counts <- function(gs, node, selected_samples = NULL) {
  if(is.null(selected_samples)) {
    selected_samples <- sampleNames(gs)
  }
  
  cat("Analyzing event counts for node:", basename(node), "\n")
  cat("Analyzing", length(selected_samples), "samples\n")
  
  event_counts <- map_dfr(selected_samples, function(sample_name) {
    gh <- gs[[sample_name]]
    
    if(!node %in% gh_get_pop_paths(gh)) {
      return(tibble(Sample = sample_name, EventCount = 0))
    }
    
    tryCatch({
      ff <- gh_pop_get_data(gh, node)
      if(is.null(ff)) {
        event_count <- 0
      } else {
        event_count <- nrow(if(inherits(ff, "cytoframe")) exprs(ff) else exprs(ff))
      }
      
      tibble(Sample = sample_name, EventCount = event_count)
    }, error = function(e) {
      tibble(Sample = sample_name, EventCount = 0)
    })
  })
  
  event_summary <- list(
    counts = event_counts,
    total_events = sum(event_counts$EventCount),
    mean_events = mean(event_counts$EventCount),
    median_events = median(event_counts$EventCount),
    min_sample = event_counts[which.min(event_counts$EventCount), ],
    max_sample = event_counts[which.max(event_counts$EventCount), ],
    samples_with_zero = sum(event_counts$EventCount == 0)
  )
  
  return(event_summary)
}

extract_single_cell_data <- function(gs, node, selected_markers, 
                                     selected_samples = NULL,
                                     keywords = c("pairing_factor", "tissue_factor"),
                                     sample_limit = NULL, 
                                     max_cells_per_sample = 5000) {
  
  if(is.null(selected_samples)) {
    selected_samples <- sampleNames(gs)
  }
  
  cat("Extracting single-cell data from node:", basename(node), "\n")
  cat("Using", length(selected_samples), "samples\n")
  
  tryCatch({
    if(exists("pData") && is.function(pData)) {
      pd <- pData(gs) %>% 
        rownames_to_column("Sample") %>%
        select(Sample, any_of(keywords))
    } else {
      pd <- tibble(Sample = selected_samples)
    }
  }, error = function(e) {
    cat("Warning: Could not load sample metadata, using sample names only\n")
    pd <- tibble(Sample = selected_samples)
  })
  
  samples_to_process <- selected_samples
  if(!is.null(sample_limit) && sample_limit < length(samples_to_process)) {
    samples_to_process <- sample(samples_to_process, sample_limit)
    cat("Randomly selected", sample_limit, "samples for analysis\n")
  }
  
  all_data <- map_dfr(samples_to_process, function(sample_name) {
    cat("Processing sample:", sample_name, "\n")
    
    gh <- gs[[sample_name]]
    
    if(!node %in% gh_get_pop_paths(gh)) {
      cat("  Node not found in sample, skipping\n")
      return(tibble())
    }
    
    tryCatch({
      ff <- gh_pop_get_data(gh, node)
      
      if(is.null(ff)) {
        cat("  No data found, skipping\n")
        return(tibble())
      }
      
      expr_data <- if(inherits(ff, "cytoframe")) exprs(ff) else exprs(ff)
      
      available_channels <- intersect(selected_markers$channels, colnames(expr_data))
      
      if(length(available_channels) == 0) {
        cat("  No selected channels found, skipping\n")
        return(tibble())
      }
      
      cell_data <- as_tibble(expr_data[, available_channels, drop = FALSE])
      
      if(nrow(cell_data) > max_cells_per_sample) {
        cell_data <- slice_sample(cell_data, n = max_cells_per_sample)
        cat("  Subsampled to", max_cells_per_sample, "cells\n")
      }
      
      cell_data <- cell_data %>%
        mutate(
          Sample = sample_name,
          CellID = paste0(sample_name, "_", row_number()),
          .before = 1
        )
      
      cat("  Extracted", nrow(cell_data), "cells\n")
      return(cell_data)
      
    }, error = function(e) {
      cat("  Error processing sample:", e$message, "\n")
      return(tibble())
    })
  })
  
  if(nrow(all_data) == 0) {
    stop("No single-cell data extracted. Check node selection and markers.")
  }
  
  marker_lookup <- selected_markers$lookup %>%
    select(colname, marker) %>%
    deframe()
  
  for(channel in names(marker_lookup)) {
    if(channel %in% names(all_data)) {
      names(all_data)[names(all_data) == channel] <- marker_lookup[channel]
    }
  }
  
  all_data <- all_data %>%
    left_join(pd, by = "Sample")
  
  cat("Final dataset:", nrow(all_data), "cells from", 
      n_distinct(all_data$Sample), "samples\n")
  cat("Markers for UMAP:", paste(marker_lookup, collapse = ", "), "\n")
  
  return(list(
    data = all_data,
    marker_names = unname(marker_lookup),
    sample_metadata = pd
  ))
}

# ============================================================================
# ENHANCED UMAP COMPUTATION (ROBUST ERROR HANDLING)
# ============================================================================

compute_umap_embedding <- function(single_cell_data, 
                                   n_neighbors = 15, 
                                   min_dist = 0.1, 
                                   n_components = 2,
                                   transform_data = "asinh") {
  
  cat("Computing UMAP embedding...\n")
  
  metadata_cols <- c("Sample", "CellID")
  if(!is.null(single_cell_data$sample_metadata)) {
    metadata_cols <- c(metadata_cols, names(single_cell_data$sample_metadata))
  }
  
  all_cols <- names(single_cell_data$data)
  numeric_cols <- all_cols[map_lgl(all_cols, ~is.numeric(single_cell_data$data[[.x]]))]
  marker_cols <- setdiff(numeric_cols, metadata_cols)
  
  cat("Using", length(marker_cols), "markers:", paste(marker_cols, collapse = ", "), "\n")
  
  umap_matrix <- single_cell_data$data %>%
    select(all_of(marker_cols)) %>%
    as.matrix()
  
  cat("Data range before transformation: [", min(umap_matrix, na.rm = TRUE), ", ", 
      max(umap_matrix, na.rm = TRUE), "]\n")
  
  if(any(is.na(umap_matrix))) {
    n_na <- sum(is.na(umap_matrix))
    cat("Warning:", n_na, "NA values found in data. These will be handled.\n")
  }
  
  if(transform_data == "asinh") {
    umap_matrix <- asinh(umap_matrix / 5)
    cat("Applied asinh transformation (cofactor = 5)\n")
    
  } else if(transform_data == "log10") {
    min_val <- min(umap_matrix, na.rm = TRUE)
    
    if(min_val <= 0) {
      shift_amount <- abs(min_val) + 1
      cat("Shifting data by", shift_amount, "to handle negative/zero values\n")
      umap_matrix <- umap_matrix + shift_amount
    } else {
      umap_matrix <- umap_matrix + 1
    }
    
    umap_matrix <- log10(umap_matrix)
    cat("Applied log10 transformation\n")
    
  } else if(transform_data == "sqrt") {
    if(any(umap_matrix < 0, na.rm = TRUE)) {
      min_val <- min(umap_matrix, na.rm = TRUE)
      shift_amount <- abs(min_val)
      cat("Shifting data by", shift_amount, "to handle negative values for sqrt\n")
      umap_matrix <- umap_matrix + shift_amount
    }
    
    umap_matrix <- sqrt(umap_matrix)
    cat("Applied square root transformation\n")
    
  } else {
    cat("No transformation applied\n")
  }
  
  # Check for NaN/Inf values after transformation
  if(any(is.nan(umap_matrix)) || any(is.infinite(umap_matrix))) {
    n_nan <- sum(is.nan(umap_matrix))
    n_inf <- sum(is.infinite(umap_matrix))
    
    if(n_nan > 0) cat("Error: Found", n_nan, "NaN values after transformation\n")
    if(n_inf > 0) cat("Error: Found", n_inf, "infinite values after transformation\n")
    
    # Replace NaN and Inf with median values
    for(i in 1:ncol(umap_matrix)) {
      col_data <- umap_matrix[, i]
      if(any(is.nan(col_data) | is.infinite(col_data))) {
        median_val <- median(col_data[is.finite(col_data)], na.rm = TRUE)
        if(is.na(median_val)) median_val <- 0
        
        umap_matrix[is.nan(col_data) | is.infinite(col_data), i] <- median_val
        cat("Replaced problematic values in", marker_cols[i], "with median:", median_val, "\n")
      }
    }
  }
  
  # Final check for remaining NA values
  if(any(is.na(umap_matrix))) {
    for(i in 1:ncol(umap_matrix)) {
      col_data <- umap_matrix[, i]
      if(any(is.na(col_data))) {
        median_val <- median(col_data, na.rm = TRUE)
        if(is.na(median_val)) median_val <- 0
        
        umap_matrix[is.na(col_data), i] <- median_val
        cat("Replaced NA values in", marker_cols[i], "with median:", median_val, "\n")
      }
    }
  }
  
  cat("Data range after transformation: [", min(umap_matrix), ", ", max(umap_matrix), "]\n")
  
  # UMAP configuration
  umap_config <- umap.defaults
  umap_config$n_neighbors <- n_neighbors
  umap_config$min_dist <- min_dist
  umap_config$n_components <- n_components
  umap_config$random_state <- 42
  
  # Compute UMAP
  cat("Running UMAP with", n_neighbors, "neighbors and", min_dist, "min_dist...\n")
  
  tryCatch({
    umap_result <- umap(umap_matrix, config = umap_config)
  }, error = function(e) {
    cat("UMAP failed with error:", e$message, "\n")
    cat("This might be due to data scaling issues. Trying with scaled data...\n")
    
    # Try with scaled data as fallback
    umap_matrix_scaled <- scale(umap_matrix)
    
    # Replace any NaN values from scaling
    if(any(is.nan(umap_matrix_scaled))) {
      umap_matrix_scaled[is.nan(umap_matrix_scaled)] <- 0
    }
    
    umap_result <- umap(umap_matrix_scaled, config = umap_config)
    cat("UMAP successful with scaled data\n")
    return(umap_result)
  }) -> umap_result
  
  # Add UMAP coordinates to original data
  umap_data <- single_cell_data$data %>%
    mutate(
      UMAP1 = umap_result$layout[,1],
      UMAP2 = umap_result$layout[,2]
    )
  
  if(n_components > 2) {
    for(i in 3:n_components) {
      umap_data[[paste0("UMAP", i)]] <- umap_result$layout[,i]
    }
  }
  
  cat("UMAP computation complete!\n")
  
  return(list(
    data = umap_data,
    embedding = umap_result,
    marker_names = marker_cols,
    config = umap_config,
    transformation = transform_data
  ))
}

# ============================================================================
# GLOBAL SESSION ENVIRONMENT
# ============================================================================

.umap_session_plots <- new.env()

# ============================================================================
# SESSION INITIALIZATION
# ============================================================================

initialize_visualization_session <- function(umap_result) {
  cat("Setting up visualization session...\n")
  
  .umap_session_plots$umap_data <- umap_result$data
  .umap_session_plots$plots <- list()
  .umap_session_plots$plot_counter <- 0
  
  marker_cols <- umap_result$marker_names
  coord_cols <- c("UMAP1", "UMAP2", "CellID")
  
  available_metadata <- setdiff(names(umap_result$data), c(marker_cols, coord_cols))
  .umap_session_plots$metadata_cols <- available_metadata
  .umap_session_plots$marker_names <- marker_cols
  
  cat("Session initialized with:\n")
  cat("- Cells:", scales::comma(nrow(umap_result$data)), "\n")
  cat("- Markers:", length(marker_cols), "\n")
  cat("- Metadata columns:", length(available_metadata), "\n")
  cat("- Available metadata:", paste(available_metadata, collapse = ", "), "\n\n")
  
  invisible(TRUE)
}

# ============================================================================
# ENHANCED INTERACTIVE VISUALIZATION WITH CONSISTENT CHOICE MENUS
# ============================================================================

# Helper function for consistent menu selection
select_from_menu <- function(options, title, allow_multiple = FALSE, show_preview = TRUE) {
  if (length(options) == 0) {
    cat("No options available\n")
    return(NULL)
  }
  
  cat("\n", title, "\n")
  cat(strrep("=", nchar(title)), "\n")
  
  # Show preview if requested and many options
  if (show_preview && length(options) > 20) {
    cat("Showing first 20 options (", length(options), "total):\n")
    display_options <- head(options, 20)
    cat("...\n")
  } else {
    display_options <- options
  }
  
  iwalk(display_options, ~cat(sprintf("%d. %s\n", .y, .x)))
  
  if (length(options) > 20) {
    cat(sprintf("... and %d more options\n", length(options) - 20))
    cat("Enter 'all' to see all options\n")
  }
  
  while (TRUE) {
    prompt <- if (allow_multiple) {
      "Enter choice(s) (numbers/ranges like 1,3,5-8 or 'all'): "
    } else {
      "Enter choice number: "
    }
    
    choice <- readline(prompt)
    
    if (choice == "") {
      cat("No selection made\n")
      return(NULL)
    }
    
    if (tolower(choice) == "all") {
      if (length(options) > 20) {
        cat("\nAll available options:\n")
        iwalk(options, ~cat(sprintf("%d. %s\n", .y, .x)))
        next
      } else {
        return(if (allow_multiple) options else options[1])
      }
    }
    
    # Parse selection
    selected_indices <- parse_selection(choice, length(options))
    
    if (length(selected_indices) == 0) {
      cat("Invalid selection. Please try again.\n")
      next
    }
    
    if (!allow_multiple && length(selected_indices) > 1) {
      cat("Multiple selections not allowed. Using first selection.\n")
      selected_indices <- selected_indices[1]
    }
    
    return(options[selected_indices])
  }
}

# Parse selection input (handles ranges, comma-separated values)
parse_selection <- function(input, max_value) {
  tryCatch({
    parts <- str_trim(str_split(input, ",")[[1]])
    indices <- numeric(0)
    
    for (part in parts) {
      if (str_detect(part, "-")) {
        range_parts <- as.numeric(str_split(part, "-")[[1]])
        if (length(range_parts) == 2 && all(!is.na(range_parts))) {
          start_idx <- range_parts[1]
          end_idx <- range_parts[2]
          if (start_idx >= 1 && end_idx <= max_value && start_idx <= end_idx) {
            indices <- c(indices, start_idx:end_idx)
          }
        }
      } else if (str_detect(part, "^\\d+$")) {
        idx <- as.numeric(part)
        if (!is.na(idx) && idx >= 1 && idx <= max_value) {
          indices <- c(indices, idx)
        }
      }
    }
    
    return(unique(indices))
  }, error = function(e) {
    return(numeric(0))
  })
}

# Enhanced variable selection with consistent menus
select_visualization_variable <- function(marker_names, metadata_cols, 
                                          purpose = "coloring", allow_multiple = FALSE) {
  
  cat("\n=== Variable Selection for", str_to_title(purpose), "===\n")
  
  all_options <- c()
  option_labels <- c()
  
  if (length(metadata_cols) > 0) {
    metadata_labels <- paste0(metadata_cols, " (metadata)")
    all_options <- c(all_options, metadata_cols)
    option_labels <- c(option_labels, metadata_labels)
  }
  
  if (length(marker_names) > 0) {
    marker_labels <- paste0(marker_names, " (marker)")
    all_options <- c(all_options, marker_names)
    option_labels <- c(option_labels, marker_labels)
  }
  
  if (length(all_options) == 0) {
    cat("No variables available for", purpose, "\n")
    return(NULL)
  }
  
  # Show categorized menu
  cat("Available variables:\n")
  
  counter <- 1
  if (length(metadata_cols) > 0) {
    cat("\nMETADATA VARIABLES:\n")
    for (i in seq_along(metadata_cols)) {
      cat(sprintf("%d. %s\n", counter, metadata_cols[i]))
      counter <- counter + 1
    }
  }
  
  if (length(marker_names) > 0) {
    cat("\nMARKER VARIABLES:\n")
    for (i in seq_along(marker_names)) {
      cat(sprintf("%d. %s\n", counter, marker_names[i]))
      counter <- counter + 1
    }
  }
  
  selected <- select_from_menu(
    all_options, 
    paste("Select variable(s) for", purpose),
    allow_multiple = allow_multiple,
    show_preview = FALSE
  )
  
  return(selected)
}

# Create and save UMAP plot function
create_and_save_umap_plot <- function(color_by = NULL, facet_by = NULL, 
                                      title = NULL, filename = NULL, 
                                      plot_type = "scatter") {
  
  umap_data <- .umap_session_plots$umap_data
  
  if (is.null(title)) {
    title_parts <- c("UMAP")
    if (!is.null(color_by)) title_parts <- c(title_parts, paste("colored by", color_by))
    if (!is.null(facet_by)) title_parts <- c(title_parts, paste("faceted by", facet_by))
    title <- paste(title_parts, collapse = " ")
  }
  
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2))
  
  if (plot_type == "density") {
    p <- p + 
      stat_density_2d_filled(alpha = 0.7) +
      scale_fill_viridis_d(option = "plasma") +
      theme_void() +
      theme(legend.position = "none")
    
  } else {
    # Regular scatter plot
    if (is.null(color_by)) {
      p <- p + geom_point(color = "steelblue", alpha = 0.6, size = 0.3)
    } else {
      if (color_by %in% .umap_session_plots$marker_names) {
        # Marker expression - continuous scale
        p <- p + 
          geom_point(aes_string(color = color_by), alpha = 0.6, size = 0.3) +
          scale_color_viridis_c(option = "plasma") +
          guides(color = guide_colorbar())
      } else {
        # Metadata - discrete scale
        p <- p + 
          geom_point(aes_string(color = color_by), alpha = 0.6, size = 0.3) +
          guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
      }
    }
  }
  
  # Add faceting if specified
  if (!is.null(facet_by)) {
    # For marker faceting, create discrete bins
    if (facet_by %in% .umap_session_plots$marker_names) {
      marker_values <- umap_data[[facet_by]]
      breaks <- quantile(marker_values, probs = c(0, 0.33, 0.67, 1), na.rm = TRUE)
      umap_data[[paste0(facet_by, "_level")]] <- cut(marker_values, 
                                                     breaks = breaks,
                                                     labels = c("Low", "Medium", "High"),
                                                     include.lowest = TRUE)
      p <- p + facet_wrap(as.formula(paste("~", paste0(facet_by, "_level"))), scales = "free")
    } else {
      p <- p + facet_wrap(as.formula(paste("~", facet_by)), scales = "free")
    }
  }
  
  # Apply theme
  p <- p +
    labs(title = title) +
    theme_minimal() +
    theme(
      legend.position = "right",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  if (!is.null(facet_by)) {
    p <- p + theme(
      strip.background = element_rect(fill = "lightgray", color = "black"),
      strip.text = element_text(face = "bold")
    )
  }
  
  print(p)
  
  # Store plot
  .umap_session_plots$plot_counter <<- .umap_session_plots$plot_counter + 1
  plot_id <- paste0("plot_", .umap_session_plots$plot_counter)
  
  plot_info <- list(
    plot = p,
    color_by = color_by,
    facet_by = facet_by,
    plot_type = plot_type,
    title = title,
    timestamp = Sys.time()
  )
  
  .umap_session_plots$plots[[plot_id]] <<- plot_info
  
  # Save if filename provided
  if (!is.null(filename) && filename != "") {
    if (!str_detect(filename, "\\.(png|pdf|jpg|jpeg|tiff|svg)$")) {
      filename <- paste0(filename, ".png")
    }
    
    tryCatch({
      ggsave(filename, plot = p, width = 10, height = 8, dpi = 300)
      cat("Plot saved as:", filename, "\n")
    }, error = function(e) {
      cat("Error saving plot:", e$message, "\n")
    })
  }
  
  cat("Plot created and stored as:", plot_id, "\n")
  return(plot_id)
}

# Enhanced multi-marker heatmap with consistent selection
create_multi_marker_heatmap <- function() {
  umap_data <- .umap_session_plots$umap_data
  marker_names <- .umap_session_plots$marker_names
  
  if (length(marker_names) == 0) {
    cat("No markers available for heatmap\n")
    return(NULL)
  }
  
  cat("\n=== Multi-Marker Heatmap Creation ===\n")
  cat("Available markers:", length(marker_names), "\n")
  
  # Menu for marker selection strategy
  selection_options <- c(
    "Select specific markers",
    "Use first 9 markers (3x3 grid)",
    "Use all markers",
    "Show marker list and return"
  )
  
  strategy <- select_from_menu(
    selection_options,
    "Choose marker selection strategy"
  )
  
  if (is.null(strategy)) return(NULL)
  
  selected_markers <- switch(
    strategy,
    "Select specific markers" = {
      select_from_menu(
        marker_names,
        "Select markers for heatmap (max 16 recommended)",
        allow_multiple = TRUE
      )
    },
    "Use first 9 markers (3x3 grid)" = {
      n_markers <- min(9, length(marker_names))
      head(marker_names, n_markers)
    },
    "Use all markers" = {
      if (length(marker_names) > 16) {
        cat("Warning: Using all", length(marker_names), "markers may create a crowded plot\n")
        confirm_options <- c("Yes, continue", "No, go back")
        confirm <- select_from_menu(confirm_options, "Continue with all markers?")
        if (is.null(confirm) || confirm == "No, go back") return(NULL)
      }
      marker_names
    },
    "Show marker list and return" = {
      cat("\nAll available markers:\n")
      iwalk(marker_names, ~cat(sprintf("%d. %s\n", .y, .x)))
      readline("Press Enter to continue...")
      return(NULL)
    }
  )
  
  if (is.null(selected_markers) || length(selected_markers) == 0) {
    cat("No markers selected\n")
    return(NULL)
  }
  
  cat("Creating heatmap with", length(selected_markers), "markers\n")
  
  # Create individual plots for each marker
  marker_plots <- map(selected_markers, function(marker) {
    ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[marker]])) +
      geom_point(size = 0.1, alpha = 0.5) +
      scale_color_viridis_c(name = "", option = "plasma") +
      labs(title = marker) +
      theme_void() +
      theme(
        plot.title = element_text(size = 10, hjust = 0.5),
        legend.key.size = unit(0.3, "cm"),
        legend.text = element_text(size = 6)
      )
  })
  
  # Determine optimal grid layout
  n_markers <- length(selected_markers)
  ncol_grid <- case_when(
    n_markers <= 4 ~ 2,
    n_markers <= 9 ~ 3,
    n_markers <= 16 ~ 4,
    TRUE ~ ceiling(sqrt(n_markers))
  )
  
  # Arrange in grid
  combined_plot <- wrap_plots(marker_plots, ncol = ncol_grid)
  
  # Store and display
  .umap_session_plots$plot_counter <<- .umap_session_plots$plot_counter + 1
  plot_id <- paste0("plot_", .umap_session_plots$plot_counter)
  
  plot_title <- paste("Multi-marker heatmap (", length(selected_markers), "markers)")
  
  plot_info <- list(
    plot = combined_plot,
    color_by = "multiple_markers",
    facet_by = NULL,
    plot_type = "heatmap",
    title = plot_title,
    markers_used = selected_markers,
    timestamp = Sys.time()
  )
  
  .umap_session_plots$plots[[plot_id]] <<- plot_info
  print(combined_plot)
  
  # Save option
  save_options <- c("Yes, save plot", "No, continue without saving")
  save_choice <- select_from_menu(save_options, "Save this heatmap?")
  
  if (!is.null(save_choice) && save_choice == "Yes, save plot") {
    filename <- readline("Enter filename: ")
    if (filename != "") {
      if (!str_detect(filename, "\\.(png|pdf|jpg|jpeg|tiff|svg)$")) {
        filename <- paste0(filename, ".png")
      }
      
      # Adjust size based on number of markers
      width <- max(8, ceiling(ncol_grid * 2.5))
      height <- max(6, ceiling((length(marker_plots) / ncol_grid) * 2))
      
      tryCatch({
        ggsave(filename, plot = combined_plot, width = width, height = height, dpi = 300)
        cat("Heatmap saved as:", filename, "\n")
      }, error = function(e) {
        cat("Error saving plot:", e$message, "\n")
      })
    }
  }
  
  cat("Plot created and stored as:", plot_id, "\n")
  return(plot_id)
}

# Enhanced density plot creation
create_density_plot <- function() {
  umap_data <- .umap_session_plots$umap_data
  metadata_cols <- .umap_session_plots$metadata_cols
  marker_names <- .umap_session_plots$marker_names
  all_options <- c(metadata_cols, marker_names)
  
  cat("\n=== Density Plot Creation ===\n")
  
  density_options <- c(
    "Overall density (no grouping)",
    "Density faceted by variable"
  )
  
  density_type <- select_from_menu(density_options, "Choose density plot type")
  
  if (is.null(density_type)) return(NULL)
  
  facet_var <- NULL
  if (density_type == "Density faceted by variable") {
    # Filter to categorical variables
    categorical_vars <- all_options[map_lgl(all_options, ~!is.numeric(umap_data[[.x]]))]
    
    if (length(categorical_vars) == 0) {
      cat("No categorical variables available for faceting\n")
      return(NULL)
    }
    
    facet_var <- select_from_menu(
      categorical_vars,
      "Select variable for faceting"
    )
    
    if (is.null(facet_var)) return(NULL)
  }
  
  # Save option
  save_options <- c("Yes, save plot", "No, continue without saving")
  save_choice <- select_from_menu(save_options, "Save this plot?")
  
  filename <- NULL
  if (!is.null(save_choice) && save_choice == "Yes, save plot") {
    filename <- readline("Enter filename: ")
    if (filename == "") filename <- NULL
  }
  
  create_and_save_umap_plot(
    facet_by = facet_var,
    plot_type = "density",
    filename = filename
  )
}

# Enhanced plot management
manage_stored_plots <- function() {
  if (length(.umap_session_plots$plots) == 0) {
    cat("No plots stored in current session\n")
    return(NULL)
  }
  
  while (TRUE) {
    cat("\n=== Stored Plot Management ===\n")
    cat("Total plots:", length(.umap_session_plots$plots), "\n\n")
    
    # Show plot list
    plot_info <- map_chr(names(.umap_session_plots$plots), function(plot_id) {
      plot_data <- .umap_session_plots$plots[[plot_id]]
      info_parts <- c(plot_data$title)
      if (!is.null(plot_data$color_by)) {
        info_parts <- c(info_parts, paste("color:", plot_data$color_by))
      }
      if (!is.null(plot_data$facet_by)) {
        info_parts <- c(info_parts, paste("facet:", plot_data$facet_by))
      }
      info_parts <- c(info_parts, format(plot_data$timestamp, "%H:%M:%S"))
      paste(info_parts, collapse = " | ")
    })
    
    management_options <- c(
      "Display a plot",
      "Save a plot to file", 
      "Export all plots",
      "Delete a plot",
      "Clear all plots",
      "Return to main menu"
    )
    
    action <- select_from_menu(management_options, "Choose action")
    
    if (is.null(action) || action == "Return to main menu") {
      break
    }
    
    if (action %in% c("Display a plot", "Save a plot to file", "Delete a plot")) {
      selected_plot_id <- select_from_menu(
        names(.umap_session_plots$plots),
        paste("Select plot to", tolower(str_remove(action, " a plot| to file")))
      )
      
      if (is.null(selected_plot_id)) next
      
      if (action == "Display a plot") {
        print(.umap_session_plots$plots[[selected_plot_id]]$plot)
        
      } else if (action == "Save a plot to file") {
        filename <- readline("Enter filename: ")
        if (filename != "") {
          if (!str_detect(filename, "\\.(png|pdf|jpg|jpeg|tiff|svg)$")) {
            filename <- paste0(filename, ".png")
          }
          tryCatch({
            ggsave(filename, plot = .umap_session_plots$plots[[selected_plot_id]]$plot)
            cat("Plot saved as:", filename, "\n")
          }, error = function(e) {
            cat("Error saving plot:", e$message, "\n")
          })
        }
        
      } else if (action == "Delete a plot") {
        confirm_options <- c("Yes, delete", "No, cancel")
        confirm <- select_from_menu(confirm_options, "Confirm deletion")
        if (!is.null(confirm) && confirm == "Yes, delete") {
          .umap_session_plots$plots[[selected_plot_id]] <<- NULL
          cat("Plot", selected_plot_id, "deleted\n")
        }
      }
      
    } else if (action == "Export all plots") {
      directory <- readline("Export directory (default: umap_plots): ")
      if (directory == "") directory <- "umap_plots"
      export_all_session_plots(directory)
      
    } else if (action == "Clear all plots") {
      confirm_options <- c("Yes, clear all", "No, cancel")
      confirm <- select_from_menu(confirm_options, "Clear all plots?")
      if (!is.null(confirm) && confirm == "Yes, clear all") {
        .umap_session_plots$plots <<- list()
        .umap_session_plots$plot_counter <<- 0
        cat("All plots cleared\n")
      }
    }
  }
}

# Main enhanced interactive visualization menu
interactive_visualization_menu <- function() {
  if (length(.umap_session_plots) == 0) {
    cat("No active visualization session. Please run analyze_flow_umap_enhanced() first.\n")
    return(NULL)
  }
  
  umap_data <- .umap_session_plots$umap_data
  marker_names <- .umap_session_plots$marker_names
  metadata_cols <- .umap_session_plots$metadata_cols
  
  while (TRUE) {
    cat("\n=== Interactive UMAP Visualization Menu ===\n")
    cat("Current session:", scales::comma(nrow(umap_data)), "cells\n")
    cat("Stored plots:", length(.umap_session_plots$plots), "\n")
    cat("Available variables:", length(c(metadata_cols, marker_names)), 
        "(", length(metadata_cols), "metadata,", length(marker_names), "markers)\n\n")
    
    main_options <- c(
      "Create basic UMAP plot",
      "Color by variable", 
      "Facet by variable",
      "Color AND facet",
      "Create density plot",
      "Multi-marker heatmap overlay",
      "Manage stored plots",
      "Session summary",
      "Exit visualization"
    )
    
    choice <- select_from_menu(main_options, "Choose visualization option")
    
    if (is.null(choice) || choice == "Exit visualization") {
      cat("Exiting visualization menu\n")
      break
    }
    
    # Execute selected action
    if (choice == "Create basic UMAP plot") {
      save_options <- c("Yes, save plot", "No, continue without saving")
      save_choice <- select_from_menu(save_options, "Save this plot?")
      
      filename <- NULL
      if (!is.null(save_choice) && save_choice == "Yes, save plot") {
        filename <- readline("Enter filename: ")
        if (filename == "") filename <- NULL
      }
      
      create_and_save_umap_plot(title = "Basic UMAP", filename = filename)
      
    } else if (choice == "Color by variable") {
      color_var <- select_visualization_variable(
        marker_names, metadata_cols, "coloring"
      )
      
      if (!is.null(color_var)) {
        var_type <- if (color_var %in% marker_names) "marker" else "metadata"
        plot_title <- paste("UMAP colored by", color_var, "(", var_type, ")")
        
        save_options <- c("Yes, save plot", "No, continue without saving")
        save_choice <- select_from_menu(save_options, "Save this plot?")
        
        filename <- NULL
        if (!is.null(save_choice) && save_choice == "Yes, save plot") {
          filename <- readline("Enter filename: ")
          if (filename == "") filename <- NULL
        }
        
        create_and_save_umap_plot(
          color_by = color_var,
          title = plot_title,
          filename = filename
        )
      }
      
    } else if (choice == "Facet by variable") {
      all_options <- c(metadata_cols, marker_names)
      # Include both categorical variables AND markers (markers will be binned automatically)
      categorical_vars <- all_options[map_lgl(all_options, ~!is.numeric(umap_data[[.x]]))]
      
      # Combine categorical metadata with all markers
      facet_options <- c(categorical_vars, marker_names)
      
      if (length(facet_options) == 0) {
        cat("No variables available for faceting\n")
        next
      }
      
      facet_var <- select_from_menu(
        facet_options,
        "Select variable for faceting (markers will be auto-binned)"
      )
      
      if (!is.null(facet_var)) {
        save_options <- c("Yes, save plot", "No, continue without saving")
        save_choice <- select_from_menu(save_options, "Save this plot?")
        
        filename <- NULL
        if (!is.null(save_choice) && save_choice == "Yes, save plot") {
          filename <- readline("Enter filename: ")
          if (filename == "") filename <- NULL
        }
        
        create_and_save_umap_plot(facet_by = facet_var, filename = filename)
      }
      
    } else if (choice == "Color AND facet") {
      all_options <- c(metadata_cols, marker_names)
      
      # Step 1: Color variable
      color_var <- select_visualization_variable(
        marker_names, metadata_cols, "coloring"
      )
      
      if (is.null(color_var)) next
      
      # Step 2: Facet variable (exclude color variable and include both categorical and markers)
      categorical_vars <- all_options[map_lgl(all_options, ~!is.numeric(umap_data[[.x]]))]
      remaining_markers <- setdiff(marker_names, color_var)
      
      # Combine remaining categorical metadata with remaining markers
      available_facet_vars <- c(setdiff(categorical_vars, color_var), remaining_markers)
      
      if (length(available_facet_vars) == 0) {
        cat("No additional variables available for faceting\n")
        next
      }
      
      facet_var <- select_from_menu(
        available_facet_vars,
        "Select variable for faceting (markers will be auto-binned)"
      )
      
      if (!is.null(facet_var)) {
        save_options <- c("Yes, save plot", "No, continue without saving")
        save_choice <- select_from_menu(save_options, "Save this plot?")
        
        filename <- NULL
        if (!is.null(save_choice) && save_choice == "Yes, save plot") {
          filename <- readline("Enter filename: ")
          if (filename == "") filename <- NULL
        }
        
        create_and_save_umap_plot(
          color_by = color_var,
          facet_by = facet_var,
          filename = filename
        )
      }
      
    } else if (choice == "Create density plot") {
      create_density_plot()
      
    } else if (choice == "Multi-marker heatmap overlay") {
      create_multi_marker_heatmap()
      
    } else if (choice == "Manage stored plots") {
      manage_stored_plots()
      
    } else if (choice == "Session summary") {
      cat("\n=== Session Summary ===\n")
      cat("Cells:", scales::comma(nrow(umap_data)), "\n")
      cat("Markers:", length(marker_names), "\n")
      cat("Metadata columns:", length(metadata_cols), "\n")
      cat("Stored plots:", length(.umap_session_plots$plots), "\n\n")
      
      if (length(marker_names) > 0) {
        cat("Markers:", paste(marker_names, collapse = ", "), "\n")
      }
      if (length(metadata_cols) > 0) {
        cat("Metadata:", paste(metadata_cols, collapse = ", "), "\n")
      }
      
      if (length(.umap_session_plots$plots) > 0) {
        cat("\nPlot history:\n")
        iwalk(.umap_session_plots$plots, function(plot_info, plot_id) {
          cat(sprintf("- %s: %s\n", plot_id, plot_info$title))
        })
      }
      
      readline("Press Enter to continue...")
    }
  }
  
  return(invisible(.umap_session_plots))
}

# ============================================================================
# ROBUST EXTERNAL METADATA IMPORT (OPTIONAL)
# ============================================================================

# ============================================================================
# FIXED EXTERNAL METADATA IMPORT WITH PROPER DUPLICATE HANDLING
# ============================================================================

import_external_metadata <- function(gs) {
  cat("\n=== External Metadata Import ===\n")
  cat("This function will help you import and merge external metadata\n")
  cat("with your flow cytometry samples.\n\n")
  
  current_samples <- sampleNames(gs)
  cat("Current samples in gating set:", length(current_samples), "\n")
  cat("Sample preview:", paste(head(current_samples, 3), collapse = ", "), "...\n\n")
  
  while(TRUE) {
    cat("Metadata Import Options:\n")
    cat("1. Load from CSV file\n")
    cat("2. Load from R object in environment\n")
    cat("3. Skip metadata import\n")
    
    choice <- readline("Choose option (1-3): ")
    
    metadata <- NULL
    
    if(choice == "1") {
      file_path <- readline("Enter CSV file path: ")
      if(file.exists(file_path)) {
        tryCatch({
          metadata <- read_csv(file_path, show_col_types = FALSE)
          metadata <- as_tibble(metadata)
          cat("Successfully loaded CSV with", nrow(metadata), "rows and", ncol(metadata), "columns\n")
        }, error = function(e) {
          cat("Error loading CSV:", e$message, "\n")
          next
        })
      } else {
        cat("File not found:", file_path, "\n")
        next
      }
      
    } else if(choice == "2") {
      cat("\n=== R Object Metadata Import ===\n")
      cat("Available objects in environment:\n")
      objects_list <- ls(envir = .GlobalEnv)
      data_objects <- objects_list[sapply(objects_list, function(x) {
        obj <- get(x, envir = .GlobalEnv)
        is.data.frame(obj) || is_tibble(obj)
      })]
      
      if(length(data_objects) == 0) {
        cat("No data frame objects found in global environment\n")
        next
      }
      
      cat("Data frame objects:\n")
      iwalk(data_objects, ~cat(sprintf("%d. %s\n", .y, .x)))
      
      obj_choice <- readline("Enter object name or number: ")
      
      selected_obj <- NULL
      if(grepl("^\\d+$", obj_choice)) {
        obj_num <- as.numeric(obj_choice)
        if(!is.na(obj_num) && obj_num >= 1 && obj_num <= length(data_objects)) {
          selected_obj <- data_objects[obj_num]
        }
      } else if(obj_choice %in% data_objects) {
        selected_obj <- obj_choice
      }
      
      if(is.null(selected_obj)) {
        cat("Invalid selection\n")
        next
      }
      
      metadata_obj <- get(selected_obj, envir = .GlobalEnv)
      metadata <- as_tibble(metadata_obj)
      cat("Using object '", selected_obj, "' with", nrow(metadata), "rows and", ncol(metadata), "columns\n")
      
    } else if(choice == "3") {
      cat("Skipping metadata import\n")
      return(NULL)
    } else {
      cat("Invalid choice\n")
      next
    }
    
    if(!is.null(metadata)) {
      # Ensure metadata is a proper data frame/tibble
      if(!is.data.frame(metadata)) {
        tryCatch({
          metadata <- as_tibble(metadata)
        }, error = function(e) {
          cat("Error converting to tibble:", e$message, "\n")
          next
        })
      }
      
      cat("\nMetadata preview:\n")
      print(head(metadata, 5))
      cat("Column names:", paste(names(metadata), collapse = ", "), "\n\n")
      
      # Sample ID Column Selection
      cat("Sample ID Column Selection:\n")
      cat("Available columns:\n")
      iwalk(names(metadata), ~cat(sprintf("%d. %s\n", .y, .x)))
      
      sample_col_choice <- readline("Enter column number/name for Sample ID: ")
      
      sample_col <- NULL
      if(grepl("^\\d+$", sample_col_choice)) {
        col_num <- as.numeric(sample_col_choice)
        if(!is.na(col_num) && col_num >= 1 && col_num <= ncol(metadata)) {
          sample_col <- names(metadata)[col_num]
        }
      } else if(sample_col_choice %in% names(metadata)) {
        sample_col <- sample_col_choice
      }
      
      if(is.null(sample_col)) {
        cat("Invalid selection\n")
        next
      }
      
      # FIXED: Removed duplicate checking logic
      # In flow cytometry, having the same sample metadata repeated is expected and normal
      metadata_samples <- metadata[[sample_col]]
      
      # Simple validation: check for data quality issues
      if(length(metadata_samples) == 0) {
        cat("Error: No sample names found in selected column\n")
        next
      }
      
      # Check for actual data issues (like empty/NA sample names)
      na_samples <- sum(is.na(metadata_samples))
      empty_samples <- sum(metadata_samples == "", na.rm = TRUE)
      
      if(na_samples > 0) {
        cat("Warning:", na_samples, "NA sample names found\n")
      }
      
      if(empty_samples > 0) {
        cat("Warning:", empty_samples, "empty sample names found\n")
      }
      
      # If there are repeated sample names with identical metadata, deduplicate to sample level
      unique_samples <- unique(metadata_samples)
      if(length(unique_samples) < length(metadata_samples)) {
        cat("Note: Found repeated sample names in metadata (", length(metadata_samples), 
            "rows ->", length(unique_samples), "unique samples)\n")
        cat("This is normal for single-cell data. Deduplicating to sample-level metadata.\n")
        
        # Convert to standard data frame to avoid Bioconductor S4 object issues
        tryCatch({
          # Force conversion to standard data types
          metadata <- as.data.frame(metadata)
          metadata <- as_tibble(metadata)
          
          # Keep first occurrence of each sample using base R approach to avoid slice() issues
          sample_indices <- match(unique_samples, metadata[[sample_col]])
          metadata <- metadata[sample_indices, ]
          
        }, error = function(e) {
          cat("Error during deduplication:", e$message, "\n")
          cat("Attempting alternative deduplication method...\n")
          
          # Alternative approach using base R
          unique_rows <- !duplicated(metadata[[sample_col]])
          metadata <- metadata[unique_rows, ]
        })
        
        metadata_samples <- metadata[[sample_col]]
        cat("Deduplicated to", nrow(metadata), "unique sample entries\n")
      }
      
      # Continue with sample matching
      matched_samples <- intersect(current_samples, metadata_samples)
      unmatched_gs <- setdiff(current_samples, metadata_samples)
      unmatched_meta <- setdiff(metadata_samples, current_samples)
      
      cat("\nSample Matching Results:\n")
      cat("Samples in both datasets:", length(matched_samples), "\n")
      cat("Samples in GatingSet only:", length(unmatched_gs), "\n")
      cat("Samples in metadata only:", length(unmatched_meta), "\n")
      
      if(length(matched_samples) == 0) {
        cat("No matching samples found! Check sample naming.\n")
        cat("GatingSet samples preview:", paste(head(current_samples, 3), collapse = ", "), "\n")
        cat("Metadata samples preview:", paste(head(metadata_samples, 3), collapse = ", "), "\n")
        
        # Show sample comparison for debugging
        cat("\nFirst 10 GatingSet samples:\n")
        iwalk(head(current_samples, 10), ~cat(sprintf("  %d. '%s'\n", .y, .x)))
        cat("\nFirst 10 metadata samples:\n") 
        iwalk(head(metadata_samples, 10), ~cat(sprintf("  %d. '%s'\n", .y, .x)))
        
        next
      }
      
      if(length(unmatched_gs) > 0) {
        cat("Unmatched GatingSet samples:", paste(head(unmatched_gs, 5), collapse = ", "))
        if(length(unmatched_gs) > 5) cat("...")
        cat("\n")
      }
      
      # Column selection
      other_cols <- setdiff(names(metadata), sample_col)
      cat("\nAvailable metadata columns to include:\n")
      iwalk(other_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
      
      col_selection <- readline("Enter column numbers (space-separated) or 'all' or 'none': ")
      
      selected_cols <- NULL
      if(tolower(col_selection) == "all") {
        selected_cols <- other_cols
      } else if(tolower(col_selection) == "none") {
        cat("No additional metadata columns selected\n")
        next
      } else if(col_selection == "") {
        selected_cols <- other_cols
        cat("No input provided, using all columns\n")
      } else {
        tryCatch({
          selections <- str_trim(str_split(col_selection, "\\s+")[[1]])
          selected_cols <- character(0)
          
          for(sel in selections) {
            if(grepl("^\\d+$", sel)) {
              col_num <- as.numeric(sel)
              if(!is.na(col_num) && col_num >= 1 && col_num <= length(other_cols)) {
                selected_cols <- c(selected_cols, other_cols[col_num])
              }
            } else if(sel %in% other_cols) {
              selected_cols <- c(selected_cols, sel)
            }
          }
          
          selected_cols <- unique(selected_cols)
          
          if(length(selected_cols) == 0) {
            cat("No valid columns selected\n")
            next
          }
          
        }, error = function(e) {
          cat("Invalid input format. Using all columns\n")
          selected_cols <- other_cols
        })
      }
      
      cat("\nSelected metadata columns:\n")
      iwalk(selected_cols, ~cat(sprintf("  %d. %s\n", .y, .x)))
      
      confirm <- readline("Confirm metadata column selection? (y/n): ")
      if(tolower(confirm) != "y") {
        cat("Column selection cancelled\n")
        next
      }
      
      # Create final metadata with proper structure
      final_metadata <- metadata %>%
        select(all_of(c(sample_col, selected_cols))) %>%
        rename(Sample = all_of(sample_col)) %>%
        dplyr::filter(Sample %in% current_samples) %>%
        as_tibble()
      
      # Final verification - ensure no duplicates at sample level
      final_sample_check <- table(final_metadata$Sample)
      final_duplicates <- names(final_sample_check)[final_sample_check > 1]
      
      if(length(final_duplicates) > 0) {
        cat("Warning: Sample-level duplicates still present. Taking first occurrence of each.\n")
        tryCatch({
          # Use base R approach to avoid slice() issues with S4 objects
          unique_sample_rows <- !duplicated(final_metadata$Sample)
          final_metadata <- final_metadata[unique_sample_rows, ]
        }, error = function(e) {
          cat("Error removing final duplicates:", e$message, "\n")
          # Force conversion and try again
          final_metadata <- as.data.frame(final_metadata)
          final_metadata <- as_tibble(final_metadata)
          unique_sample_rows <- !duplicated(final_metadata$Sample)
          final_metadata <- final_metadata[unique_sample_rows, ]
        })
      }
      
      cat("Final metadata prepared successfully:\n")
      cat("- Samples:", nrow(final_metadata), "\n")
      cat("- Metadata columns:", ncol(final_metadata)-1, "\n")
      cat("- Column names:", paste(names(final_metadata), collapse = ", "), "\n")
      cat("- Sample-level metadata confirmed\n")
      
      return(final_metadata)
    }
  }
}

# ============================================================================
# FIXED SECTION IN MAIN FUNCTION FOR METADATA MERGING
# ============================================================================

merge_external_metadata_safely <- function(single_cell_data, external_metadata) {
  cat("Merging external metadata with single-cell data...\n")
  cat("External metadata dimensions:", nrow(external_metadata), "rows,", ncol(external_metadata), "columns\n")
  
  # Ensure Sample column exists
  if(!"Sample" %in% names(external_metadata)) {
    stop("No 'Sample' column found in external metadata")
  }
  
  # Check for column conflicts
  current_meta_cols <- setdiff(names(single_cell_data$data), 
                               c("CellID", single_cell_data$marker_names))
  ext_meta_cols <- setdiff(names(external_metadata), "Sample")
  conflicting_cols <- intersect(current_meta_cols, ext_meta_cols)
  
  if(length(conflicting_cols) > 0) {
    cat("Column name conflicts found:", paste(conflicting_cols, collapse = ", "), "\n")
    cat("External metadata columns will be prefixed with 'ext_'\n")
    
    for(col in conflicting_cols) {
      names(external_metadata)[names(external_metadata) == col] <- paste0("ext_", col)
    }
  }
  
  original_nrow <- nrow(single_cell_data$data)
  
  # FIXED: Check for duplicates ONLY in external metadata (sample-level duplicates)
  # This checks if the same sample has different metadata (which would be an error)
  cat("Checking external metadata for sample-level inconsistencies...\n")
  
  # Group by Sample and check if all metadata is identical for each sample
  metadata_consistency_check <- external_metadata %>%
    group_by(Sample) %>%
    summarise(
      n_rows = n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_rows > 1)
  
  if(nrow(metadata_consistency_check) > 0) {
    cat("Found samples with multiple metadata rows in external data:\n")
    print(head(metadata_consistency_check, 10))
    
    cat("\nChecking if metadata is consistent across duplicate sample entries...\n")
    
    # Check if duplicated samples have identical metadata
    problem_samples <- c()
    for(sample_id in metadata_consistency_check$Sample) {
      sample_rows <- external_metadata %>% dplyr::filter(Sample == sample_id)
      
      # Check if all rows are identical (excluding Sample column)
      metadata_cols <- setdiff(names(sample_rows), "Sample")
      if(length(metadata_cols) > 0) {
        unique_metadata_combinations <- sample_rows %>%
          select(all_of(metadata_cols)) %>%
          distinct() %>%
          nrow()
        
        if(unique_metadata_combinations > 1) {
          problem_samples <- c(problem_samples, sample_id)
        }
      }
    }
    
    if(length(problem_samples) > 0) {
      cat("ERROR: Found samples with inconsistent metadata:\n")
      for(sample_id in head(problem_samples, 5)) {
        cat("Sample:", sample_id, "\n")
        problem_rows <- external_metadata %>% dplyr::filter(Sample == sample_id)
        print(problem_rows)
        cat("\n")
      }
      stop("Cannot merge metadata with inconsistent values for the same sample. Please clean your external metadata.")
    } else {
      cat("Duplicate sample entries have consistent metadata. Using first occurrence of each.\n")
      external_metadata <- external_metadata %>%
        group_by(Sample) %>%
        slice(1) %>%
        ungroup()
      cat("Deduplicated external metadata to", nrow(external_metadata), "unique samples\n")
    }
  } else {
    cat("No duplicate samples found in external metadata - each sample has unique metadata\n")
  }
  
  # Check sample overlap between datasets
  sc_samples <- unique(single_cell_data$data$Sample)
  meta_samples <- external_metadata$Sample
  
  cat("Sample matching summary:\n")
  cat("- Single-cell data samples:", length(sc_samples), "\n")
  cat("- External metadata samples:", length(meta_samples), "\n")
  cat("- Overlapping samples:", length(intersect(sc_samples, meta_samples)), "\n")
  cat("- SC samples without metadata:", length(setdiff(sc_samples, meta_samples)), "\n")
  cat("- Metadata samples not in SC data:", length(setdiff(meta_samples, sc_samples)), "\n")
  
  # Perform the join - this is a many-to-one relationship (many cells to one sample metadata)
  tryCatch({
    merged_data <- single_cell_data$data %>%
      left_join(external_metadata, by = "Sample", relationship = "many-to-one")
    
    if(nrow(merged_data) != original_nrow) {
      stop("Row count changed during merge (", original_nrow, "->", nrow(merged_data), 
           "). This suggests a join issue.")
    }
    
    cat("Merge successful: row count maintained at", nrow(merged_data), "\n")
    
    # Check how many cells got metadata
    cells_with_metadata <- merged_data %>%
      select(all_of(setdiff(names(external_metadata), "Sample"))) %>%
      complete.cases() %>%
      sum()
    
    cat("Cells with complete metadata:", scales::comma(cells_with_metadata), 
        "out of", scales::comma(nrow(merged_data)), "\n")
    
    cat("External metadata successfully merged!\n")
    cat("Added columns:", paste(setdiff(names(external_metadata), "Sample"), collapse = ", "), "\n")
    
    # Update the single-cell data object
    single_cell_data$data <- merged_data
    
    return(single_cell_data)
    
  }, error = function(e) {
    cat("Error during join operation:", conditionMessage(e), "\n")
    
    # Provide detailed debugging information
    cat("\nDebugging information:\n")
    cat("Single-cell data samples (unique):", length(unique(single_cell_data$data$Sample)), "\n")
    cat("External metadata samples:", nrow(external_metadata), "\n")
    cat("Matching samples:", length(intersect(unique(single_cell_data$data$Sample), external_metadata$Sample)), "\n")
    
    stop("Join operation failed. See debugging information above.")
  })
}
# ============================================================================
# MAIN ENHANCED FUNCTION WITH SAMPLE EXCLUSION
# ============================================================================

analyze_flow_umap_enhanced <- function(gs, keywords = c("pairing_factor", "tissue_factor")) {
  
  cat("=== Enhanced Interactive UMAP Analysis for Flow Cytometry ===\n")
  cat("Features:\n")
  cat("- Sample exclusion (unstained controls, etc.)\n")
  cat("- Event counting and population assessment\n") 
  cat("- Optional external metadata import\n")
  cat("- Persistent visualization with plot saving\n")
  cat("- Session-based plot management\n\n")
  
  # Step 1: Sample exclusion
  selected_samples <- interactive_sample_exclusion(gs)
  
  # Step 2: Select node
  selected_node <- select_single_node_for_umap(gs)
  cat("\nSelected node:", basename(selected_node), "\n")
  
  # Step 3: Select markers
  selected_markers <- select_markers_for_umap(gs)
  if(is.character(selected_markers) && selected_markers == "BACK") {
    cat("Analysis cancelled\n")
    return(NULL)
  }
  
  # Step 4: Enhanced data extraction with event counting
  cat("\n=== Data Extraction Parameters ===\n")
  event_summary <- get_population_event_counts(gs, selected_node, selected_samples)
  
  cat("Population Event Summary (selected samples only):\n")
  cat("Total events:", scales::comma(event_summary$total_events), "\n")
  cat("Mean events per sample:", round(event_summary$mean_events), "\n")
  cat("Samples with zero events:", event_summary$samples_with_zero, "\n")
  
  sample_limit_input <- readline("Max samples to analyze (Enter for all): ")
  sample_limit <- NULL
  if(sample_limit_input != "" && grepl("^\\d+$", sample_limit_input)) {
    sample_limit <- as.numeric(sample_limit_input)
  }
  
  max_cells_input <- readline("Max cells per sample (default 5000): ")
  max_cells <- 5000
  if(max_cells_input != "" && grepl("^\\d+$", max_cells_input)) {
    max_cells <- as.numeric(max_cells_input)
  }
  
  single_cell_data <- extract_single_cell_data(
    gs = gs,
    node = selected_node,
    selected_markers = selected_markers,
    selected_samples = selected_samples,
    keywords = keywords,
    sample_limit = sample_limit,
    max_cells_per_sample = max_cells
  )
  
  # Step 5: Import external metadata (optional) - FIXED VERSION
  cat("\n=== Metadata Import ===\n")
  import_choice <- readline("Import external metadata? (y/n, default n): ")
  
  external_metadata <- NULL
  if(tolower(import_choice) == "y") {
    external_metadata <- import_external_metadata(gs)
    
    if(!is.null(external_metadata)) {
      # FIXED: More robust type checking
      cat("Debug: external_metadata class:", paste(class(external_metadata), collapse = ", "), "\n")
      cat("Debug: is.data.frame():", is.data.frame(external_metadata), "\n")
      cat("Debug: inherits data.frame:", inherits(external_metadata, "data.frame"), "\n")
      
      # Use inherits() instead of is.data.frame() for more robust checking
      is_dataframe_like <- inherits(external_metadata, "data.frame") || 
        inherits(external_metadata, "tbl_df") || 
        inherits(external_metadata, "tbl")
      
      if(!is_dataframe_like) {
        cat("Warning: External metadata is not a data frame-like object\n")
        
        # Only try to extract from list if it's actually not data.frame-like
        if(is.list(external_metadata)) {
          cat("Attempting to extract data frame from list...\n")
          
          # Try to find a data frame in the list
          df_elements <- external_metadata[sapply(external_metadata, function(x) {
            inherits(x, "data.frame") || inherits(x, "tbl_df") || inherits(x, "tbl")
          })]
          
          if(length(df_elements) > 0) {
            external_metadata <- df_elements[[1]]
            cat("Extracted data frame from list\n")
          } else {
            cat("Error: No data frame found in returned list. Skipping external metadata.\n")
            external_metadata <- NULL
          }
        } else {
          cat("Error: Metadata is not a list or data frame. Skipping external metadata.\n")
          external_metadata <- NULL
        }
      }
      
      # Final validation and conversion if needed
      if(!is.null(external_metadata)) {
        # Ensure it's a proper tibble for consistency
        if(!inherits(external_metadata, "tbl_df")) {
          tryCatch({
            external_metadata <- as_tibble(external_metadata)
            cat("Converted to tibble\n")
          }, error = function(e) {
            cat("Error: Could not convert external metadata to tibble:", e$message, "\n")
            external_metadata <- NULL
          })
        } else {
          cat("Metadata is already a proper tibble\n")
        }
      }
      
      # FIXED: Use the new safe merging function
      if(!is.null(external_metadata) && inherits(external_metadata, "data.frame")) {
        single_cell_data <- merge_external_metadata_safely(single_cell_data, external_metadata)
      }
    }
  }
  
  # Step 6: UMAP parameters
  cat("\n=== UMAP Parameters ===\n")
  
  n_neighbors_input <- readline("Number of neighbors (default 15): ")
  n_neighbors <- 15
  if(n_neighbors_input != "" && grepl("^\\d+$", n_neighbors_input)) {
    n_neighbors <- as.numeric(n_neighbors_input)
  }
  
  min_dist_input <- readline("Minimum distance (default 0.1): ")
  min_dist <- 0.1
  if(min_dist_input != "" && grepl("^[0-9.]+$", min_dist_input)) {
    min_dist <- as.numeric(min_dist_input)
  }
  
  cat("Data transformation options:\n")
  cat("1. asinh (recommended for flow cytometry)\n")
  cat("2. log10\n")
  cat("3. sqrt\n")
  cat("4. none\n")
  transform_choice <- readline("Choose transformation (1-4, default 1): ")
  
  transform_data <- switch(
    if(transform_choice == "") "1" else transform_choice,
    "1" = "asinh",
    "2" = "log10", 
    "3" = "sqrt",
    "4" = "none",
    "asinh"
  )
  
  # Step 7: Compute UMAP
  umap_result <- compute_umap_embedding(
    single_cell_data = single_cell_data,
    n_neighbors = n_neighbors,
    min_dist = min_dist,
    transform_data = transform_data
  )
  
  # Step 8: Initialize enhanced visualization session
  cat("\n=== Initializing Visualization Session ===\n")
  initialize_visualization_session(umap_result)
  
  # Step 9: Interactive visualization menu
  cat("\nStarting interactive visualization...\n")
  cat("You can now create, save, and manage multiple plots in this session.\n")
  
  visualization_results <- interactive_visualization_menu()
  
  # Return comprehensive results
  results <- list(
    node = selected_node,
    markers = selected_markers,
    single_cell_data = single_cell_data,
    external_metadata = external_metadata,
    selected_samples = selected_samples,
    excluded_samples = setdiff(sampleNames(gs), selected_samples),
    event_summary = event_summary,
    umap_result = umap_result,
    visualization_session = visualization_results,
    final_data = umap_result$data
  )
  
  cat("\n=== Enhanced Analysis Complete ===\n")
  cat("Results contain:\n")
  cat("- Original node:", basename(selected_node), "\n")
  cat("- Markers analyzed:", length(selected_markers$markers), "\n")
  cat("- Samples included:", length(selected_samples), "\n")
  cat("- Samples excluded:", length(results$excluded_samples), "\n")
  cat("- Cells analyzed:", scales::comma(nrow(results$final_data)), "\n")
  cat("- External metadata:", !is.null(external_metadata), "\n")
  cat("- Plots created:", length(visualization_results$plots), "\n")
  cat("- Session plots stored in global variable for continued access\n")
  
  return(results)
}

# ============================================================================
# UTILITY FUNCTIONS FOR ENHANCED FEATURES
# ============================================================================

export_all_session_plots <- function(directory = "umap_plots", 
                                     width = 8, 
                                     height = 6, 
                                     dpi = 300,
                                     format = "png") {
  
  if(length(.umap_session_plots$plots) == 0) {
    cat("No plots to export\n")
    return(invisible(NULL))
  }
  
  if(!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    cat("Created directory:", directory, "\n")
  }
  
  exported_files <- map_chr(names(.umap_session_plots$plots), function(plot_id) {
    plot_info <- .umap_session_plots$plots[[plot_id]]
    
    clean_title <- str_replace_all(plot_info$title, "[^A-Za-z0-9_-]", "_")
    clean_title <- str_replace_all(clean_title, "_+", "_")
    clean_title <- str_trim(clean_title, side = "both")
    
    filename <- paste0(plot_id, "_", clean_title, ".", format)
    filepath <- file.path(directory, filename)
    
    tryCatch({
      ggsave(filepath, plot = plot_info$plot, 
             width = width, height = height, dpi = dpi)
      cat("Exported:", filename, "\n")
      return(filepath)
    }, error = function(e) {
      cat("Error exporting", filename, ":", e$message, "\n")
      return(NA_character_)
    })
  })
  
  successful_exports <- sum(!is.na(exported_files))
  cat("\nExport complete:", successful_exports, "plots exported to", directory, "\n")
  
  return(invisible(exported_files[!is.na(exported_files)]))
}

restart_visualization_session <- function(umap_results) {
  cat("Restarting visualization session with existing UMAP results...\n")
  
  if(is.null(umap_results$umap_result)) {
    stop("Invalid UMAP results object")
  }
  
  initialize_visualization_session(umap_results$umap_result)
  interactive_visualization_menu()
  
  return(invisible(.umap_session_plots))
}

create_quick_umap_plot <- function(umap_results, color_by, plot_title = NULL) {
  
  umap_data <- umap_results$final_data
  
  if(!color_by %in% names(umap_data)) {
    stop("Column '", color_by, "' not found in UMAP data")
  }
  
  if(is.null(plot_title)) {
    plot_title <- paste("UMAP colored by", color_by)
  }
  
  p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = .data[[color_by]])) +
    geom_point(size = 0.3, alpha = 0.6) +
    labs(
      title = plot_title,
      subtitle = paste("Total cells:", scales::comma(nrow(umap_data))),
      color = color_by
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  if(is.numeric(umap_data[[color_by]])) {
    p <- p + scale_color_viridis_c(option = "plasma")
  } else {
    n_levels <- n_distinct(umap_data[[color_by]], na.rm = TRUE)
    if(n_levels <= 12) {
      p <- p + scale_color_brewer(type = "qual", palette = "Set3")
    } else {
      p <- p + scale_color_viridis_d()
    }
  }
  
  return(p)
}

export_umap_data <- function(umap_results, filename = "umap_data.csv") {
  
  export_data <- umap_results$final_data %>%
    select(-CellID)
  
  write_csv(export_data, filename)
  cat("UMAP data exported to:", filename, "\n")
  cat("Columns exported:", paste(names(export_data), collapse = ", "), "\n")
  
  return(invisible(export_data))
}

summarize_umap_by_group <- function(umap_results, group_column) {
  
  umap_data <- umap_results$final_data
  marker_names <- umap_results$umap_result$marker_names
  
  if(!group_column %in% names(umap_data)) {
    stop("Column '", group_column, "' not found in UMAP data")
  }
  
  marker_summary <- umap_data %>%
    select(all_of(c(group_column, marker_names))) %>%
    pivot_longer(cols = all_of(marker_names), names_to = "Marker", values_to = "Expression") %>%
    group_by(.data[[group_column]], Marker) %>%
    summarise(
      n_cells = n(),
      mean_expr = mean(Expression, na.rm = TRUE),
      median_expr = median(Expression, na.rm = TRUE),
      sd_expr = sd(Expression, na.rm = TRUE),
      q25 = quantile(Expression, 0.25, na.rm = TRUE),
      q75 = quantile(Expression, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  coord_summary <- umap_data %>%
    group_by(.data[[group_column]]) %>%
    summarise(
      n_cells = n(),
      mean_umap1 = mean(UMAP1),
      mean_umap2 = mean(UMAP2),
      sd_umap1 = sd(UMAP1),
      sd_umap2 = sd(UMAP2),
      .groups = "drop"
    )
  
  return(list(
    marker_summary = marker_summary,
    coordinate_summary = coord_summary
  ))
}

# ============================================================================
# TODO ENHANCED CLUSTERING ANALYSIS EXTENSION FOR UMAP FLOW CYTOMETRY
# Integrates with existing UMAP visualization framework (IN-Devleopment)
# ============================================================================

# ============================================================================
# GLOBAL CLUSTERING SESSION ENVIRONMENT
# ============================================================================

.clustering_session <- new.env()

# ============================================================================
# UTILITY FUNCTIONS TO HANDLE S4/TIBBLE COUNTING ISSUES
# ============================================================================

safe_tibble_operations <- function(data) {
  # Convert any S4 or complex objects to standard data frame, then back to tibble
  # This resolves issues with FlowSOM results and dplyr operations
  if (inherits(data, c("tbl_df", "tbl", "data.frame"))) {
    return(as_tibble(as.data.frame(data)))
  } else {
    return(as_tibble(data))
  }
}

safe_count <- function(data, group_col, count_name = "n") {
  # Use base R table instead of dplyr count to avoid S4/vector issues
  if (is.character(group_col)) {
    group_values <- data[[group_col]]
  } else {
    group_values <- pull(data, !!group_col)
  }
  
  counts <- table(group_values)
  
  result <- tibble(
    !!group_col := names(counts),
    !!count_name := as.numeric(counts)
  ) %>%
    arrange(desc(.data[[count_name]]))
  
  return(result)
}

check_clustering_packages <- function() {
  required_packages <- c("FlowSOM", "ConsensusClusterPlus", "cluster", "dbscan")
  available <- sapply(required_packages, function(pkg) {
    requireNamespace(pkg, quietly = TRUE)
  })
  
  if (!all(available)) {
    missing <- required_packages[!available]
    cat("Missing packages:", paste(missing, collapse = ", "), "\n")
    cat("Install with: BiocManager::install(c(", paste0("'", missing, "'", collapse = ", "), "))\n")
    return(FALSE)
  }
  
  return(TRUE)
}

# ============================================================================
# INTERACTIVE CLUSTERING METHOD SELECTION
# ============================================================================

select_clustering_method <- function() {
  clustering_methods <- c(
    "FlowSOM (Self-Organizing Maps) - Recommended for flow cytometry",
    "FlowSOM + ConsensusClusterPlus - Most robust option",
    "K-means clustering",
    "Hierarchical clustering", 
    "DBSCAN density-based clustering",
    "Multiple methods comparison"
  )
  
  cat("\n=== Clustering Method Selection ===\n")
  cat("Choose clustering approach for your UMAP data:\n\n")
  
  iwalk(clustering_methods, ~cat(sprintf("%d. %s\n", .y, .x)))
  
  while (TRUE) {
    choice <- readline("\nEnter choice number (1-6): ")
    
    if (grepl("^[1-6]$", choice)) {
      method_choice <- as.numeric(choice)
      selected_method <- switch(method_choice,
                                "flowsom",
                                "flowsom_consensus", 
                                "kmeans",
                                "hierarchical",
                                "dbscan",
                                "comparison"
      )
      
      cat("Selected:", str_split(clustering_methods[method_choice], " - ")[[1]][1], "\n")
      return(selected_method)
    }
    
    cat("Invalid choice. Please enter 1-6.\n")
  }
}

# ============================================================================
# FLOWSOM CLUSTERING IMPLEMENTATION
# ============================================================================

perform_flowsom_clustering <- function(umap_data, marker_names, 
                                       use_umap_coords = TRUE,
                                       xdim = 10, ydim = 10, 
                                       n_metaclusters = NULL) {
  
  cat("\n=== FlowSOM Clustering ===\n")
  
  if (!requireNamespace("FlowSOM", quietly = TRUE)) {
    stop("FlowSOM package required. Install with: BiocManager::install('FlowSOM')")
  }
  
  # Clean data structure first
  umap_data <- safe_tibble_operations(umap_data)
  
  # Data preparation
  if (use_umap_coords) {
    clustering_data <- umap_data %>%
      select(UMAP1, UMAP2, all_of(marker_names)) %>%
      as.matrix()
    cat("Using UMAP coordinates + original markers for clustering\n")
  } else {
    clustering_data <- umap_data %>%
      select(all_of(marker_names)) %>%
      as.matrix()
    cat("Using original markers only for clustering\n")
  }
  
  cat("Clustering data dimensions:", nrow(clustering_data), "cells x", ncol(clustering_data), "features\n")
  
  # Handle missing values
  if (any(is.na(clustering_data))) {
    n_na <- sum(is.na(clustering_data))
    cat("Warning: Found", n_na, "NA values. Replacing with column medians.\n")
    
    for (i in 1:ncol(clustering_data)) {
      col_data <- clustering_data[, i]
      if (any(is.na(col_data))) {
        median_val <- median(col_data, na.rm = TRUE)
        clustering_data[is.na(col_data), i] <- median_val
      }
    }
  }
  
  # FlowSOM parameters
  if (is.null(n_metaclusters)) {
    # Estimate optimal number of metaclusters
    n_cells <- nrow(clustering_data)
    n_metaclusters <- case_when(
      n_cells < 1000 ~ 5,
      n_cells < 5000 ~ 8,
      n_cells < 20000 ~ 12,
      n_cells < 50000 ~ 15,
      TRUE ~ 20
    )
    
    cat("Auto-selected", n_metaclusters, "metaclusters based on cell count\n")
    
    # Allow user to override
    override_input <- readline(paste0("Override with custom number (current: ", n_metaclusters, ", Enter to keep): "))
    if (override_input != "" && grepl("^\\d+$", override_input)) {
      n_metaclusters <- as.numeric(override_input)
      cat("Using", n_metaclusters, "metaclusters\n")
    }
  }
  
  cat("FlowSOM grid size:", xdim, "x", ydim, "\n")
  cat("Number of metaclusters:", n_metaclusters, "\n")
  
  # Run FlowSOM
  cat("Running FlowSOM clustering...\n")
  
  tryCatch({
    # Create FlowSOM object
    flowsom_result <- FlowSOM::FlowSOM(
      input = clustering_data,
      compensate = FALSE,
      transform = FALSE,
      scale = TRUE,
      colsToUse = 1:ncol(clustering_data),
      nClus = n_metaclusters,
      xdim = xdim,
      ydim = ydim,
      seed = 42
    )
    
    # Extract cluster assignments
    som_clusters <- FlowSOM::GetClusters(flowsom_result)
    metaclusters <- FlowSOM::GetMetaclusters(flowsom_result)
    
    # Add cluster information to original data
    clustered_data <- umap_data %>%
      mutate(
        SOM_Cluster = as.factor(som_clusters),
        MetaCluster = as.factor(metaclusters),
        ClusterID = paste0("Cluster_", metaclusters)
      )
    
    cat("FlowSOM clustering complete!\n")
    cat("SOM nodes:", xdim * ydim, "\n")
    cat("Metaclusters identified:", n_distinct(metaclusters), "\n")
    
    # FIXED: More robust cluster size summary creation
    cat("Creating cluster summary...\n")
    
    cluster_summary <- safe_count(clustered_data, "MetaCluster", "n_cells") %>%
      mutate(
        percentage = round(n_cells / sum(n_cells) * 100, 1),
        cumulative_pct = round(cumsum(n_cells) / sum(n_cells) * 100, 1)
      )
    
    cat("\nCluster size distribution:\n")
    print(as.data.frame(cluster_summary))
    
    return(list(
      method = "FlowSOM",
      clustered_data = clustered_data,
      flowsom_object = flowsom_result,
      cluster_summary = cluster_summary,
      n_clusters = n_distinct(metaclusters),
      parameters = list(
        xdim = xdim,
        ydim = ydim, 
        n_metaclusters = n_metaclusters,
        use_umap_coords = use_umap_coords
      )
    ))
    
  }, error = function(e) {
    stop("FlowSOM clustering failed: ", e$message)
  })
}

# ============================================================================
# FLOWSOM + CONSENSUSCLUSTERPLUS IMPLEMENTATION
# ============================================================================

perform_consensus_clustering <- function(umap_data, marker_names,
                                         use_umap_coords = TRUE,
                                         max_clusters = 20,
                                         reps = 50,
                                         sample_proportion = 0.8) {
  
  cat("\n=== FlowSOM + ConsensusClusterPlus Clustering ===\n")
  
  if (!requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
    stop("ConsensusClusterPlus package required. Install with: BiocManager::install('ConsensusClusterPlus')")
  }
  
  # Clean data structure first
  umap_data <- safe_tibble_operations(umap_data)
  
  # Data preparation  
  if (use_umap_coords) {
    clustering_data <- umap_data %>%
      select(UMAP1, UMAP2, all_of(marker_names)) %>%
      as.matrix()
  } else {
    clustering_data <- umap_data %>%
      select(all_of(marker_names)) %>%
      as.matrix()
  }
  
  # Handle large datasets by subsampling for consensus
  if (nrow(clustering_data) > 10000) {
    subsample_size <- 10000
    cat("Large dataset detected. Subsampling to", subsample_size, "cells for consensus clustering\n")
    
    subsample_indices <- sample(nrow(clustering_data), subsample_size)
    consensus_data <- clustering_data[subsample_indices, ]
  } else {
    consensus_data <- clustering_data
    subsample_indices <- 1:nrow(clustering_data)
  }
  
  cat("Consensus clustering on", nrow(consensus_data), "cells\n")
  cat("Testing up to", max_clusters, "clusters with", reps, "repetitions\n")
  
  # Run consensus clustering
  cat("Running ConsensusClusterPlus (this may take a while)...\n")
  
  tryCatch({
    consensus_result <- ConsensusClusterPlus::ConsensusClusterPlus(
      d = t(consensus_data),
      maxK = max_clusters,
      reps = reps,
      pItem = sample_proportion,
      pFeature = 1,
      clusterAlg = "km",
      distance = "euclidean",
      seed = 42,
      plot = NULL,  # Suppress automatic plots
      verbose = FALSE
    )
    
    # Analyze consensus results to find optimal k
    consensus_stats <- map_dfr(2:max_clusters, function(k) {
      consensus_matrix <- consensus_result[[k]]$consensusMatrix
      consensus_cdf <- ConsensusClusterPlus::calcICL(consensus_result, k)
      
      tibble(
        k = k,
        consensus_cdf = consensus_cdf[["clusterConsensus"]],
        area_under_cdf = sum(consensus_cdf[["cdf"]][, "y"])
      )
    })
    
    # Find optimal k using elbow method on area under CDF
    optimal_k_idx <- which.max(diff(consensus_stats$area_under_cdf)) + 1
    optimal_k <- consensus_stats$k[optimal_k_idx]
    
    cat("Optimal number of clusters (by consensus):", optimal_k, "\n")
    
    # Allow user to override optimal k
    override_input <- readline(paste0("Use different number of clusters? (current: ", optimal_k, ", Enter to keep): "))
    if (override_input != "" && grepl("^\\d+$", override_input)) {
      selected_k <- as.numeric(override_input)
      if (selected_k >= 2 && selected_k <= max_clusters) {
        optimal_k <- selected_k
        cat("Using", optimal_k, "clusters\n")
      } else {
        cat("Invalid cluster number, using optimal:", optimal_k, "\n")
      }
    }
    
    # Get final cluster assignments for subsampled data
    consensus_clusters <- consensus_result[[optimal_k]]$consensusClass
    
    # Apply clustering to full dataset using k-means with optimal k
    cat("Applying", optimal_k, "clusters to full dataset...\n")
    
    set.seed(42)
    full_kmeans <- kmeans(clustering_data, centers = optimal_k, nstart = 25, iter.max = 100)
    
    # Add cluster information to original data
    clustered_data <- umap_data %>%
      mutate(
        ConsensusCluster = as.factor(full_kmeans$cluster),
        ClusterID = paste0("Cluster_", full_kmeans$cluster)
      )
    
    # FIXED: Cluster size summary
    cluster_summary <- safe_count(clustered_data, "ConsensusCluster", "n_cells") %>%
      mutate(
        percentage = round(n_cells / sum(n_cells) * 100, 1),
        cumulative_pct = round(cumsum(n_cells) / sum(n_cells) * 100, 1)
      )
    
    cat("\nCluster size distribution:\n")
    print(as.data.frame(cluster_summary))
    
    return(list(
      method = "ConsensusClusterPlus",
      clustered_data = clustered_data,
      consensus_result = consensus_result,
      consensus_stats = consensus_stats,
      kmeans_result = full_kmeans,
      cluster_summary = cluster_summary,
      n_clusters = optimal_k,
      optimal_k = optimal_k,
      parameters = list(
        max_clusters = max_clusters,
        reps = reps,
        sample_proportion = sample_proportion,
        use_umap_coords = use_umap_coords
      )
    ))
    
  }, error = function(e) {
    stop("Consensus clustering failed: ", e$message)
  })
}

# ============================================================================
# ALTERNATIVE CLUSTERING METHODS
# ============================================================================

perform_kmeans_clustering <- function(umap_data, marker_names, 
                                      use_umap_coords = TRUE,
                                      max_clusters = 15) {
  
  cat("\n=== K-means Clustering ===\n")
  
  # Clean data structure first
  umap_data <- safe_tibble_operations(umap_data)
  
  # Data preparation
  if (use_umap_coords) {
    clustering_data <- umap_data %>%
      select(UMAP1, UMAP2, all_of(marker_names)) %>%
      as.matrix()
  } else {
    clustering_data <- umap_data %>%
      select(all_of(marker_names)) %>%
      as.matrix()
  }
  
  # Scale the data
  clustering_data_scaled <- scale(clustering_data)
  
  # Find optimal k using elbow method
  cat("Finding optimal number of clusters (testing 1 to", max_clusters, ")...\n")
  
  wss_values <- map_dbl(1:max_clusters, function(k) {
    if (k == 1) {
      return(sum(scale(clustering_data, center = TRUE, scale = FALSE)^2))
    }
    
    set.seed(42)
    kmeans_result <- kmeans(clustering_data_scaled, centers = k, nstart = 25)
    return(kmeans_result$tot.withinss)
  })
  
  # Find elbow point
  elbow_k <- 2
  if (length(wss_values) > 3) {
    diffs <- diff(wss_values)
    second_diffs <- diff(diffs)
    elbow_k <- which.max(second_diffs) + 2
  }
  
  cat("Optimal k by elbow method:", elbow_k, "\n")
  
  # Allow user to override
  override_input <- readline(paste0("Use different k? (current: ", elbow_k, ", Enter to keep): "))
  if (override_input != "" && grepl("^\\d+$", override_input)) {
    selected_k <- as.numeric(override_input)
    if (selected_k >= 1 && selected_k <= max_clusters) {
      elbow_k <- selected_k
      cat("Using k =", elbow_k, "\n")
    }
  }
  
  # Perform final k-means clustering
  set.seed(42)
  final_kmeans <- kmeans(clustering_data_scaled, centers = elbow_k, nstart = 25, iter.max = 100)
  
  # Add cluster information to original data
  clustered_data <- umap_data %>%
    mutate(
      KMeansCluster = as.factor(final_kmeans$cluster),
      ClusterID = paste0("Cluster_", final_kmeans$cluster)
    )
  
  # FIXED: Cluster size summary
  cluster_summary <- safe_count(clustered_data, "KMeansCluster", "n_cells") %>%
    mutate(percentage = round(n_cells / sum(n_cells) * 100, 1))
  
  cat("\nK-means clustering complete!\n")
  cat("Final k:", elbow_k, "\n")
  print(as.data.frame(cluster_summary))
  
  return(list(
    method = "K-means",
    clustered_data = clustered_data,
    kmeans_result = final_kmeans,
    wss_values = wss_values,
    cluster_summary = cluster_summary,
    n_clusters = elbow_k,
    optimal_k = elbow_k
  ))
}

perform_dbscan_clustering <- function(umap_data, marker_names, 
                                      use_umap_coords = TRUE,
                                      eps = NULL, min_pts = 5) {
  
  cat("\n=== DBSCAN Clustering ===\n")
  
  if (!requireNamespace("dbscan", quietly = TRUE)) {
    cat("dbscan package not available. Install with: install.packages('dbscan')\n")
    cat("Using alternative density-based clustering...\n")
    
    # Fallback to hierarchical clustering with density estimation
    return(perform_hierarchical_clustering(umap_data, marker_names, use_umap_coords))
  }
  
  # Clean data structure first
  umap_data <- safe_tibble_operations(umap_data)
  
  # Data preparation
  if (use_umap_coords) {
    clustering_data <- umap_data %>%
      select(UMAP1, UMAP2) %>%  # DBSCAN works best on UMAP coordinates
      as.matrix()
    cat("Using UMAP coordinates for DBSCAN\n")
  } else {
    clustering_data <- umap_data %>%
      select(all_of(marker_names)) %>%
      as.matrix()
    cat("Using original markers for DBSCAN\n")
  }
  
  # Scale the data
  clustering_data_scaled <- scale(clustering_data)
  
  # Estimate eps if not provided
  if (is.null(eps)) {
    # Use k-distance plot to estimate eps
    k_distances <- dbscan::kNNdist(clustering_data_scaled, k = min_pts)
    eps <- quantile(sort(k_distances), 0.95)  # Use 95th percentile
    cat("Auto-estimated eps:", round(eps, 3), "\n")
    
    # Allow user to override
    override_input <- readline(paste0("Use different eps? (current: ", round(eps, 3), ", Enter to keep): "))
    if (override_input != "" && grepl("^[0-9.]+$", override_input)) {
      eps <- as.numeric(override_input)
      cat("Using eps =", eps, "\n")
    }
  }
  
  cat("DBSCAN parameters: eps =", eps, ", min_pts =", min_pts, "\n")
  
  # Perform DBSCAN clustering
  dbscan_result <- dbscan::dbscan(clustering_data_scaled, eps = eps, minPts = min_pts)
  
  # Add cluster information to original data
  clustered_data <- umap_data %>%
    mutate(
      DBSCANCluster = as.factor(ifelse(dbscan_result$cluster == 0, "Noise", 
                                       paste0("Cluster_", dbscan_result$cluster))),
      ClusterID = as.character(DBSCANCluster)
    )
  
  n_clusters <- max(dbscan_result$cluster)
  n_noise <- sum(dbscan_result$cluster == 0)
  
  # FIXED: Cluster size summary
  cluster_summary <- safe_count(clustered_data, "DBSCANCluster", "n_cells") %>%
    mutate(percentage = round(n_cells / sum(n_cells) * 100, 1))
  
  cat("\nDBSCAN clustering complete!\n")
  cat("Clusters found:", n_clusters, "\n")
  cat("Noise points:", n_noise, "(", round(n_noise/nrow(clustered_data)*100, 1), "%)\n")
  print(as.data.frame(cluster_summary))
  
  return(list(
    method = "DBSCAN",
    clustered_data = clustered_data,
    dbscan_result = dbscan_result,
    cluster_summary = cluster_summary,
    n_clusters = n_clusters,
    n_noise = n_noise,
    parameters = list(eps = eps, min_pts = min_pts)
  ))
}

perform_hierarchical_clustering <- function(umap_data, marker_names, 
                                            use_umap_coords = TRUE,
                                            max_clusters = 15,
                                            linkage_method = "ward.D2") {
  
  cat("\n=== Hierarchical Clustering ===\n")
  
  # Clean data structure first
  umap_data <- safe_tibble_operations(umap_data)
  
  # Data preparation
  if (use_umap_coords) {
    clustering_data <- umap_data %>%
      select(UMAP1, UMAP2, all_of(marker_names)) %>%
      as.matrix()
  } else {
    clustering_data <- umap_data %>%
      select(all_of(marker_names)) %>%
      as.matrix()
  }
  
  # Scale the data
  clustering_data_scaled <- scale(clustering_data)
  
  # For large datasets, subsample for dendrogram computation
  if (nrow(clustering_data_scaled) > 5000) {
    subsample_size <- 5000
    cat("Large dataset detected. Using", subsample_size, "cells for dendrogram\n")
    
    subsample_indices <- sample(nrow(clustering_data_scaled), subsample_size)
    hclust_data <- clustering_data_scaled[subsample_indices, ]
  } else {
    hclust_data <- clustering_data_scaled
  }
  
  cat("Computing hierarchical clustering (method:", linkage_method, ")...\n")
  
  # Compute distance matrix and hierarchical clustering
  dist_matrix <- dist(hclust_data, method = "euclidean")
  hclust_result <- hclust(dist_matrix, method = linkage_method)
  
  # Find optimal number of clusters using silhouette analysis
  cat("Finding optimal number of clusters...\n")
  
  sil_scores <- map_dbl(2:min(max_clusters, nrow(hclust_data)-1), function(k) {
    clusters <- cutree(hclust_result, k = k)
    sil_result <- cluster::silhouette(clusters, dist_matrix)
    mean(sil_result[, "sil_width"])
  })
  
  optimal_k <- which.max(sil_scores) + 1
  cat("Optimal k by silhouette analysis:", optimal_k, "\n")
  
  # Allow user to override
  override_input <- readline(paste0("Use different k? (current: ", optimal_k, ", Enter to keep): "))
  if (override_input != "" && grepl("^\\d+$", override_input)) {
    selected_k <- as.numeric(override_input)
    if (selected_k >= 2 && selected_k <= max_clusters) {
      optimal_k <- selected_k
      cat("Using k =", optimal_k, "\n")
    }
  }
  
  # Cut dendrogram to get clusters for full dataset
  if (nrow(clustering_data_scaled) > 5000) {
    # Apply clustering to full dataset using same tree structure
    full_dist <- dist(clustering_data_scaled, method = "euclidean")
    full_hclust <- hclust(full_dist, method = linkage_method)
    clusters <- cutree(full_hclust, k = optimal_k)
  } else {
    clusters <- cutree(hclust_result, k = optimal_k)
  }
  
  # Add cluster information to original data
  clustered_data <- umap_data %>%
    mutate(
      HierarchicalCluster = as.factor(clusters),
      ClusterID = paste0("Cluster_", clusters)
    )
  
  # FIXED: Cluster size summary
  cluster_summary <- safe_count(clustered_data, "HierarchicalCluster", "n_cells") %>%
    mutate(percentage = round(n_cells / sum(n_cells) * 100, 1))
  
  cat("\nHierarchical clustering complete!\n")
  cat("Final k:", optimal_k, "\n")
  print(as.data.frame(cluster_summary))
  
  return(list(
    method = "Hierarchical",
    clustered_data = clustered_data,
    hclust_result = if (nrow(clustering_data_scaled) <= 5000) hclust_result else full_hclust,
    sil_scores = sil_scores,
    cluster_summary = cluster_summary,
    n_clusters = optimal_k,
    optimal_k = optimal_k,
    parameters = list(linkage_method = linkage_method)
  ))
}

# ============================================================================
# CLUSTER CHARACTERIZATION AND MFI ANALYSIS
# ============================================================================

analyze_cluster_mfi <- function(clustered_data, marker_names, 
                                cluster_column = "ClusterID") {
  
  cat("\n=== Cluster MFI Analysis ===\n")
  
  if (!cluster_column %in% names(clustered_data)) {
    available_cluster_cols <- names(clustered_data)[str_detect(names(clustered_data), "Cluster")]
    cat("Cluster column not found. Available cluster columns:\n")
    iwalk(available_cluster_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
    
    choice <- readline("Enter column number or name: ")
    
    if (grepl("^\\d+$", choice)) {
      col_num <- as.numeric(choice)
      if (!is.na(col_num) && col_num >= 1 && col_num <= length(available_cluster_cols)) {
        cluster_column <- available_cluster_cols[col_num]
      }
    } else if (choice %in% available_cluster_cols) {
      cluster_column <- choice
    } else {
      stop("Invalid cluster column selection")
    }
  }
  
  cat("Analyzing MFI by", cluster_column, "\n")
  
  # Clean data structure first
  clustered_data <- safe_tibble_operations(clustered_data)
  
  # FIXED: Convert to standard data frame first and handle operations more robustly
  analysis_data <- clustered_data %>%
    select(all_of(c(cluster_column, marker_names))) %>%
    as.data.frame()  # Convert to standard data frame to avoid tibble/S4 issues
  
  # Convert back to tibble for dplyr operations
  analysis_data <- as_tibble(analysis_data)
  
  # Calculate MFI statistics for each cluster
  cat("Computing MFI statistics...\n")
  
  mfi_summary <- analysis_data %>%
    pivot_longer(cols = all_of(marker_names), names_to = "Marker", values_to = "Expression") %>%
    group_by(.data[[cluster_column]], Marker) %>%
    summarise(
      n_cells = n(),
      mean_expr = mean(Expression, na.rm = TRUE),
      median_expr = median(Expression, na.rm = TRUE),
      sd_expr = sd(Expression, na.rm = TRUE),
      q25 = quantile(Expression, 0.25, na.rm = TRUE),
      q75 = quantile(Expression, 0.75, na.rm = TRUE),
      min_expr = min(Expression, na.rm = TRUE),
      max_expr = max(Expression, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(.data[[cluster_column]], Marker)
  
  # Calculate relative expression (z-scores across clusters for each marker)
  cat("Computing z-scores...\n")
  
  mfi_zscore <- mfi_summary %>%
    group_by(Marker) %>%
    mutate(
      mean_zscore = scale(mean_expr)[,1],
      median_zscore = scale(median_expr)[,1]
    ) %>%
    ungroup()
  
  # Create cluster signatures (top markers per cluster)
  cat("Creating cluster signatures...\n")
  
  cluster_signatures <- mfi_zscore %>%
    group_by(.data[[cluster_column]]) %>%
    arrange(desc(mean_zscore)) %>%
    slice_head(n = 5) %>%
    summarise(
      top_markers = paste(Marker, collapse = ", "),
      signature_score = round(mean(mean_zscore), 2),
      .groups = "drop"
    )
  
  # FIXED: More robust cluster overview creation
  cat("Creating cluster overview...\n")
  
  cluster_overview <- safe_count(analysis_data, cluster_column, "n_cells") %>%
    mutate(percentage = round(n_cells / sum(n_cells) * 100, 1)) %>%
    left_join(cluster_signatures, by = cluster_column)
  
  cat("\nCluster MFI analysis complete!\n")
  cat("Clusters analyzed:", n_distinct(analysis_data[[cluster_column]]), "\n")
  cat("Markers analyzed:", length(marker_names), "\n")
  
  cat("\nCluster overview:\n")
  # Force proper printing
  overview_df <- as.data.frame(cluster_overview)
  print(overview_df)
  
  return(list(
    mfi_summary = mfi_summary,
    mfi_zscore = mfi_zscore,
    cluster_signatures = cluster_signatures,
    cluster_overview = cluster_overview,
    cluster_column = cluster_column
  ))
}

# ============================================================================
# CLUSTER VISUALIZATION FUNCTIONS
# ============================================================================

create_cluster_umap_plot <- function(clustered_data, cluster_column = "ClusterID",
                                     title = NULL, show_percentages = TRUE) {
  
  if (is.null(title)) {
    title <- paste("UMAP with", str_replace(cluster_column, "Cluster", ""), "Clusters")
  }
  
  # Clean data structure first
  clustered_data <- safe_tibble_operations(clustered_data)
  
  # FIXED: Use safe_count instead of dplyr count()
  cluster_stats <- safe_count(clustered_data, cluster_column, "n_cells") %>%
    mutate(
      percentage = round(n_cells / sum(n_cells) * 100, 1),
      label = if (show_percentages) {
        paste0(.data[[cluster_column]], "\n(", percentage, "%)")
      } else {
        as.character(.data[[cluster_column]])
      }
    )
  
  # Add labels to data - convert to data frame first to avoid S4/tibble issues
  plot_data <- as.data.frame(clustered_data) %>%
    as_tibble() %>%
    left_join(cluster_stats %>% select(all_of(cluster_column), label), 
              by = cluster_column)
  
  # Create color palette
  n_clusters <- n_distinct(clustered_data[[cluster_column]])
  
  if (n_clusters <= 12) {
    colors <- RColorBrewer::brewer.pal(min(n_clusters, 11), "Set3")
    if (n_clusters == 12) colors <- c(colors, "#FF1493")  # Add pink for 12th cluster
  } else {
    colors <- rainbow(n_clusters, alpha = 0.8)
  }
  
  p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2, color = .data[[cluster_column]])) +
    geom_point(size = 0.5, alpha = 0.7) +
    scale_color_manual(values = colors, name = "Cluster") +
    labs(
      title = title,
      subtitle = paste("Total cells:", scales::comma(nrow(clustered_data)), 
                       "| Clusters:", n_clusters)
    ) +
    theme_minimal() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(), 
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
  
  return(p)
}

create_cluster_heatmap <- function(mfi_analysis, cluster_column = "ClusterID",
                                   value_type = "mean_zscore", 
                                   title = NULL) {
  
  if (is.null(title)) {
    title <- paste("Cluster Expression Heatmap (", str_replace(value_type, "_", " "), ")", sep = "")
  }
  
  # Prepare data for heatmap
  heatmap_data <- mfi_analysis$mfi_zscore %>%
    select(all_of(c(cluster_column, "Marker", value_type))) %>%
    pivot_wider(names_from = Marker, values_from = all_of(value_type)) %>%
    column_to_rownames(cluster_column) %>%
    as.matrix()
  
  # Convert to long format for ggplot
  heatmap_long <- heatmap_data %>%
    as_tibble(rownames = "Cluster") %>%
    pivot_longer(-Cluster, names_to = "Marker", values_to = "Expression") %>%
    mutate(
      Cluster = factor(Cluster, levels = rownames(heatmap_data)),
      Marker = factor(Marker, levels = colnames(heatmap_data))
    )
  
  # Create heatmap
  p <- ggplot(heatmap_long, aes(x = Marker, y = Cluster, fill = Expression)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = 0, name = str_replace(value_type, "_", "\n")
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      panel.grid = element_blank()
    )
  
  return(p)
}

create_mfi_violin_plots <- function(clustered_data, marker_names, 
                                    cluster_column = "ClusterID",
                                    max_markers = 9) {
  
  # Clean data structure first
  clustered_data <- safe_tibble_operations(clustered_data)
  
  # Select top markers or limit number
  if (length(marker_names) > max_markers) {
    cat("Too many markers for violin plots. Selecting first", max_markers, "\n")
    selected_markers <- head(marker_names, max_markers)
  } else {
    selected_markers <- marker_names
  }
  
  # Prepare data
  violin_data <- clustered_data %>%
    select(all_of(c(cluster_column, selected_markers))) %>%
    pivot_longer(cols = all_of(selected_markers), 
                 names_to = "Marker", values_to = "Expression")
  
  # Create violin plots
  violin_plots <- map(selected_markers, function(marker) {
    marker_data <- violin_data %>% dplyr::filter(Marker == marker)
    
    ggplot(marker_data, aes(x = .data[[cluster_column]], y = Expression, 
                            fill = .data[[cluster_column]])) +
      geom_violin(alpha = 0.7, scale = "width") +
      geom_boxplot(width = 0.1, fill = "white", alpha = 0.8) +
      labs(
        title = marker,
        x = "Cluster",
        y = "Expression"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 10, hjust = 0.5)
      )
  })
  
  # Arrange plots
  ncol_grid <- ceiling(sqrt(length(violin_plots)))
  combined_plot <- wrap_plots(violin_plots, ncol = ncol_grid)
  
  return(combined_plot)
}

# ============================================================================
# CLUSTER COMPARISON AND DIFFERENTIAL EXPRESSION
# ============================================================================

perform_cluster_differential_analysis <- function(mfi_analysis, 
                                                  reference_cluster = NULL,
                                                  min_fold_change = 1.5,
                                                  cluster_column = "ClusterID") {
  
  cat("\n=== Cluster Differential Expression Analysis ===\n")
  
  mfi_data <- mfi_analysis$mfi_summary
  
  # Select reference cluster if not provided
  if (is.null(reference_cluster)) {
    cluster_sizes <- mfi_analysis$cluster_overview %>%
      arrange(desc(n_cells))
    
    cat("Available clusters:\n")
    iwalk(cluster_sizes[[cluster_column]], ~cat(sprintf("%d. %s (%d cells)\n", 
                                                        .y, .x, 
                                                        cluster_sizes$n_cells[.y])))
    
    ref_choice <- readline("Select reference cluster number (Enter for largest): ")
    
    if (ref_choice == "") {
      reference_cluster <- cluster_sizes[[cluster_column]][1]  # Largest cluster
    } else if (grepl("^\\d+$", ref_choice)) {
      ref_num <- as.numeric(ref_choice)
      if (!is.na(ref_num) && ref_num >= 1 && ref_num <= nrow(cluster_sizes)) {
        reference_cluster <- cluster_sizes[[cluster_column]][ref_num]
      } else {
        reference_cluster <- cluster_sizes[[cluster_column]][1]
      }
    } else {
      reference_cluster <- cluster_sizes[[cluster_column]][1]
    }
  }
  
  cat("Using reference cluster:", reference_cluster, "\n")
  
  # Get reference cluster expression values
  reference_expr <- mfi_data %>%
    dplyr::filter(.data[[cluster_column]] == reference_cluster) %>%
    select(Marker, mean_expr) %>%
    rename(ref_mean = mean_expr)
  
  # Calculate fold changes and differences for all clusters
  diff_expr <- mfi_data %>%
    left_join(reference_expr, by = "Marker") %>%
    mutate(
      fold_change = mean_expr / ref_mean,
      log2_fold_change = log2(fold_change),
      abs_diff = abs(mean_expr - ref_mean),
      is_upregulated = fold_change >= min_fold_change,
      is_downregulated = fold_change <= (1/min_fold_change),
      is_differential = is_upregulated | is_downregulated
    ) %>%
    arrange(.data[[cluster_column]], desc(abs(log2_fold_change)))
  
  # FIXED: Summarize differential markers per cluster with proper variable handling
  cluster_diff_summary <- diff_expr %>%
    dplyr::filter(.data[[cluster_column]] != reference_cluster) %>%
    group_by(.data[[cluster_column]]) %>%
    summarise(
      n_upregulated = sum(is_upregulated, na.rm = TRUE),
      n_downregulated = sum(is_downregulated, na.rm = TRUE),
      n_differential = sum(is_differential, na.rm = TRUE),
      # Fix the variable access issue by using different approach
      top_upregulated = {
        up_markers <- Marker[is_upregulated]
        paste(head(up_markers, 3), collapse = ", ")
      },
      top_downregulated = {
        down_markers <- Marker[is_downregulated] 
        paste(head(down_markers, 3), collapse = ", ")
      },
      .groups = "drop"
    ) %>%
    arrange(desc(n_differential))
  
  cat("\nDifferential expression summary (vs", reference_cluster, "):\n")
  print(as.data.frame(cluster_diff_summary))
  
  return(list(
    differential_expression = diff_expr,
    cluster_summary = cluster_diff_summary,
    reference_cluster = reference_cluster,
    min_fold_change = min_fold_change
  ))
}

# ============================================================================
# CLUSTER DATA EXTRACTION AND EXPORT
# ============================================================================

extract_cluster_data <- function(clustered_data, target_cluster, 
                                 cluster_column = "ClusterID") {
  
  cat("\n=== Cluster Data Extraction ===\n")
  
  # Clean data structure first
  clustered_data <- safe_tibble_operations(clustered_data)
  
  if (!target_cluster %in% clustered_data[[cluster_column]]) {
    available_clusters <- unique(clustered_data[[cluster_column]])
    cat("Cluster not found. Available clusters:\n")
    iwalk(available_clusters, ~cat(sprintf("%d. %s\n", .y, .x)))
    
    cluster_choice <- readline("Enter cluster number or name: ")
    
    if (grepl("^\\d+$", cluster_choice)) {
      cluster_num <- as.numeric(cluster_choice)
      if (!is.na(cluster_num) && cluster_num >= 1 && cluster_num <= length(available_clusters)) {
        target_cluster <- available_clusters[cluster_num]
      } else {
        stop("Invalid cluster selection")
      }
    } else if (cluster_choice %in% available_clusters) {
      target_cluster <- cluster_choice
    } else {
      stop("Invalid cluster selection")
    }
  }
  
  # Extract cluster data
  cluster_data <- clustered_data %>%
    dplyr::filter(.data[[cluster_column]] == target_cluster)
  
  cat("Extracted", nrow(cluster_data), "cells from", target_cluster, "\n")
  
  # Summary statistics for the cluster
  numeric_cols <- names(cluster_data)[map_lgl(cluster_data, is.numeric)]
  marker_cols <- setdiff(numeric_cols, c("UMAP1", "UMAP2"))
  
  if (length(marker_cols) > 0) {
    cluster_summary <- cluster_data %>%
      select(all_of(marker_cols)) %>%
      summarise(across(everything(), list(
        mean = ~mean(.x, na.rm = TRUE),
        median = ~median(.x, na.rm = TRUE),
        sd = ~sd(.x, na.rm = TRUE)
      ))) %>%
      pivot_longer(everything(), names_to = "stat", values_to = "value") %>%
      separate(stat, into = c("marker", "statistic"), sep = "_(?=[^_]*$)") %>%
      pivot_wider(names_from = statistic, values_from = value) %>%
      arrange(desc(mean))
    
    cat("\nCluster expression summary (top markers by mean):\n")
    print(head(as.data.frame(cluster_summary), 10))
  }
  
  return(list(
    cluster_data = cluster_data,
    cluster_name = target_cluster,
    n_cells = nrow(cluster_data),
    summary_stats = if(exists("cluster_summary")) cluster_summary else NULL
  ))
}

export_cluster_results <- function(clustering_result, mfi_analysis = NULL,
                                   directory = "cluster_analysis", 
                                   include_plots = TRUE) {
  
  cat("\n=== Exporting Cluster Results ===\n")
  
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    cat("Created directory:", directory, "\n")
  }
  
  # Export clustered data
  clustered_data_file <- file.path(directory, "clustered_data.csv")
  write_csv(clustering_result$clustered_data, clustered_data_file)
  cat("Exported clustered data to:", clustered_data_file, "\n")
  
  # Export cluster summary
  summary_file <- file.path(directory, "cluster_summary.csv")
  write_csv(clustering_result$cluster_summary, summary_file)
  cat("Exported cluster summary to:", summary_file, "\n")
  
  # Export MFI analysis if available
  if (!is.null(mfi_analysis)) {
    mfi_file <- file.path(directory, "cluster_mfi_analysis.csv")
    write_csv(mfi_analysis$mfi_summary, mfi_file)
    
    mfi_zscore_file <- file.path(directory, "cluster_mfi_zscore.csv")
    write_csv(mfi_analysis$mfi_zscore, mfi_zscore_file)
    
    signatures_file <- file.path(directory, "cluster_signatures.csv")
    write_csv(mfi_analysis$cluster_signatures, signatures_file)
    
    cat("Exported MFI analysis files\n")
  }
  
  # Export plots if requested
  if (include_plots && !is.null(mfi_analysis)) {
    tryCatch({
      # UMAP cluster plot
      cluster_plot <- create_cluster_umap_plot(clustering_result$clustered_data)
      ggsave(file.path(directory, "cluster_umap_plot.png"), 
             plot = cluster_plot, width = 10, height = 8, dpi = 300)
      
      # MFI heatmap
      heatmap_plot <- create_cluster_heatmap(mfi_analysis)
      ggsave(file.path(directory, "cluster_heatmap.png"), 
             plot = heatmap_plot, width = 12, height = 8, dpi = 300)
      
      cat("Exported visualization plots\n")
      
    }, error = function(e) {
      cat("Warning: Could not export plots:", e$message, "\n")
    })
  }
  
  # Create analysis summary report
  report_file <- file.path(directory, "analysis_report.txt")
  
  report_content <- paste0(
    "Cluster Analysis Report\n",
    "======================\n\n",
    "Analysis Date: ", Sys.time(), "\n",
    "Clustering Method: ", clustering_result$method, "\n",
    "Total Cells: ", scales::comma(nrow(clustering_result$clustered_data)), "\n",
    "Number of Clusters: ", clustering_result$n_clusters, "\n\n",
    
    "Cluster Size Distribution:\n",
    paste(capture.output(print(clustering_result$cluster_summary)), collapse = "\n"),
    "\n\n",
    
    if (!is.null(mfi_analysis)) {
      paste0(
        "Cluster Signatures (Top Markers):\n",
        paste(capture.output(print(mfi_analysis$cluster_signatures)), collapse = "\n"),
        "\n\n"
      )
    } else "",
    
    "Files Exported:\n",
    "- clustered_data.csv: Full dataset with cluster assignments\n",
    "- cluster_summary.csv: Cluster size and percentage summary\n",
    if (!is.null(mfi_analysis)) {
      paste0(
        "- cluster_mfi_analysis.csv: Detailed MFI statistics per cluster\n",
        "- cluster_mfi_zscore.csv: Z-scored MFI values for comparison\n",
        "- cluster_signatures.csv: Top markers defining each cluster\n"
      )
    } else "",
    if (include_plots) "- PNG plots for visualization\n" else "",
    "\nEnd of Report\n"
  )
  
  writeLines(report_content, report_file)
  cat("Exported analysis report to:", report_file, "\n")
  
  cat("\nExport complete! Files saved to:", directory, "\n")
  
  return(invisible(directory))
}

# ============================================================================
# MAIN CLUSTERING WORKFLOW FUNCTION
# ============================================================================

perform_clustering_analysis <- function(umap_results, method = NULL) {
  
  cat("\n=== Starting Clustering Analysis ===\n")
  
  if (!check_clustering_packages()) {
    stop("Required packages not available. Please install missing packages.")
  }
  
  # Extract necessary data from UMAP results
  if (is.null(umap_results$final_data) || is.null(umap_results$umap_result$marker_names)) {
    stop("Invalid UMAP results. Please run UMAP analysis first.")
  }
  
  umap_data <- safe_tibble_operations(umap_results$final_data)
  marker_names <- umap_results$umap_result$marker_names
  
  cat("Input data:", scales::comma(nrow(umap_data)), "cells\n")
  cat("Available markers:", length(marker_names), "\n")
  cat("Markers:", paste(head(marker_names, 5), collapse = ", "), 
      if(length(marker_names) > 5) "..." else "", "\n")
  
  # Select clustering method if not provided
  if (is.null(method)) {
    method <- select_clustering_method()
  }
  
  cat("Selected clustering method:", method, "\n")
  
  # Store in session environment
  .clustering_session$umap_data <- umap_data
  .clustering_session$marker_names <- marker_names
  .clustering_session$method <- method
  
  # Perform clustering based on selected method
  clustering_result <- switch(method,
                              "flowsom" = perform_flowsom_clustering(umap_data, marker_names),
                              "flowsom_consensus" = perform_consensus_clustering(umap_data, marker_names),
                              "kmeans" = perform_kmeans_clustering(umap_data, marker_names),
                              "hierarchical" = perform_hierarchical_clustering(umap_data, marker_names),
                              "dbscan" = perform_dbscan_clustering(umap_data, marker_names),
                              "comparison" = {
                                cat("Multiple methods comparison not yet implemented. Using FlowSOM.\n")
                                perform_flowsom_clustering(umap_data, marker_names)
                              }
  )
  
  if (is.null(clustering_result)) {
    stop("Clustering failed")
  }
  
  # Store clustering results
  .clustering_session$clustering_result <- clustering_result
  
  # Perform MFI analysis
  cat("\nPerforming MFI analysis...\n")
  mfi_analysis <- analyze_cluster_mfi(clustering_result$clustered_data, marker_names)
  .clustering_session$mfi_analysis <- mfi_analysis
  
  # Interactive post-clustering analysis menu
  interactive_clustering_menu(clustering_result, mfi_analysis, marker_names)
  
  # Return comprehensive results
  final_results <- list(
    clustering_result = clustering_result,
    mfi_analysis = mfi_analysis,
    method = method,
    umap_data = umap_data,
    marker_names = marker_names
  )
  
  return(final_results)
}

# ============================================================================
# INTERACTIVE CLUSTERING ANALYSIS MENU
# ============================================================================

interactive_clustering_menu <- function(clustering_result, mfi_analysis, marker_names) {
  
  while (TRUE) {
    cat("\n=== Interactive Clustering Analysis Menu ===\n")
    cat("Method:", clustering_result$method, "\n")
    cat("Clusters:", clustering_result$n_clusters, "\n")
    cat("Cells:", scales::comma(nrow(clustering_result$clustered_data)), "\n\n")
    
    menu_options <- c(
      "View cluster UMAP plot",
      "View cluster MFI heatmap", 
      "View marker expression violin plots",
      "Perform differential expression analysis",
      "Extract specific cluster data",
      "Export all clustering results",
      "Re-run clustering with different parameters",
      "Return to main analysis"
    )
    
    iwalk(menu_options, ~cat(sprintf("%d. %s\n", .y, .x)))
    
    choice <- readline("\nEnter choice (1-8): ")
    
    if (!grepl("^[1-8]$", choice)) {
      cat("Invalid choice. Please enter 1-8.\n")
      next
    }
    
    choice_num <- as.numeric(choice)
    
    if (choice_num == 1) {
      # UMAP plot
      cluster_plot <- create_cluster_umap_plot(clustering_result$clustered_data)
      print(cluster_plot)
      
      save_choice <- readline("Save this plot? (y/n): ")
      if (tolower(save_choice) == "y") {
        filename <- readline("Enter filename (without extension): ")
        if (filename != "") {
          ggsave(paste0(filename, ".png"), plot = cluster_plot, 
                 width = 10, height = 8, dpi = 300)
          cat("Plot saved as:", paste0(filename, ".png"), "\n")
        }
      }
      
    } else if (choice_num == 2) {
      # MFI heatmap
      heatmap_plot <- create_cluster_heatmap(mfi_analysis)
      print(heatmap_plot)
      
      save_choice <- readline("Save this plot? (y/n): ")
      if (tolower(save_choice) == "y") {
        filename <- readline("Enter filename (without extension): ")
        if (filename != "") {
          ggsave(paste0(filename, ".png"), plot = heatmap_plot,
                 width = 12, height = 8, dpi = 300)
          cat("Plot saved as:", paste0(filename, ".png"), "\n")
        }
      }
      
    } else if (choice_num == 3) {
      # Violin plots
      violin_plot <- create_mfi_violin_plots(clustering_result$clustered_data, marker_names)
      print(violin_plot)
      
      save_choice <- readline("Save this plot? (y/n): ")
      if (tolower(save_choice) == "y") {
        filename <- readline("Enter filename (without extension): ")
        if (filename != "") {
          ggsave(paste0(filename, ".png"), plot = violin_plot,
                 width = 12, height = 10, dpi = 300)
          cat("Plot saved as:", paste0(filename, ".png"), "\n")
        }
      }
      
    } else if (choice_num == 4) {
      # Differential expression
      diff_analysis <- perform_cluster_differential_analysis(mfi_analysis)
      .clustering_session$diff_analysis <- diff_analysis
      
      readline("Press Enter to continue...")
      
    } else if (choice_num == 5) {
      # Extract cluster data
      available_clusters <- unique(clustering_result$clustered_data$ClusterID)
      cat("Available clusters:\n")
      iwalk(available_clusters, ~cat(sprintf("%d. %s\n", .y, .x)))
      
      cluster_choice <- readline("Enter cluster number or name: ")
      
      target_cluster <- NULL
      if (grepl("^\\d+$", cluster_choice)) {
        cluster_num <- as.numeric(cluster_choice)
        if (!is.na(cluster_num) && cluster_num >= 1 && cluster_num <= length(available_clusters)) {
          target_cluster <- available_clusters[cluster_num]
        }
      } else if (cluster_choice %in% available_clusters) {
        target_cluster <- cluster_choice
      }
      
      if (!is.null(target_cluster)) {
        extracted_data <- extract_cluster_data(clustering_result$clustered_data, target_cluster)
        .clustering_session$extracted_cluster <- extracted_data
        
        export_choice <- readline("Export this cluster data to CSV? (y/n): ")
        if (tolower(export_choice) == "y") {
          filename <- readline("Enter filename (without .csv extension): ")
          if (filename == "") filename <- paste0("cluster_", target_cluster, "_data")
          
          write_csv(extracted_data$cluster_data, paste0(filename, ".csv"))
          cat("Cluster data exported to:", paste0(filename, ".csv"), "\n")
        }
      } else {
        cat("Invalid cluster selection\n")
      }
      
    } else if (choice_num == 6) {
      # Export all results
      export_dir <- readline("Export directory (default: cluster_analysis): ")
      if (export_dir == "") export_dir <- "cluster_analysis"
      
      export_cluster_results(clustering_result, mfi_analysis, export_dir)
      
    } else if (choice_num == 7) {
      # Re-run clustering
      cat("Re-running clustering with different parameters...\n")
      new_method <- select_clustering_method()
      
      new_clustering_result <- switch(new_method,
                                      "flowsom" = perform_flowsom_clustering(.clustering_session$umap_data, marker_names),
                                      "flowsom_consensus" = perform_consensus_clustering(.clustering_session$umap_data, marker_names),
                                      "kmeans" = perform_kmeans_clustering(.clustering_session$umap_data, marker_names),
                                      "hierarchical" = perform_hierarchical_clustering(.clustering_session$umap_data, marker_names),
                                      "dbscan" = perform_dbscan_clustering(.clustering_session$umap_data, marker_names)
      )
      
      if (!is.null(new_clustering_result)) {
        clustering_result <- new_clustering_result
        .clustering_session$clustering_result <- clustering_result
        
        mfi_analysis <- analyze_cluster_mfi(clustering_result$clustered_data, marker_names)
        .clustering_session$mfi_analysis <- mfi_analysis
        
        cat("Clustering re-run complete!\n")
      }
      
    } else if (choice_num == 8) {
      # Return to main
      cat("Returning to main analysis...\n")
      break
    }
  }
  
  return(invisible(.clustering_session))
}

# ============================================================================
# INTEGRATION FUNCTION WITH EXISTING UMAP WORKFLOW
# ============================================================================

add_clustering_to_umap_session <- function() {
  
  if (length(.umap_session_plots) == 0) {
    cat("No active UMAP session found. Please run UMAP analysis first.\n")
    return(NULL)
  }
  
  cat("\n=== Adding Clustering Analysis to UMAP Session ===\n")
  
  # Create mock UMAP results from session data
  umap_results <- list(
    final_data = .umap_session_plots$umap_data,
    umap_result = list(marker_names = .umap_session_plots$marker_names)
  )
  
  # Perform clustering analysis
  clustering_results <- perform_clustering_analysis(umap_results)
  
  # Add cluster information back to UMAP session
  if (!is.null(clustering_results$clustering_result)) {
    .umap_session_plots$clustered_data <- clustering_results$clustering_result$clustered_data
    .umap_session_plots$clustering_results <- clustering_results
    
    cat("\nClustering results added to UMAP session!\n")
    cat("You can now create cluster-colored UMAP plots in the visualization menu.\n")
  }
  
  return(clustering_results)
}

# ============================================================================
# UTILITY AND HELPER FUNCTIONS
# ============================================================================

get_clustering_session_summary <- function() {
  if (length(.clustering_session) == 0) {
    cat("No active clustering session\n")
    return(NULL)
  }
  
  cat("=== Clustering Session Summary ===\n")
  
  if (!is.null(.clustering_session$clustering_result)) {
    result <- .clustering_session$clustering_result
    cat("Method:", result$method, "\n")
    cat("Clusters:", result$n_clusters, "\n")
    cat("Cells:", scales::comma(nrow(result$clustered_data)), "\n")
    
    if (!is.null(.clustering_session$mfi_analysis)) {
      cat("MFI analysis: Available\n")
    }
    
    if (!is.null(.clustering_session$diff_analysis)) {
      cat("Differential analysis: Available\n")
    }
    
    if (!is.null(.clustering_session$extracted_cluster)) {
      cat("Extracted cluster:", .clustering_session$extracted_cluster$cluster_name, "\n")
    }
  }
  
  return(invisible(.clustering_session))
}

restart_clustering_session <- function(clustering_results) {
  cat("Restarting clustering session with existing results...\n")
  
  if (is.null(clustering_results$clustering_result)) {
    stop("Invalid clustering results object")
  }
  
  .clustering_session$clustering_result <- clustering_results$clustering_result
  .clustering_session$mfi_analysis <- clustering_results$mfi_analysis
  .clustering_session$method <- clustering_results$method
  .clustering_session$umap_data <- clustering_results$umap_data
  .clustering_session$marker_names <- clustering_results$marker_names
  
  interactive_clustering_menu(
    clustering_results$clustering_result,
    clustering_results$mfi_analysis,
    clustering_results$marker_names
  )
  
  return(invisible(.clustering_session))
}
