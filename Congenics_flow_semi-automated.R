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

# ============================================================================
# SECTION 1: CORE SETUP AND DATA LOADING (test123)
# ============================================================================

setup_flowjo_workspace <- function(xml_path, fcs_path, keywords = c("pairing_factor", "tissue_factor")) {
  ws <- open_flowjo_xml(xml_path)
  gs <- flowjo_to_gatingset(
    ws, 
    keywords = keywords,
    keywords.source = "FCS", 
    path = fcs_path, 
    extend_val = -10000
  )
  return(gs)
}

# ============================================================================
# SECTION 2: INTERACTIVE KEYWORD SELECTION FOR FLOWJO WORKSPACE SETUP
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
    filter(non_empty_values == total_files, unique_values > 1) %>%
    arrange(desc(unique_values))
  
  constant_keywords <- summary_df %>% 
    filter(non_empty_values == total_files, unique_values == 1)
  
  mostly_empty <- summary_df %>% 
    filter(non_empty_values < total_files)
  
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
    filter(keyword == keyword_name)
  
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
        filter(
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

# ============================================================================
# CORRECTED WORKFLOW - FACTOR DEFINITION BEFORE DATA EXTRACTION
# ============================================================================

# First, fix the extract functions to not assume standardized columns exist
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
        Count = result$count,
        ParentCount = result$parent_count
      )
    })
  })
  
  if(nrow(results) == 0) {
    stop("No counts extracted. Check that nodes exist in the gating set.")
  }
  
  # Get ALL available sample metadata - don't filter by keywords yet
  pd <- pData(gs) %>% 
    rownames_to_column("Sample")
  
  # Join with all available metadata
  results <- results %>%
    left_join(pd, by = "Sample") %>%
    rename(Subpop = Count)
  
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
      filter(str_detect(colname, "(?i)comp")) %>%
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
  
  cat("=== Flow Cytometry Analysis Pipeline (Corrected) ===\n")
  
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
                                             ...) {
  
  cat("=== Enhanced Flow Analysis - Corrected Workflow ===\n")
  
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
    
    # Apply the SAME factor definitions to MFI data
    if(!is.null(results$mfi) && nrow(results$mfi) > 0 && exists("factor_info")) {
      cat("Applying same factor definitions to MFI data...\n")
      
      # Apply the same transformations that were applied to counts data
      if(factor_info$original_tissue_col %in% names(results$mfi)) {
        results$mfi <- results$mfi %>%
          mutate(tissue_factor = .data[[factor_info$original_tissue_col]])
      }
      
      if(factor_info$pairing_enabled && factor_info$original_pairing_col %in% names(results$mfi)) {
        # Handle the pairing factor creation logic
        if(factor_info$original_pairing_col == "pairing_factor") {
          # This means pairing factor was created from well ID
          if("pairing_factor" %in% names(results$mfi)) {
            results$mfi <- results$mfi %>%
              mutate(pairing_factor = str_extract(pairing_factor, "^[A-H]"))
          }
        } else {
          results$mfi <- results$mfi %>%
            mutate(pairing_factor = .data[[factor_info$original_pairing_col]])
        }
      } else {
        results$mfi <- results$mfi %>%
          mutate(pairing_factor = "no_pairing")
      }
      
      cat("✅ Applied factor definitions to MFI data\n")
    }
    
    # Store factor information in results
    results$factor_info <- factor_info
    
    cat("\n✅ Factor definition complete!\n")
    cat("All subsequent analyses will use:\n")
    cat("• tissue_factor for grouping\n")
    cat("• pairing_factor for paired analyses\n")
  }
  
  return(results)
}

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
    cat("
=== Parent Node Selection ===
")
    cat("Selected nodes:
")
    iwalk(nodes, ~cat(sprintf("  %d. %s
", .y, basename(.x))))
    
    cat("
Available methods:
")
    cat("1. Auto (hierarchical parent)
")
    cat("2. Same parent for all nodes
")
    cat("3. Individual parent for each node
")
    cat("4. Congenic-specific (CD45.1, CD45.2, CD90.1, CD90.2, CD45.1.2)
")
    cat("5. Back to node selection
")
    
    method <- readline("Choose method (1-5): ")
    
    if (method == "5") {
      return("BACK")
    }
    
    result <- switch(method,
                     "1" = {
                       # Auto hierarchical
                       cat("
Using automatic hierarchical parent detection...
")
                       mapping <- set_names(map_chr(nodes, get_hierarchical_parent), nodes)
                       list(mapping = mapping, method = "auto")
                     },
                     
                     "2" = {
                       # Same parent for all
                       while (TRUE) {
                         cat("
Select common parent for all nodes:
")
                         all_pops <- get_all_populations(gs)
                         parent <- select.list(all_pops, multiple = FALSE, title = "Select common parent (Cancel to go back)")
                         
                         if (length(parent) == 0) {
                           cat("No parent selected.
")
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
                           cat(sprintf("
--- Node %d/%d: %s ---
", i, length(nodes), basename(node)))
                           
                           # Filter to likely parents (hierarchically related)
                           likely_parents <- all_pops[
                             map_lgl(all_pops, ~ str_detect(node, fixed(.x))) &
                               nchar(all_pops) < nchar(node) &
                               all_pops != node
                           ]
                           
                           if (length(likely_parents) > 0) {
                             cat("Suggested parents:
")
                             iwalk(likely_parents, ~cat(sprintf("%d. %s
", .y, .x)))
                             cat(sprintf("%d. Choose from all populations
", length(likely_parents) + 1))
                             cat(sprintf("%d. Back to previous node
", length(likely_parents) + 2))
                             
                             choice <- readline("Select parent: ")
                             
                             if (choice == as.character(length(likely_parents) + 2)) {
                               if (i == 1) {
                                 cat("Already at first node. Going back to method selection.
")
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
                                 cat("Invalid choice. Please try again.
")
                               }
                             }
                           } else {
                             cat("No hierarchically related parents found.
")
                             cat("1. Choose from all populations
")
                             cat("2. Back to previous node
")
                             
                             choice <- readline("Select option: ")
                             if (choice == "2") {
                               if (i == 1) {
                                 cat("Already at first node. Going back to method selection.
")
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
                       cat("
Using congenic-specific mapping...
")
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
                       cat("Invalid choice. Please select 1-5.
")
                       NULL
                     }
    )
    
    if (!is.null(result)) {
      if (is.character(result) && length(result) == 1 && result == "BACK_TO_METHOD") next
      if (is.list(result)) return(result$mapping)
      if (is.character(result) && length(result) > 1) return(result)  # For congenic mapping
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
    filter(str_detect(colname, "(?i)comp")) %>%
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
               comp_lookup <- lookup %>% filter(colname %in% comp_channels)
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
      filter(str_detect(colname, "(?i)comp")) %>%
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
    
    # Show unique values for potential tissue factors
    cat("Column preview with unique values:\n")
    for(i in seq_along(available_cols)) {
      col_name <- available_cols[i]
      unique_vals <- unique(df[[col_name]])
      
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
    
    # Show tissue factor summary
    tissue_summary <- df %>%
      count(.data[[selected_tissue_factor]], name = "n_samples") %>%
      arrange(desc(n_samples))
    
    cat(sprintf("\n✅ Selected tissue factor: %s\n", selected_tissue_factor))
    cat("Groups found:\n")
    print(tissue_summary)
    
    confirm <- readline("Confirm this tissue factor? (y/n): ")
    if(tolower(confirm) == "y") {
      tissue_factor_col <- selected_tissue_factor
    } else {
      cat("Selection cancelled. Choose again.\n")
    }
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
      # Use existing column
      cat("\nColumns that might contain pairing information:\n")
      
      # Show columns with reasonable number of unique values for pairing
      pairing_candidates <- character(0)
      counter <- 1
      
      for(col_name in available_cols) {
        unique_vals <- unique(df[[col_name]])
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
          pairing_factor_col <- pairing_candidates[choice_num]
        }
      } else if(pair_choice %in% available_cols) {
        pairing_factor_col <- pair_choice
      }
      
    } else if(pairing_option == "2") {
      # Extract from WELL ID
      cat("\nWELL ID Pairing Options:\n")
      cat("a. Row letter pairing (A1,A2,A3 = Mouse 1; B1,B2,B3 = Mouse 2)\n")
      cat("b. Custom pattern extraction\n")
      
      # Look for well ID columns
      wellid_cols <- available_cols[grepl("well|WellID|Well_ID|pairing_factor", available_cols, ignore.case = TRUE)]
      
      if(length(wellid_cols) > 0) {
        cat("\nFound potential well ID columns:\n")
        iwalk(wellid_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
        
        # Show examples
        for(col in wellid_cols) {
          sample_vals <- head(unique(df[[col]]), 5)
          cat(sprintf("  %s examples: %s\n", col, paste(sample_vals, collapse = ", ")))
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
        
        if(!is.null(selected_wellid)) {
          cat("\nPairing method:\n")
          cat("a. Row letter (A1,A2→A; B1,B2→B)\n")
          cat("b. Custom pattern\n")
          
          method_choice <- readline("Choose method (a/b): ")
          
          if(tolower(method_choice) == "a") {
            # Create row letter pairing
            pairing_factor_col <- "pairing_factor_created"  # Signal that we'll create this
            df <- df %>%
              mutate(pairing_factor_created = str_extract(.data[[selected_wellid]], "^[A-H]"))
            
            # Show pairing preview
            pairing_preview <- df %>%
              count(pairing_factor_created, name = "n_samples") %>%
              arrange(pairing_factor_created)
            
            cat("\nPairing groups created from row letters:\n")
            print(pairing_preview)
            
          } else if(tolower(method_choice) == "b") {
            pattern <- readline("Enter regex pattern to extract pairing ID: ")
            if(pattern != "") {
              pairing_factor_col <- "pairing_factor_created"
              df <- df %>%
                mutate(pairing_factor_created = str_extract(.data[[selected_wellid]], pattern))
              
              # Show results
              pairing_preview <- df %>%
                count(pairing_factor_created, name = "n_samples") %>%
                arrange(pairing_factor_created)
              
              cat("\nPairing groups created from pattern:\n")
              print(pairing_preview)
            }
          }
        }
        
      } else {
        cat("No well ID columns found. Please choose option 1 or 3.\n")
        next
      }
      
    } else if(pairing_option == "3") {
      # Custom pairing
      cat("\nCustom Pairing Creation:\n")
      cat("This will create pairing groups based on a pattern or rule you define.\n")
      
      cat("Available columns:\n")
      iwalk(available_cols, ~cat(sprintf("%d. %s\n", .y, .x)))
      
      base_col_choice <- readline("Enter base column number or name: ")
      
      # Parse base column choice
      base_col <- NULL
      if(grepl("^\\d+$", base_col_choice)) {
        choice_num <- as.numeric(base_col_choice)
        if(!is.na(choice_num) && choice_num >= 1 && choice_num <= length(available_cols)) {
          base_col <- available_cols[choice_num]
        }
      } else if(base_col_choice %in% available_cols) {
        base_col <- base_col_choice
      }
      
      if(!is.null(base_col)) {
        cat("\nExamples from", base_col, ":\n")
        sample_vals <- head(unique(df[[base_col]]), 10)
        iwalk(sample_vals, ~cat(sprintf("%d. %s\n", .y, .x)))
        
        cat("\nTransformation options:\n")
        cat("1. Extract pattern with regex\n")
        cat("2. Remove suffix/prefix\n")
        cat("3. Use first N characters\n")
        cat("4. Use as-is\n")
        
        transform_choice <- readline("Choose transformation (1-4): ")
        
        if(transform_choice == "1") {
          pattern <- readline("Enter regex pattern to extract: ")
          if(pattern != "") {
            pairing_factor_col <- "pairing_factor_created"
            df <- df %>%
              mutate(pairing_factor_created = str_extract(.data[[base_col]], pattern))
          }
        } else if(transform_choice == "2") {
          cat("Remove:\n1. Prefix\n2. Suffix\n")
          remove_choice <- readline("Choose (1/2): ")
          pattern_to_remove <- readline("Enter pattern to remove: ")
          
          if(pattern_to_remove != "") {
            pairing_factor_col <- "pairing_factor_created"
            if(remove_choice == "1") {
              df <- df %>%
                mutate(pairing_factor_created = str_remove(.data[[base_col]], paste0("^", pattern_to_remove)))
            } else {
              df <- df %>%
                mutate(pairing_factor_created = str_remove(.data[[base_col]], paste0(pattern_to_remove, "$")))
            }
          }
        } else if(transform_choice == "3") {
          n_chars <- readline("Enter number of characters to keep: ")
          if(grepl("^\\d+$", n_chars)) {
            n <- as.numeric(n_chars)
            pairing_factor_col <- "pairing_factor_created"
            df <- df %>%
              mutate(pairing_factor_created = str_sub(.data[[base_col]], 1, n))
          }
        } else if(transform_choice == "4") {
          pairing_factor_col <- base_col
        }
        
        # Show results if we created a new column
        if(pairing_factor_col == "pairing_factor_created") {
          pairing_preview <- df %>%
            count(pairing_factor_created, name = "n_samples") %>%
            arrange(pairing_factor_created)
          
          cat("\nPairing groups created:\n")
          print(pairing_preview)
        }
      }
      
    } else if(pairing_option == "4") {
      # No pairing
      pairing_factor_col <- "none"
      cat("✅ No pairing selected. Only unpaired analyses will be available.\n")
    }
    
    if(is.null(pairing_factor_col)) {
      cat("❌ Invalid selection. Please try again.\n")
      next
    }
    
    # Confirm pairing factor (unless "none")
    if(pairing_factor_col != "none") {
      # Show pairing summary
      if(pairing_factor_col %in% names(df)) {
        pairing_summary <- df %>%
          count(.data[[pairing_factor_col]], name = "n_samples") %>%
          arrange(desc(n_samples))
        
        cat(sprintf("\n✅ Selected pairing factor: %s\n", pairing_factor_col))
        cat("Pairing groups:\n")
        print(head(pairing_summary, 10))
        if(nrow(pairing_summary) > 10) {
          cat("... and", nrow(pairing_summary) - 10, "more groups\n")
        }
      }
      
      confirm <- readline("Confirm this pairing factor? (y/n): ")
      if(tolower(confirm) != "y") {
        pairing_factor_col <- NULL
        cat("Selection cancelled. Choose again.\n")
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
    mutate(tissue_factor = .data[[tissue_factor_col]])
  
  cat("✅ Created 'tissue_factor' column based on:", tissue_factor_col, "\n")
  
  # Create pairing_factor column if not "none"
  if(pairing_factor_col != "none") {
    df <- df %>%
      mutate(pairing_factor = .data[[pairing_factor_col]])
    cat("✅ Created 'pairing_factor' column based on:", pairing_factor_col, "\n")
  } else {
    df <- df %>%
      mutate(pairing_factor = "no_pairing")
    cat("✅ Created 'pairing_factor' column with value: no_pairing\n")
  }
  
  # Final summary
  cat("\n=== FACTOR DEFINITION COMPLETE ===\n")
  cat("Analysis factors created:\n")
  
  tissue_final <- df %>% count(tissue_factor, name = "n") %>% arrange(desc(n))
  cat("\nTissue Factor Groups:\n")
  print(tissue_final)
  
  if(pairing_factor_col != "none") {
    pairing_final <- df %>% count(pairing_factor, name = "n") %>% arrange(desc(n))
    cat("\nPairing Factor Groups:\n")
    print(head(pairing_final, 15))
    if(nrow(pairing_final) > 15) {
      cat("... and", nrow(pairing_final) - 15, "more pairing groups\n")
    }
  } else {
    cat("\nPairing Factor: No pairing (unpaired analyses only)\n")
  }
  
  return(list(
    data = df,
    tissue_factor = "tissue_factor",
    pairing_factor = "pairing_factor",
    original_tissue_col = tissue_factor_col,
    original_pairing_col = if(pairing_factor_col != "none") pairing_factor_col else NULL,
    pairing_enabled = pairing_factor_col != "none"
  ))
}

#' Enhanced metadata assignment with factor definition
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
data_clean_custom <- function(data) {
  # Function to handle unstained sample detection and removal
  handle_unstained_samples <- function(df, data_type = "data") {
    if (!"Sample" %in% names(df)) {
      return(df)
    }
    
    unstained_patterns <- c("unstained", "no stain", "no_stain")
    pattern <- paste(unstained_patterns, collapse = "|")
    
    unstained_rows <- df %>%
      mutate(row_id = row_number()) %>%
      filter(str_detect(tolower(Sample), pattern)) %>%
      pull(row_id)
    
    if (length(unstained_rows) == 0) {
      return(df)
    }
    
    # Only show interactive prompt for the first dataset
    if (!exists(".unstained_decision", envir = .GlobalEnv)) {
      cat("\n=== Unstained Samples Detected in", data_type, "===\n")
      cat(paste("Found", length(unstained_rows), "unstained sample(s)\n"))
      
      unstained_data <- df %>% slice(unstained_rows)
      if ("NodeShort" %in% names(df)) {
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
      .unstained_decision <<- if (user_choice == 1) "remove" else "keep"
      
      if (.unstained_decision == "remove") {
        cat("Will remove unstained samples from all datasets\n")
      } else {
        cat("Will keep unstained samples in all datasets\n")
      }
    }
    
    # Apply the stored decision
    if (.unstained_decision == "remove") {
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
    
    if ("NodeShort" %in% names(df)) {
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
  if (exists(".unstained_decision", envir = .GlobalEnv)) {
    rm(.unstained_decision, envir = .GlobalEnv)
  }
  
  # Check if input is a list
  if (is.list(data) && !is.data.frame(data)) {
    # Apply cleaning function to each element in the list
    cleaned_data <- imap(data, ~ {
      if (is.data.frame(.x)) {
        data_type <- if (.y == "counts") "counts data" else if (.y == "mfi") "MFI data" else .y
        clean_single_df(.x, data_type)
      } else {
        .x
      }
    })
    
    # Clean up the global decision variable
    if (exists(".unstained_decision", envir = .GlobalEnv)) {
      rm(.unstained_decision, envir = .GlobalEnv)
    }
    
    return(cleaned_data)
  } else if (is.data.frame(data)) {
    result <- clean_single_df(data, "single dataset")
    
    # Clean up the global decision variable
    if (exists(".unstained_decision", envir = .GlobalEnv)) {
      rm(.unstained_decision, envir = .GlobalEnv)
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
        filter(n_distinct(congenics) == 2) %>%
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
        filter(n_distinct(congenics) == 2) %>%
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

# Main interactive function
create_paired_comparison_plots <- function(data) {
  # Load required libraries
  require(dplyr)
  require(purrr)
  require(ggplot2)
  require(rstatix)
  require(ggpubr)
  require(scales)
  
  # Validate inputs
  if (missing(data) || !is.data.frame(data)) {
    stop("Please provide a valid data frame")
  }
  
  required_cols <- c("NodeShort", "tissue_factor", "pairing_factor", "congenics", "Freq")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
  }
  
  # Interactive test selection
  test_choice <- menu(c("Paired t-test", "Wilcoxon signed-rank test"),
                      title = "Choose statistical test:")
  test_type <- c("t_test", "wilcox_test")[test_choice]
  
  # Interactive faceting selection
  facet_choice <- menu(c("tissue_factor (tissue)", "NodeShort (subpopulation of interest)", "No faceting"),
                       title = "Choose how to facet the plots:")
  facet_var <- c("tissue_factor", "NodeShort", "none")[facet_choice]
  
  # Print information about the analysis
  test_description <- switch(test_type, 
                             't_test' = 'paired t-test', 
                             'wilcox_test' = 'Wilcoxon signed-rank test')
  facet_description <- switch(facet_var,
                              'tissue_factor' = 'tissue',
                              'NodeShort' = 'subpopulation of interest',
                              'none' = 'no faceting (separate plots)')
  
  # FIXED: Choose split variable based on what you want separate plots for
  if (facet_var == "tissue_factor") {
    # If faceting BY tissue, you want separate plots FOR each NodeShort
    split_var <- "NodeShort"
    plot_description <- "cell population(s)"
  } else if (facet_var == "NodeShort") {
    # If faceting BY NodeShort, you want separate plots FOR each tissue
    split_var <- "tissue_factor"  
    plot_description <- "tissue(s)"
  } else {
    # No faceting - create separate plots for each NodeShort-tissue combination
    split_var <- c("NodeShort", "tissue_factor")
    plot_description <- "NodeShort-tissue combination(s)"
  }
  
  cat("Creating paired comparison plots using", test_description, "\n")
  cat("Faceting by", facet_description, "\n")
  
  # Calculate number of plots that will be created
  if (facet_var == "none") {
    n_plots <- data %>% 
      distinct(NodeShort, tissue_factor) %>% 
      nrow()
    cat("Processing", n_plots, plot_description, "\n\n")
  } else {
    n_plots <- n_distinct(data[[split_var]])
    cat("Processing", n_plots, plot_description, "\n\n")
  }
  
  # FIXED: Create plots based on the correct splitting logic
  if (facet_var == "tissue_factor") {
    # Split by NodeShort to get separate plots for each cell population
    # Each plot will be faceted by tissue (tissue_factor)
    plots <- data %>%
      group_split(NodeShort, .keep = TRUE) %>%
      map(~ create_subgroup_plot(.x, test_type = test_type, facet_var = "tissue_factor")) %>%
      set_names(map_chr(data %>% distinct(NodeShort) %>% pull(NodeShort), as.character))
  } else if (facet_var == "NodeShort") {
    # Split by tissue_factor to get separate plots for each tissue
    # Each plot will be faceted by cell population (NodeShort)  
    plots <- data %>%
      group_split(tissue_factor, .keep = TRUE) %>%
      map(~ create_subgroup_plot(.x, test_type = test_type, facet_var = "NodeShort")) %>%
      set_names(map_chr(data %>% distinct(tissue_factor) %>% pull(tissue_factor), as.character))
  } else {
    # No faceting - create separate plots for each NodeShort-tissue combination
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
    filter(!is.na(congenics)) %>%
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
      filter(congenics == congenic) %>%
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
        filter(!is.na(congenics)) %>%
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
        filter(!is.na(congenics)) %>%
        count(tissue_factor, congenics) %>%
        arrange(tissue_factor, congenics)
      
      cat("Available tissues:", paste(available_groups, collapse = ", "), "\n")
      print(distribution)
      
      marker_dist <- mfi_data %>%
        filter(!is.na(congenics)) %>%
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
    filter(!is.na(congenics)) %>%
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
      filter(marker == .env$marker, !is.na(congenics)) %>%
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
        filter(!is.na(congenics)) %>%
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
    cat("Common pairing variables: Sample, NodeShort, Subject_ID, etc.\n")
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
      filter(!is.na(congenics), !is.na(tissue_factor)) %>%
      count(tissue_factor, congenics) %>%
      pivot_wider(names_from = congenics, values_from = n, values_fill = 0)
    
    print(sample_dist)
  }
}

show_data_distribution <- function(mfi_data, selected_congenics, selected_markers) {
  cat("\n=== DATA DISTRIBUTION ANALYSIS ===\n")
  
  filtered_data <- mfi_data %>%
    filter(congenics %in% selected_congenics,
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
  normal_pct <- norm_summary %>% filter(normal_likely == "Yes") %>% pull(percentage)
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
    marker_data <- filtered_data %>% filter(marker == !!marker)
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
    filter(congenics %in% selected_congenics,
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
    
    tissue_data <- filtered_data %>% filter(tissue_factor == current_tissue)
    
    if (nrow(tissue_data) == 0) {
      warning(paste("No data for tissue:", current_tissue))
      next
    }
    
    for (current_marker in selected_markers) {
      marker_tissue_data <- tissue_data %>% filter(marker == current_marker)
      
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
  group1_data <- marker_data %>% filter(congenics == groups[1]) %>% pull(MFI)
  group2_data <- marker_data %>% filter(congenics == groups[2]) %>% pull(MFI)
  
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
      filter(!is.na(.data[[groups[1]]]), !is.na(.data[[groups[2]]]))
    
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
  group1_data <- marker_data %>% filter(congenics == groups[1]) %>% pull(MFI)
  group2_data <- marker_data %>% filter(congenics == groups[2]) %>% pull(MFI)
  
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
    filter(!is.na(.data[[groups[1]]]), !is.na(.data[[groups[2]]]))
  
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
    filter(!is.na(MFI))
  
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
    filter(!is.na(MFI))
  
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
    filter(!is.na(MFI)) %>%
    group_by(.data[[stats_config$pairing_var]]) %>%
    filter(n_distinct(congenics) == length(stats_config$groups)) %>%
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
    filter(!is.na(MFI)) %>%
    pivot_wider(names_from = congenics, values_from = MFI) %>%
    filter(if_all(all_of(stats_config$groups), ~ !is.na(.)))
  
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
    filter(!is.na(congenics))
  
  if (!is.null(selected_congenics)) {
    mfi_clean <- mfi_clean %>%
      filter(congenics %in% selected_congenics)
    cat("Filtered to selected congenics:", paste(selected_congenics, collapse = ", "), "\n")
  }
  
  if (!is.null(selected_markers)) {
    mfi_clean <- mfi_clean %>%
      filter(marker %in% selected_markers)
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
      filter(tissue_factor == group)
    
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
    filter(tissue == current_tissue)
  
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

create_mfi_heatmaps_interactive_enhanced <- function(mfi_data, ...) {
  cat("=== Enhanced Interactive MFI Heatmap Creation with Statistics ===\n")
  
  while (TRUE) {
    # Step 1: Select congenics
    selected_congenics <- select_congenics_interactive(mfi_data)
    if (is.null(selected_congenics)) {
      cat("Congenic selection cancelled.\n")
      return(NULL)
    }
    
    # Step 2: Select grouping option
    grouping_option <- select_grouping_option(mfi_data)
    
    # Step 3: Select markers
    filtered_data <- mfi_data %>%
      filter(!is.na(congenics), congenics %in% selected_congenics)
    
    selected_markers <- select_markers_interactive(filtered_data)
    if (is.null(selected_markers)) {
      cat("Marker selection cancelled.\n")
      return(NULL)
    }
    
    # Step 4: Select statistical testing
    stats_config <- select_statistical_test_interactive(mfi_data, selected_congenics, selected_markers)
    
    # Step 5: Get scaling method
    scaling_params <- select_scaling_method_interactive()
    
    if (is.null(scaling_params)) {
      cat("Scaling method selection cancelled.\n")
      return(NULL)
    }
    
    if (is.character(scaling_params) && scaling_params == "DIAGNOSTIC") {
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
    if (stats_config$perform_stats) {
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
    
    if (length(heatmaps) > 0) {
      cat(sprintf("\nDisplaying heatmap: %s\n", names(heatmaps)[1]))
      draw(heatmaps[[1]])
      return(heatmaps)
    } else {
      cat("No heatmaps were created. Please check your selections.\n")
      return(NULL)
    }
  }
}

# ===== UTILITY FUNCTIONS =====

# Export statistics results to data frame
export_stats_results <- function(stats_results, file_path = NULL) {
  if (is.null(stats_results) || is.null(stats_results$summary)) {
    warning("No statistical results to export")
    return(NULL)
  }
  
  export_df <- stats_results$summary %>%
    mutate(
      test_type = stats_results$config$test_type,
      alpha_level = stats_results$config$alpha,
      significance_display = stats_results$config$sig_display
    ) %>%
    select(tissue, marker, test_type, p_value, p_formatted, stars, significant, 
           alpha_level, significance_display, test_name) %>%
    arrange(tissue, marker)
  
  if (!is.null(file_path)) {
    write_csv(export_df, file_path)
    cat("Statistical results exported to:", file_path, "\n")
    cat("Results include", nrow(export_df), "tissue-marker combinations\n")
  }
  
  return(export_df)
}

# Quick statistical testing function for existing heatmaps
add_stats_to_existing_heatmap <- function(mfi_data, selected_congenics, selected_markers) {
  cat("=== Adding Statistical Analysis to Existing Selection ===\n")
  
  stats_config <- select_statistical_test_interactive(mfi_data, selected_congenics, selected_markers)
  
  if (stats_config$perform_stats) {
    filtered_data <- mfi_data %>%
      filter(!is.na(congenics), 
             congenics %in% selected_congenics,
             marker %in% selected_markers)
    
    stats_results <- perform_statistical_analysis(filtered_data, selected_congenics, 
                                                  selected_markers, stats_config)
    
    if (!is.null(stats_results)) {
      display_stats_summary(stats_results)
      return(stats_results)
    }
  }
  
  return(NULL)
}

# ===== USAGE DOCUMENTATION =====

cat("
=== ENHANCED MFI HEATMAP CODE WITH STATISTICAL TESTING - CLEANED VERSION ===

MAIN FUNCTION:
create_mfi_heatmaps_interactive_enhanced(mfi_data)
- Complete interactive workflow with statistics
- Streamlined for optimal performance

UTILITY FUNCTIONS:
1. add_stats_to_existing_heatmap(mfi_data, congenics, markers)
   - Add statistics to specific selections
   
2. export_stats_results(stats_results, 'file.csv')
   - Export statistical results to CSV

EXAMPLE USAGE:
# Interactive mode (recommended)
heatmaps <- create_mfi_heatmaps_interactive_enhanced(mfi_data)

# Add stats to existing selection
stats <- add_stats_to_existing_heatmap(mfi_data, 
                                      c('CD45.1', 'CD45.2'), 
                                      c('CD3', 'CD4', 'CD8'))

# Export results
export_df <- export_stats_results(stats, 'my_stats_results.csv')

FEATURES:
✓ Interactive congenic, marker, and test selection
✓ Data distribution analysis and test recommendations  
✓ Multiple scaling methods (raw, Z-score, log, etc.)
✓ Statistical tests: t-test, ANOVA, non-parametric options
✓ Significance annotations on heatmaps
✓ Comprehensive results summaries
✓ Export functionality
✓ Tidy code style with reduced redundancy

")
#===============================================================================
# Engraftment plot of congenic markers
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
                                    interactive_stats = TRUE) {
  
  # Check if required columns exist
  required_cols <- c(wellid_col, group_col, marker_col, freq_col)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", "), 
               "\nAvailable columns:", paste(names(data), collapse = ", ")))
  }
  
  # Interactive statistical method selection
  if (interactive_stats && add_stats) {
    cat("\n=== Statistical Analysis Options ===\n")
    cat("1: Paired t-test (recommended for engraftment data)\n")
    cat("2: Unpaired t-test\n") 
    cat("3: ANOVA with post-hoc comparisons\n")
    cat("4: No statistical testing\n")
    
    stat_choice <- as.integer(readline(prompt = "Select statistical method (enter number): "))
    
    if (is.na(stat_choice) || stat_choice < 1 || stat_choice > 4) {
      stop("Invalid selection. Please run the function again and choose a valid number.")
    }
    
    if (stat_choice == 1) {
      stat_method <- "paired_t_test"
      cat("Selected: Paired t-test\n")
    } else if (stat_choice == 2) {
      stat_method <- "t_test"
      cat("Selected: Unpaired t-test\n")
    } else if (stat_choice == 3) {
      stat_method <- "anova"
      cat("Selected: ANOVA with post-hoc\n")
    } else {
      add_stats <- FALSE
      cat("Selected: No statistical testing\n")
    }
    
    # If doing comparisons, ask about control group
    if (add_stats && stat_choice %in% c(1, 2, 3)) {
      cat("\n=== Comparison Options ===\n")
      cat("1: All pairwise comparisons\n")
      cat("2: Compare all groups to a control group\n")
      
      comp_choice <- as.integer(readline(prompt = "Select comparison type (enter number): "))
      
      if (is.na(comp_choice) || comp_choice < 1 || comp_choice > 2) {
        comparison_type <- "pairwise"
        cat("Invalid selection - defaulting to pairwise comparisons\n")
      } else if (comp_choice == 1) {
        comparison_type <- "pairwise"
        cat("Selected: All pairwise comparisons\n")
      } else {
        comparison_type <- "vs_control"
        
        # Select control group
        available_groups <- unique(data[[group_col]])
        cat("\nAvailable groups for control:\n")
        for (i in seq_along(available_groups)) {
          cat(paste(i, ":", available_groups[i], "\n"))
        }
        
        control_choice <- as.integer(readline(prompt = "Select control group (enter number): "))
        
        if (is.na(control_choice) || control_choice < 1 || control_choice > length(available_groups)) {
          control_group <- available_groups[1]
          cat(paste("Invalid selection - defaulting to:", control_group, "\n"))
        } else {
          control_group <- available_groups[control_choice]
          cat(paste("Selected control group:", control_group, "\n"))
        }
      }
    }
  } else {
    # Non-interactive mode
    comparison_type <- "pairwise"
    control_group <- NULL
  }
  # Interactive selection of normalization tissue
  normalization_tissue <- NULL
  use_normalization <- FALSE
  
  if (interactive) {
    # Show available tissue types
    available_tissues <- unique(data[[group_col]])
    cat("\nAvailable tissues for normalization:\n")
    for (i in seq_along(available_tissues)) {
      cat(paste(i, ":", available_tissues[i], "\n"))
    }
    cat(paste(length(available_tissues) + 1, ": No normalization (use raw frequencies)\n"))
    
    # Get user choice
    choice <- as.integer(readline(prompt = "Select normalization tissue (enter number): "))
    
    if (is.na(choice) || choice < 1 || choice > length(available_tissues) + 1) {
      stop("Invalid selection. Please run the function again and choose a valid number.")
    }
    
    if (choice <= length(available_tissues)) {
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
    if (!is.null(spleen_group) && spleen_group %in% unique(data[[group_col]])) {
      use_normalization <- TRUE
      normalization_tissue <- spleen_group
    } else {
      use_normalization <- FALSE
      normalization_tissue <- NULL
    }
  }
  
  # Process the data to calculate engraftment ratios
  congenic_engraftment <- data %>%
    select(all_of(c(wellid_col, group_col, marker_col, freq_col))) %>%
    # Extract the first letter from pairing_factor to create matching groups
    mutate(
      well_letter = substr(.data[[wellid_col]], 1, 1)
    )
  
  # Apply normalization if selected
  if (use_normalization) {
    congenic_engraftment <- congenic_engraftment %>%
      # For each marker and well letter, get the normalization tissue frequency
      group_by(across(all_of(c(marker_col, "well_letter")))) %>%
      mutate(
        norm_freq = .data[[freq_col]][.data[[group_col]] == normalization_tissue],
        normalized_freq = .data[[freq_col]] / norm_freq
      ) %>%
      ungroup() %>%
      select(all_of(c(wellid_col, group_col, marker_col)), normalized_freq) %>%
      pivot_wider(names_from = all_of(marker_col), values_from = normalized_freq)
    
    # Update plot labels for normalized data
    y_label <- paste0("Log2(ratio KO:WT normalized to ", normalization_tissue, ")")
    caption_text <- paste0("*Normalized to paired ", normalization_tissue, " controls")
    
  } else {
    congenic_engraftment <- congenic_engraftment %>%
      ungroup() %>%
      select(all_of(c(wellid_col, group_col, marker_col, freq_col))) %>%
      pivot_wider(names_from = all_of(marker_col), values_from = all_of(freq_col))
    
    # Update plot labels for non-normalized data
    y_label <- "Log2(ratio KO:WT - raw frequencies)"
    caption_text <- "*Using raw frequencies (no normalization)"
  }
  
  # Calculate engraftment ratio
  congenic_engraftment <- congenic_engraftment %>%
    mutate(
      engraftment_ratio = log2(.data[[ko_marker]] / .data[[wt_marker]])
    ) %>%
    # Filter out infinite values
    filter(is.finite(engraftment_ratio)) %>%
    select(all_of(c(wellid_col, group_col)), engraftment_ratio)
  
  # Check if we have data after processing
  if (nrow(congenic_engraftment) == 0) {
    stop("No valid engraftment ratios calculated. Check your marker names and data.")
  }
  
  # Calculate summary statistics
  summary_stats <- congenic_engraftment %>%
    group_by(across(all_of(group_col))) %>%
    summarise(
      mean_val = mean(engraftment_ratio, na.rm = TRUE),
      sd_val = sd(engraftment_ratio, na.rm = TRUE),
      n = n(),
      sem_val = sd_val / sqrt(n),
      .groups = 'drop'
    )
  
  # Rename the grouping column for consistency
  names(summary_stats)[1] <- group_col
  
  # Create factor levels with normalization group first, then order by mean if requested
  all_groups <- unique(summary_stats[[group_col]])
  
  if (use_normalization && normalization_tissue %in% all_groups) {
    # Put normalization tissue first
    other_groups <- setdiff(all_groups, normalization_tissue)
    
    if (order_by_mean) {
      # Order other groups by mean (excluding normalization tissue)
      other_stats <- summary_stats %>% 
        filter(.data[[group_col]] != normalization_tissue)
      other_groups_ordered <- other_stats[[group_col]][order(-other_stats$mean_val)]
      factor_levels <- c(normalization_tissue, other_groups_ordered)
    } else {
      factor_levels <- c(normalization_tissue, other_groups)
    }
  } else {
    # No normalization, just order by mean if requested
    if (order_by_mean) {
      factor_levels <- summary_stats[[group_col]][order(-summary_stats$mean_val)]
    } else {
      factor_levels <- all_groups
    }
  }
  
  # Apply factor levels to both summary stats and original data
  summary_stats <- summary_stats %>%
    mutate(!!sym(group_col) := factor(.data[[group_col]], levels = factor_levels))
  
  congenic_engraftment <- congenic_engraftment %>%
    mutate(!!sym(group_col) := factor(.data[[group_col]], levels = factor_levels))
  
  # Calculate error bar positions based on bar direction
  summary_stats <- summary_stats %>%
    mutate(
      error_ymin = ifelse(mean_val >= 0, mean_val, mean_val - sem_val),
      error_ymax = ifelse(mean_val >= 0, mean_val + sem_val, mean_val)
    )
  
  # Perform statistical tests if requested
  stat_results <- NULL
  if (add_stats && nrow(summary_stats) > 1) {
    
    # Perform pairwise comparisons
    # Create formula dynamically
    formula_str <- paste("engraftment_ratio ~", group_col)
    test_formula <- as.formula(formula_str)
    
    if (stat_method == "paired_t_test") {
      # For paired t-test, use the existing congenic_engraftment data
      # but add the well_letter for pairing
      congenic_for_stats <- data %>%
        select(all_of(c(wellid_col, group_col, marker_col, freq_col))) %>%
        mutate(well_letter = substr(.data[[wellid_col]], 1, 1))
      
      # Apply the same processing as the main data
      if (use_normalization) {
        congenic_for_stats <- congenic_for_stats %>%
          group_by(across(all_of(c(marker_col, "well_letter")))) %>%
          mutate(
            norm_freq = .data[[freq_col]][.data[[group_col]] == normalization_tissue],
            normalized_freq = .data[[freq_col]] / norm_freq
          ) %>%
          ungroup() %>%
          select(all_of(c(wellid_col, group_col, marker_col, "well_letter")), normalized_freq) %>%
          pivot_wider(names_from = all_of(marker_col), values_from = normalized_freq)
      } else {
        congenic_for_stats <- congenic_for_stats %>%
          select(all_of(c(wellid_col, group_col, marker_col, freq_col, "well_letter"))) %>%
          pivot_wider(names_from = all_of(marker_col), values_from = all_of(freq_col))
      }
      
      congenic_for_stats <- congenic_for_stats %>%
        mutate(
          engraftment_ratio = log2(.data[[ko_marker]] / .data[[wt_marker]])
        ) %>%
        filter(is.finite(engraftment_ratio)) %>%
        select(all_of(c(wellid_col, group_col, "well_letter")), engraftment_ratio) %>%
        mutate(!!sym(group_col) := factor(.data[[group_col]], levels = factor_levels))
      
      # Determine comparisons based on user choice
      all_groups <- levels(congenic_for_stats[[group_col]])
      
      if (comparison_type == "vs_control" && !is.null(control_group)) {
        # Only compare each group to control
        other_groups <- setdiff(all_groups, control_group)
        comparisons <- lapply(other_groups, function(g) c(control_group, g))
      } else {
        # All pairwise comparisons
        comparisons <- combn(all_groups, 2, simplify = FALSE)
      }
      
      # Perform paired t-test for each comparison
      stat_list <- list()
      for (i in seq_along(comparisons)) {
        group1_name <- comparisons[[i]][1]
        group2_name <- comparisons[[i]][2]
        
        # Get data for both groups, ensuring we have paired data
        paired_data <- congenic_for_stats %>%
          filter(.data[[group_col]] %in% c(group1_name, group2_name)) %>%
          select(all_of(c(group_col, "well_letter")), engraftment_ratio) %>%
          pivot_wider(names_from = all_of(group_col), values_from = engraftment_ratio) %>%
          filter(complete.cases(.))
        
        if (nrow(paired_data) >= 3) {  # Need at least 3 pairs for meaningful test
          # Perform paired t-test
          test_result <- t.test(paired_data[[group1_name]], 
                                paired_data[[group2_name]], 
                                paired = TRUE)
          
          stat_list[[i]] <- data.frame(
            group1 = group1_name,
            group2 = group2_name,
            n1 = nrow(paired_data),
            n2 = nrow(paired_data),
            statistic = test_result$statistic,
            p = test_result$p.value,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if (length(stat_list) > 0) {
        stat_results <- bind_rows(stat_list) %>%
          mutate(
            p.adj = p.adjust(p, method = "bonferroni"),
            p.adj.signif = case_when(
              p.adj <= 0.001 ~ "***",
              p.adj <= 0.01 ~ "**", 
              p.adj <= 0.05 ~ "*",
              TRUE ~ "ns"
            )
          )
      } else {
        stat_results <- NULL
      }
      
    } else if (stat_method == "t_test") {
      # Determine comparisons
      all_groups <- levels(congenic_engraftment[[group_col]])
      
      if (comparison_type == "vs_control" && !is.null(control_group)) {
        # Filter to only control vs others
        stat_results <- congenic_engraftment %>%
          pairwise_t_test(test_formula, 
                          ref.group = control_group,
                          p.adjust.method = "bonferroni") %>%
          add_significance("p.adj")
      } else {
        # All pairwise
        stat_results <- congenic_engraftment %>%
          pairwise_t_test(test_formula, 
                          p.adjust.method = "bonferroni") %>%
          add_significance("p.adj")
      }
      
    } else if (stat_method == "wilcox_test") {
      # Determine comparisons
      if (comparison_type == "vs_control" && !is.null(control_group)) {
        stat_results <- congenic_engraftment %>%
          pairwise_wilcox_test(test_formula, 
                               ref.group = control_group,
                               p.adjust.method = "bonferroni") %>%
          add_significance("p.adj")
      } else {
        stat_results <- congenic_engraftment %>%
          pairwise_wilcox_test(test_formula, 
                               p.adjust.method = "bonferroni") %>%
          add_significance("p.adj")
      }
      
    } else if (stat_method == "anova") {
      # Perform ANOVA first
      anova_result <- congenic_engraftment %>%
        anova_test(test_formula)
      
      cat("\n=== ANOVA Results ===\n")
      print(anova_result)
      
      # If significant, proceed with post-hoc
      if (anova_result$p < 0.05) {
        if (comparison_type == "vs_control" && !is.null(control_group)) {
          # Dunnett's test for multiple comparisons with control
          stat_results <- congenic_engraftment %>%
            dunn_test(test_formula, p.adjust.method = "bonferroni") %>%
            filter(group1 == control_group | group2 == control_group) %>%
            add_significance("p.adj")
        } else {
          # Tukey's HSD for all pairwise comparisons
          stat_results <- congenic_engraftment %>%
            tukey_hsd(test_formula) %>%
            add_significance("p.adj")
        }
      } else {
        stat_results <- NULL
        cat("ANOVA not significant - no post-hoc testing performed\n")
      }
    }
    
    # Filter out non-significant results (p > 0.05)
    if (!is.null(stat_results)) {
      stat_results <- stat_results %>%
        filter(p.adj <= 0.05) %>%
        mutate(
          group1 = factor(group1, levels = factor_levels),
          group2 = factor(group2, levels = factor_levels)
        )
    }
    
    # Add y positions for stat annotations with proper spacing
    if (!is.null(stat_results) && nrow(stat_results) > 0) {
      y_max <- max(summary_stats$error_ymax, na.rm = TRUE)
      y_min <- min(summary_stats$error_ymin, na.rm = TRUE)
      y_range <- y_max - y_min
      
      # Sort by x-distance to minimize overlap
      stat_results <- stat_results %>%
        mutate(
          x1_pos = as.numeric(group1),
          x2_pos = as.numeric(group2),
          x_distance = abs(x2_pos - x1_pos)
        ) %>%
        arrange(x_distance, x1_pos)
      
      # Calculate y positions with adequate spacing and proper text placement
      base_height <- y_max + y_range * 0.08
      spacing <- y_range * 0.15
      
      stat_results <- stat_results %>%
        mutate(
          layer = row_number(),
          y.position = base_height + (layer - 1) * spacing,
          p_label = case_when(
            p.adj < 0.001 ~ paste0("p=", format(p.adj, scientific = TRUE, digits = 2)),
            p.adj < 0.01 ~ paste0("p=", format(round(p.adj, 4), nsmall = 4)),
            TRUE ~ paste0("p=", format(round(p.adj, 3), nsmall = 3))
          )
        )
    }
  }
  
  # Create base plot using aes() instead of deprecated aes_string()
  p <- ggplot(summary_stats, aes(x = .data[[group_col]], y = mean_val)) +
    geom_col(aes(fill = .data[[group_col]]), alpha = 0.8, 
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = error_ymin, ymax = error_ymax), 
                  width = 0.3, color = "black", linewidth = 0.3) +
    # Always add individual data points
    # Always add individual data points with jitter
    geom_jitter(data = congenic_engraftment, 
                aes(x = .data[[group_col]], y = engraftment_ratio), 
                size = 2, alpha = 0.7, width = 0.2, height = 0,
                inherit.aes = FALSE) +
    labs(
      title = plot_title,
      x = x_label,
      y = y_label,
      caption = caption_text
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5)
  
  # Add statistical annotations if available
  if (add_stats && !is.null(stat_results) && nrow(stat_results) > 0) {
    # Add comparison brackets and p-values with proper spacing
    for (i in 1:nrow(stat_results)) {
      row <- stat_results[i, ]
      
      # Get x positions for the groups
      x1 <- as.numeric(row$group1)
      x2 <- as.numeric(row$group2)
      
      # Calculate bracket and text positioning
      bracket_height <- row$y.position
      text_height <- bracket_height + max(summary_stats$error_ymax) * 0.04
      bracket_offset <- max(summary_stats$error_ymax) * 0.03
      
      # Add bracket lines
      p <- p + 
        annotate("segment", 
                 x = x1, xend = x1, 
                 y = bracket_height - bracket_offset, yend = bracket_height,
                 color = "black", linewidth = 0.4) +
        annotate("segment", 
                 x = x1, xend = x2, 
                 y = bracket_height, yend = bracket_height,
                 color = "black", linewidth = 0.4) +
        annotate("segment", 
                 x = x2, xend = x2, 
                 y = bracket_height, yend = bracket_height - bracket_offset,
                 color = "black", linewidth = 0.4) +
        annotate("text", 
                 x = (x1 + x2) / 2, y = text_height,
                 label = row$p_label, size = 2.8, hjust = 0.5, vjust = 0)
    }
  }
  
  # Return plot and optionally statistical results
  result <- list(plot = p)
  
  if (add_stats && !is.null(stat_results)) {
    result$statistics <- stat_results
    result$summary_stats <- summary_stats
  }
  
  # If only plot requested, return just the plot
  if (!add_stats) {
    return(p)
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
            filter(n_distinct(congenics) == 2) %>%
            pairwise_t_test(as.formula(paste("Freq ~", "congenics")), 
                            paired = TRUE, 
                            p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else if(test_type == "wilcox_test" && pairing_enabled) {
          stat_test_results <- df %>%
            group_by(NodeShort) %>%
            filter(n_distinct(congenics) == 2) %>%
            pairwise_wilcox_test(as.formula(paste("Freq ~", "congenics")), 
                                 paired = TRUE,
                                 p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else {
          # Unpaired tests
          stat_test_results <- df %>%
            group_by(NodeShort) %>%
            filter(n_distinct(congenics) == 2) %>%
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
            filter(n_distinct(congenics) == 2) %>%
            pairwise_t_test(as.formula(paste("Freq ~", "congenics")), 
                            paired = TRUE,
                            p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else if(test_type == "wilcox_test" && pairing_enabled) {
          stat_test_results <- df %>%
            group_by(.data[[tissue_col]]) %>%
            filter(n_distinct(congenics) == 2) %>%
            pairwise_wilcox_test(as.formula(paste("Freq ~", "congenics")), 
                                 paired = TRUE,
                                 p.adjust.method = "bonferroni") %>%
            add_xy_position(x = "congenics") %>%
            mutate(p.format = scales::pvalue(p.adj, accuracy = 0.001, add_p = TRUE))
        } else {
          # Unpaired tests
          stat_test_results <- df %>%
            group_by(.data[[tissue_col]]) %>%
            filter(n_distinct(congenics) == 2) %>%
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
#' Updated engraftment plot using standardized factors
create_engraftment_plot_enhanced <- function(data, 
                                             tissue_col = "tissue_factor",
                                             pairing_col = "tissue_factor",
                                             marker_col = "congenics",
                                             freq_col = "Freq",
                                             ko_marker = "CD45.1.2",
                                             wt_marker = "CD45.1",
                                             normalization_tissue = "Spleen",
                                             fill_colors = c("Spleen" = "white", "SG" = "black", "IEL" = "grey60"),
                                             plot_title = "Engraftment Ratio by Tissue",
                                             x_label = "Tissue",
                                             y_label = "Log2(ratio KO:WT)",
                                             caption_text = "*Normalized to paired controls",
                                             order_by_mean = TRUE,
                                             interactive = TRUE,
                                             add_stats = TRUE,
                                             stat_method = "paired_t_test",
                                             interactive_stats = TRUE) {
  
  # Check if required columns exist
  required_cols <- c(pairing_col, tissue_col, marker_col, freq_col)
  missing_cols <- required_cols[!required_cols %in% names(data)]
  
  if(length(missing_cols) > 0) {
    stop(paste("Missing columns:", paste(missing_cols, collapse = ", "), 
               "\nAvailable columns:", paste(names(data), collapse = ", ")))
  }
  
  # Check if pairing is enabled
  pairing_enabled <- !all(data[[pairing_col]] %in% c("no_pairing", "none", NA))
  
  # Interactive statistical method selection
  if(interactive_stats && add_stats) {
    cat("\n=== Statistical Analysis Options ===\n")
    
    stat_options <- if(pairing_enabled) {
      c("Paired t-test (recommended for engraftment data)", 
        "Unpaired t-test", 
        "ANOVA with post-hoc comparisons", 
        "No statistical testing")
    } else {
      c("Unpaired t-test", 
        "ANOVA with post-hoc comparisons", 
        "No statistical testing")
    }
    
    iwalk(stat_options, ~cat(sprintf("%d: %s\n", .y, .x)))
    
    stat_choice <- as.integer(readline(prompt = "Select statistical method (enter number): "))
    
    if(is.na(stat_choice) || stat_choice < 1 || stat_choice > length(stat_options)) {
      stop("Invalid selection. Please run the function again and choose a valid number.")
    }
    
    if(pairing_enabled) {
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
    } else {
      if(stat_choice == 1) {
        stat_method <- "t_test"
        cat("Selected: Unpaired t-test\n")
      } else if(stat_choice == 2) {
        stat_method <- "anova"
        cat("Selected: ANOVA with post-hoc\n")
      } else {
        add_stats <- FALSE
        cat("Selected: No statistical testing\n")
      }
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
        available_groups <- unique(data[[tissue_col]])
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
    }
  } else {
    # Non-interactive mode
    comparison_type <- "pairwise"
    control_group <- NULL
  }
  
  # Interactive selection of normalization tissue
  use_normalization <- FALSE
  
  if(interactive) {
    # Show available tissue types
    available_tissues <- unique(data[[tissue_col]])
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
    # Non-interactive mode - use provided normalization_tissue if it exists
    if(!is.null(normalization_tissue) && normalization_tissue %in% unique(data[[tissue_col]])) {
      use_normalization <- TRUE
    } else {
      use_normalization <- FALSE
      normalization_tissue <- NULL
    }
  }
  
  # Process the data to calculate engraftment ratios
  # Create a generic pairing identifier
  if(pairing_enabled) {
    pairing_id_col <- pairing_col
  } else {
    # If no pairing, create a unique ID for each row
    data <- data %>%
      mutate(temp_pairing_id = row_number())
    pairing_id_col <- "temp_pairing_id"
  }
  
  congenic_engraftment <- data %>%
    select(all_of(c(pairing_id_col, tissue_col, marker_col, freq_col)))
  
  # Apply normalization if selected
  if(use_normalization && pairing_enabled) {
    congenic_engraftment <- congenic_engraftment %>%
      # For each marker and pairing group, get the normalization tissue frequency
      group_by(across(all_of(c(marker_col, pairing_id_col)))) %>%
      mutate(
        norm_freq = .data[[freq_col]][.data[[tissue_col]] == normalization_tissue],
        normalized_freq = .data[[freq_col]] / norm_freq
      ) %>%
      ungroup() %>%
      select(all_of(c(pairing_id_col, tissue_col, marker_col)), normalized_freq) %>%
      pivot_wider(names_from = all_of(marker_col), values_from = normalized_freq)
    
    # Update plot labels for normalized data
    y_label <- paste0("Log2(ratio KO:WT normalized to ", normalization_tissue, ")")
    caption_text <- paste0("*Normalized to paired ", normalization_tissue, " controls")
    
  } else {
    congenic_engraftment <- congenic_engraftment %>%
      select(all_of(c(pairing_id_col, tissue_col, marker_col, freq_col))) %>%
      pivot_wider(names_from = all_of(marker_col), values_from = all_of(freq_col))
    
    # Update plot labels for non-normalized data
    y_label <- "Log2(ratio KO:WT - raw frequencies)"
    caption_text <- "*Using raw frequencies (no normalization)"
  }
  
  # Calculate engraftment ratio
  congenic_engraftment <- congenic_engraftment %>%
    mutate(
      engraftment_ratio = log2(.data[[ko_marker]] / .data[[wt_marker]])
    ) %>%
    # Filter out infinite values
    filter(is.finite(engraftment_ratio)) %>%
    select(all_of(c(pairing_id_col, tissue_col)), engraftment_ratio)
  
  # Check if we have data after processing
  if(nrow(congenic_engraftment) == 0) {
    stop("No valid engraftment ratios calculated. Check your marker names and data.")
  }
  
  # Calculate summary statistics
  summary_stats <- congenic_engraftment %>%
    group_by(across(all_of(tissue_col))) %>%
    summarise(
      mean_val = mean(engraftment_ratio, na.rm = TRUE),
      sd_val = sd(engraftment_ratio, na.rm = TRUE),
      n = n(),
      sem_val = sd_val / sqrt(n),
      .groups = 'drop'
    )
  
  # Create factor levels with normalization group first, then order by mean if requested
  all_groups <- unique(summary_stats[[tissue_col]])
  
  if(use_normalization && normalization_tissue %in% all_groups) {
    # Put normalization tissue first
    other_groups <- setdiff(all_groups, normalization_tissue)
    
    if(order_by_mean) {
      # Order other groups by mean (excluding normalization tissue)
      other_stats <- summary_stats %>% 
        filter(.data[[tissue_col]] != normalization_tissue)
      other_groups_ordered <- other_stats[[tissue_col]][order(-other_stats$mean_val)]
      factor_levels <- c(normalization_tissue, other_groups_ordered)
    } else {
      factor_levels <- c(normalization_tissue, other_groups)
    }
  } else {
    # No normalization, just order by mean if requested
    if(order_by_mean) {
      factor_levels <- summary_stats[[tissue_col]][order(-summary_stats$mean_val)]
    } else {
      factor_levels <- all_groups
    }
  }
  
  # Apply factor levels to both summary stats and original data
  summary_stats <- summary_stats %>%
    mutate(!!sym(tissue_col) := factor(.data[[tissue_col]], levels = factor_levels))
  
  congenic_engraftment <- congenic_engraftment %>%
    mutate(!!sym(tissue_col) := factor(.data[[tissue_col]], levels = factor_levels))
  
  # Calculate error bar positions based on bar direction
  summary_stats <- summary_stats %>%
    mutate(
      error_ymin = ifelse(mean_val >= 0, mean_val, mean_val - sem_val),
      error_ymax = ifelse(mean_val >= 0, mean_val + sem_val, mean_val)
    )
  
  # Perform statistical tests if requested
  stat_results <- NULL
  if(add_stats && nrow(summary_stats) > 1) {
    
    # Create formula dynamically
    formula_str <- paste("engraftment_ratio ~", tissue_col)
    test_formula <- as.formula(formula_str)
    
    if(stat_method == "paired_t_test" && pairing_enabled) {
      # For paired t-test, use the existing congenic_engraftment data
      # Determine comparisons based on user choice
      all_groups <- levels(congenic_engraftment[[tissue_col]])
      
      if(comparison_type == "vs_control" && !is.null(control_group)) {
        # Only compare each group to control
        other_groups <- setdiff(all_groups, control_group)
        comparisons <- lapply(other_groups, function(g) c(control_group, g))
      } else {
        # All pairwise comparisons
        comparisons <- combn(all_groups, 2, simplify = FALSE)
      }
      
      # Perform paired t-test for each comparison
      stat_list <- list()
      for(i in seq_along(comparisons)) {
        group1_name <- comparisons[[i]][1]
        group2_name <- comparisons[[i]][2]
        
        # Get data for both groups, ensuring we have paired data
        paired_data <- congenic_engraftment %>%
          filter(.data[[tissue_col]] %in% c(group1_name, group2_name)) %>%
          select(all_of(c(tissue_col, pairing_id_col)), engraftment_ratio) %>%
          pivot_wider(names_from = all_of(tissue_col), values_from = engraftment_ratio) %>%
          filter(complete.cases(.))
        
        if(nrow(paired_data) >= 3) {  # Need at least 3 pairs for meaningful test
          # Perform paired t-test
          test_result <- t.test(paired_data[[group1_name]], 
                                paired_data[[group2_name]], 
                                paired = TRUE)
          
          stat_list[[i]] <- data.frame(
            group1 = group1_name,
            group2 = group2_name,
            n1 = nrow(paired_data),
            n2 = nrow(paired_data),
            statistic = test_result$statistic,
            p = test_result$p.value,
            stringsAsFactors = FALSE
          )
        }
      }
      
      if(length(stat_list) > 0) {
        stat_results <- bind_rows(stat_list) %>%
          mutate(
            p.adj = p.adjust(p, method = "bonferroni"),
            p.adj.signif = case_when(
              p.adj <= 0.001 ~ "***",
              p.adj <= 0.01 ~ "**", 
              p.adj <= 0.05 ~ "*",
              TRUE ~ "ns"
            )
          )
      } else {
        stat_results <- NULL
      }
      
    } else if(stat_method == "t_test") {
      # Unpaired t-test
      if(comparison_type == "vs_control" && !is.null(control_group)) {
        stat_results <- congenic_engraftment %>%
          pairwise_t_test(test_formula, 
                          ref.group = control_group,
                          p.adjust.method = "bonferroni") %>%
          add_significance("p.adj")
      } else {
        stat_results <- congenic_engraftment %>%
          pairwise_t_test(test_formula, 
                          p.adjust.method = "bonferroni") %>%
          add_significance("p.adj")
      }
      
    } else if(stat_method == "anova") {
      # Perform ANOVA first
      anova_result <- congenic_engraftment %>%
        anova_test(test_formula)
      
      cat("\n=== ANOVA Results ===\n")
      print(anova_result)
      
      # If significant, proceed with post-hoc
      if(anova_result$p < 0.05) {
        if(comparison_type == "vs_control" && !is.null(control_group)) {
          stat_results <- congenic_engraftment %>%
            dunn_test(test_formula, p.adjust.method = "bonferroni") %>%
            filter(group1 == control_group | group2 == control_group) %>%
            add_significance("p.adj")
        } else {
          stat_results <- congenic_engraftment %>%
            tukey_hsd(test_formula) %>%
            add_significance("p.adj")
        }
      } else {
        stat_results <- NULL
        cat("ANOVA not significant - no post-hoc testing performed\n")
      }
    }
    
    # Filter out non-significant results (p > 0.05)
    if(!is.null(stat_results)) {
      stat_results <- stat_results %>%
        filter(p.adj <= 0.05) %>%
        mutate(
          group1 = factor(group1, levels = factor_levels),
          group2 = factor(group2, levels = factor_levels)
        )
    }
    
    # Add y positions for stat annotations with proper spacing
    if(!is.null(stat_results) && nrow(stat_results) > 0) {
      y_max <- max(summary_stats$error_ymax, na.rm = TRUE)
      y_min <- min(summary_stats$error_ymin, na.rm = TRUE)
      y_range <- y_max - y_min
      
      # Sort by x-distance to minimize overlap
      stat_results <- stat_results %>%
        mutate(
          x1_pos = as.numeric(group1),
          x2_pos = as.numeric(group2),
          x_distance = abs(x2_pos - x1_pos)
        ) %>%
        arrange(x_distance, x1_pos)
      
      # Calculate y positions with adequate spacing and proper text placement
      base_height <- y_max + y_range * 0.08
      spacing <- y_range * 0.15
      
      stat_results <- stat_results %>%
        mutate(
          layer = row_number(),
          y.position = base_height + (layer - 1) * spacing,
          p_label = case_when(
            p.adj < 0.001 ~ paste0("p=", format(p.adj, scientific = TRUE, digits = 2)),
            p.adj < 0.01 ~ paste0("p=", format(round(p.adj, 4), nsmall = 4)),
            TRUE ~ paste0("p=", format(round(p.adj, 3), nsmall = 3))
          )
        )
    }
  }
  
  # Create base plot using aes() instead of deprecated aes_string()
  p <- ggplot(summary_stats, aes(x = .data[[tissue_col]], y = mean_val)) +
    geom_col(aes(fill = .data[[tissue_col]]), alpha = 0.8, 
             color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = error_ymin, ymax = error_ymax), 
                  width = 0.3, color = "black", linewidth = 0.3) +
    # Always add individual data points with jitter
    geom_jitter(data = congenic_engraftment, 
                aes(x = .data[[tissue_col]], y = engraftment_ratio), 
                size = 2, alpha = 0.7, width = 0.2, height = 0,
                inherit.aes = FALSE) +
    labs(
      title = plot_title,
      x = x_label,
      y = y_label,
      caption = caption_text
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.3),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", alpha = 0.5)
  
  # Add statistical annotations if available
  if(add_stats && !is.null(stat_results) && nrow(stat_results) > 0) {
    # Add comparison brackets and p-values with proper spacing
    for(i in 1:nrow(stat_results)) {
      row <- stat_results[i, ]
      
      # Get x positions for the groups
      x1 <- as.numeric(row$group1)
      x2 <- as.numeric(row$group2)
      
      # Calculate bracket and text positioning
      bracket_height <- row$y.position
      text_height <- bracket_height + max(summary_stats$error_ymax) * 0.04
      bracket_offset <- max(summary_stats$error_ymax) * 0.03
      
      # Add bracket lines
      p <- p + 
        annotate("segment", 
                 x = x1, xend = x1, 
                 y = bracket_height - bracket_offset, yend = bracket_height,
                 color = "black", linewidth = 0.4) +
        annotate("segment", 
                 x = x1, xend = x2, 
                 y = bracket_height, yend = bracket_height,
                 color = "black", linewidth = 0.4) +
        annotate("segment", 
                 x = x2, xend = x2, 
                 y = bracket_height, yend = bracket_height - bracket_offset,
                 color = "black", linewidth = 0.4) +
        annotate("text", 
                 x = (x1 + x2) / 2, y = text_height,
                 label = row$p_label, size = 2.8, hjust = 0.5, vjust = 0)
    }
  }
  
  # Return plot and optionally statistical results
  result <- list(plot = p)
  
  if(add_stats && !is.null(stat_results)) {
    result$statistics <- stat_results
    result$summary_stats <- summary_stats
  }
  
  # If only plot requested, return just the plot
  if(!add_stats) {
    return(p)
  } else {
    return(result)
  }
}

#' Updated MFI heatmap functions using standardized factors
create_mfi_heatmaps_enhanced <- function(mfi_data,
                                         tissue_col = "tissue_factor",
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
  
  # Updated required columns to use tissue_factor instead of tissue_factor
  required_cols <- c("marker", "MFI", tissue_col, "congenics", "NodeShort")
  missing_cols <- setdiff(required_cols, names(mfi_data))
  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ', ')))
  }
  
  # Rest of the function remains the same, but replace tissue_factor with tissue_col
  mfi_clean <- mfi_data %>%
    filter(!is.na(congenics))
  
  if(!is.null(selected_congenics)) {
    mfi_clean <- mfi_clean %>%
      filter(congenics %in% selected_congenics)
    cat("Filtered to selected congenics:", paste(selected_congenics, collapse = ", "), "\n")
  }
  
  if(!is.null(selected_markers)) {
    mfi_clean <- mfi_clean %>%
      filter(marker %in% selected_markers)
    cat("Filtered to selected markers:", paste(selected_markers, collapse = ", "), "\n")
  }
  
  if(nrow(mfi_clean) == 0) {
    stop("No data remaining after filtering")
  }
  
  # Update references to use tissue_col
  if(grouping_option$type == "separate") {
    group_names <- unique(mfi_clean[[tissue_col]])
    cat("Creating separate heatmaps for", length(group_names), "tissue(s):", paste(group_names, collapse = ", "), "\n")
    
    heatmap_list <- list()
    
    for(group in group_names) {
      cat("Processing tissue:", group, "\n")
      
      group_data <- mfi_clean %>%
        filter(.data[[tissue_col]] == group)
      
      if(nrow(group_data) == 0) {
        warning(paste("No data found for tissue:", group))
        next
      }
      
      heatmap_matrix <- group_data %>%
        group_by(marker, congenics) %>%
        summarise(agg_MFI = mean(MFI, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = congenics, values_from = agg_MFI) %>%
        column_to_rownames("marker") %>%
        as.matrix()
      
      if(nrow(heatmap_matrix) == 0 || ncol(heatmap_matrix) == 0) {
        warning(paste("Empty matrix for tissue:", group))
        next
      }
      
      # Create individual heatmap (simplified version)
      scaled_result <- apply_scaling_method(heatmap_matrix, scale_method, "Mean MFI", 
                                            log_base, percentile_range)
      
      ht <- Heatmap(
        scaled_result$matrix,
        name = scaled_result$legend_title,
        col = scaled_result$color_function,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_cols,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        column_title = group,
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        row_title = "Markers",
        row_title_gp = gpar(fontsize = 12)
      )
      
      heatmap_list[[group]] <- ht
      
      cat(sprintf("  - Created heatmap with %d markers and %d congenics\n", 
                  nrow(heatmap_matrix), ncol(heatmap_matrix)))
    }
    
    return(heatmap_list)
  } else {
    # Combined heatmap logic would go here with similar tissue_col replacements
    stop("Combined heatmap not implemented in this simplified version")
  }
}

#' Updated data cleaning function using standardized factors
data_clean_custom_enhanced <- function(data) {
  
  # Function to handle unstained sample detection and removal
  handle_unstained_samples <- function(df) {
    # Check if Sample column exists
    if(!"Sample" %in% names(df)) {
      return(df)
    }
    
    # Define unstained patterns (case insensitive)
    unstained_patterns <- c("unstained", "no stain", "no_stain")
    
    # Create regex pattern for case-insensitive matching
    pattern <- paste(unstained_patterns, collapse = "|")
    
    # Find rows with unstained samples
    unstained_rows <- df %>%
      mutate(row_id = row_number()) %>%
      filter(str_detect(tolower(Sample), pattern)) %>%
      pull(row_id)
    
    # If no unstained samples found, return original dataframe
    if(length(unstained_rows) == 0) {
      cat("No unstained samples detected.\n")
      return(df)
    }
    
    # Display detected unstained samples
    cat("\n=== Unstained Samples Detected ===\n")
    cat(paste("Found", length(unstained_rows), "unstained sample(s):\n"))
    
    unstained_data <- df %>%
      slice(unstained_rows)
    
    # Select relevant columns for display
    if("NodeShort" %in% names(df)) {
      unstained_data <- unstained_data %>%
        select(Sample, NodeShort)
    } else {
      unstained_data <- unstained_data %>%
        select(Sample)
    }
    
    print(unstained_data)
    
    # Interactive prompt for user decision
    cat("\n=== Action Required ===\n")
    action_options <- c(
      "Remove unstained samples from dataset",
      "Keep unstained samples in dataset"
    )
    
    user_choice <- menu(action_options, title = "What would you like to do with these unstained samples?")
    
    if(user_choice == 0) {
      cat("Action cancelled. Keeping unstained samples in dataset.\n")
      return(df)
    } else if(user_choice == 1) {
      # Remove unstained samples
      df_cleaned <- df %>%
        slice(-unstained_rows)
      
      cat(paste("Removed", length(unstained_rows), "unstained sample(s).\n"))
      cat(paste("Dataset reduced from", nrow(df), "to", nrow(df_cleaned), "rows.\n"))
      
      return(df_cleaned)
    } else {
      # Keep unstained samples
      cat("Keeping unstained samples in dataset.\n")
      return(df)
    }
  }
  
  # Function to clean a single data frame
  clean_single_df <- function(df) {
    # Clean column names - remove "$" symbols
    names(df) <- str_remove_all(names(df), "\\$")
    
    # Clean NodeShort column if it exists
    if("NodeShort" %in% names(df)) {
      df <- df %>%
        mutate(
          NodeShort = NodeShort %>%
            str_remove(".*:") %>%        # Remove everything before and including ":"
            str_remove_all(" ") %>%      # Remove all spaces
            str_remove_all(",") %>%      # Remove all commas
            str_remove_all(":")          # Remove any remaining colons
        )
    }
    
    # Handle unstained samples
    df <- handle_unstained_samples(df)
    
    # Validate standardized factor columns
    if("tissue_factor" %in% names(df)) {
      cat("✅ Found tissue_factor column\n")
      tissue_summary <- df %>% count(tissue_factor, name = "n") %>% arrange(desc(n))
      cat("Tissue groups:", nrow(tissue_summary), "\n")
    } else {
      cat("ℹ️  No tissue_factor column found\n")
    }
    
    if("tissue_factor" %in% names(df)) {
      cat("✅ Found tissue_factor column\n")
      pairing_summary <- df %>% count(tissue_factor, name = "n") %>% arrange(desc(n))
      cat("Pairing groups:", nrow(pairing_summary), "\n")
    } else {
      cat("ℹ️  No tissue_factor column found\n")
    }
    
    return(df)
  }
  
  # Check if input is a list
  if(is.list(data) && !is.data.frame(data)) {
    # Apply cleaning function to each element in the list
    cleaned_data <- map(data, ~ {
      if(is.data.frame(.x)) {
        clean_single_df(.x)
      } else {
        .x  # Return unchanged if not a data frame
      }
    })
    return(cleaned_data)
  } else if(is.data.frame(data)) {
    # Apply cleaning function to single data frame
    return(clean_single_df(data))
  } else {
    # Return unchanged if neither list nor data frame
    warning("Input is neither a list nor a data frame. Returning unchanged.")
    return(data)
  }
}

#=============================================================================
# UMAP Interactive tool
# ============================================================================

# ============================================================================
# ENHANCED INTERACTIVE UMAP ANALYSIS FOR FLOW CYTOMETRY DATA
# Fixed version with sample exclusion and streamlined workflow
# ============================================================================

library(tidyverse)
library(umap)
library(ggplot2)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(scales)

# ============================================================================
# MISSING DEPENDENCY FUNCTIONS - REPLACE WITH YOUR ACTUAL IMPLEMENTATIONS
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
      return(exclude_samples_by_pattern(all_samples))
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
      # Multi-marker heatmap
      cat("\n=== Multi-Marker Heatmap Creation ===\n")
      cat("Available markers:", length(marker_names), "\n")
      
      if(length(marker_names) == 0) {
        cat("No markers available for heatmap\n")
        next
      }
      
      cat("Options:\n")
      cat("1. Select specific markers\n")
      cat("2. Use first 9 markers (3x3 grid)\n")
      cat("3. Use all markers (may be crowded)\n")
      cat("4. Show marker list first\n")
      
      heatmap_choice <- readline("Choose option (1-4): ")
      
      selected_markers <- NULL
      
      if(heatmap_choice == "1") {
        cat("\nAvailable markers:\n")
        iwalk(marker_names, ~cat(sprintf("%d. %s\n", .y, .x)))
        
        cat("\nEnter marker numbers (space-separated, max 12 recommended): ")
        marker_selection <- readline()
        
        if(marker_selection == "") {
          cat("No markers selected\n")
          next
        }
        
        # Parse space-separated marker numbers
        tryCatch({
          selections <- str_trim(str_split(marker_selection, "\\s+")[[1]])
          marker_indices <- numeric(0)
          
          for(sel in selections) {
            if(grepl("^\\d+$", sel)) {
              marker_num <- as.numeric(sel)
              if(!is.na(marker_num) && marker_num >= 1 && marker_num <= length(marker_names)) {
                marker_indices <- c(marker_indices, marker_num)
              }
            }
          }
          
          if(length(marker_indices) == 0) {
            cat("No valid markers selected\n")
            next
          }
          
          selected_markers <- marker_names[marker_indices]
          
        }, error = function(e) {
          cat("Invalid input format\n")
          next
        })
        
      } else if(heatmap_choice == "2") {
        n_markers <- min(9, length(marker_names))
        selected_markers <- marker_names[1:n_markers]
        cat("Using first", n_markers, "markers\n")
        
      } else if(heatmap_choice == "3") {
        if(length(marker_names) > 16) {
          cat("Warning: Using all", length(marker_names), "markers may create a crowded plot\n")
          confirm <- readline("Continue? (y/n): ")
          if(tolower(confirm) != "y") next
        }
        selected_markers <- marker_names
        
      } else if(heatmap_choice == "4") {
        cat("\nAll available markers:\n")
        iwalk(marker_names, ~cat(sprintf("%d. %s\n", .y, .x)))
        readline("Press Enter to continue...")
        next
        
      } else {
        cat("Invalid choice\n")
        next
      }
      
      if(length(selected_markers) == 0) {
        cat("No markers selected\n")
        next
      }
      
      cat("Creating heatmap with", length(selected_markers), "markers:", 
          paste(head(selected_markers, 5), collapse = ", "))
      if(length(selected_markers) > 5) cat("...")
      cat("\n")
      
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
      
      # Determine grid layout
      n_markers <- length(selected_markers)
      if(n_markers <= 4) {
        ncol_grid <- 2
      } else if(n_markers <= 9) {
        ncol_grid <- 3
      } else if(n_markers <= 16) {
        ncol_grid <- 4
      } else {
        ncol_grid <- ceiling(sqrt(n_markers))
      }
      
      # Arrange in grid
      combined_plot <- wrap_plots(marker_plots, ncol = ncol_grid)
      
      # Store and display
      .umap_session_plots$plot_counter <<- .umap_session_plots$plot_counter + 1
      plot_id <- paste0("plot_", .umap_session_plots$plot_counter)
      
      plot_title <- paste("Multi-marker heatmap (", length(selected_markers), "markers )")
      
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
      
      # Option to save
      save_choice <- readline("Save this heatmap? (y/n): ")
      if(tolower(save_choice) == "y") {
        filename <- readline("Enter filename: ")
        if(filename != "") {
          if(!str_detect(filename, "\\.(png|pdf|jpg|jpeg|tiff|svg)$")) {
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
      cat("Session now contains", length(.umap_session_plots$plots), "plots\n")
      
    } else if(choice == "11") {
      # This would be used if we had current exclusions
      cat("Proceeding with current sample list\n")
      return(all_samples)
    } else {
      cat("Invalid choice. Please select 1-6.\n")
    }
  }
}

exclude_samples_by_pattern <- function(all_samples) {
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
      
      # Fixed: Handle both numeric and text input safely
      selected_pattern <- NULL
      
      # Try numeric first
      if(grepl("^\\d+$", pattern_choice)) {
        pattern_num <- as.numeric(pattern_choice)
        if(!is.na(pattern_num) && pattern_num >= 1 && pattern_num <= length(common_patterns)) {
          selected_pattern <- common_patterns[[pattern_num]]
        }
      } else if(pattern_choice != "") {
        # Use as custom pattern
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
  
  # Use select.list for cross-platform compatibility
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
  
  # Fixed: Safe column selection
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
  
  # Show unique values in selected column
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
    filter(.data[[selected_col]] %in% values_to_exclude) %>%
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
                           
                           # Fixed: Safe numeric conversion
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
                             
                             # Fixed: Safe numeric conversion
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
  
  # Identify potential markers (exclude scatter, time, etc.)
  exclude_patterns <- c("FSC", "SSC", "Time", "Event", "Original", "Width", "Height")
  exclude_regex <- paste0("(?i)", paste(exclude_patterns, collapse = "|"))
  
  potential_markers <- lookup %>%
    filter(!str_detect(colname, exclude_regex)) %>%
    filter(!str_detect(marker, exclude_regex))
  
  # Identify compensated channels
  comp_channels <- potential_markers %>%
    filter(str_detect(colname, "(?i)comp")) %>%
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
                 lookup = potential_markers %>% filter(colname %in% selected_channels)
               ))
             }
           },
           
           "2" = {
             if(length(comp_channels) == 0) {
               cat("No compensated channels found. Try option 1 or 3.\n")
               next
             }
             
             selected_lookup <- potential_markers %>% filter(colname %in% comp_channels)
             
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
    
    # Check if node exists
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
  
  # Add summary statistics
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
  
  # Use selected samples if provided, otherwise use all
  if(is.null(selected_samples)) {
    selected_samples <- sampleNames(gs)
  }
  
  cat("Extracting single-cell data from node:", basename(node), "\n")
  cat("Using", length(selected_samples), "samples\n")
  
  # Get sample metadata - with error handling
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
  
  # Determine samples to process
  samples_to_process <- selected_samples
  if(!is.null(sample_limit) && sample_limit < length(samples_to_process)) {
    samples_to_process <- sample(samples_to_process, sample_limit)
    cat("Randomly selected", sample_limit, "samples for analysis\n")
  }
  
  # Extract data from each sample
  all_data <- map_dfr(samples_to_process, function(sample_name) {
    cat("Processing sample:", sample_name, "\n")
    
    gh <- gs[[sample_name]]
    
    # Check if node exists in this sample
    if(!node %in% gh_get_pop_paths(gh)) {
      cat("  Node not found in sample, skipping\n")
      return(tibble())
    }
    
    # Get single-cell data
    tryCatch({
      ff <- gh_pop_get_data(gh, node)
      
      if(is.null(ff)) {
        cat("  No data found, skipping\n")
        return(tibble())
      }
      
      # Extract expression matrix
      expr_data <- if(inherits(ff, "cytoframe")) exprs(ff) else exprs(ff)
      
      # Check if we have the required channels
      available_channels <- intersect(selected_markers$channels, colnames(expr_data))
      
      if(length(available_channels) == 0) {
        cat("  No selected channels found, skipping\n")
        return(tibble())
      }
      
      # Convert to tibble and add sample info
      cell_data <- as_tibble(expr_data[, available_channels, drop = FALSE])
      
      # Subsample if too many cells
      if(nrow(cell_data) > max_cells_per_sample) {
        cell_data <- slice_sample(cell_data, n = max_cells_per_sample)
        cat("  Subsampled to", max_cells_per_sample, "cells\n")
      }
      
      # Add sample identifier and cell ID
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
  
  # Replace channel names with marker names
  marker_lookup <- selected_markers$lookup %>%
    select(colname, marker) %>%
    deframe()
  
  # Rename columns (keeping Sample and CellID)
  for(channel in names(marker_lookup)) {
    if(channel %in% names(all_data)) {
      names(all_data)[names(all_data) == channel] <- marker_lookup[channel]
    }
  }
  
  # Add sample metadata
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
  
  # Get marker columns (exclude Sample, CellID, and metadata)
  metadata_cols <- c("Sample", "CellID")
  if(!is.null(single_cell_data$sample_metadata)) {
    metadata_cols <- c(metadata_cols, names(single_cell_data$sample_metadata))
  }
  
  # Also exclude any non-numeric columns that might have been added via external metadata
  all_cols <- names(single_cell_data$data)
  numeric_cols <- all_cols[map_lgl(all_cols, ~is.numeric(single_cell_data$data[[.x]]))]
  marker_cols <- setdiff(numeric_cols, metadata_cols)
  
  cat("Using", length(marker_cols), "markers:", paste(marker_cols, collapse = ", "), "\n")
  
  # Prepare data matrix
  umap_matrix <- single_cell_data$data %>%
    select(all_of(marker_cols)) %>%
    as.matrix()
  
  # Check for problematic values before transformation
  cat("Data range before transformation: [", min(umap_matrix, na.rm = TRUE), ", ", 
      max(umap_matrix, na.rm = TRUE), "]\n")
  
  if(any(is.na(umap_matrix))) {
    n_na <- sum(is.na(umap_matrix))
    cat("Warning:", n_na, "NA values found in data. These will be handled.\n")
  }
  
  # Apply transformation with proper handling of edge cases
  if(transform_data == "asinh") {
    umap_matrix <- asinh(umap_matrix / 5)  # Typical cofactor for flow cytometry
    cat("Applied asinh transformation (cofactor = 5)\n")
    
  } else if(transform_data == "log10") {
    # Handle negative and zero values properly
    min_val <- min(umap_matrix, na.rm = TRUE)
    
    if(min_val <= 0) {
      # Shift data to ensure all values are positive
      shift_amount <- abs(min_val) + 1
      cat("Shifting data by", shift_amount, "to handle negative/zero values\n")
      umap_matrix <- umap_matrix + shift_amount
    } else {
      # Add small constant to avoid log(0) 
      umap_matrix <- umap_matrix + 1
    }
    
    umap_matrix <- log10(umap_matrix)
    cat("Applied log10 transformation\n")
    
  } else if(transform_data == "sqrt") {
    # Handle negative values for sqrt
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
        if(is.na(median_val)) median_val <- 0  # fallback if all values are problematic
        
        umap_matrix[is.nan(col_data) | is.infinite(col_data), i] <- median_val
        cat("Replaced problematic values in", marker_cols[i], "with median:", median_val, "\n")
      }
    }
  }
  
  # Final check for remaining NA values
  if(any(is.na(umap_matrix))) {
    # Replace remaining NAs with column medians
    for(i in 1:ncol(umap_matrix)) {
      col_data <- umap_matrix[, i]
      if(any(is.na(col_data))) {
        median_val <- median(col_data, na.rm = TRUE)
        if(is.na(median_val)) median_val <- 0  # fallback
        
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
    # Remove spaces and split by commas
    parts <- str_trim(str_split(input, ",")[[1]])
    indices <- numeric(0)
    
    for (part in parts) {
      if (str_detect(part, "-")) {
        # Handle ranges like "1-5"
        range_parts <- as.numeric(str_split(part, "-")[[1]])
        if (length(range_parts) == 2 && all(!is.na(range_parts))) {
          start_idx <- range_parts[1]
          end_idx <- range_parts[2]
          if (start_idx >= 1 && end_idx <= max_value && start_idx <= end_idx) {
            indices <- c(indices, start_idx:end_idx)
          }
        }
      } else if (str_detect(part, "^\\d+$")) {
        # Handle single numbers
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
  
  # Combine and categorize options
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
        legend.key.size = unit(0.6, "cm"), # Legend key dot size
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
      categorical_vars <- all_options[map_lgl(all_options, ~!is.numeric(umap_data[[.x]]))]
      
      if (length(categorical_vars) == 0) {
        cat("No categorical variables available for faceting\n")
        next
      }
      
      facet_var <- select_from_menu(
        categorical_vars,
        "Select variable for faceting"
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
      
      # Step 2: Facet variable (exclude color variable if categorical)
      categorical_vars <- all_options[map_lgl(all_options, ~!is.numeric(umap_data[[.x]]))]
      available_facet_vars <- setdiff(categorical_vars, color_var)
      
      if (length(available_facet_vars) == 0) {
        cat("No additional categorical variables available for faceting\n")
        next
      }
      
      facet_var <- select_from_menu(
        available_facet_vars,
        "Select variable for faceting"
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

import_external_metadata <- function(gs) {
  cat("\n=== External Metadata Import ===\n")
  cat("This function will help you import and merge external metadata\n")
  cat("with your flow cytometry samples.\n\n")
  
  # Get current sample names for matching
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
      # CSV import
      file_path <- readline("Enter CSV file path: ")
      if(file.exists(file_path)) {
        tryCatch({
          metadata <- read_csv(file_path, show_col_types = FALSE)
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
      # R object import
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
      
      # Fixed: Safe object selection
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
      
      metadata <- get(selected_obj, envir = .GlobalEnv)
      cat("Using object '", selected_obj, "' with", nrow(metadata), "rows and", ncol(metadata), "columns\n")
      
    } else if(choice == "3") {
      cat("Skipping metadata import\n")
      return(NULL)
    } else {
      cat("Invalid choice\n")
      next
    }
    
    # If we have metadata, process it
    if(!is.null(metadata)) {
      # Show preview
      cat("\nMetadata preview:\n")
      print(head(metadata, 5))
      cat("Column names:", paste(names(metadata), collapse = ", "), "\n\n")
      
      # Identify sample ID column
      cat("Sample ID Column Selection:\n")
      cat("Available columns:\n")
      iwalk(names(metadata), ~cat(sprintf("%d. %s\n", .y, .x)))
      
      sample_col_choice <- readline("Enter column number/name for Sample ID: ")
      
      # Fixed: Safe column selection
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
      
      # Check sample matching
      metadata_samples <- metadata[[sample_col]]
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
        next
      }
      
      if(length(unmatched_gs) > 0) {
        cat("Unmatched GatingSet samples:", paste(head(unmatched_gs, 5), collapse = ", "))
        if(length(unmatched_gs) > 5) cat("...")
        cat("\n")
      }
      
      # Select metadata columns to include
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
        # Parse space-separated column numbers/names
        tryCatch({
          # Split by spaces and clean up
          selections <- str_trim(str_split(col_selection, "\\s+")[[1]])
          selected_cols <- character(0)
          
          for(sel in selections) {
            if(grepl("^\\d+$", sel)) {
              # Numeric selection
              col_num <- as.numeric(sel)
              if(!is.na(col_num) && col_num >= 1 && col_num <= length(other_cols)) {
                selected_cols <- c(selected_cols, other_cols[col_num])
              }
            } else if(sel %in% other_cols) {
              # Column name selection
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
      
      # Prepare final metadata
      final_metadata <- metadata %>%
        select(all_of(c(sample_col, selected_cols))) %>%
        rename(Sample = all_of(sample_col)) %>%
        filter(Sample %in% current_samples)
      
      cat("Final metadata prepared with", nrow(final_metadata), 
          "samples and", ncol(final_metadata)-1, "metadata columns\n")
      cat("Final columns:", paste(names(final_metadata), collapse = ", "), "\n")
      
      return(final_metadata)
    }
  }
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
  
  # Fixed: Safe numeric input handling
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
  
  # Step 5: Import external metadata (optional)
  cat("\n=== Metadata Import ===\n")
  import_choice <- readline("Import external metadata? (y/n, default n): ")
  
  external_metadata <- NULL
  if(tolower(import_choice) == "y") {
    external_metadata <- import_external_metadata(gs)
    
    if(!is.null(external_metadata)) {
      # Merge with existing data
      cat("Merging external metadata with single-cell data...\n")
      
      # Get current metadata columns to avoid conflicts
      current_meta_cols <- setdiff(names(single_cell_data$data), 
                                   c("CellID", single_cell_data$marker_names))
      
      # Handle column name conflicts
      ext_meta_cols <- setdiff(names(external_metadata), "Sample")
      conflicting_cols <- intersect(current_meta_cols, ext_meta_cols)
      
      if(length(conflicting_cols) > 0) {
        cat("Column name conflicts found:", paste(conflicting_cols, collapse = ", "), "\n")
        cat("External metadata columns will be prefixed with 'ext_'\n")
        
        for(col in conflicting_cols) {
          names(external_metadata)[names(external_metadata) == col] <- paste0("ext_", col)
        }
      }
      
      # Merge data with many-to-many relationship handling
      original_nrow <- nrow(single_cell_data$data)
      
      # Check for duplicates before merging
      duplicate_samples_ext <- external_metadata %>%
        count(Sample) %>%
        filter(n > 1)
      
      if(nrow(duplicate_samples_ext) > 0) {
        cat("Warning: Duplicate Sample IDs found in external metadata:\n")
        print(duplicate_samples_ext)
        cat("Using first occurrence of each Sample ID\n")
        
        external_metadata <- external_metadata %>%
          group_by(Sample) %>%
          slice(1) %>%
          ungroup()
      }
      
      single_cell_data$data <- single_cell_data$data %>%
        left_join(external_metadata, by = "Sample", relationship = "many-to-one")
      
      if(nrow(single_cell_data$data) != original_nrow) {
        cat("Warning: Row count changed during merge (", original_nrow, "->", nrow(single_cell_data$data), ")\n")
        cat("This suggests duplicate Sample IDs. Please check your data.\n")
      } else {
        cat("Merge successful: row count maintained\n")
      }
      
      cat("External metadata successfully merged!\n")
      cat("Added columns:", paste(setdiff(names(external_metadata), "Sample"), collapse = ", "), "\n")
    }
  }
  
  # Step 6: UMAP parameters
  cat("\n=== UMAP Parameters ===\n")
  
  # Fixed: Safe numeric input handling for UMAP parameters
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
    "asinh"  # default
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

# Function to export all session plots at once
export_all_session_plots <- function(directory = "umap_plots", 
                                     width = 8, 
                                     height = 6, 
                                     dpi = 300,
                                     format = "png") {
  
  if(length(.umap_session_plots$plots) == 0) {
    cat("No plots to export\n")
    return(invisible(NULL))
  }
  
  # Create directory if it doesn't exist
  if(!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
    cat("Created directory:", directory, "\n")
  }
  
  # Export each plot
  exported_files <- map_chr(names(.umap_session_plots$plots), function(plot_id) {
    plot_info <- .umap_session_plots$plots[[plot_id]]
    
    # Create clean filename
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

# Function to restart visualization session with existing data
restart_visualization_session <- function(umap_results) {
  cat("Restarting visualization session with existing UMAP results...\n")
  
  if(is.null(umap_results$umap_result)) {
    stop("Invalid UMAP results object")
  }
  
  initialize_visualization_session(umap_results$umap_result)
  interactive_visualization_menu()
  
  return(invisible(.umap_session_plots))
}

# Function to create quick plots from results
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
  
  # Choose appropriate color scale
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

# Function to export UMAP data
export_umap_data <- function(umap_results, filename = "umap_data.csv") {
  
  export_data <- umap_results$final_data %>%
    select(-CellID)  # Remove cell ID for cleaner export
  
  write_csv(export_data, filename)
  cat("UMAP data exported to:", filename, "\n")
  cat("Columns exported:", paste(names(export_data), collapse = ", "), "\n")
  
  return(invisible(export_data))
}

# Function to get summary statistics by group
summarize_umap_by_group <- function(umap_results, group_column) {
  
  umap_data <- umap_results$final_data
  marker_names <- umap_results$umap_result$marker_names
  
  if(!group_column %in% names(umap_data)) {
    stop("Column '", group_column, "' not found in UMAP data")
  }
  
  # Summary statistics for markers by group
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
  
  # UMAP coordinate summary
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
# EXAMPLE USAGE AND DOCUMENTATION
# ============================================================================

cat("Enhanced Interactive UMAP Flow Cytometry Analysis Functions Loaded!\n")
cat("Main function: analyze_flow_umap_enhanced(gs)\n")
cat("\nNew features:\n")
cat("- Interactive sample exclusion (remove unstained, controls, etc.)\n")
cat("- Event counting with population assessment\n")
cat("- Optional external metadata import from CSV\n")
cat("- Persistent plot session with save/reload capability\n")
cat("- Interactive visualization menu with multiple plot types\n")
cat("- Session management and export utilities\n")
cat("- Fixed all numeric input handling to prevent NA errors\n")

# Example workflow:
# 
# # 1. Setup your gating set (using your existing function)
# gs <- setup_flowjo_workspace(
#   xml_path = here("data/your_workspace.wsp"),
#   fcs_path = here("data/fcs_files/")
# )
# 
# # 2. Run enhanced interactive UMAP analysis
#results <- analyze_flow_umap_enhanced(gs)
# 
# # 3. Continue visualization session later (plots persist)
# restart_visualization_session(results)
# 
# # 4. Create additional plots from results
# sample_plot <- create_quick_umap_plot(results, "tissue_factor")
# print(sample_plot)
# 
# # 5. Export data and session
# export_umap_data(results, "my_analysis.csv")
# export_all_session_plots("plots", format = "png")

#=============================================================================
# EXAMPLE USAGE
# ============================================================================
# 
# # # Setup
# gs <- setup_flowjo_workspace(
#   xml_path = here("data/20-Jun-2025.wsp"),
#   fcs_path = here("data/fcs_files/")
#  )
# # 
# # # Run interactive analysis to get subpop-specifc DFs
# congenics_results <-analyze_flow_data_auto(gs)
# 
# #Clean-Up Data Node naming data 
# congenics_results <- data_clean_custom(congenics_results)
# 
# # Add genotype information (Only works if the code can detect congenic markers)
# congenics_results_with_genotype <- assign_genotypes_menu(congenics_results1$counts)

#===============================================================================
