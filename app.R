library(shiny)
library(DT)
library(jsonlite)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)

# Base data directories for flat structure
BASE_DIRS <- list(
  bgee_summary = "./data/bgee_summary_stats/",
  fasta_compare = "./data/fasta_compare_stats/",
  fastp_json = "./data/fastp_json/",
  metadata = "./data/"
)

BASE_DATA_DIR <- "./data"

# Load and parse sample metadata
load_sample_metadata <- function() {
  if (!dir.exists(BASE_DATA_DIR)) {
    warning("Data directory not found at ", BASE_DATA_DIR)
    return(data.frame())
  }
  
  # Look for CSV files in the data directory
  csv_files <- list.files(BASE_DATA_DIR, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) {
    warning("No CSV files found in ", BASE_DATA_DIR)
    return(data.frame())
  }
  
  # Try each CSV file to find one that looks like metadata
  for (csv_file in csv_files) {
    tryCatch({
      metadata <- read.csv(csv_file, stringsAsFactors = FALSE)
      
      # Check for different possible column name variations
      has_sample_id <- any(c("sample_id", "Sample ID", "Sample.ID", "sampleid") %in% names(metadata))
      has_process_id <- any(c("process_id", "Process ID", "Process.ID", "processid") %in% names(metadata))
      
      if (has_sample_id && has_process_id) {
        # Standardize column names
        names(metadata)[names(metadata) %in% c("Sample ID", "Sample.ID", "sampleid")] <- "sample_id"
        names(metadata)[names(metadata) %in% c("Process ID", "Process.ID", "processid")] <- "process_id"
        
        message("Found metadata file: ", basename(csv_file))
        message("Loaded ", nrow(metadata), " sample records")
        return(metadata)
      }
    }, error = function(e) {
      # Skip files that can't be read
      message("Could not read ", basename(csv_file), ": ", e$message)
      next
    })
  }
  
  warning("No suitable metadata CSV file found in ", BASE_DATA_DIR)
  warning("Expected columns: 'Sample ID' and 'Process ID' (or variations)")
  return(data.frame())
}

# Extract plate ID from Sample ID with proper validation
extract_plate_id_from_sample <- function(sample_ids) {
  # Use more specific regex that preserves leading zeros
  plate_ids <- character(length(sample_ids))
  
  for (i in seq_along(sample_ids)) {
    sample_id <- sample_ids[i]
    
    # Check if it matches BGE format with proper structure
    if (grepl("^BGE_\\d{5}_[A-Z]\\d{2}$", sample_id, ignore.case = TRUE)) {
      # Extract the plate part (BGE_XXXXX) preserving leading zeros
      plate_match <- regmatches(sample_id, regexpr("^BGE_\\d{5}", sample_id))
      
      # Validate minimum plate number (BGE_00001 or higher)
      plate_number <- as.numeric(gsub("BGE_", "", plate_match))
      if (!is.na(plate_number) && plate_number >= 1) {
        plate_ids[i] <- plate_match
      } else {
        plate_ids[i] <- NA_character_
      }
    } else {
      plate_ids[i] <- NA_character_
    }
  }
  
  return(plate_ids)
}

# Get all available process IDs across all datasets
get_available_process_ids <- function() {
  all_process_ids <- character(0)
  
  # Get process IDs from BGEE summary stats files
  if (dir.exists(BASE_DIRS$bgee_summary)) {
    csv_files <- list.files(BASE_DIRS$bgee_summary, pattern = "\\.csv$", full.names = TRUE)
    for (file in csv_files) {
      tryCatch({
        data <- read.csv(file, stringsAsFactors = FALSE)
        # Handle different possible column names
        if ("Process.ID" %in% names(data)) {
          process_ids <- data$Process.ID
        } else if ("Process ID" %in% names(data)) {
          process_ids <- data$`Process ID`
        } else if ("process_id" %in% names(data)) {
          process_ids <- data$process_id
        } else if ("ID" %in% names(data)) {
          process_ids <- data$ID
        } else {
          next
        }
        all_process_ids <- c(all_process_ids, process_ids)
      }, error = function(e) {
        message(paste("Error reading", file, ":", e$message))
      })
    }
  }
  
  return(unique(all_process_ids[!is.na(all_process_ids) & all_process_ids != ""]))
}


# Parse fastp JSON files
parse_json_file <- function(json_path, process_id) {
  tryCatch({
    data <- fromJSON(json_path)
    
    # Helper function to safely extract numeric values
    safe_numeric <- function(value) {
      if (is.null(value) || length(value) == 0) {
        return(NA_real_)
      }
      as.double(value)
    }
    
    result <- data.frame(
      process_id = process_id,
      # Before filtering
      before_total_reads = safe_numeric(data$summary$before_filtering$total_reads),
      before_total_bases = safe_numeric(data$summary$before_filtering$total_bases),
      before_q20_bases = safe_numeric(data$summary$before_filtering$q20_bases),
      before_q30_bases = safe_numeric(data$summary$before_filtering$q30_bases),
      before_q20_rate = safe_numeric(data$summary$before_filtering$q20_rate),
      before_q30_rate = safe_numeric(data$summary$before_filtering$q30_rate),
      before_read1_mean_length = safe_numeric(data$summary$before_filtering$read1_mean_length),
      before_read2_mean_length = safe_numeric(data$summary$before_filtering$read2_mean_length),
      before_gc_content = safe_numeric(data$summary$before_filtering$gc_content),
      # After filtering
      after_total_reads = safe_numeric(data$summary$after_filtering$total_reads),
      after_total_bases = safe_numeric(data$summary$after_filtering$total_bases),
      after_q20_bases = safe_numeric(data$summary$after_filtering$q20_bases),
      after_q30_bases = safe_numeric(data$summary$after_filtering$q30_bases),
      after_q20_rate = safe_numeric(data$summary$after_filtering$q20_rate),
      after_q30_rate = safe_numeric(data$summary$after_filtering$q30_rate),
      after_read1_mean_length = safe_numeric(data$summary$after_filtering$read1_mean_length),
      after_read2_mean_length = safe_numeric(data$summary$after_filtering$read2_mean_length),
      after_gc_content = safe_numeric(data$summary$after_filtering$gc_content),
      # Filtering results
      passed_filter_reads = safe_numeric(data$filtering_result$passed_filter_reads),
      low_quality_reads = safe_numeric(data$filtering_result$low_quality_reads),
      too_many_N_reads = safe_numeric(data$filtering_result$too_many_N_reads),
      too_short_reads = safe_numeric(data$filtering_result$too_short_reads),
      too_long_reads = safe_numeric(data$filtering_result$too_long_reads),
      # Duplication
      duplication_rate = safe_numeric(data$duplication$rate),
      stringsAsFactors = FALSE
    )
    
    return(result)
  }, error = function(e) {
    message(paste("Error parsing JSON for", process_id, ":", e$message))
    return(NULL)
  })
}
# Enhanced data loading with proper validation
load_plate_data <- function(plate_id, target_process_ids) {
  result <- list(
    summary_stats = data.frame(),
    fasta_compare = data.frame(),
    json_data = data.frame(),
    status = "",
    debug_info = list()
  )
  
  # Track actual records found for target process IDs
  records_found <- list(
    bgee_records = 0,
    fasta_records = 0, 
    json_records = 0
  )
  
  process_ids_found <- list(
    bgee_ids = character(),
    fasta_ids = character(),
    json_ids = character()
  )
  
  message(paste("Loading data for plate:", plate_id))
  message(paste("Target Process IDs:", paste(target_process_ids, collapse = ", ")))
  
  # Load BGEE summary stats
  if (dir.exists(BASE_DIRS$bgee_summary)) {
    csv_files <- list.files(BASE_DIRS$bgee_summary, pattern = "\\.csv$", full.names = TRUE)
    all_summary_data <- data.frame()
    
    message(paste("Found", length(csv_files), "CSV files in bgee_summary directory"))
    
    for (file in csv_files) {
      tryCatch({
        data <- read.csv(file, stringsAsFactors = FALSE)
        
        # The BGEE files use "ID" column based on your sample
        if ("ID" %in% names(data)) {
          names(data)[names(data) == "ID"] <- "process_id"
        } else if ("Process.ID" %in% names(data)) {
          names(data)[names(data) == "Process.ID"] <- "process_id"
        } else if ("Process ID" %in% names(data)) {
          names(data)[names(data) == "Process ID"] <- "process_id"
        }
        
        if ("process_id" %in% names(data)) {
          all_summary_data <- rbind(all_summary_data, data)
        } else {
          message(paste("No ID/process_id column found in", basename(file), "- columns:", paste(names(data), collapse = ", ")))
        }
      }, error = function(e) {
        message(paste("Error reading", file, ":", e$message))
      })
    }
    
    # Filter for target process IDs and count actual records
    if (nrow(all_summary_data) > 0 && "process_id" %in% names(all_summary_data)) {
      result$summary_stats <- all_summary_data[all_summary_data$process_id %in% target_process_ids, ]
      records_found$bgee_records <- nrow(result$summary_stats)
      
      if (records_found$bgee_records > 0) {
        process_ids_found$bgee_ids <- unique(result$summary_stats$process_id)
        message(paste("BGEE Summary: Found", records_found$bgee_records, "records for", length(process_ids_found$bgee_ids), "process IDs"))
      } else {
        message("BGEE Summary: No records found for target process IDs")
      }
    }
  } else {
    message("BGEE summary directory does not exist")
  }
  
  # Load FASTA compare data
  if (dir.exists(BASE_DIRS$fasta_compare)) {
    csv_files <- list.files(BASE_DIRS$fasta_compare, pattern = "\\.csv$", full.names = TRUE)
    all_fasta_data <- data.frame()
    
    message(paste("Found", length(csv_files), "CSV files in fasta_compare directory"))
    
    for (file in csv_files) {
      tryCatch({
        data <- read.csv(file, stringsAsFactors = FALSE)
        
        # FASTA compare files use "process_id" column based on your sample
        if ("process_id" %in% names(data)) {
          all_fasta_data <- rbind(all_fasta_data, data)
        } else if ("Process.ID" %in% names(data)) {
          names(data)[names(data) == "Process.ID"] <- "process_id"
          all_fasta_data <- rbind(all_fasta_data, data)
        } else if ("Process ID" %in% names(data)) {
          names(data)[names(data) == "Process ID"] <- "process_id"
          all_fasta_data <- rbind(all_fasta_data, data)
        } else if ("ID" %in% names(data)) {
          names(data)[names(data) == "ID"] <- "process_id"
          all_fasta_data <- rbind(all_fasta_data, data)
        } else {
          message(paste("No process_id column found in", basename(file), "- columns:", paste(names(data), collapse = ", ")))
        }
      }, error = function(e) {
        message(paste("Error reading", file, ":", e$message))
      })
    }
    
    # Filter for target process IDs and count actual records
    if (nrow(all_fasta_data) > 0 && "process_id" %in% names(all_fasta_data)) {
      result$fasta_compare <- all_fasta_data[all_fasta_data$process_id %in% target_process_ids, ]
      records_found$fasta_records <- nrow(result$fasta_compare)
      
      if (records_found$fasta_records > 0) {
        process_ids_found$fasta_ids <- unique(result$fasta_compare$process_id)
        message(paste("FASTA Compare: Found", records_found$fasta_records, "records for", length(process_ids_found$fasta_ids), "process IDs"))
      } else {
        message("FASTA Compare: No records found for target process IDs")
      }
    }
  } else {
    message("FASTA compare directory does not exist")
  }
  
  # Load JSON data with corrected filename parsing
  if (dir.exists(BASE_DIRS$fastp_json)) {
    json_files <- list.files(BASE_DIRS$fastp_json, pattern = "\\.json$", full.names = TRUE)
    json_data_list <- list()
    
    message(paste("Found", length(json_files), "JSON files in fastp_json directory"))
    
    for (json_file in json_files) {
      # Extract process ID from filename: [ProcessID]_fastp_report.json
      filename <- basename(json_file)
      
      # Remove _fastp_report.json suffix to get process ID
      if (grepl("_fastp_report\\.json$", filename)) {
        process_id <- gsub("_fastp_report\\.json$", "", filename)
      } else {
        # Fallback: remove .json extension
        process_id <- gsub("\\.json$", "", filename)
      }
      
      # Only process if this process_id is in our target list
      if (process_id %in% target_process_ids) {
        parsed_data <- parse_json_file(json_file, process_id)
        if (!is.null(parsed_data)) {
          json_data_list[[length(json_data_list) + 1]] <- parsed_data
          process_ids_found$json_ids <- c(process_ids_found$json_ids, process_id)
          message(paste("Successfully parsed JSON for process ID:", process_id))
        } else {
          message(paste("Failed to parse JSON for process ID:", process_id))
        }
      }
    }
    
    if (length(json_data_list) > 0) {
      result$json_data <- bind_rows(json_data_list)
      records_found$json_records <- nrow(result$json_data)
      message(paste("JSON Data: Successfully loaded", records_found$json_records, "records"))
    } else {
      message("JSON Data: No JSON files found for target process IDs")
    }
  } else {
    message("Fastp JSON directory does not exist")
  }
  
  # Generate status based on ACTUAL RECORDS FOUND
  total_records_found <- records_found$bgee_records + records_found$fasta_records + records_found$json_records
  
  if (total_records_found == 0) {
    # No actual data records found
    result$status <- paste("ERROR: No data found for plate", plate_id, 
                           "- Searched", length(target_process_ids), "Process IDs but found 0 records")
  } else {
    # Build success message with actual record counts
    loaded_items <- c()
    
    if (records_found$bgee_records > 0) {
      loaded_items <- c(loaded_items, paste("BGEE:", records_found$bgee_records, "records"))
    }
    
    if (records_found$fasta_records > 0) {
      loaded_items <- c(loaded_items, paste("FASTA:", records_found$fasta_records, "records"))
    }
    
    if (records_found$json_records > 0) {
      loaded_items <- c(loaded_items, paste("JSON:", records_found$json_records, "records"))
    }
    
    result$status <- paste("Successfully loaded:", paste(loaded_items, collapse = ", "))
    
    # Add process ID summary
    unique_process_ids_found <- unique(c(process_ids_found$bgee_ids, process_ids_found$fasta_ids, process_ids_found$json_ids))
    result$status <- paste(result$status, "- Found data for", length(unique_process_ids_found), "of", length(target_process_ids), "Process IDs")
    
    # Add warning if some data types are completely missing
    missing_types <- c()
    if (records_found$bgee_records == 0) missing_types <- c(missing_types, "BGEE Summary")
    if (records_found$fasta_records == 0) missing_types <- c(missing_types, "FASTA Compare")
    if (records_found$json_records == 0) missing_types <- c(missing_types, "Fastp JSON")
    
    if (length(missing_types) > 0) {
      result$status <- paste(result$status, "- WARNING: Missing data types:", paste(missing_types, collapse = ", "))
    }
  }
  
  # Store debug info
  result$debug_info <- list(
    records_found = records_found,
    process_ids_found = process_ids_found,
    target_process_ids = target_process_ids
  )
  
  return(result)
}

# Enhanced build_plate_index with validation
build_plate_index <- function() {
  # Load metadata
  metadata <- load_sample_metadata()
  
  if (nrow(metadata) == 0) {
    message("No metadata available - cannot build plate index")
    return(data.frame(
      plate_id = character(),
      process_ids = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  message(paste("Building plate index from", nrow(metadata), "sample records"))
  
  # Extract plate IDs from Sample IDs
  metadata$plate_id <- extract_plate_id_from_sample(metadata$sample_id)
  
  # Show extraction results for debugging
  valid_plates <- sum(!is.na(metadata$plate_id))
  invalid_samples <- sum(is.na(metadata$plate_id))
  
  message(paste("Plate ID extraction results:"))
  message(paste("- Valid plates extracted:", valid_plates))
  message(paste("- Invalid/excluded samples:", invalid_samples))
  
  if (invalid_samples > 0) {
    invalid_sample_ids <- metadata$sample_id[is.na(metadata$plate_id)]
    message(paste("- Invalid sample IDs:", paste(head(invalid_sample_ids, 5), collapse = ", "), 
                  if(length(invalid_sample_ids) > 5) "..." else ""))
  }
  
  # Remove rows where plate_id extraction failed
  metadata <- metadata[!is.na(metadata$plate_id) & metadata$plate_id != "", ]
  
  if (nrow(metadata) == 0) {
    message("No valid plate IDs found in sample metadata")
    return(data.frame(
      plate_id = character(),
      process_ids = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Group process IDs by plate
  plate_data <- metadata %>%
    group_by(plate_id) %>%
    summarise(
      process_ids = paste(process_id, collapse = ","),
      sample_count = n(),
      .groups = 'drop'
    ) %>%
    as.data.frame()
  
  # Show plate summary
  message(paste("Built plate index with", nrow(plate_data), "plates:"))
  for (i in 1:min(nrow(plate_data), 5)) {
    plate_info <- plate_data[i, ]
    process_count <- length(unlist(strsplit(plate_info$process_ids, ",")))
    message(paste("-", plate_info$plate_id, ":", process_count, "process IDs"))
  }
  if (nrow(plate_data) > 5) {
    message("- ... and more")
  }
  
  return(plate_data %>% select(plate_id, process_ids))
}



# Smart number formatting function
smart_format_number <- function(x) {
  if (is.na(x) || !is.numeric(x)) return(as.character(x))
  
  rounded <- round(x, 2)
  
  if (rounded == round(rounded, 0)) {
    return(as.character(as.integer(rounded)))
  } else {
    return(sprintf("%.2f", rounded))
  }
}

# Define UI
ui <- fluidPage(
  useShinyjs(),
  tags$head(
    tags$style(HTML("
      .dataTables_wrapper { font-size: 14px; }
      .download-btn { margin: 10px 0; }
      .search-section {
        background-color: #f8f9fa;
        padding: 20px;
        border-radius: 8px;
        margin-bottom: 20px;
      }
      .status-message {
        margin-top: 15px;
        padding: 10px;
        border-radius: 4px;
      }
      .status-success {
        background-color: #d4edda;
        color: #155724;
        border: 1px solid #c3e6cb;
      }
      .status-error {
        background-color: #f8d7da;
        color: #721c24;
        border: 1px solid #f5c6cb;
      }
      .progress-container {
        background-color: white;
        padding: 20px;
        border-radius: 8px;
        border: 1px solid #dee2e6;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .bold-link-info {
        background-color: #e3f2fd;
        padding: 12px 15px;
        border-radius: 6px;
        margin-bottom: 15px;
        border-left: 4px solid #1976d2;
        font-size: 0.9em;
        color: #1565c0;
      }
      .bold-link-info a {
        color: #1976d2;
        text-decoration: none;
        font-weight: 500;
      }
      .bold-link-info a:hover {
        text-decoration: underline;
      }
      .proportions-control {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        border-left: 4px solid #007bff;
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('updateProgress', function(data) {
        $('#progress_bar').css('width', data.percent + '%');
        $('#progress_text').text(data.message + ' (' + Math.round(data.percent) + '%)');
      });
    "))
  ),
  
  titlePanel("Barcoding Results Dashboard"),
  
  tabsetPanel(
    # Search Results Tab
    tabPanel("Search Results",
             br(),
             div(class = "search-section",
                 h3("Search by Plate ID"),
                 p("Select a plate to view results for all Process IDs in that plate:"),
                 fluidRow(
                   column(6,
                          selectInput("plate_search", 
                                      "Available Plates:",
                                      choices = c("Loading plates..." = ""),
                                      width = "100%")
                   ),
                   column(6,
                          br(),
                          actionButton("load_plate", "Load Plate Results", 
                                       class = "btn-primary btn-lg",
                                       style = "margin-top: 5px;")
                   )
                 ),
                 
                 div(id = "loading_section", style = "display: none; margin-top: 20px;",
                     div(class = "progress-container",
                         h4("Loading Plate Data...", style = "color: #495057; margin-bottom: 15px;"),
                         div(class = "progress", style = "height: 25px; background-color: #e9ecef; border-radius: 12px; overflow: hidden;",
                             div(id = "progress_bar", class = "progress-bar progress-bar-striped progress-bar-animated", 
                                 style = "width: 0%; background-color: #007bff; transition: width 0.3s ease;")
                         ),
                         div(id = "progress_text", style = "text-align: center; margin-top: 10px; font-weight: 500; color: #495057;",
                             "Initializing...")
                     )
                 ),
                 div(id = "plate_status", class = "status-message")
             ),
             br(),
             h4("About the Data"),
             p("Once you select and load a plate, you can navigate to the other tabs to view:"),
             tags$ul(
               tags$li(strong("BGEE Summary Statistics:"), " Combined statistics from the Barcode Gene Extractor & Evaluator (BGEE) pipeline."),
               tags$li(strong("Fastp Metrics:"), " Before and after read trimming statistics from Fastp (parsed from json files)."),
               tags$li(strong("FASTA Compare Results:"), " Barcode consensus sequence comparison and quality metrics.")
             )
    ),
    
    # Barcoding Outcome Tab
    tabPanel("Barcoding Outcome",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               # Static description text area
               div(
                 h3("🚦 Barcoding Outcome Overview"),
                 p("Overview of BIN-compliant barcode extraction results for all Process IDs in the loaded plate."),
                 style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;"
               ),
               
               # BOLD Systems link information
               div(class = "bold-link-info",
                   p("🔗 ", strong("BOLD Systems Integration:"), " Process IDs in the table below are clickable links that will take you directly to the corresponding record in the ",
                     a("BOLD Systems database", href = "https://portal.boldsystems.org", target = "_blank"), 
                     ". Click on any Process ID to view detailed specimen and sequence information.",
                     style = "margin: 0;")
               ),
               
               # Color coding explanation
               div(
                 h4("Result Categories", style = "color: #495057; margin-bottom: 15px;"),
                 div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin-bottom: 20px;",
                     div(style = "background-color: #d4edda; padding: 15px; border-radius: 8px; border-left: 4px solid #28a745;",
                         div(style = "font-weight: bold; color: #155724; margin-bottom: 5px;", "🟢 Green - Success"),
                         p("BIN-compliant barcode successfully extracted (selected_barcode_fasta = TRUE)", 
                           style = "margin: 0; color: #155724; font-size: 0.9em;")
                     ),
                     div(style = "background-color: #fff3cd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffc107;",
                         div(style = "font-weight: bold; color: #856404; margin-bottom: 5px;", "🟡 Amber - Partial"),
                         p("Barcode rank 4+", 
                           style = "margin: 0; color: #856404; font-size: 0.9em;")
                     ),
                     div(style = "background-color: #f8d7da; padding: 15px; border-radius: 8px; border-left: 4px solid #dc3545;",
                         div(style = "font-weight: bold; color: #721c24; margin-bottom: 5px;", "🔴 Red - Failed"),
                         p("No BIN-compliant barcode extracted (no selected_barcode_fasta = TRUE)", 
                           style = "margin: 0; color: #721c24; font-size: 0.9em;")
                     )
                 ),
                 style = "background-color: white; padding: 20px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;"
               ),
               
               # Barcode rank explanation
               div(
                 h4("Barcode Rank Explanation", style = "color: #495057; margin-bottom: 15px;"),
                 tags$ul(
                   tags$li(strong("Rank 1:"), " No ambiguous bases, longest stretch ≥ 650"),
                   tags$li(strong("Rank 2:"), " No ambiguous bases, longest stretch ≥ 500"),
                   tags$li(strong("Rank 3:"), " No ambiguous bases, 300 ≤ longest stretch ≤ 499"),
                   tags$li(strong("Rank 4:"), " No ambiguous bases, 1 ≤ longest stretch ≤ 299"),
                   tags$li(strong("Rank 5:"), " Has ambiguous bases"),
                   style = "margin: 0; color: #495057;"
                 ),
                 p("Lower ranks (1-3) indicate higher quality barcode sequences suitable for BOLD.", 
                   style = "margin-top: 15px; margin-bottom: 0; color: #6c757d; font-style: italic;"),
                 style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; border-left: 4px solid #007bff; margin-bottom: 20px;"
               ),
               
               # Download button
               downloadButton("downloadOutcomeTable", "Download Outcome Table", class = "download-btn btn-primary"),
               
               # Outcome table
               DTOutput("outcomeTable"),
               
               # Summary statistics
               br(),
               div(id = "outcome_summary_container",
                   style = "margin-top: 20px;")
             )
    ),
    
    # Fastp Metrics Tab
    tabPanel("Fastp Metrics",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("About Fastp Metrics"),
                 HTML("Fastp metrics are generated from JSON output files of the <a href='https://github.com/OpenGene/fastp' target='_blank'>fastp preprocessing tool</a>. These files contain detailed statistics about read quality before and after filtering."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               downloadButton("downloadJsonData", "Download Filtered JSON Data", class = "download-btn btn-primary"),
               DTOutput("jsonDataTable"),
               
               hr(),
               h4("Process ID Overview"),
               fluidRow(
                 column(4, selectInput("json_bar_ycol", "Y Axis:", choices = NULL)),
                 column(4, div(class = "proportions-control",
                               checkboxInput("show_proportions", "Show as proportions (%)", value = FALSE),
                               style = "margin-top: 5px;")),
                 column(4, div(style = "height: 34px;"))
               ),
               plotlyOutput("jsonBarPlot", height = "400px"),
               
               hr(),
               h4("Interactive Scatter Plot"),
               fluidRow(
                 column(6, selectInput("json_xcol", "X Axis:", choices = NULL)),
                 column(6, selectInput("json_ycol", "Y Axis:", choices = NULL))
               ),
               plotlyOutput("jsonDataPlot", height = "500px")
             )
    ),
    
    # BGEE Summary Statistics Tab
    tabPanel("BGEE Summary Statistics",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("About BGEE Summary Statistics"),
                 HTML("<a href='https://github.com/bge-barcoding/BGEE' target='_blank'>BGEE (Barcode Gene Extractor & Evaluator)</a> summary statistics are generated from the combined output of the MGE and fasta_cleaner steps."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               downloadButton("downloadSummaryStats", "Download Filtered Summary Stats", class = "download-btn btn-primary"),
               DTOutput("summaryStatsTable"),
               
               hr(),
               h4("Interactive Scatter Plot"),
               fluidRow(
                 column(6, selectInput("ycol", "Y Axis:", choices = NULL)),
                 column(6, selectInput("xcol", "X Axis:", choices = NULL))
               ),
               plotlyOutput("summaryStatsPlot", height = "500px")
             )
    ),
    
    # FASTA Compare Results Tab
    tabPanel("FASTA Compare Results",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("About FASTA Compare Results"),
                 HTML("<a href='https://github.com/bge-barcoding/fasta_compare' target='_blank'>FASTA Compare</a> is the final post-processing step of the BGEE workflow."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               downloadButton("downloadFastaCompare", "Download Filtered FASTA Compare", class = "download-btn btn-primary"),
               DTOutput("fastaCompareTable"),
               
               hr(),
               h4("Interactive Plot"),
               fluidRow(
                 column(6, selectInput("fasta_xcol", "X Axis:", choices = NULL)),
                 column(6, selectInput("fasta_ycol", "Y Axis:", choices = NULL))
               ),
               plotlyOutput("fastaComparePlot", height = "500px")
             )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  values <- reactiveValues(
    summary_stats = data.frame(),
    fasta_compare = data.frame(),
    json_data = data.frame(),
    dataset_loaded = FALSE,
    plate_index = data.frame(),
    current_plate = "",
    current_process_ids = character()
  )
  
  # Initialize plate index
  observe({
    values$plate_index <- build_plate_index()
    
    if (nrow(values$plate_index) > 0) {
      available_plates <- sort(unique(values$plate_index$plate_id))
      choices <- setNames(available_plates, available_plates)
      choices <- c("Select a plate..." = "", choices)
    } else {
      choices <- c("No plates found" = "")
    }
    
    updateSelectInput(session, "plate_search", choices = choices)
  })
  
  # Create output for conditionalPanel
  output$dataset_loaded <- reactive({
    values$dataset_loaded
  })
  outputOptions(output, "dataset_loaded", suspendWhenHidden = FALSE)
  
  # Load plate data when button is clicked
  observeEvent(input$load_plate, {
    if (input$plate_search == "" || is.null(input$plate_search)) {
      showNotification("Please select a plate first.", type = "warning")
      return()
    }
    
    # Get process IDs for this plate
    plate_row <- values$plate_index[values$plate_index$plate_id == input$plate_search, ]
    if (nrow(plate_row) == 0) {
      showNotification("No data found for this plate.", type = "error")
      return()
    }
    
    target_process_ids <- unlist(strsplit(plate_row$process_ids, ","))
    
    # Disable button and clear status
    updateActionButton(session, "load_plate", label = "Loading...", icon = NULL)
    shinyjs::disable("load_plate")
    removeUI("#plate_status div")
    
    tryCatch({
      # Load data first (without showing progress)
      result <- load_plate_data(input$plate_search, target_process_ids)
      
      # Check if any actual data was found
      total_records <- sum(unlist(result$debug_info$records_found))
      missing_data_types <- c()
      if (result$debug_info$records_found$bgee_records == 0) missing_data_types <- c(missing_data_types, "BGEE")
      if (result$debug_info$records_found$fasta_records == 0) missing_data_types <- c(missing_data_types, "FASTA")
      if (result$debug_info$records_found$json_records == 0) missing_data_types <- c(missing_data_types, "JSON")
      
      if (total_records == 0 || length(missing_data_types) > 0) {
        # ERROR CASE: No data found - skip loading animation entirely
        values$summary_stats <- data.frame()
        values$fasta_compare <- data.frame()
        values$json_data <- data.frame()
        values$current_plate <- ""
        values$current_process_ids <- character()
        values$dataset_loaded <- FALSE
        
        error_msg <- if (total_records == 0) {
          paste("ERROR: No data files found for plate", input$plate_search, 
                "- Searched", length(target_process_ids), "Process IDs but found 0 records")
        } else {
          paste("ERROR: Missing data types for plate", input$plate_search, 
                "- Missing:", paste(missing_data_types, collapse = ", "))
        }
        
        insertUI("#plate_status", "beforeEnd",
                 div(class = "status-message status-error",
                     p(error_msg, style = "margin: 0;")))
        
        showNotification("Failed to load plate - no data found.", type = "error")
        
      } else {
        # SUCCESS CASE: Data found - show loading animation and process
        shinyjs::show("loading_section")
        
        # Simulate progress updates
        session$sendCustomMessage("updateProgress", list(percent = 25, message = "Processing BGEE data"))
        Sys.sleep(0.5)
        session$sendCustomMessage("updateProgress", list(percent = 50, message = "Processing FASTA compare data"))
        Sys.sleep(0.5)
        session$sendCustomMessage("updateProgress", list(percent = 75, message = "Processing Fastp JSON data"))
        Sys.sleep(0.5)
        session$sendCustomMessage("updateProgress", list(percent = 100, message = "Complete"))
        
        # Hide loading and update values
        shinyjs::hide("loading_section")
        
        values$summary_stats <- result$summary_stats
        values$fasta_compare <- result$fasta_compare
        values$json_data <- result$json_data
        values$current_plate <- input$plate_search
        values$current_process_ids <- target_process_ids
        values$dataset_loaded <- TRUE
        
        insertUI("#plate_status", "beforeEnd",
                 div(class = "status-message status-success",
                     p(paste("Successfully loaded plate:", input$plate_search, 
                             "- Process IDs:", length(target_process_ids)), style = "margin: 0;")))
        
        showNotification(paste("Plate", input$plate_search, "loaded successfully!"), type = "message")
      }
      
    }, error = function(e) {
      shinyjs::hide("loading_section")
      
      insertUI("#plate_status", "beforeEnd",
               div(class = "status-message status-error",
                   p(paste("Error loading plate:", e$message), style = "margin: 0;")))
      
      showNotification("Failed to load plate. Please try again.", type = "error")
    })
    
    # Re-enable button
    updateActionButton(session, "load_plate", label = "Load Plate Results", icon = NULL)
    shinyjs::enable("load_plate")
  })
  
  # Prepare outcome data
  prepare_outcome_data <- reactive({
    req(values$dataset_loaded, nrow(values$fasta_compare) > 0)
    
    if (!"process_id" %in% names(values$fasta_compare)) {
      return(data.frame())
    }
    
    # Group by process_id and get outcome information
    outcome_data <- values$fasta_compare %>%
      group_by(process_id) %>%
      summarise(
        has_selected_barcode = any(tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
        
        selected_barcode_longest_stretch = ifelse(
          any(tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          barcode_longest_stretch[tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t")][1],
          NA
        ),
        
        selected_barcode_rank = ifelse(
          any(tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          barcode_rank[tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t")][1],
          NA
        ),
        
        best_barcode_rank = ifelse(
          any(tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          barcode_rank[tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t")][1],
          min(barcode_rank, na.rm = TRUE)
        ),
        
        fallback_barcode_longest_stretch = ifelse(
          any(tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          barcode_longest_stretch[tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t")][1],
          barcode_longest_stretch[barcode_rank == min(barcode_rank, na.rm = TRUE)][1]
        ),
        
        .groups = 'drop'
      ) %>%
      mutate(
        best_barcode_rank = ifelse(is.infinite(best_barcode_rank), NA, best_barcode_rank),
        selected_barcode_longest_stretch = ifelse(is.na(selected_barcode_longest_stretch), 0, selected_barcode_longest_stretch),
        fallback_barcode_longest_stretch = ifelse(is.na(fallback_barcode_longest_stretch), 0, fallback_barcode_longest_stretch)
      )
    
    # Determine outcome category
    outcome_data <- outcome_data %>%
      mutate(
        outcome_category = case_when(
          has_selected_barcode & (is.na(selected_barcode_rank) | selected_barcode_rank <= 3) ~ "success",
          (!is.na(selected_barcode_rank) & selected_barcode_rank >= 4) | 
            (!is.na(best_barcode_rank) & best_barcode_rank >= 4) ~ "partial",
          !has_selected_barcode & (!is.na(best_barcode_rank) & best_barcode_rank < 4) ~ "failed",
          TRUE ~ "failed"
        ),
        
        outcome_status = case_when(
          outcome_category == "success" ~ "✅ Success",
          outcome_category == "partial" ~ "⚠️ Partial",
          TRUE ~ "❌ Failed"
        ),
        
        display_barcode_stretch = case_when(
          has_selected_barcode ~ as.character(selected_barcode_longest_stretch),
          !has_selected_barcode ~ paste0("(", fallback_barcode_longest_stretch, ")"),
          TRUE ~ NA_character_
        ),
        
        display_barcode_rank = case_when(
          has_selected_barcode ~ as.character(selected_barcode_rank),
          !has_selected_barcode ~ paste0("(", best_barcode_rank, ")"),
          TRUE ~ NA_character_
        ),
        
        # Create BOLD link for Process ID
        process_id_link = paste0('<a href="https://portal.boldsystems.org/record/', process_id, '" target="_blank">', process_id, '</a>')
      ) %>%
      select(
        `Process ID` = process_id_link,
        `Status` = outcome_status,
        `Barcode Longest Stretch` = display_barcode_stretch,
        `Best Barcode Rank` = display_barcode_rank,
        outcome_category
      )
    
    return(outcome_data)
  })
  
  # Render outcome table
  output$outcomeTable <- renderDT({
    outcome_data <- prepare_outcome_data()
    
    if (nrow(outcome_data) == 0) {
      return(datatable(data.frame("No data available" = "Please check that FASTA compare data is loaded correctly.")))
    }
    
    display_data <- outcome_data %>% select(-outcome_category)
    numeric_columns <- which(sapply(display_data, is.numeric))
    
    datatable(
      display_data,
      escape = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-center', targets = c(1)),
          list(className = 'dt-right', targets = numeric_columns - 1)
        ),
        rowCallback = JS(
          "function(row, data, index) {",
          "  var status = data[1];",
          "  if (status.includes('Success')) {",
          "    $(row).css('background-color', '#d4edda');",
          "  } else if (status.includes('Partial')) {",
          "    $(row).css('background-color', '#fff3cd');", 
          "  } else if (status.includes('Failed')) {",
          "    $(row).css('background-color', '#f8d7da');",
          "  }",
          "}"
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    ) %>%
      {if(length(numeric_columns) > 0) formatRound(., columns = numeric_columns, digits = 2) else .}
  })
  
  # Render outcome summary
  output$outcome_summary <- renderUI({
    outcome_data <- prepare_outcome_data()
    
    if (nrow(outcome_data) == 0) {
      return(div())
    }
    
    total_processes <- nrow(outcome_data)
    success_count <- sum(outcome_data$outcome_category == "success", na.rm = TRUE)
    partial_count <- sum(outcome_data$outcome_category == "partial", na.rm = TRUE)
    failed_count <- sum(outcome_data$outcome_category == "failed", na.rm = TRUE)
    
    success_rate <- round((success_count / total_processes) * 100, 1)
    
    div(
      h4("Summary Statistics", style = "color: #495057; margin-bottom: 15px;"),
      div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 15px;",
          div(class = "outcome-stat", style = "background-color: white; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #495057;", total_processes),
              div(style = "color: #6c757d; margin-top: 5px;", "Total Process IDs")
          ),
          div(class = "outcome-stat", style = "background-color: #d4edda; padding: 15px; border-radius: 8px; border: 1px solid #c3e6cb; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #155724;", success_count),
              div(style = "color: #155724; margin-top: 5px;", "Successful")
          ),
          div(class = "outcome-stat", style = "background-color: #fff3cd; padding: 15px; border-radius: 8px; border: 1px solid #ffeaa7; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #856404;", partial_count),
              div(style = "color: #856404; margin-top: 5px;", "Partial")
          ),
          div(class = "outcome-stat", style = "background-color: #f8d7da; padding: 15px; border-radius: 8px; border: 1px solid #f5c6cb; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #721c24;", failed_count),
              div(style = "color: #721c24; margin-top: 5px;", "Failed")
          ),
          div(class = "outcome-stat", style = "background-color: #e3f2fd; padding: 15px; border-radius: 8px; border: 1px solid #bbdefb; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #1976d2;", paste0(success_rate, "%")),
              div(style = "color: #1976d2; margin-top: 5px;", "Success Rate")
          )
      ),
      style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-top: 20px;"
    )
  })
  
  # Update the outcome summary when data changes
  observe({
    if (values$dataset_loaded) {
      removeUI("#outcome_summary_container > div")
      insertUI("#outcome_summary_container", "beforeEnd", 
               uiOutput("outcome_summary"))
    }
  })
  
  # Render tables
  output$summaryStatsTable <- renderDT({
    req(values$dataset_loaded, nrow(values$summary_stats) > 0)
    
    numeric_cols <- which(sapply(values$summary_stats, is.numeric))
    
    datatable(
      values$summary_stats,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-right', targets = numeric_cols - 1),
          list(
            targets = numeric_cols - 1,
            render = JS("
              function(data, type, row) {
                if (type === 'display' && data !== null && data !== undefined) {
                  var num = Number(data);
                  if (!isNaN(num)) {
                    if (num % 1 === 0) {
                      return num.toString();
                    } else {
                      return num.toFixed(2);
                    }
                  }
                }
                return data;
              }
            ")
          )
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    )
  })
  
  output$fastaCompareTable <- renderDT({
    req(values$dataset_loaded, nrow(values$fasta_compare) > 0)
    
    numeric_cols <- which(sapply(values$fasta_compare, is.numeric))
    
    datatable(
      values$fasta_compare,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-right', targets = numeric_cols - 1),
          list(
            targets = numeric_cols - 1,
            render = JS("
              function(data, type, row) {
                if (type === 'display' && data !== null && data !== undefined) {
                  var num = Number(data);
                  if (!isNaN(num)) {
                    if (num % 1 === 0) {
                      return num.toString();
                    } else {
                      return num.toFixed(2);
                    }
                  }
                }
                return data;
              }
            ")
          )
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    )
  })
  
  output$jsonDataTable <- renderDT({
    req(values$dataset_loaded, nrow(values$json_data) > 0)
    
    # Filter out unwanted columns for display
    columns_to_remove <- c(
      "before_q20_bases", "before_q30_bases", "before_q20_rate", "before_q30_rate", 
      "before_read1_mean_length", "before_read2_mean_length", "before_gc_content",
      "after_q20_bases", "after_q30_bases", "after_q20_rate", "after_q30_rate", 
      "too_long_reads"
    )
    
    display_data <- values$json_data %>%
      select(-any_of(columns_to_remove))
    
    numeric_cols <- which(sapply(display_data, is.numeric))
    
    datatable(
      display_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-right', targets = numeric_cols - 1),
          list(
            targets = numeric_cols - 1,
            render = JS("
              function(data, type, row) {
                if (type === 'display' && data !== null && data !== undefined) {
                  var num = Number(data);
                  if (!isNaN(num)) {
                    if (num % 1 === 0) {
                      return num.toString();
                    } else {
                      return num.toFixed(2);
                    }
                  }
                }
                return data;
              }
            ")
          )
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    )
  })
  
  # Update column selectors
  observe({
    if (values$dataset_loaded && nrow(values$summary_stats) > 0) {
      numeric_cols <- names(values$summary_stats)[sapply(values$summary_stats, is.numeric)]
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
        updateSelectInput(session, "xcol", choices = numeric_cols, selected = numeric_cols[1])
      }
    }
  })
  
  observe({
    if (values$dataset_loaded && nrow(values$fasta_compare) > 0) {
      numeric_cols <- names(values$fasta_compare)[sapply(values$fasta_compare, is.numeric)]
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "fasta_xcol", choices = numeric_cols, selected = numeric_cols[1])
        updateSelectInput(session, "fasta_ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
      }
    }
  })
  
  observe({
    if (values$dataset_loaded && nrow(values$json_data) > 0) {
      columns_to_remove <- c(
        "before_q20_bases", "before_q30_bases", "before_q20_rate", "before_q30_rate", 
        "before_read1_mean_length", "before_read2_mean_length", "before_gc_content",
        "after_q20_bases", "after_q30_bases", "after_q20_rate", "after_q30_rate", 
        "too_long_reads"
      )
      
      available_cols <- names(values$json_data)[!names(values$json_data) %in% columns_to_remove]
      numeric_cols <- available_cols[sapply(values$json_data[available_cols], is.numeric)]
      
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "json_xcol", choices = numeric_cols, selected = numeric_cols[1])
        updateSelectInput(session, "json_ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
        
        bar_choices <- available_cols[available_cols != "process_id"]
        updateSelectInput(session, "json_bar_ycol", choices = bar_choices, selected = bar_choices[1])
      }
    }
  })
  
  # Render plots
  output$summaryStatsPlot <- renderPlotly({
    req(input$xcol, input$ycol, values$dataset_loaded, nrow(values$summary_stats) > 0)
    
    plot_ly(
      data = values$summary_stats,
      x = ~.data[[input$xcol]],
      y = ~.data[[input$ycol]],
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 10, color = 'rgba(51, 122, 183, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)'))
    ) %>%
      layout(
        xaxis = list(title = input$xcol),
        yaxis = list(title = input$ycol),
        margin = list(t = 30)
      )
  })
  
  output$fastaComparePlot <- renderPlotly({
    req(input$fasta_xcol, input$fasta_ycol, values$dataset_loaded, nrow(values$fasta_compare) > 0)
    
    plot_ly(
      data = values$fasta_compare,
      x = ~.data[[input$fasta_xcol]],
      y = ~.data[[input$fasta_ycol]],
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 10, color = 'rgba(0, 180, 135, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)'))
    ) %>%
      layout(
        xaxis = list(title = input$fasta_xcol),
        yaxis = list(title = input$fasta_ycol),
        margin = list(t = 30)
      )
  })
  
  output$jsonDataPlot <- renderPlotly({
    req(input$json_xcol, input$json_ycol, values$dataset_loaded, nrow(values$json_data) > 0)
    
    plot_ly(
      data = values$json_data,
      x = ~.data[[input$json_xcol]],
      y = ~.data[[input$json_ycol]],
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 10, color = 'rgba(255, 99, 71, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)'))
    ) %>%
      layout(
        xaxis = list(title = input$json_xcol),
        yaxis = list(title = input$json_ycol),
        margin = list(t = 30)
      )
  })
  
  output$jsonBarPlot <- renderPlotly({
    req(input$json_bar_ycol, values$dataset_loaded, nrow(values$json_data) > 0)
    
    if (!"process_id" %in% names(values$json_data) || !input$json_bar_ycol %in% names(values$json_data)) {
      return(NULL)
    }
    
    plot_data <- values$json_data %>%
      select(process_id, !!sym(input$json_bar_ycol)) %>%
      distinct() %>%
      arrange(process_id) %>%
      filter(!is.na(!!sym(input$json_bar_ycol)))
    
    if (nrow(plot_data) == 0) {
      return(plot_ly() %>% 
               layout(title = "No data available for selected column"))
    }
    
    y_values <- plot_data[[input$json_bar_ycol]]
    y_title <- input$json_bar_ycol
    
    if (isTRUE(input$show_proportions)) {
      max_val <- max(y_values, na.rm = TRUE)
      if (max_val > 0 && is.finite(max_val)) {
        y_values <- (y_values / max_val) * 100
        y_title <- paste(input$json_bar_ycol, "(%)")
      }
    }
    
    plot_ly(
      data = plot_data,
      x = ~factor(process_id),
      y = y_values,
      type = 'bar',
      marker = list(color = 'rgba(255, 99, 71, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)'))
    ) %>%
      layout(
        xaxis = list(title = "Process ID", showticklabels = FALSE),
        yaxis = list(title = y_title),
        margin = list(t = 30, b = 50)
      )
  })
  
  # Download handlers
  output$downloadOutcomeTable <- downloadHandler(
    filename = function() {
      paste0("barcoding_outcome_", values$current_plate, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      outcome_data <- prepare_outcome_data()
      if (nrow(outcome_data) > 0) {
        export_data <- outcome_data %>% 
          select(-outcome_category) %>%
          mutate(`Process ID` = gsub('<.*?>', '', `Process ID`))
        write.csv(export_data, file, row.names = FALSE)
      } else {
        write.csv(data.frame("No data available"), file, row.names = FALSE)
      }
    }
  )
  
  output$downloadSummaryStats <- downloadHandler(
    filename = function() {
      paste0("summary_stats_", values$current_plate, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$summary_stats, file, row.names = FALSE)
    }
  )
  
  output$downloadFastaCompare <- downloadHandler(
    filename = function() {
      paste0("fasta_compare_", values$current_plate, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$fasta_compare, file, row.names = FALSE)
    }
  )
  
  output$downloadJsonData <- downloadHandler(
    filename = function() {
      paste0("json_data_", values$current_plate, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$json_data, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
