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
  barcode_validation = "./data/barcode_validation/",
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
  
  # Also get process IDs from barcode validation files
  if (dir.exists(BASE_DIRS$barcode_validation)) {
    tsv_files <- list.files(BASE_DIRS$barcode_validation, pattern = "\\.tsv$", full.names = TRUE)
    for (file in tsv_files) {
      tryCatch({
        data <- read.delim(file, stringsAsFactors = FALSE, sep = "\t")
        if ("group_id" %in% names(data)) {
          process_ids <- data$group_id
        } else if ("process_id" %in% names(data)) {
          process_ids <- data$process_id
        } else if ("Process.ID" %in% names(data)) {
          process_ids <- data$Process.ID
        } else if ("Process ID" %in% names(data)) {
          process_ids <- data$`Process ID`
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
  
  # Also get process IDs from JSON files
  if (dir.exists(BASE_DIRS$fastp_json)) {
    json_files <- list.files(BASE_DIRS$fastp_json, pattern = "\\.json$", full.names = TRUE)
    for (json_file in json_files) {
      filename <- basename(json_file)
      if (grepl("_fastp_report\\.json$", filename)) {
        process_id <- gsub("_fastp_report\\.json$", "", filename)
      } else {
        process_id <- gsub("\\.json$", "", filename)
      }
      all_process_ids <- c(all_process_ids, process_id)
    }
  }
  
  return(unique(all_process_ids[!is.na(all_process_ids) & all_process_ids != ""]))
}

# Get process IDs by project code
get_process_ids_by_project <- function(project_code) {
  all_process_ids <- get_available_process_ids()
  
  # Filter process IDs containing the project code
  matching_ids <- all_process_ids[grepl(project_code, all_process_ids, ignore.case = TRUE)]
  
  message(paste("Found", length(matching_ids), "Process IDs containing project code:", project_code))
  if (length(matching_ids) > 0) {
    message(paste("Sample Process IDs:", paste(head(matching_ids, 5), collapse = ", "), 
                  if(length(matching_ids) > 5) "..." else ""))
  }
  
  return(matching_ids)
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

# Enhanced data loading with proper validation (renamed from load_plate_data)
load_data <- function(identifier, target_process_ids, identifier_type = "plate") {
  result <- list(
    summary_stats = data.frame(),
    barcode_validation = data.frame(),
    json_data = data.frame(),
    status = "",
    debug_info = list()
  )
  
  # Track actual records found for target process IDs
  records_found <- list(
    bgee_records = 0,
    validation_records = 0, 
    json_records = 0
  )
  
  process_ids_found <- list(
    bgee_ids = character(),
    validation_ids = character(),
    json_ids = character()
  )
  
  message(paste("Loading data for", identifier_type, ":", identifier))
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
  
  # Load Barcode Validation data (TSV files)
  if (dir.exists(BASE_DIRS$barcode_validation)) {
    tsv_files <- list.files(BASE_DIRS$barcode_validation, pattern = "\\.tsv$", full.names = TRUE)
    all_validation_data <- data.frame()
    
    message(paste("Found", length(tsv_files), "TSV files in barcode_validation directory"))
    
    for (file in tsv_files) {
      tryCatch({
        data <- read.delim(file, stringsAsFactors = FALSE, sep = "\t")
        
        # Barcode validation files use "group_id" column
        if ("group_id" %in% names(data)) {
          names(data)[names(data) == "group_id"] <- "process_id"
          all_validation_data <- rbind(all_validation_data, data)
        } else if ("process_id" %in% names(data)) {
          all_validation_data <- rbind(all_validation_data, data)
        } else if ("Process.ID" %in% names(data)) {
          names(data)[names(data) == "Process.ID"] <- "process_id"
          all_validation_data <- rbind(all_validation_data, data)
        } else if ("Process ID" %in% names(data)) {
          names(data)[names(data) == "Process ID"] <- "process_id"
          all_validation_data <- rbind(all_validation_data, data)
        } else if ("ID" %in% names(data)) {
          names(data)[names(data) == "ID"] <- "process_id"
          all_validation_data <- rbind(all_validation_data, data)
        } else {
          message(paste("No group_id/process_id column found in", basename(file), "- columns:", paste(names(data), collapse = ", ")))
        }
      }, error = function(e) {
        message(paste("Error reading", file, ":", e$message))
      })
    }
    
    # Filter for target process IDs and count actual records
    if (nrow(all_validation_data) > 0 && "process_id" %in% names(all_validation_data)) {
      result$barcode_validation <- all_validation_data[all_validation_data$process_id %in% target_process_ids, ]
      records_found$validation_records <- nrow(result$barcode_validation)
      
      if (records_found$validation_records > 0) {
        process_ids_found$validation_ids <- unique(result$barcode_validation$process_id)
        message(paste("Barcode Validation: Found", records_found$validation_records, "records for", length(process_ids_found$validation_ids), "process IDs"))
      } else {
        message("Barcode Validation: No records found for target process IDs")
      }
    }
  } else {
    message("Barcode Validation directory does not exist")
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
  total_records_found <- records_found$bgee_records + records_found$validation_records + records_found$json_records
  
  if (total_records_found == 0) {
    # No actual data records found
    result$status <- paste("ERROR: No data found for", identifier_type, identifier, 
                           "- Searched", length(target_process_ids), "Process IDs but found 0 records")
  } else {
    # Build success message with actual record counts
    loaded_items <- c()
    
    if (records_found$bgee_records > 0) {
      loaded_items <- c(loaded_items, paste("BGEE:", records_found$bgee_records, "records"))
    }
    
    if (records_found$validation_records > 0) {
      loaded_items <- c(loaded_items, paste("Validation:", records_found$validation_records, "records"))
    }
    
    if (records_found$json_records > 0) {
      loaded_items <- c(loaded_items, paste("JSON:", records_found$json_records, "records"))
    }
    
    result$status <- paste("Successfully loaded:", paste(loaded_items, collapse = ", "))
    
    # Add process ID summary
    unique_process_ids_found <- unique(c(process_ids_found$bgee_ids, process_ids_found$validation_ids, process_ids_found$json_ids))
    result$status <- paste(result$status, "- Found data for", length(unique_process_ids_found), "of", length(target_process_ids), "Process IDs")
    
    # Add warning if some data types are completely missing
    missing_types <- c()
    if (records_found$bgee_records == 0) missing_types <- c(missing_types, "BGEE Summary")
    if (records_found$validation_records == 0) missing_types <- c(missing_types, "Barcode Validation")
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

# Function to check taxonomic match
check_taxonomic_match <- function(bold_family, taxonomic_id) {
  # Handle missing/empty/NA data
  if (is.na(bold_family) || is.na(taxonomic_id) || 
      bold_family == "" || taxonomic_id == "") {
    return("NO")
  }
  
  # Split taxonomic_id by commas and trim whitespace
  taxonomic_parts <- trimws(unlist(strsplit(as.character(taxonomic_id), ",")))
  
  # Check if bold_family appears as an exact word (case-insensitive)
  bold_family_clean <- trimws(as.character(bold_family))
  
  # Case-insensitive exact match
  if (tolower(bold_family_clean) %in% tolower(taxonomic_parts)) {
    return("YES")
  } else {
    return("NO")
  }
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
      .proportions-control {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 15px;
        border-left: 4px solid #007bff;
      }
      .workflow-image {
        max-width: 100%;
        height: auto;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin: 20px 0;
      }
      .fasta-content {
        background-color: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        padding: 15px;
        font-family: 'Courier New', monospace;
        font-size: 12px;
        max-height: 500px;
        overflow-y: auto;
        white-space: pre-wrap;
        word-break: break-all;
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('updateProgress', function(data) {
        $('#progress_bar').css('width', data.percent + '%');
        $('#progress_text').text(data.message + ' (' + Math.round(data.percent) + '%)');
      });
      
      Shiny.addCustomMessageHandler('copyToClipboard', function(text) {
        navigator.clipboard.writeText(text).then(function() {
          // Show success message
          Shiny.setInputValue('clipboard_success', Math.random());
        }, function(err) {
          // Show error message
          Shiny.setInputValue('clipboard_error', Math.random());
        });
      });
    "))
  ),
  
  titlePanel("Barcoding Results Dashboard"),
  
  tabsetPanel(
    # Search Results Tab
    tabPanel("Search Results",
             br(),
             div(class = "search-section",
                 h3("Search by Plate ID or Project"),
                 p("Select either a plate to view results for all Process IDs in that plate, OR select a project to view results for all Process IDs containing that project code:"),
                 p("Choose one option - either plate or project", style = "font-style: italic; color: #6c757d;"),
                 fluidRow(
                   column(5,
                          h5("Search by Plate"),
                          selectInput("plate_search", 
                                      "Available Plates:",
                                      choices = c("Loading plates... Please wait..." = ""),
                                      width = "100%")
                   ),
                   column(2,
                          div(style = "text-align: center; padding-top: 40px;",
                              strong("OR"))
                   ),
                   column(5,
                          h5("Search by Project"),
                          selectInput("project_search", 
                                      "Available Projects:",
                                      choices = c("Select a project..." = "", 
                                                  "BSNHM" = "BSNHM", "NHMCG" = "NHMCG", "BGETR" = "BGETR", 
                                                  "BSUIO" = "BSUIO", "BGLIB" = "BGLIB", "BSNTN" = "BSNTN", 
                                                  "BGEGR" = "BGEGR", "DTAUT" = "DTAUT", "HMAUT" = "HMAUT", 
                                                  "BGENL" = "BGENL", "DTKNU" = "DTKNU", "BBIOP" = "BBIOP", 
                                                  "BHNHM" = "BHNHM", "UNIFI" = "UNIFI", "DTULO" = "DTULO", 
                                                  "MEAMP" = "MEAMP", "MUSBA" = "MUSBA", "BGSNL" = "BGSNL", 
                                                  "BGSNH" = "BGSNH", "BGEPL" = "BGEPL", "EBGEP" = "EBGEP", 
                                                  "BSCRO" = "BSCRO", "BIOSC" = "BIOSC", "INVBG" = "INVBG", 
                                                  "BCEMI" = "BCEMI", "ILECA" = "ILECA", "ALPFU" = "ALPFU"),
                                      width = "100%")
                   )
                 ),
                 br(),
                 fluidRow(
                   column(12,
                          div(style = "text-align: center;",
                              actionButton("load_data", "Load Results", 
                                           class = "btn-primary btn-lg"))
                   )
                 ),
                 
                 div(id = "loading_section", style = "display: none; margin-top: 20px;",
                     div(class = "progress-container",
                         h4("Loading Data...", style = "color: #495057; margin-bottom: 15px;"),
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
             p("Once you select and load a plate or project, you can navigate to the other tabs to view:"),
             tags$ul(
               tags$li(strong("BGEE Summary Statistics:"), " Combined statistics from the Barcode Gene Extractor & Evaluator (BGEE) pipeline."),
               tags$li(strong("Fastp Metrics:"), " Before and after read trimming statistics from Fastp (parsed from json files)."),
               tags$li(strong("Barcode Validation:"), " Structural and taxonomic validation of barcode consensus sequences output by BGEE pipeline.")
             ),
             
             br(),
             h4("Barcode Gene Extractor & Evaluator (BGEE) Workflow"),
             div(
               p("The BGEE (Barcode Gene Extractor & Evaluator) workflow is a comprehensive pipeline for recovery of high-quality barcode sequences from raw genome skim sequencing data derived from museum specimens. The workflow includes several quality control steps, consensus sequence generation, and validation processes to ensure reliable barcode extraction."),
               p("The pipeline supports two main processing modes: 'concat' mode for concatenating filtered reads, and 'merge' mode for merging paired-end reads. Both pathways converge at the multi-parameter barcode recovery step (MitoGeneExtractor), followed by consensus cleaning (Fasta_cleaner), barcode consensus selection (Fasta_compare), and final barcode validation (barcode_validator)."),
               style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
             ),
             div(
               img(src = "images/workflow.png", 
                   class = "workflow-image",
                   alt = "BGEE Workflow Diagram showing the complete pipeline from raw reads through barcode validation"),
               style = "text-align: center;"
             )
    ),
    
    # Barcoding Outcome Tab
    tabPanel("Barcoding Outcome",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               # Static description text area
               div(
                 h3("🚦 Barcoding Outcome Overview"),
                 p("Overview of barcode extraction results for all Process IDs in the loaded plate or project, and whether they are BOLD BIN-compliant (pass) or not (fail).",
                   "🔗 ", strong("BOLD Systems Integration:"), " Process IDs in the table below are clickable links that will take you directly to the corresponding record in the ",
                   a("BOLD Systems database", href = "https://portal.boldsystems.org", target = "_blank"), 
                   ". Click on any Process ID to view detailed specimen and sequence information.",
                 ),
                 style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;"
               ),
               
               # Summary
               div(id = "outcome_summary_container",
                   style = "margin-bottom: 20px;"),
               
               # Traffick light colour explanation
               div(
                 h4("Result Categories", style = "color: #495057; margin-bottom: 15px;"),
                 div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin-bottom: 20px;",
                     div(style = "background-color: #d4edda; padding: 15px; border-radius: 8px; border-left: 4px solid #28a745;",
                         div(style = "font-weight: bold; color: #155724; margin-bottom: 5px;", "🟢 Green - Pass"),
                         p("Structural validation successful (≥500bp in length, 0 ambiguous bases, no stop codons, Reading frame 1/2/3) AND Taxonomic validation successful (BOLD family exactly matches a term in taxonomic ID).", 
                           style = "margin: 0; color: #155724; font-size: 0.9em;")
                     ),
                     div(style = "background-color: #fff3cd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffc107;",
                         div(style = "font-weight: bold; color: #856404; margin-bottom: 5px;", "🟡 Amber - Partial"),
                         p("Either structural OR taxonomic validation successful, but not both. Still considered a fail overall.", 
                           style = "margin: 0; color: #856404; font-size: 0.9em;")
                     ),
                     div(style = "background-color: #f8d7da; padding: 15px; border-radius: 8px; border-left: 4px solid #dc3545;",
                         div(style = "font-weight: bold; color: #721c24; margin-bottom: 5px;", "🔴 Red - Fail"),
                         p("Neither structural nor taxonomic validation successful", 
                           style = "margin: 0; color: #721c24; font-size: 0.9em;")
                     )
                 ),
                 style = "background-color: white; padding: 20px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;"
               ),
               
               # Download button
               downloadButton("downloadOutcomeTable", "Download Outcome Table", class = "download-btn btn-primary"),
               
               # Outcome table
               DTOutput("outcomeTable")
             )
    ),
    
    # Fastp Metrics Tab
    tabPanel("Fastp Metrics",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
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
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
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
    
    # Barcode Validation Results Tab
    tabPanel("Barcode Validation",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("About Barcode Validation"),
                 HTML("<a href='https://github.com/bge-barcoding/barcode_validator' target='_blank'>Barcode Validation</a> is the final validation step of the BGEE workflow, and consists of structural and taxonmic validation of each barcode consensus sequence."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               downloadButton("downloadBarcodeValidation", "Download Filtered Barcode Validation", class = "download-btn btn-primary"),
               DTOutput("barcodeValidationTable")
             )
    ),
    
    # Barcodes Tab
    tabPanel("Barcodes",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("Barcode Sequences (FASTA Format)"),
                 p("View and download barcode consensus sequences in FASTA format. Sequences are derived from the Barcode Validation results."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               # Filter options
               div(
                 h5("Filter Options"),
                 fluidRow(
                   column(3, checkboxInput("show_all_sequences", "All Sequences", value = TRUE)),
                   column(3, checkboxInput("show_pass_sequences", "Pass Only", value = FALSE)),
                   column(3, checkboxInput("show_amber_sequences", "Partial Only", value = FALSE)),
                   column(3, checkboxInput("show_fail_sequences", "Fail Only", value = FALSE))
                 ),
                 style = "background-color: white; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;"
               ),
               
               # Action buttons
               fluidRow(
                 column(6, downloadButton("downloadFasta", "Download FASTA File", class = "download-btn btn-primary")),
                 column(6, actionButton("copyFasta", "Copy All to Clipboard", class = "btn-secondary", style = "margin: 10px 0;"))
               ),
               
               # FASTA content display
               div(
                 h5("FASTA Content"),
                 verbatimTextOutput("fastaContent"),
                 style = "margin-top: 20px;"
               )
             )
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Add resource path for serving images from data directory
  addResourcePath("images", "./data")
  
  values <- reactiveValues(
    summary_stats = data.frame(),
    barcode_validation = data.frame(),
    json_data = data.frame(),
    dataset_loaded = FALSE,
    plate_index = data.frame(),
    current_identifier = "",
    current_identifier_type = "",
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
  
  # Clear the other dropdown when one is selected
  observe({
    if (!is.null(input$plate_search) && input$plate_search != "") {
      updateSelectInput(session, "project_search", selected = "")
    }
  })
  
  observe({
    if (!is.null(input$project_search) && input$project_search != "") {
      updateSelectInput(session, "plate_search", selected = "")
    }
  })
  
  # Create output for conditionalPanel
  output$dataset_loaded <- reactive({
    values$dataset_loaded
  })
  outputOptions(output, "dataset_loaded", suspendWhenHidden = FALSE)
  
  # Load data when button is clicked
  observeEvent(input$load_data, {
    # Check which selection method was used
    plate_selected <- !is.null(input$plate_search) && input$plate_search != ""
    project_selected <- !is.null(input$project_search) && input$project_search != ""
    
    # Validation logic
    if (!plate_selected && !project_selected) {
      showNotification("Please select either a plate or a project first.", type = "warning")
      return()
    }
    
    if (plate_selected && project_selected) {
      showNotification("Please select either a plate OR a project, not both.", type = "warning")
      return()
    }
    
    # Determine identifier and process IDs based on selection
    if (plate_selected) {
      # Plate selection logic
      plate_row <- values$plate_index[values$plate_index$plate_id == input$plate_search, ]
      if (nrow(plate_row) == 0) {
        showNotification("No data found for this plate.", type = "error")
        return()
      }
      
      identifier <- input$plate_search
      identifier_type <- "plate"
      target_process_ids <- unlist(strsplit(plate_row$process_ids, ","))
      
    } else {
      # Project selection logic
      identifier <- input$project_search
      identifier_type <- "project"
      target_process_ids <- get_process_ids_by_project(input$project_search)
      
      if (length(target_process_ids) == 0) {
        showNotification(paste("No Process IDs found containing project code:", input$project_search), type = "error")
        return()
      }
    }
    
    # Disable button and clear status
    updateActionButton(session, "load_data", label = "Loading...", icon = NULL)
    shinyjs::disable("load_data")
    removeUI("#plate_status div")
    
    tryCatch({
      # Load data first (without showing progress)
      result <- load_data(identifier, target_process_ids, identifier_type)
      
      # Check if any actual data was found
      total_records <- sum(unlist(result$debug_info$records_found))
      missing_data_types <- c()
      if (result$debug_info$records_found$bgee_records == 0) missing_data_types <- c(missing_data_types, "BGEE")
      if (result$debug_info$records_found$validation_records == 0) missing_data_types <- c(missing_data_types, "Validation")
      if (result$debug_info$records_found$json_records == 0) missing_data_types <- c(missing_data_types, "JSON")
      
      if (total_records == 0 || length(missing_data_types) > 0) {
        # ERROR CASE: No data found - skip loading animation entirely
        values$summary_stats <- data.frame()
        values$barcode_validation <- data.frame()
        values$json_data <- data.frame()
        values$current_identifier <- ""
        values$current_identifier_type <- ""
        values$current_process_ids <- character()
        values$dataset_loaded <- FALSE
        
        error_msg <- if (total_records == 0) {
          paste("ERROR: No data files found for", identifier_type, identifier, 
                "- Searched", length(target_process_ids), "Process IDs but found 0 records")
        } else {
          paste("ERROR: Missing data types for", identifier_type, identifier, 
                "- Missing:", paste(missing_data_types, collapse = ", "))
        }
        
        insertUI("#plate_status", "beforeEnd",
                 div(class = "status-message status-error",
                     p(error_msg, style = "margin: 0;")))
        
        showNotification("Failed to load data - no data found.", type = "error")
        
      } else {
        # SUCCESS CASE: Data found - show loading animation and process
        shinyjs::show("loading_section")
        
        # Simulate progress updates
        session$sendCustomMessage("updateProgress", list(percent = 25, message = "Processing BGEE data"))
        Sys.sleep(0.5)
        session$sendCustomMessage("updateProgress", list(percent = 50, message = "Processing Barcode Validation data"))
        Sys.sleep(0.5)
        session$sendCustomMessage("updateProgress", list(percent = 75, message = "Processing Fastp JSON data"))
        Sys.sleep(0.5)
        session$sendCustomMessage("updateProgress", list(percent = 100, message = "Complete"))
        
        # Hide loading and update values
        shinyjs::hide("loading_section")
        
        values$summary_stats <- result$summary_stats
        values$barcode_validation <- result$barcode_validation
        values$json_data <- result$json_data
        values$current_identifier <- identifier
        values$current_identifier_type <- identifier_type
        values$current_process_ids <- target_process_ids
        values$dataset_loaded <- TRUE
        
        insertUI("#plate_status", "beforeEnd",
                 div(class = "status-message status-success",
                     p(paste("Successfully loaded", identifier_type, ":", identifier, 
                             "- Process IDs:", length(target_process_ids)), style = "margin: 0;")))
        
        showNotification(paste(toupper(substring(identifier_type, 1, 1)), substring(identifier_type, 2), identifier, "loaded successfully!"), type = "message")
      }
      
    }, error = function(e) {
      shinyjs::hide("loading_section")
      
      insertUI("#plate_status", "beforeEnd",
               div(class = "status-message status-error",
                   p(paste("Error loading data:", e$message), style = "margin: 0;")))
      
      showNotification("Failed to load data. Please try again.", type = "error")
    })
    
    # Re-enable button
    updateActionButton(session, "load_data", label = "Load Results", icon = NULL)
    shinyjs::enable("load_data")
  })
  
  # Prepare outcome data based on barcode validation
  prepare_outcome_data <- reactive({
    req(values$dataset_loaded, nrow(values$barcode_validation) > 0)
    
    if (!"process_id" %in% names(values$barcode_validation)) {
      return(data.frame())
    }
    
    # Group by process_id and determine validation outcome
    outcome_data <- values$barcode_validation %>%
      group_by(process_id) %>%
      summarise(
        # Check structural validation criteria - UPDATED CRITERIA
        has_sufficient_length = any(nuc_basecount >= 500, na.rm = TRUE),
        has_no_ambiguous = any(ambig_basecount == 0, na.rm = TRUE),  # Changed from <= 6 to == 0
        has_no_stop_codons = any(stop_codons == 0, na.rm = TRUE),
        has_reading_frame = any(reading_frame %in% c(1, 2, 3), na.rm = TRUE),  # Changed from == 1 to 1,2,3
        
        # Get best values for display
        best_length = ifelse(any(!is.na(nuc_basecount)), max(nuc_basecount, na.rm = TRUE), NA),
        best_ambig = ifelse(any(!is.na(ambig_basecount)), min(ambig_basecount, na.rm = TRUE), NA),
        best_stop_codons = ifelse(any(!is.na(stop_codons)), min(stop_codons, na.rm = TRUE), NA),
        best_frame = ifelse(any(!is.na(reading_frame)), max(reading_frame, na.rm = TRUE), NA),
        best_identification = ifelse(any(!is.na(identification) & identification != ""), 
                                     first(identification[!is.na(identification) & identification != ""]), NA_character_),
        best_taxonomic_id = ifelse(any(!is.na(obs_taxon) & obs_taxon != ""), 
                                   first(obs_taxon[!is.na(obs_taxon) & obs_taxon != ""]), NA_character_),
        
        .groups = 'drop'
      ) %>%
      mutate(
        # Check taxonomic match using new function
        taxonomic_match = mapply(check_taxonomic_match, best_identification, best_taxonomic_id),
        taxonomic_pass = taxonomic_match == "YES",
        
        # Determine overall validation status with amber category
        structural_pass = has_sufficient_length & has_no_ambiguous & has_no_stop_codons & has_reading_frame,
        overall_pass = structural_pass & taxonomic_pass,
        partial_pass = (structural_pass & !taxonomic_pass) | (!structural_pass & taxonomic_pass),
        
        outcome_category = case_when(
          overall_pass ~ "pass",
          partial_pass ~ "amber", 
          TRUE ~ "fail"
        ),
        outcome_status = case_when(
          overall_pass ~ "✅ Pass",
          partial_pass ~ "🟡 Partial",
          TRUE ~ "❌ Fail"
        ),
        
        # Display values
        display_length = ifelse(!is.na(best_length), as.character(best_length), "N/A"),
        display_ambig = ifelse(!is.na(best_ambig), as.character(best_ambig), "N/A"),
        display_stop = ifelse(!is.na(best_stop_codons), as.character(best_stop_codons), "N/A"),
        display_frame = ifelse(!is.na(best_frame), as.character(best_frame), "N/A"),
        display_family = ifelse(!is.na(best_identification), best_identification, "N/A"),
        display_taxonomic_id = ifelse(!is.na(best_taxonomic_id), best_taxonomic_id, "N/A"),
        
        # Create BOLD link for Process ID
        process_id_link = paste0('<a href="https://portal.boldsystems.org/record/', process_id, '" target="_blank">', process_id, '</a>')
      ) %>%
      select(
        `Process ID` = process_id_link,
        `Status` = outcome_status,
        `Barcode Length` = display_length,
        `N Count` = display_ambig,
        `Stop Codons` = display_stop,
        `Reading Frame` = display_frame,
        `BOLD Family` = display_family,
        `Taxonomic ID` = display_taxonomic_id,
        `Taxonomic match` = taxonomic_match,
        outcome_category
      )
    
    return(outcome_data)
  })
  
  # Render outcome table
  output$outcomeTable <- renderDT({
    outcome_data <- prepare_outcome_data()
    
    if (nrow(outcome_data) == 0) {
      return(datatable(data.frame("No data available" = "Please check that Barcode Validation data is loaded correctly.")))
    }
    
    display_data <- outcome_data %>% select(-outcome_category)
    
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
          list(className = 'dt-center', targets = c(1, 8))  # Center Status and Taxonomic match columns
        ),
        rowCallback = JS(
          "function(row, data, index) {",
          "  var status = data[1];",
          "  if (status.includes('Pass')) {",
          "    $(row).css('background-color', '#d4edda');",
          "  } else if (status.includes('Partial')) {",
          "    $(row).css('background-color', '#fff3cd');",
          "  } else if (status.includes('Fail')) {",
          "    $(row).css('background-color', '#f8d7da');",
          "  }",
          "}"
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    )
  })
  
  # Render outcome summary
  output$outcome_summary <- renderUI({
    outcome_data <- prepare_outcome_data()
    
    if (nrow(outcome_data) == 0) {
      return(div())
    }
    
    total_processes <- nrow(outcome_data)
    pass_count <- sum(outcome_data$outcome_category == "pass", na.rm = TRUE)
    amber_count <- sum(outcome_data$outcome_category == "amber", na.rm = TRUE)
    fail_count <- sum(outcome_data$outcome_category == "fail", na.rm = TRUE)
    
    pass_rate <- round((pass_count / total_processes) * 100, 1)
    
    div(
      h4("Summary Statistics", style = "color: #495057; margin-bottom: 15px;"),
      div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 15px;",
          div(class = "outcome-stat", style = "background-color: white; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #495057;", total_processes),
              div(style = "color: #6c757d; margin-top: 5px;", "Total Process IDs")
          ),
          div(class = "outcome-stat", style = "background-color: #d4edda; padding: 15px; border-radius: 8px; border: 1px solid #c3e6cb; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #155724;", pass_count),
              div(style = "color: #155724; margin-top: 5px;", "Passed")
          ),
          div(class = "outcome-stat", style = "background-color: #fff3cd; padding: 15px; border-radius: 8px; border: 1px solid #ffeaa7; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #856404;", amber_count),
              div(style = "color: #856404; margin-top: 5px;", "Partial")
          ),
          div(class = "outcome-stat", style = "background-color: #f8d7da; padding: 15px; border-radius: 8px; border: 1px solid #f5c6cb; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #721c24;", fail_count),
              div(style = "color: #721c24; margin-top: 5px;", "Failed")
          ),
          div(class = "outcome-stat", style = "background-color: #e3f2fd; padding: 15px; border-radius: 8px; border: 1px solid #bbdefb; text-align: center;",
              div(style = "font-size: 1.5em; font-weight: bold; color: #1976d2;", paste0(pass_rate, "%")),
              div(style = "color: #1976d2; margin-top: 5px;", "Pass Rate")
          )
      ),
      style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px;"
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
  
  # Handle checkbox logic for FASTA filtering
  observe({
    if (input$show_all_sequences) {
      updateCheckboxInput(session, "show_pass_sequences", value = FALSE)
      updateCheckboxInput(session, "show_amber_sequences", value = FALSE)
      updateCheckboxInput(session, "show_fail_sequences", value = FALSE)
    }
  })
  
  observe({
    if (input$show_pass_sequences || input$show_amber_sequences || input$show_fail_sequences) {
      updateCheckboxInput(session, "show_all_sequences", value = FALSE)
    }
  })
  
  # Prepare FASTA data based on filters
  prepare_fasta_data <- reactive({
    req(values$dataset_loaded, nrow(values$barcode_validation) > 0)
    
    if (!"process_id" %in% names(values$barcode_validation) || !"nuc" %in% names(values$barcode_validation)) {
      return("")
    }
    
    # Get outcome data to determine pass/fail status
    outcome_data <- prepare_outcome_data()
    
    # Merge barcode validation with outcome status
    validation_with_status <- values$barcode_validation %>%
      left_join(outcome_data %>% 
                  mutate(process_id = gsub('<.*?>', '', `Process ID`)) %>%
                  select(process_id, outcome_category), 
                by = "process_id")
    
    # Filter based on selected checkboxes
    if (input$show_all_sequences) {
      filtered_data <- validation_with_status
    } else {
      show_statuses <- c()
      if (input$show_pass_sequences) show_statuses <- c(show_statuses, "pass")
      if (input$show_amber_sequences) show_statuses <- c(show_statuses, "amber")
      if (input$show_fail_sequences) show_statuses <- c(show_statuses, "fail")
      
      if (length(show_statuses) == 0) {
        return("No filter options selected.")
      }
      
      filtered_data <- validation_with_status %>%
        filter(outcome_category %in% show_statuses)
    }
    
    if (nrow(filtered_data) == 0) {
      return("No sequences match the selected filters.")
    }
    
    # Create FASTA format
    fasta_lines <- c()
    for (i in 1:nrow(filtered_data)) {
      process_id <- filtered_data$process_id[i]
      sequence <- filtered_data$nuc[i]
      
      # Add header
      fasta_lines <- c(fasta_lines, paste0(">", process_id))
      
      # Add sequence (or empty line if NA/missing)
      if (is.na(sequence) || sequence == "" || sequence == "NULL") {
        fasta_lines <- c(fasta_lines, "")
      } else {
        fasta_lines <- c(fasta_lines, as.character(sequence))
      }
    }
    
    return(paste(fasta_lines, collapse = "\n"))
  })
  
  # Render FASTA content
  output$fastaContent <- renderText({
    prepare_fasta_data()
  })
  
  # Copy to clipboard functionality
  observeEvent(input$copyFasta, {
    fasta_content <- prepare_fasta_data()
    session$sendCustomMessage("copyToClipboard", fasta_content)
  })
  
  # Show clipboard success/error messages
  observeEvent(input$clipboard_success, {
    showNotification("FASTA content copied to clipboard!", type = "message")
  })
  
  observeEvent(input$clipboard_error, {
    showNotification("Failed to copy to clipboard. Please try manual selection.", type = "warning")
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
  
  output$barcodeValidationTable <- renderDT({
    req(values$dataset_loaded, nrow(values$barcode_validation) > 0)
    
    # Select and rename the specific columns requested
    display_data <- values$barcode_validation %>%
      select(
        `Process ID` = process_id,
        `N count` = ambig_basecount,
        `Barcode length` = nuc_basecount,
        `Frame` = reading_frame,
        `Stop codons` = stop_codons,
        `BOLD species` = species,
        `BOLD family` = identification,
        `Taxonomic ID` = obs_taxon
      )
    
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
      paste0("barcoding_outcome_", values$current_identifier, "_", Sys.Date(), ".csv")
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
      paste0("summary_stats_", values$current_identifier, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$summary_stats, file, row.names = FALSE)
    }
  )
  
  output$downloadBarcodeValidation <- downloadHandler(
    filename = function() {
      paste0("barcode_validation_", values$current_identifier, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$barcode_validation, file, row.names = FALSE)
    }
  )
  
  output$downloadJsonData <- downloadHandler(
    filename = function() {
      paste0("json_data_", values$current_identifier, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$json_data, file, row.names = FALSE)
    }
  )
  
  output$downloadFasta <- downloadHandler(
    filename = function() {
      filter_suffix <- if (input$show_all_sequences) {
        "all"
      } else {
        filters <- c()
        if (input$show_pass_sequences) filters <- c(filters, "pass")
        if (input$show_amber_sequences) filters <- c(filters, "amber")
        if (input$show_fail_sequences) filters <- c(filters, "fail")
        paste(filters, collapse = "_")
      }
      paste0("barcodes_", values$current_identifier, "_", filter_suffix, "_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      fasta_content <- prepare_fasta_data()
      writeLines(fasta_content, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)