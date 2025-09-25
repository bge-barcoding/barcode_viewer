library(shiny)
library(DT)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)

# Base data directories for flat structure
BASE_DIRS <- list(
  BeeGees_summary = "./data/BeeGees_summary_stats/",
  barcode_validation = "./data/barcode_validation/",
  barcode_validation_merged = "./data/barcode_validation/merged/",
  barcodes = "./data/barcodes/"
)

BASE_DATA_DIR <- "./data"

# Utility function to standardise column names across datasets
standardise_column_names <- function(data, target_name, possible_names) {
  if (nrow(data) == 0 || ncol(data) == 0) {
    return(data)
  }
  
  # Find which of the possible names exists in the data
  matching_col <- intersect(possible_names, names(data))
  
  if (length(matching_col) > 0) {
    # Rename the first matching column to the target name
    names(data)[names(data) == matching_col[1]] <- target_name
  }
  
  return(data)
}

# File reading with error handling
safe_read_file <- function(file, reader_function, processor_function = NULL, file_type = "file") {
  tryCatch({
    # Read the file using the provided reader function
    data <- reader_function(file)
    
    # Apply processing function if provided
    if (!is.null(processor_function)) {
      data <- processor_function(data)
    }
    
    return(data)
  }, error = function(e) {
    message(paste("Error reading", file_type, basename(file), ":", e$message))
    return(NULL)
  })
}

# Specific reader functions
read_csv_safe <- function(file) {
  read.csv(file, stringsAsFactors = FALSE)
}

read_tsv_safe <- function(file) {
  read.delim(file, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
}

# Utility function to filter data by process IDs with consistent cleaning
filter_by_process_ids <- function(data, target_process_ids, id_column = "process_id") {
  if (nrow(data) == 0 || length(target_process_ids) == 0) {
    return(data)
  }
  
  if (!id_column %in% names(data)) {
    warning(paste("Column", id_column, "not found in data"))
    return(data.frame())
  }
  
  # Clean both the data column and target IDs consistently
  target_process_ids_clean <- trimws(as.character(target_process_ids))
  data_ids_clean <- trimws(as.character(data[[id_column]]))
  
  # Filter and return
  filtered_data <- data[data_ids_clean %in% target_process_ids_clean, ]
  return(filtered_data)
}

# Function to wrap sequence text every 100 characters for display
wrap_sequence <- function(sequence, wrap_length = 100) {
  if (is.na(sequence) || sequence == "" || sequence == "NULL") {
    return("")
  }
  
  seq_char <- as.character(sequence)
  if (nchar(seq_char) <= wrap_length) {
    return(seq_char)
  }
  
  # Split sequence into chunks of wrap_length
  chunks <- substring(seq_char, 
                      seq(1, nchar(seq_char), wrap_length),
                      seq(wrap_length, nchar(seq_char), wrap_length))
  
  # Join with regular line breaks for display
  paste(chunks, collapse = "\n")
}

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
    # Define processor function for metadata
    process_metadata <- function(data) {
      data <- standardise_column_names(data, "sample_id", c("sample_id", "Sample ID", "Sample.ID", "sampleid"))
      data <- standardise_column_names(data, "process_id", c("process_id", "Process ID", "Process.ID", "processid"))
      return(data)
    }
    
    metadata <- safe_read_file(csv_file, read_csv_safe, process_metadata, "CSV")
    
    if (!is.null(metadata) && "sample_id" %in% names(metadata) && "process_id" %in% names(metadata)) {
      message("Found metadata file: ", basename(csv_file))
      message("Loaded ", nrow(metadata), " sample records")
      return(metadata)
    }
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
  
  # Get process IDs from BeeGees summary stats files
  if (dir.exists(BASE_DIRS$BeeGees_summary)) {
    csv_files <- list.files(BASE_DIRS$BeeGees_summary, pattern = "\\.csv$", full.names = TRUE)
    
    process_BeeGees_data <- function(data) {
      standardise_column_names(data, "process_id", c("Process.ID", "Process ID", "process_id", "ID"))
    }
    
    for (file in csv_files) {
      data <- safe_read_file(file, read_csv_safe, process_BeeGees_data, "CSV")
      if (!is.null(data) && "process_id" %in% names(data)) {
        all_process_ids <- c(all_process_ids, data$process_id)
      }
    }
  }
  
  # Also get process IDs from barcode validation files (both structval and taxval)
  if (dir.exists(BASE_DIRS$barcode_validation)) {
    tsv_files <- list.files(BASE_DIRS$barcode_validation, pattern = "\\.tsv$", full.names = TRUE)
    
    process_validation_data <- function(data) {
      standardise_column_names(data, "process_id", c("ID", "group_id", "process_id", "Process.ID", "Process ID"))
    }
    
    for (file in tsv_files) {
      data <- safe_read_file(file, read_tsv_safe, process_validation_data, "TSV")
      if (!is.null(data) && "process_id" %in% names(data)) {
        all_process_ids <- c(all_process_ids, data$process_id)
      }
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

# Load FASTA sequences and extract Process IDs
load_fasta_data <- function() {
  if (!dir.exists(BASE_DIRS$barcodes)) {
    message("Barcodes directory does not exist")
    return(list(sequences = data.frame(), process_ids = character()))
  }
  
  fasta_files <- list.files(BASE_DIRS$barcodes, pattern = "\\.(fasta|fas|fa)$", full.names = TRUE, ignore.case = TRUE)
  
  if (length(fasta_files) == 0) {
    message("No FASTA files found in barcodes directory")
    return(list(sequences = data.frame(), process_ids = character()))
  }
  
  message(paste("Found", length(fasta_files), "FASTA files in barcodes directory"))
  
  all_sequences <- data.frame()
  all_process_ids <- character()
  
  for (file in fasta_files) {
    tryCatch({
      # Read FASTA file
      lines <- readLines(file, warn = FALSE)
      
      # Parse FASTA format
      header_indices <- which(grepl("^>", lines))
      
      if (length(header_indices) > 0) {
        for (i in seq_along(header_indices)) {
          header_line <- lines[header_indices[i]]
          
          # Extract Process ID from header (everything before first underscore after >)
          header_clean <- gsub("^>", "", header_line)
          process_id <- strsplit(header_clean, "_")[[1]][1]
          
          if (!is.na(process_id) && process_id != "") {
            # Get sequence (lines between this header and next header or end of file)
            start_line <- header_indices[i] + 1
            end_line <- if (i < length(header_indices)) {
              header_indices[i + 1] - 1
            } else {
              length(lines)
            }
            
            if (start_line <= end_line) {
              sequence <- paste(lines[start_line:end_line], collapse = "")
              
              # Add to results
              all_sequences <- rbind(all_sequences, data.frame(
                process_id = process_id,
                header = header_clean,
                sequence = sequence,
                stringsAsFactors = FALSE
              ))
              
              all_process_ids <- c(all_process_ids, process_id)
            }
          }
        }
      }
    }, error = function(e) {
      message(paste("Error reading FASTA file", basename(file), ":", e$message))
    })
  }
  
  unique_process_ids <- unique(all_process_ids)
  message(paste("Loaded", nrow(all_sequences), "sequences for", length(unique_process_ids), "unique Process IDs from FASTA files"))
  
  # Sort sequences alphabetically by process_id
  if (nrow(all_sequences) > 0) {
    all_sequences <- all_sequences[order(all_sequences$process_id), ]
  }
  
  return(list(sequences = all_sequences, process_ids = unique_process_ids))
}

# Check if a value is considered null/empty
is_null_value <- function(value) {
  null_values <- c("None", "none", "Null", "null", "NA", "na", "", "NULL")
  return(is.na(value) || trimws(as.character(value)) %in% null_values)
}

# Merge structval and taxval data (similar to Python script logic)
merge_validation_data <- function(structval_data, taxval_data, key_column = "ID") {
  if (nrow(structval_data) == 0 && nrow(taxval_data) == 0) {
    return(data.frame())
  }
  
  # If one dataset is empty, return the other (sorted)
  if (nrow(structval_data) == 0) {
    if ("ID" %in% names(taxval_data)) {
      names(taxval_data)[names(taxval_data) == "ID"] <- "process_id"
    }
    return(taxval_data[order(taxval_data$process_id), ])
  }
  if (nrow(taxval_data) == 0) {
    if ("ID" %in% names(structval_data)) {
      names(structval_data)[names(structval_data) == "ID"] <- "process_id"
    }
    return(structval_data[order(structval_data$process_id), ])
  }
  
  message("Merging structval and taxval data...")
  message(paste("Structval shape:", nrow(structval_data), "x", ncol(structval_data)))
  message(paste("Taxval shape:", nrow(taxval_data), "x", ncol(taxval_data)))
  
  # Get all unique columns
  all_columns <- unique(c(names(structval_data), names(taxval_data)))
  
  # Find common and unique columns
  common_columns <- intersect(names(structval_data), names(taxval_data))
  structval_only <- setdiff(names(structval_data), names(taxval_data))
  taxval_only <- setdiff(names(taxval_data), names(structval_data))
  
  message(paste("Common columns:", length(common_columns)))
  message(paste("Structval only:", paste(structval_only, collapse = ", ")))
  message(paste("Taxval only:", paste(taxval_only, collapse = ", ")))
  
  # Standardise key column name in both datasets
  structval_data <- standardise_column_names(structval_data, "process_id", c(key_column))
  taxval_data <- standardise_column_names(taxval_data, "process_id", c(key_column))
  
  # Perform outer merge
  merged_data <- full_join(taxval_data, structval_data, by = "process_id", suffix = c("_taxval", "_structval"))
  
  # Initialise result dataframe
  result_data <- data.frame(process_id = merged_data$process_id, stringsAsFactors = FALSE)
  
  # Process each column with merge logic
  dataset_combinations <- 0
  null_replacements <- 0
  
  for (col in all_columns) {
    if (col == "process_id" || col == key_column) {
      next
    }
    
    taxval_col <- paste0(col, "_taxval")
    structval_col <- paste0(col, "_structval")
    
    # Get column values (handle cases where column doesn't exist in one dataset)
    taxval_values <- if (taxval_col %in% names(merged_data)) {
      merged_data[[taxval_col]]
    } else if (col %in% names(merged_data)) {
      merged_data[[col]]
    } else {
      rep(NA, nrow(merged_data))
    }
    
    structval_values <- if (structval_col %in% names(merged_data)) {
      merged_data[[structval_col]]
    } else if (col %in% names(merged_data) && !paste0(col, "_taxval") %in% names(merged_data)) {
      merged_data[[col]]
    } else {
      rep(NA, nrow(merged_data))
    }
    
    # Apply merge rules
    result_values <- character(nrow(merged_data))
    
    for (i in 1:nrow(merged_data)) {
      taxval_val <- taxval_values[i]
      structval_val <- structval_values[i]
      
      taxval_is_null <- is_null_value(taxval_val)
      structval_is_null <- is_null_value(structval_val)
      
      # Special handling for 'dataset' column - always combine
      if (col == "dataset") {
        paths <- c()
        if (!structval_is_null) paths <- c(paths, as.character(structval_val))
        if (!taxval_is_null) paths <- c(paths, as.character(taxval_val))
        
        if (length(paths) == 2) {
          result_values[i] <- paste(paths, collapse = "; ")
          dataset_combinations <- dataset_combinations + 1
        } else if (length(paths) == 1) {
          result_values[i] <- paths[1]
        } else {
          result_values[i] <- ""
        }
      } else {
        # Standard null replacement logic
        if (taxval_is_null && structval_is_null) {
          result_values[i] <- ""
        } else if (taxval_is_null && !structval_is_null) {
          result_values[i] <- as.character(structval_val)
          null_replacements <- null_replacements + 1
        } else if (!taxval_is_null && structval_is_null) {
          result_values[i] <- as.character(taxval_val)
          null_replacements <- null_replacements + 1
        } else if (isTRUE(all.equal(taxval_val, structval_val))) {
          result_values[i] <- as.character(taxval_val)
        } else {
          # Conflict - default to structval
          result_values[i] <- as.character(structval_val)
        }
      }
    }
    
    result_data[[col]] <- result_values
  }
  
  message(paste("Dataset combinations:", dataset_combinations))
  message(paste("Null replacements:", null_replacements))
  message(paste("Merged result shape:", nrow(result_data), "x", ncol(result_data)))
  
  # Sort the result by process_id alphabetically
  result_data <- result_data[order(result_data$process_id), ]
  
  return(result_data)
}

# Function to clean and convert potentially numeric columns
clean_numeric_columns <- function(data) {
  if (nrow(data) == 0) return(data)
  
  null_values <- c("None", "none", "Null", "null", "NA", "na", "", "NULL")
  
  for (col_name in names(data)) {
    col_data <- data[[col_name]]
    
    # Skip if already numeric
    if (is.numeric(col_data)) next
    
    # Convert null representations to NA
    col_data_clean <- ifelse(trimws(as.character(col_data)) %in% null_values | is.na(col_data), 
                             NA, as.character(col_data))
    
    # Check if remaining non-NA values can be converted to numeric
    non_na_values <- col_data_clean[!is.na(col_data_clean)]
    if (length(non_na_values) > 0) {
      # Try to convert to numeric
      numeric_test <- suppressWarnings(as.numeric(non_na_values))
      
      # If all non-NA values successfully convert to numeric, convert the whole column
      if (all(!is.na(numeric_test))) {
        data[[col_name]] <- as.numeric(col_data_clean)
        message(paste("Converted column", col_name, "from character to numeric"))
      }
    }
  }
  
  return(data)
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

# Enhanced data loading with corrected validation logic
load_data <- function(identifier, target_process_ids, identifier_type = "plate") {
  result <- list(
    summary_stats = data.frame(),
    barcode_validation = data.frame(),
    fasta_data = list(),
    taxval_process_ids = character(),
    status = "",
    debug_info = list()
  )
  
  # Track actual records found for target process IDs
  records_found <- list(
    BeeGees_records = 0,
    validation_records = 0,
    fasta_sequences = 0
  )
  
  process_ids_found <- list(
    BeeGees_ids = character(),
    validation_ids = character(),
    fasta_ids = character(),
    taxval_ids = character()
  )
  
  message(paste("Loading data for", identifier_type, ":", identifier))
  message(paste("Target Process IDs count:", length(target_process_ids)))
  message(paste("Sample target IDs:", paste(head(target_process_ids, 3), collapse = ", ")))
  
  # Load BeeGees summary stats with numeric conversion
  if (dir.exists(BASE_DIRS$BeeGees_summary)) {
    csv_files <- list.files(BASE_DIRS$BeeGees_summary, pattern = "\\.csv$", full.names = TRUE)
    all_summary_data <- data.frame()
    
    message(paste("Found", length(csv_files), "CSV files in BeeGees_summary directory"))
    
    for (file in csv_files) {
      tryCatch({
        data <- read.csv(file, stringsAsFactors = FALSE)
        
        # The BeeGees files use "ID" column based on sample
        data <- standardise_column_names(data, "process_id", c("ID", "Process.ID", "Process ID", "process_id"))
        
        if ("process_id" %in% names(data)) {
          # Apply numeric conversion BEFORE combining data
          data <- clean_numeric_columns(data)
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
      result$summary_stats <- filter_by_process_ids(all_summary_data, target_process_ids)
      records_found$BeeGees_records <- nrow(result$summary_stats)
      
      if (records_found$BeeGees_records > 0) {
        process_ids_found$BeeGees_ids <- unique(result$summary_stats$process_id)
        message(paste("BeeGees Summary: Found", records_found$BeeGees_records, "records for", length(process_ids_found$BeeGees_ids), "process IDs"))
        
        # Debug: Show which columns are now numeric
        numeric_cols <- names(result$summary_stats)[sapply(result$summary_stats, is.numeric)]
        message(paste("Numeric columns after conversion:", paste(numeric_cols, collapse = ", ")))
      } else {
        message("BeeGees Summary: No records found for target process IDs")
      }
    }
  } else {
    message("BeeGees summary directory does not exist")
  }
  
  # Load validation data (structval and taxval for logic, merged for display)
  if (dir.exists(BASE_DIRS$barcode_validation)) {
    structval_files <- list.files(BASE_DIRS$barcode_validation, pattern = "structval.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)
    taxval_files <- list.files(BASE_DIRS$barcode_validation, pattern = "taxval.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)
    
    all_taxval_data <- data.frame()
    
    message(paste("Found", length(structval_files), "structval TSV files"))
    message(paste("Found", length(taxval_files), "taxval TSV files"))
    
    # Load taxval files ONLY for validation logic (we only need these to determine structural pass)
    for (file in taxval_files) {
      tryCatch({
        data <- read.delim(file, stringsAsFactors = FALSE, sep = "\t", header = TRUE, quote = "")
        if ("ID" %in% names(data)) {
          all_taxval_data <- rbind(all_taxval_data, data)
        }
      }, error = function(e) {
        message(paste("Error reading taxval file", basename(file), ":", e$message))
      })
    }
    
    # Extract Process IDs that exist in taxval (these passed structural validation)
    if (nrow(all_taxval_data) > 0 && "ID" %in% names(all_taxval_data)) {
      result$taxval_process_ids <- unique(trimws(as.character(all_taxval_data$ID)))
      process_ids_found$taxval_ids <- result$taxval_process_ids
      message(paste("STRUCTURAL VALIDATION: Found", length(result$taxval_process_ids), "Process IDs in taxval files (passed structural validation)"))
    } else {
      message("WARNING: No taxval data found - cannot determine structural validation status")
      result$taxval_process_ids <- character()
    }
  } else {
    message("Barcode Validation directory does not exist")
  }
  
  # Load pre-merged validation data for display
  if (dir.exists(BASE_DIRS$barcode_validation_merged)) {
    merged_files <- list.files(BASE_DIRS$barcode_validation_merged, pattern = "\\.tsv$", full.names = TRUE)
    all_merged_data <- data.frame()
    
    message(paste("Found", length(merged_files), "pre-merged TSV files"))
    
    process_merged_data <- function(data) {
      standardise_column_names(data, "process_id", c("process_id", "ID"))
    }
    
    for (file in merged_files) {
      data <- safe_read_file(file, read_tsv_safe, process_merged_data, "merged TSV")
      if (!is.null(data) && "process_id" %in% names(data)) {
        all_merged_data <- rbind(all_merged_data, data)
      }
    }
    
    # Filter and sort merged data for target process IDs
    if (nrow(all_merged_data) > 0 && "process_id" %in% names(all_merged_data)) {
      result$barcode_validation <- filter_by_process_ids(all_merged_data, target_process_ids)
      
      # Sort by process_id alphabetically
      if (nrow(result$barcode_validation) > 0) {
        result$barcode_validation <- result$barcode_validation[order(result$barcode_validation$process_id), ]
      }
      
      records_found$validation_records <- nrow(result$barcode_validation)
      
      if (records_found$validation_records > 0) {
        process_ids_found$validation_ids <- unique(result$barcode_validation$process_id)
        message(paste("Barcode Validation (Pre-merged): Found", records_found$validation_records, "records for", length(process_ids_found$validation_ids), "process IDs"))
      }
    }
  } else {
    message("Merged Barcode Validation directory does not exist")
  }
  
  # Load FASTA data
  fasta_result <- load_fasta_data()
  result$fasta_data <- fasta_result
  
  if (length(fasta_result$process_ids) > 0) {
    # Filter FASTA data for target process IDs
    target_process_ids_clean <- trimws(as.character(target_process_ids))
    fasta_ids_clean <- trimws(as.character(fasta_result$process_ids))
    matching_fasta_ids <- fasta_result$process_ids[fasta_ids_clean %in% target_process_ids_clean]
    records_found$fasta_sequences <- length(matching_fasta_ids)
    process_ids_found$fasta_ids <- unique(matching_fasta_ids)
    
    if (records_found$fasta_sequences > 0) {
      message(paste("TAXONOMIC VALIDATION: Found", records_found$fasta_sequences, "sequences for target process IDs (passed taxonomic validation)"))
    }
  }
  
  # Generate status based on ACTUAL RECORDS FOUND
  total_records_found <- sum(unlist(records_found))
  
  if (total_records_found == 0) {
    result$status <- paste("ERROR: No data found for", identifier_type, identifier, 
                           "- Searched", length(target_process_ids), "Process IDs but found 0 records")
  } else {
    loaded_items <- c()
    
    if (records_found$BeeGees_records > 0) {
      loaded_items <- c(loaded_items, paste("BeeGees:", records_found$BeeGees_records, "records"))
    }
    
    if (records_found$validation_records > 0) {
      loaded_items <- c(loaded_items, paste("Validation:", records_found$validation_records, "records"))
    }
    
    if (records_found$fasta_sequences > 0) {
      loaded_items <- c(loaded_items, paste("FASTA:", records_found$fasta_sequences, "sequences"))
    }
    
    result$status <- paste("Successfully loaded:", paste(loaded_items, collapse = ", "))
    
    # Add process ID summary
    unique_process_ids_found <- unique(c(process_ids_found$BeeGees_ids, process_ids_found$validation_ids, process_ids_found$fasta_ids))
    result$status <- paste(result$status, "- Found data for", length(unique_process_ids_found), "of", length(target_process_ids), "Process IDs")
    
    # Add warning if some data types are completely missing
    missing_types <- c()
    if (records_found$BeeGees_records == 0) missing_types <- c(missing_types, "BeeGees Summary")
    if (records_found$validation_records == 0) missing_types <- c(missing_types, "Barcode Validation")
    if (records_found$fasta_sequences == 0) missing_types <- c(missing_types, "FASTA Sequences")
    
    if (length(missing_types) > 0) {
      result$status <- paste(result$status, "- WARNING: Missing data types:", paste(missing_types, collapse = ", "))
    }
    
    # Add validation status summary
    target_clean <- trimws(as.character(target_process_ids))
    taxval_clean <- trimws(as.character(result$taxval_process_ids))
    structural_pass_count <- length(target_clean[target_clean %in% taxval_clean])
    taxonomic_pass_count <- length(process_ids_found$fasta_ids[process_ids_found$fasta_ids %in% target_clean])
    
    message(paste("VALIDATION SUMMARY:"))
    message(paste("- Target Process IDs:", length(target_clean)))
    message(paste("- Passed Structural Validation:", structural_pass_count, "(present in taxval)"))
    message(paste("- Passed Taxonomic Validation:", taxonomic_pass_count, "(present in FASTA)"))
    message(paste("- Failed Both Validations (expected):", length(target_clean) - structural_pass_count))
  }
  
  # Store debug info
  result$debug_info <- list(
    records_found = records_found,
    process_ids_found = process_ids_found,
    target_process_ids = target_process_ids
  )
  
  return(result)
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
                 p("Please allow a minute for the data to load behind the scenes after clicking the 'Load Results' button. A loading bar will appear momentarily."),
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
                             "Initialising...")
                     )
                 ),
                 div(id = "plate_status", class = "status-message")
             ),
             br(),
             h4("About the Data"),
             p("Once you select and load a plate or project, you can navigate to the other tabs to view:"),
             tags$ul(
               tags$li(strong("BeeGees Summary Statistics:"), " Combined statistics from the Barcode gene Extraction and Evaluation from Genome Skims (BeeGees) pipeline."),
               tags$li(strong("Barcode Validation:"), " Merged structural and taxonomic validation results from the BeeGees pipeline."),
               tags$li(strong("Visualise Data:"), " Interactive plots for exploring the loaded data."),
               tags$li(strong("Barcodes:"), " FASTA sequences for validated barcodes.")
             ),
             
             br(),
             h4("Barcode gene Extraction and Evaluation from Genome Skims (BeeGees) pipeline"),
             div(
               p("The BeeGees (Barcode gene Extraction and Evaluation from Genome Skims) pipeline is a comprehensive workflow for recovery of high-quality barcode sequences from raw genome skim sequencing data derived from museum specimens. The workflow includes several quality control steps, barcode consensus sequence generation, cleaning, and validation processes to ensure reliable extraction of high-quality barcodes."),
               p("The pipeline supports two main processing modes: 'concat' mode for concatenating paired-end reads, and 'merge' mode for merging paired-end reads. Both pathways converge at the multi-parameter barcode recovery step (MitoGeneExtractor), followed by consensus cleaning (Fasta_cleaner), barcode consensus structural validation, taxonomic validation, and finally, seleciton of validated barcode sequences (barcode_validator)."),
               style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
             ),
             
             br(),
             h4("About the Barcoding Results Dashboard (barcode_viewer)"),
             div(
               p(strong("Author:"), " Dr. Daniel A.J. Parsons (NHMUK) - D.parsons@nhm.ac.uk"),
               p(strong("Contributors:"), " Dr. Benajmin W. Price (NHMUK) - B.price@nhm.ac.uk"),
               p(strong("Funding:"), " Development was supported by Biodiversity Genomics Europe (Grant no.101059492) which is funded by Horizon Europe under the Biodiversity, Circular Economy and Environment call (REA.B.3); co-funded by the Swiss State Secretariat for Education, Research and Innovation (SERI) under contract numbers 22.00173 and 24.00054; and by UK Research and Innovation (UKRI) under the Department for Business, Energy and Industrial Strategy’s Horizon Europe Guarantee Scheme."),
               p(strong("Computational resources:"), " The authors acknowledge Research Computing at the James Hutton Institute for providing computational resources and technical support for the “UK’s Crop Diversity Bioinformatics HPC” (BBSRC grants BB/S019669/1 and BB/X019683/1), use of which has contributed to the results reported within this dashboard"),
               p(strong("Version:"), " v1.0.2 (September 2025)"),
               style = "background-color: #f1f3f5; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
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
                 h3("Barcoding Outcome Overview"),
                 p("Overview of barcode recovery outcomes for all Process IDs in the loaded plate or project, showing structural and taxonomic validation status."),
                 p("Process IDs in the table below will take you directly to the corresponding record in the ",
                   a("BOLD Systems database", href = "https://portal.boldsystems.org", target = "_blank"))
                 ,
                 style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;"
               ),
               
               # Summary
               div(id = "outcome_summary_container",
                   style = "margin-bottom: 20px;"),
               
               # Result categories explanation
               div(
                 h4("Result Categories", style = "color: #495057; margin-bottom: 15px;"),
                 div(style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin-bottom: 20px;",
                     div(style = "background-color: #d4edda; padding: 15px; border-radius: 8px; border-left: 4px solid #28a745;",
                         div(style = "font-weight: bold; color: #155724; margin-bottom: 5px;", "PASS"),
                         p("Both structural validation AND taxonomic validation were successful.", 
                           style = "margin: 0; color: #155724; font-size: 0.9em;")
                     ),
                     div(style = "background-color: #fff3cd; padding: 15px; border-radius: 8px; border-left: 4px solid #ffc107;",
                         div(style = "font-weight: bold; color: #856404; margin-bottom: 5px;", "PARTIAL"),
                         p("Structural validation successful but taxonomic validation failed.", 
                           style = "margin: 0; color: #856404; font-size: 0.9em;")
                     ),
                     div(style = "background-color: #f8d7da; padding: 15px; border-radius: 8px; border-left: 4px solid #dc3545;",
                         div(style = "font-weight: bold; color: #721c24; margin-bottom: 5px;", "FAIL"),
                         p("Both structural AND taxonomic validation failed.", 
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
    
    # BeeGees Summary Statistics Tab (without scatter plot)
    tabPanel("BeeGees Summary Statistics",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("About BeeGees Summary Statistics"),
                 HTML("<a href='https://github.com/bge-barcoding/BeeGees' target='_blank'>BeeGees (Barcode gene Extraction and Evaluation from Genome Skims)</a> summary statistics are generated from the combined output of the Gene Fetch, MGE and fasta_cleaner steps of the pipeline."),
                 p("These statistics help evaluate barcode quality across samples."),
                 p("All summary statistics can be downloaded with the 'Download ALL Summary Stats' button below. 'Copy', 'CSV', 'Excel' buttons only copy/download the first 'page' of the table"),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               downloadButton("downloadSummaryStats", "Download Summary Stats", class = "download-btn btn-primary"),
               DTOutput("summaryStatsTable")
             )
    ),
    
    # Barcode Validation Tab (without scatter plot)
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
                 HTML("<a href='https://github.com/bge-barcoding/barcode_validator' target='_blank'>Barcode Validation</a> is the final step of the BeeGees workflow, consisting of merged structural and taxonomic validation results for each barcode consensus sequence."),
                 p("These statistics help evaluate barcode quality and recovery across samples."),
                 p("All statistics can be downloaded with the 'Download ALL Barcode Validation Data' button below. 'Copy', 'CSV', 'Excel' buttons only copy/download the first 'page' of the table"),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               downloadButton("downloadBarcodeValidation", "Download Barcode Validation Data", class = "download-btn btn-primary"),
               DTOutput("barcodeValidationTable")
             )
    ),
    
    # Visualise Data Tab - UPDATED SECTION
    tabPanel("Visualise Data",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a plate or project from the 'Search Results' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               div(
                 h4("Interactive Data Visualisation"),
                 p("Explore the loaded data with interactive plots. Choose different variables to create custom visualisations."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               # BeeGees Summary Statistics Scatter Plot - UPDATED ORDER
               h4("BeeGees Summary Statistics - Scatter Plot"),
               conditionalPanel(
                 condition = "output.has_summary_data == true",
                 fluidRow(
                   column(6, selectInput("xcol", "X Axis:", choices = NULL)),  # X Axis now comes first
                   column(6, selectInput("ycol", "Y Axis:", choices = NULL))   # Y Axis now comes second
                 ),
                 plotlyOutput("summaryStatsPlot", height = "500px")
               ),
               conditionalPanel(
                 condition = "output.has_summary_data == false",
                 div(p("No BeeGees Summary Statistics data available for visualisation."),
                     style = "color: #6c757d; font-style: italic; margin: 20px 0;")
               ),
               
               hr(),
               
               # Barcode Validation Bar Chart
               h4("Barcode Validation - Bar Chart"),
               conditionalPanel(
                 condition = "output.has_validation_data == true",
                 fluidRow(
                   column(6, selectInput("validation_xcol", "X Axis:", choices = NULL)),
                   column(6, selectInput("validation_ycol", "Y Axis (Numeric):", choices = NULL))
                 ),
                 div(
                   p("Note: The chart is scrollable if there are many entries. Use your mouse wheel or drag to navigate."),
                   style = "color: #6c757d; font-size: 0.9em; margin-bottom: 15px;"
                 ),
                 plotlyOutput("barcodeValidationPlot", height = "500px")
               ),
               conditionalPanel(
                 condition = "output.has_validation_data == false",
                 div(p("No Barcode Validation data available for visualisation."),
                     style = "color: #6c757d; font-style: italic; margin: 20px 0;")
               )
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
                 p("View and download barcode consensus sequences that PASSED validation, in FASTA format."),
                 p("If you are interested in FAILED barcode consensus sequences, these can be retrieved from the 'sequence' column of the Barcode Validation tab."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
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
    fasta_data = list(),
    taxval_process_ids = character(), # Store which process IDs passed structural validation
    dataset_loaded = FALSE,
    plate_index = data.frame(),
    current_identifier = "",
    current_identifier_type = "",
    current_process_ids = character()
  )
  
  # Initialise plate index
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
  
  # Create outputs for conditionalPanels
  output$dataset_loaded <- reactive({
    values$dataset_loaded
  })
  outputOptions(output, "dataset_loaded", suspendWhenHidden = FALSE)
  
  output$has_summary_data <- reactive({
    values$dataset_loaded && nrow(values$summary_stats) > 0
  })
  outputOptions(output, "has_summary_data", suspendWhenHidden = FALSE)
  
  output$has_validation_data <- reactive({
    values$dataset_loaded && nrow(values$barcode_validation) > 0
  })
  outputOptions(output, "has_validation_data", suspendWhenHidden = FALSE)
  
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
      # Load data
      result <- load_data(identifier, target_process_ids, identifier_type)
      
      # Check if any actual data was found
      total_records <- sum(unlist(result$debug_info$records_found))
      missing_data_types <- c()
      if (result$debug_info$records_found$BeeGees_records == 0) missing_data_types <- c(missing_data_types, "BeeGees")
      if (result$debug_info$records_found$validation_records == 0) missing_data_types <- c(missing_data_types, "Validation")
      if (result$debug_info$records_found$fasta_sequences == 0) missing_data_types <- c(missing_data_types, "FASTA")
      
      if (total_records == 0) {
        # ERROR CASE: No data found
        values$summary_stats <- data.frame()
        values$barcode_validation <- data.frame()
        values$fasta_data <- list()
        values$taxval_process_ids <- character()
        values$current_identifier <- ""
        values$current_identifier_type <- ""
        values$current_process_ids <- character()
        values$dataset_loaded <- FALSE
        
        error_msg <- paste("ERROR: No data files found for", identifier_type, identifier, 
                           "- Searched", length(target_process_ids), "Process IDs but found 0 records")
        
        insertUI("#plate_status", "beforeEnd",
                 div(class = "status-message status-error",
                     p(error_msg, style = "margin: 0;")))
        
        showNotification("Failed to load data - no data found.", type = "error")
        
      } else {
        # SUCCESS CASE: Data found
        shinyjs::show("loading_section")
        
        # Simulate progress updates
        session$sendCustomMessage("updateProgress", list(percent = 33, message = "Processing BeeGees data"))
        Sys.sleep(0.3)
        session$sendCustomMessage("updateProgress", list(percent = 66, message = "Merging Validation data"))
        Sys.sleep(0.3)
        session$sendCustomMessage("updateProgress", list(percent = 100, message = "Loading FASTA sequences"))
        Sys.sleep(0.3)
        session$sendCustomMessage("updateProgress", list(percent = 100, message = "Complete"))
        
        # Hide loading and update values
        shinyjs::hide("loading_section")
        
        values$summary_stats <- result$summary_stats
        values$barcode_validation <- result$barcode_validation
        values$fasta_data <- result$fasta_data
        values$taxval_process_ids <- result$taxval_process_ids
        values$current_identifier <- identifier
        values$current_identifier_type <- identifier_type
        values$current_process_ids <- target_process_ids
        values$dataset_loaded <- TRUE
        
        success_msg <- paste("Successfully loaded", identifier_type, ":", identifier, 
                             "- Process IDs:", length(target_process_ids))
        if (length(missing_data_types) > 0) {
          success_msg <- paste(success_msg, "- WARNING: Missing data types:", paste(missing_data_types, collapse = ", "))
        }
        
        insertUI("#plate_status", "beforeEnd",
                 div(class = "status-message status-success",
                     p(success_msg, style = "margin: 0;")))
        
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
  
  # Prepare outcome data with validation logic
  prepare_outcome_data <- reactive({
    req(values$dataset_loaded)
    
    # Get all unique process IDs from target list
    all_process_ids <- unique(values$current_process_ids)
    
    if (length(all_process_ids) == 0) {
      return(data.frame())
    }
    
    # Structural validation PASS = Process ID exists in taxval data
    # Taxonomic validation PASS = Process ID exists in FASTA data
    structural_pass_ids <- values$taxval_process_ids
    taxonomic_pass_ids <- if (length(values$fasta_data$process_ids) > 0) {
      unique(values$fasta_data$process_ids)
    } else {
      character()
    }
    
    # Clean and compare process IDs
    all_process_ids_clean <- trimws(as.character(all_process_ids))
    structural_pass_ids_clean <- trimws(as.character(structural_pass_ids))
    taxonomic_pass_ids_clean <- trimws(as.character(taxonomic_pass_ids))
    
    # Create outcome data for all process IDs
    outcome_data <- data.frame(
      process_id = all_process_ids_clean,
      structural_pass = all_process_ids_clean %in% structural_pass_ids_clean,
      taxonomic_pass = all_process_ids_clean %in% taxonomic_pass_ids_clean,
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        # Determine overall status
        overall_status = case_when(
          structural_pass & taxonomic_pass ~ "PASS",
          structural_pass & !taxonomic_pass ~ "PARTIAL",
          TRUE ~ "FAIL"
        ),
        
        # Create display values
        structural_status = ifelse(structural_pass, "✓ PASS", "✗ FAIL"),
        taxonomic_status = ifelse(taxonomic_pass, "✓ PASS", "✗ FAIL"),
        
        # Create BOLD link for Process ID
        process_id_link = paste0('<a href="https://portal.boldsystems.org/record/', process_id, '" target="_blank">', process_id, '</a>')
      ) %>%
      select(
        `Process ID` = process_id_link,
        `Structural Validation` = structural_status,
        `Taxonomic Validation` = taxonomic_status,
        `Overall Status` = overall_status,
        overall_status_code = overall_status
      )
    
    return(outcome_data)
  })
  
  # Function to filter and reorder barcode validation columns with sequence wrapping
  prepare_filtered_barcode_validation <- reactive({
    req(values$dataset_loaded, nrow(values$barcode_validation) > 0)
    
    # Define the desired column order - sequence after reading_frame
    desired_columns <- c("Filename", "process_id", "fasta_header", "mge_params", 
                         "marker_code", "nuc_basecount", "nuc_full_basecount", 
                         "ambig_original", "ambig_basecount", 
                         "ambig_full_basecount", "stop_codons", "reading_frame",
                         "sequence", # MOVED HERE - after reading_frame
                         "identification", "species", "identification_method", 
                         "identification_rank", "obs_taxon", "BOLD_submissions", "error")
    
    # Get available columns from the data
    available_columns <- names(values$barcode_validation)
    
    # Find which desired columns are actually available
    columns_to_select <- desired_columns[desired_columns %in% available_columns]
    
    if (length(columns_to_select) == 0) {
      # If none of the desired columns are available, return the original data
      return(values$barcode_validation)
    }
    
    # Select and reorder the columns
    filtered_data <- values$barcode_validation %>%
      select(all_of(columns_to_select))
    
    # Apply sequence wrapping to the 'sequence' column if it exists
    if ("sequence" %in% names(filtered_data)) {
      filtered_data <- filtered_data %>%
        mutate(sequence = sapply(sequence, wrap_sequence, USE.NAMES = FALSE))
    }
    
    return(filtered_data)
  })
  
  # Render outcome table
  output$outcomeTable <- renderDT({
    outcome_data <- prepare_outcome_data()
    
    if (nrow(outcome_data) == 0) {
      return(datatable(data.frame("No data available" = "Please check that validation data is loaded correctly.")))
    }
    
    display_data <- outcome_data %>% select(-overall_status_code)
    
    datatable(
      display_data,
      escape = FALSE,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list(
          list(extend = 'copy', exportOptions = list(modifier = list(search = 'applied'))),
          list(extend = 'csv', exportOptions = list(modifier = list(search = 'applied'))),
          list(extend = 'excel', exportOptions = list(modifier = list(search = 'applied')))
        ),
        columnDefs = list(
          list(className = 'dt-center', targets = c(1, 2, 3))
        ),
        rowCallback = JS(
          "function(row, data, index) {",
          "  var status = data[3];",
          "  if (status === 'PASS') {",
          "    $(row).css('background-color', '#d4edda');",
          "  } else if (status === 'PARTIAL') {",
          "    $(row).css('background-color', '#fff3cd');",
          "  } else if (status === 'FAIL') {",
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
    pass_count <- sum(outcome_data$overall_status_code == "PASS", na.rm = TRUE)
    partial_count <- sum(outcome_data$overall_status_code == "PARTIAL", na.rm = TRUE)
    fail_count <- sum(outcome_data$overall_status_code == "FAIL", na.rm = TRUE)
    
    # Pass rate only counts full PASS, not PARTIAL
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
              div(style = "font-size: 1.5em; font-weight: bold; color: #856404;", partial_count),
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
  
  # Prepare FASTA data
  prepare_fasta_data <- reactive({
    req(values$dataset_loaded, length(values$fasta_data$process_ids) > 0)
    
    # Filter FASTA sequences for target process IDs only
    filtered_sequences <- filter_by_process_ids(values$fasta_data$sequences, values$current_process_ids)
    
    if (nrow(filtered_sequences) == 0) {
      return("No FASTA sequences found for the selected Process IDs.")
    }
    
    # Create FASTA format
    fasta_lines <- c()
    for (i in 1:nrow(filtered_sequences)) {
      header <- filtered_sequences$header[i]
      sequence <- filtered_sequences$sequence[i]
      
      # Add header
      fasta_lines <- c(fasta_lines, paste0(">", header))
      
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
  
  # Render tables with updated export options
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
        buttons = list(
          list(extend = 'copy', exportOptions = list(modifier = list(search = 'applied'))),
          list(extend = 'csv', exportOptions = list(modifier = list(search = 'applied'))),
          list(extend = 'excel', exportOptions = list(modifier = list(search = 'applied')))
        ),
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
  
  # Render barcode validation table with filtered columns and sequence wrapping
  output$barcodeValidationTable <- renderDT({
    req(values$dataset_loaded, nrow(values$barcode_validation) > 0)
    
    # Use the filtered and reordered data with sequence wrapping
    display_data <- prepare_filtered_barcode_validation()
    
    numeric_cols <- which(sapply(display_data, is.numeric))
    
    datatable(
      display_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list(
          list(extend = 'copy', exportOptions = list(modifier = list(search = 'applied'))),
          list(extend = 'csv', exportOptions = list(modifier = list(search = 'applied'))),
          list(extend = 'excel', exportOptions = list(modifier = list(search = 'applied')))
        ),
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
  
  # UPDATED: Column selectors for summary stats - now includes ALL numeric columns
  observe({
    if (values$dataset_loaded && nrow(values$summary_stats) > 0) {
      numeric_cols <- names(values$summary_stats)[sapply(values$summary_stats, is.numeric)]
      
      # Debug: Print available columns to console
      message("DEBUG: Available columns in summary_stats:")
      message(paste("All columns:", paste(names(values$summary_stats), collapse = ", ")))
      message(paste("Numeric columns:", paste(numeric_cols, collapse = ", ")))
      
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "xcol", choices = numeric_cols, selected = numeric_cols[1])
        updateSelectInput(session, "ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
      }
    }
  })
  
  # Update column selectors for barcode validation
  observe({
    if (values$dataset_loaded && nrow(values$barcode_validation) > 0) {
      filtered_data <- prepare_filtered_barcode_validation()
      numeric_cols <- names(filtered_data)[sapply(filtered_data, is.numeric)]
      
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "validation_ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
      }
      
      # For X-axis, use specific columns as requested
      x_axis_options <- c("process_id", "mge_params", "identification", "species")
      available_x_options <- x_axis_options[x_axis_options %in% names(filtered_data)]
      
      if (length(available_x_options) > 0) {
        updateSelectInput(session, "validation_xcol", choices = available_x_options, selected = available_x_options[1])
      } else {
        # Fallback to fasta_header if none of the preferred options are available
        if ("fasta_header" %in% names(filtered_data)) {
          updateSelectInput(session, "validation_xcol", choices = "fasta_header", selected = "fasta_header")
        }
      }
    }
  })
  
  # Render summary stats plot (filter NAs)
  output$summaryStatsPlot <- renderPlotly({
    req(input$xcol, input$ycol, values$dataset_loaded, nrow(values$summary_stats) > 0)
    
    plot_data <- values$summary_stats
    
    # Simple check that both columns exist and are numeric
    if (!input$xcol %in% names(plot_data) || !input$ycol %in% names(plot_data)) {
      return(plot_ly() %>%
               add_text(x = 0.5, y = 0.5, text = "Selected columns not available", 
                        textposition = "middle center") %>%
               layout(
                 xaxis = list(title = input$xcol, range = c(0, 1), showticklabels = FALSE),
                 yaxis = list(title = input$ycol, range = c(0, 1), showticklabels = FALSE),
                 showlegend = FALSE
               ))
    }
    
    # Filter out rows with NA values in the selected columns
    plot_data <- plot_data %>%
      filter(!is.na(.data[[input$xcol]]) & !is.na(.data[[input$ycol]]))
    
    if (nrow(plot_data) == 0) {
      return(plot_ly() %>%
               add_text(x = 0.5, y = 0.5, text = "No valid data points available", 
                        textposition = "middle center") %>%
               layout(
                 xaxis = list(title = input$xcol, range = c(0, 1), showticklabels = FALSE),
                 yaxis = list(title = input$ycol, range = c(0, 1), showticklabels = FALSE),
                 showlegend = FALSE
               ))
    }
    
    plot_ly(
      data = plot_data,
      x = ~.data[[input$xcol]],
      y = ~.data[[input$ycol]],
      type = 'scatter',
      mode = 'markers',
      marker = list(size = 10, color = 'rgba(51, 122, 183, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
      text = ~paste("Process ID:", process_id, "<br>Fasta Header:", fasta_header),
      textposition = "none",
      hovertemplate = "%{text}<br>%{xaxis.title.text}: %{x}<br>%{yaxis.title.text}: %{y}<extra></extra>"
    ) %>%
      layout(
        xaxis = list(title = input$xcol),
        yaxis = list(title = input$ycol),
        margin = list(t = 30)
      )
  })
  
  # Render barcode validation bar chart
  output$barcodeValidationPlot <- renderPlotly({
    req(input$validation_xcol, input$validation_ycol, values$dataset_loaded, nrow(values$barcode_validation) > 0)
    
    filtered_data <- prepare_filtered_barcode_validation()
    
    # Remove HTML tags from the data for plotting
    plot_data <- filtered_data
    if (input$validation_xcol == "fasta_header" && "fasta_header" %in% names(plot_data)) {
      plot_data$fasta_header <- gsub("<.*?>", "", plot_data$fasta_header)
    }
    
    # Create bar chart with config for scrolling/zooming
    p <- plot_ly(
      data = plot_data,
      x = ~.data[[input$validation_xcol]],
      y = ~.data[[input$validation_ycol]],
      type = 'bar',
      marker = list(color = 'rgba(220, 53, 69, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
      text = ~paste("Process ID:", process_id, "<br>Fasta Header:", fasta_header),
      textposition = "none",
      hovertemplate = "%{text}<br>%{xaxis.title.text}: %{x}<br>%{yaxis.title.text}: %{y}<extra></extra>"
    ) %>%
      layout(
        xaxis = list(
          title = input$validation_xcol,
          tickangle = -45
        ),
        yaxis = list(title = input$validation_ycol),
        margin = list(t = 30, b = 100),
        dragmode = "pan"
      ) %>%
      config(scrollZoom = TRUE)
    
    return(p)
  })
  
  # Download handlers with updated logic for filtered data
  output$downloadOutcomeTable <- downloadHandler(
    filename = function() {
      paste0("barcoding_outcome_", values$current_identifier, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      outcome_data <- prepare_outcome_data()
      if (nrow(outcome_data) > 0) {
        export_data <- outcome_data %>% 
          select(-overall_status_code) %>%
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
  
  # Download handler for filtered barcode validation data (unwrap sequences for export)
  output$downloadBarcodeValidation <- downloadHandler(
    filename = function() {
      paste0("barcode_validation_filtered_", values$current_identifier, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      filtered_data <- prepare_filtered_barcode_validation()
      
      # Remove line breaks from sequence column for export
      if ("sequence" %in% names(filtered_data)) {
        export_data <- filtered_data %>%
          mutate(sequence = gsub("\n", "", sequence))
      } else {
        export_data <- filtered_data
      }
      
      write.csv(export_data, file, row.names = FALSE)
    }
  )
  
  output$downloadFasta <- downloadHandler(
    filename = function() {
      paste0("barcodes_", values$current_identifier, "_", Sys.Date(), ".fasta")
    },
    content = function(file) {
      fasta_content <- prepare_fasta_data()
      writeLines(fasta_content, file)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)