#install.packages(c("shiny", "DT", "jsonlite", "dplyr", "tidyr"))
#install.packages("openxlsx")
#install.packages("plotly")

library(shiny)
library(DT)
library(jsonlite)
library(dplyr)
library(tidyr)
library(plotly)
library(shinyjs)

# Base data directory
BASE_DATA_DIR <- "./data"
BGE_PLATES_PATH <- file.path(BASE_DATA_DIR, "bge_plates_ids.csv")

# Load BGE plates data
load_bge_plates <- function() {
  if (file.exists(BGE_PLATES_PATH)) {
    tryCatch({
      bge_data <- read.csv(BGE_PLATES_PATH, stringsAsFactors = FALSE)
      # Clean column names in case of spaces/formatting issues
      names(bge_data) <- trimws(names(bge_data))
      return(bge_data)
    }, error = function(e) {
      message("Error loading BGE plates file:", e$message)
      return(data.frame())
    })
  } else {
    message("BGE plates file not found:", BGE_PLATES_PATH)
    return(data.frame())
  }
}

# Search BGE plates data
search_bge_plates <- function(bge_data, plate_search = "", sample_search = "", process_search = "") {
  if (nrow(bge_data) == 0) return(data.frame())
  
  # Start with all data
  filtered_data <- bge_data
  
  # Apply filters based on non-empty search terms
  if (plate_search != "" && !is.null(plate_search)) {
    filtered_data <- filtered_data[grepl(plate_search, filtered_data$Plate_Number, ignore.case = TRUE), ]
  }
  
  if (sample_search != "" && !is.null(sample_search)) {
    filtered_data <- filtered_data[grepl(sample_search, filtered_data$Sample_ID, ignore.case = TRUE), ]
  }
  
  if (process_search != "" && !is.null(process_search)) {
    filtered_data <- filtered_data[grepl(process_search, filtered_data$Process_ID, ignore.case = TRUE), ]
  }
  
  return(filtered_data)
}

get_available_datasets <- function() {
  if (!dir.exists(BASE_DATA_DIR)) {
    return(character(0))
  }
  
  dirs <- list.dirs(BASE_DATA_DIR, recursive = FALSE, full.names = FALSE)
  return(dirs[dirs != ""])
}

# Construct file paths for a given dataset
get_dataset_paths <- function(dataset_name) {
  base_path <- file.path(BASE_DATA_DIR, dataset_name)
  
  list(
    sum_stats = file.path(base_path, paste0(dataset_name, "_combined_stats.csv")),
    fcomp = file.path(base_path, paste0(dataset_name, "_fasta_compare.csv")),
    json_dir = file.path(base_path, paste0(dataset_name, "_fastp_json"))
  )
}

# Parse fastp JSON files
parse_json_file <- function(json_path, process_id) {
  tryCatch({
    data <- fromJSON(json_path)
    
    # Define null-coalescing operator if not available
    `%||%` <- function(x, y) if (is.null(x)) y else x
    
    result <- data.frame(
      process_id = process_id,
      # Before filtering
      before_total_reads = data$summary$before_filtering$total_reads %||% NA,
      before_total_bases = data$summary$before_filtering$total_bases %||% NA,
      before_q20_bases = data$summary$before_filtering$q20_bases %||% NA,
      before_q30_bases = data$summary$before_filtering$q30_bases %||% NA,
      before_q20_rate = data$summary$before_filtering$q20_rate %||% NA,
      before_q30_rate = data$summary$before_filtering$q30_rate %||% NA,
      before_read1_mean_length = data$summary$before_filtering$read1_mean_length %||% NA,
      before_read2_mean_length = data$summary$before_filtering$read2_mean_length %||% NA,
      before_gc_content = data$summary$before_filtering$gc_content %||% NA,
      # After filtering
      after_total_reads = data$summary$after_filtering$total_reads %||% NA,
      after_total_bases = data$summary$after_filtering$total_bases %||% NA,
      after_q20_bases = data$summary$after_filtering$q20_bases %||% NA,
      after_q30_bases = data$summary$after_filtering$q30_bases %||% NA,
      after_q20_rate = data$summary$after_filtering$q20_rate %||% NA,
      after_q30_rate = data$summary$after_filtering$q30_rate %||% NA,
      after_read1_mean_length = data$summary$after_filtering$read1_mean_length %||% NA,
      after_read2_mean_length = data$summary$after_filtering$read2_mean_length %||% NA,
      after_gc_content = data$summary$after_filtering$gc_content %||% NA,
      # Filtering results
      passed_filter_reads = data$filtering_result$passed_filter_reads %||% NA,
      low_quality_reads = data$filtering_result$low_quality_reads %||% NA,
      too_many_N_reads = data$filtering_result$too_many_N_reads %||% NA,
      too_short_reads = data$filtering_result$too_short_reads %||% NA,
      too_long_reads = data$filtering_result$too_long_reads %||% NA,
      # Duplication
      duplication_rate = data$duplication$rate %||% NA,
      stringsAsFactors = FALSE
    )
    return(result)
  }, error = function(e) {
    message(paste("Error parsing JSON for", process_id, ":", e$message))
    return(NULL)
  })
}

# Load dataset with progress updates
load_dataset_with_progress <- function(dataset_name, progress_callback = NULL) {
  if (is.null(dataset_name) || dataset_name == "") {
    return(list(
      summary_stats = data.frame(),
      fasta_compare = data.frame(),
      json_data = data.frame(),
      status = "No dataset selected"
    ))
  }
  
  # Initialise progress
  if (!is.null(progress_callback)) progress_callback(0, "Initializing...")
  
  paths <- get_dataset_paths(dataset_name)
  result <- list(
    summary_stats = data.frame(),
    fasta_compare = data.frame(),
    json_data = data.frame(),
    status = ""
  )
  
  # Load summary stats (33% progress)
  if (!is.null(progress_callback)) progress_callback(33, "Loading summary statistics...")
  Sys.sleep(0.5) # Brief pause to show progress
  
  if (file.exists(paths$sum_stats)) {
    tryCatch({
      result$summary_stats <- read.csv(paths$sum_stats, stringsAsFactors = FALSE)
      # Rename ID column to process_id if it exists
      if ("ID" %in% names(result$summary_stats)) {
        names(result$summary_stats)[names(result$summary_stats) == "ID"] <- "process_id"
      }
      message(paste("Loaded", nrow(result$summary_stats), "rows from summary stats"))
    }, error = function(e) {
      result$status <- paste(result$status, "Error loading summary stats:", e$message, "\n")
    })
  } else {
    result$status <- paste(result$status, "Summary stats file not found:", paths$sum_stats, "\n")
  }
  
  # Load FASTA compare (66% progress)
  if (!is.null(progress_callback)) progress_callback(66, "Loading FASTA compare data...")
  Sys.sleep(0.5) # Brief pause to show progress
  
  if (file.exists(paths$fcomp)) {
    tryCatch({
      result$fasta_compare <- read.csv(paths$fcomp, stringsAsFactors = FALSE)
      message(paste("Loaded", nrow(result$fasta_compare), "rows from FASTA compare"))
    }, error = function(e) {
      result$status <- paste(result$status, "Error loading FASTA compare:", e$message, "\n")
    })
  } else {
    result$status <- paste(result$status, "FASTA compare file not found:", paths$fcomp, "\n")
  }
  
  # Load JSON files (100% progress)
  if (!is.null(progress_callback)) progress_callback(90, "Loading JSON quality metrics...")
  Sys.sleep(0.5) # Brief pause to show progress
  
  if (dir.exists(paths$json_dir)) {
    tryCatch({
      json_files <- list.files(paths$json_dir, pattern = "\\.json$", full.names = TRUE)
      json_data_list <- list()
      
      total_json_files <- length(json_files)
      for (i in seq_along(json_files)) {
        json_file <- json_files[i]
        # Update progress for JSON processing
        json_progress <- 90 + (i / total_json_files) * 8 # 90-98%
        if (!is.null(progress_callback)) {
          progress_callback(json_progress, paste("Processing JSON file", i, "of", total_json_files))
        }
        
        # Extract process ID from filename
        process_id <- gsub("_fastp_report\\.json$", "", basename(json_file))
        
        parsed_data <- parse_json_file(json_file, process_id)
        if (!is.null(parsed_data)) {
          json_data_list[[length(json_data_list) + 1]] <- parsed_data
        }
      }
      
      if (length(json_data_list) > 0) {
        result$json_data <- bind_rows(json_data_list)
        message(paste("Loaded", nrow(result$json_data), "JSON records"))
      }
    }, error = function(e) {
      result$status <- paste(result$status, "Error loading JSON files:", e$message, "\n")
    })
  } else {
    result$status <- paste(result$status, "JSON directory not found:", paths$json_dir, "\n")
  }
  
  # Final step
  if (!is.null(progress_callback)) progress_callback(100, "Loading complete!")
  Sys.sleep(0.3) # Brief pause to show completion
  
  if (result$status == "") {
    result$status <- paste("Successfully loaded dataset:", dataset_name)
  }
  
  return(result)
}

# Define UI
ui <- fluidPage(
  useShinyjs(), # Enable shinyjs
  tags$head(
    tags$style(HTML("
      .dataTables_wrapper {
        font-size: 14px;
      }
      .download-btn {
        margin: 10px 0;
      }
      .dataset-selection {
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
      .sample-finder-section {
        background-color: #f8f9fa;
        padding: 20px;
        border-radius: 8px;
        margin-bottom: 20px;
      }
      .plate-card {
        background-color: white;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 15px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.05);
        transition: box-shadow 0.2s ease;
      }
      .plate-card:hover {
        box-shadow: 0 4px 8px rgba(0,0,0,0.1);
      }
      .plate-header {
        font-size: 1.2em;
        font-weight: bold;
        color: #495057;
        margin-bottom: 15px;
        padding-bottom: 10px;
        border-bottom: 2px solid #e9ecef;
      }
      .process-list {
        list-style: none;
        padding: 0;
        margin: 0;
      }
      .process-item {
        background-color: #f8f9fa;
        padding: 8px 12px;
        margin-bottom: 5px;
        border-radius: 4px;
        border-left: 3px solid #007bff;
        font-family: 'Courier New', monospace;
        font-size: 0.9em;
      }
      .process-item.highlighted {
        background-color: #fff3cd;
        border-left-color: #ffc107;
      }
      .sample-id {
        color: #6c757d;
        font-size: 0.85em;
      }
      .search-results-header {
        background-color: #e7f3ff;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 20px;
        border-left: 4px solid #007bff;
      }
      .no-results {
        text-align: center;
        color: #6c757d;
        font-style: italic;
        padding: 40px;
        background-color: #f8f9fa;
        border-radius: 8px;
        border: 2px dashed #dee2e6;
      }
      .process-search-section {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        margin-bottom: 20px;
      }
      .outcome-result {
        background-color: white;
        border-radius: 12px;
        padding: 25px;
        margin: 20px 0;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        border-left: 6px solid;
        transition: all 0.3s ease;
      }
      .outcome-result:hover {
        box-shadow: 0 6px 20px rgba(0,0,0,0.15);
        transform: translateY(-2px);
      }
      .outcome-success {
        border-left-color: #28a745;
        background: linear-gradient(135deg, #ffffff 0%, #f8fff9 100%);
      }
      .outcome-error {
        border-left-color: #dc3545;
        background: linear-gradient(135deg, #ffffff 0%, #fff8f8 100%);
      }
      .outcome-header {
        font-size: 1.3em;
        font-weight: bold;
        margin-bottom: 15px;
        display: flex;
        align-items: center;
        gap: 12px;
      }
      .outcome-icon {
        font-size: 1.8em;
        padding: 8px;
        border-radius: 50%;
        min-width: 45px;
        text-align: center;
        color: white;
        font-weight: bold;
      }
      .outcome-icon.success {
        background-color: #28a745;
      }
      .outcome-icon.error {
        background-color: #dc3545;
      }
      .outcome-details {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 8px;
        margin-top: 15px;
        border: 1px solid #e9ecef;
      }
      .outcome-stats {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
        gap: 15px;
        margin-top: 15px;
      }
      .outcome-stat {
        background-color: white;
        padding: 12px;
        border-radius: 6px;
        border: 1px solid #dee2e6;
        text-align: center;
      }
      .outcome-stat-value {
        font-size: 1.4em;
        font-weight: bold;
        color: #495057;
      }
      .outcome-stat-label {
        font-size: 0.9em;
        color: #6c757d;
        margin-top: 5px;
      }
    ")),
    # JavaScript for progress updates
    tags$script(HTML("
      Shiny.addCustomMessageHandler('updateProgress', function(data) {
        $('#progress_bar').css('width', data.percent + '%');
        $('#progress_text').text(data.message + ' (' + Math.round(data.percent) + '%)');
      });
    "))
  ),
  
  titlePanel("Barcoding Results Dashboard"),
  
  # Set up tabs for each section
  tabsetPanel(
    # Sample Finder Tab (1st)
    tabPanel("Sample Finder",
             br(),
             div(class = "sample-finder-section",
                 h3("🔍 Sample Finder"),
                 p("Search for samples across all datasets to find the information you need:"),
                 
                 # Search inputs
                 fluidRow(
                   column(4,
                          textInput("plate_search", 
                                    "Search by Plate Number:",
                                    placeholder = "e.g., BGE_00001")
                   ),
                   column(4,
                          textInput("sample_search", 
                                    "Search by Sample ID:",
                                    placeholder = "e.g., BGE_00001_A01")
                   ),
                   column(4,
                          textInput("process_search", 
                                    "Search by Process ID:",
                                    placeholder = "e.g., BSNHM002-24")
                   )
                 ),
                 
                 div(class = "search-info",
                     p("💡 ", strong("Tip:"), " Enter any search term above to find related samples. Use this information to identify which dataset to load in the next tab."),
                     style = "background-color: #e3f2fd; padding: 10px; border-radius: 5px; margin: 15px 0; color: #1565c0;"
                 )
             ),
             
             # Search results
             div(id = "search_results_container"),
             
             # Instructions
             br(),
             div(class = "instructions-panel",
                 h4("How to Use Sample Finder"),
                 tags$ul(
                   tags$li("Enter a ", strong("Plate Number"), " to see all samples from that plate"),
                   tags$li("Enter a ", strong("Sample ID"), " to find its Process ID and Plate"),
                   tags$li("Enter a ", strong("Process ID"), " to find its Sample ID and Plate"),
                   tags$li("Use the results to identify which dataset contains your samples"),
                   tags$li("Navigate to 'Dataset Selection' tab to load the appropriate dataset")
                 ),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; border-left: 4px solid #28a745;"
             )
    ),
    
    # Dataset Selection Tab (2nd)
    tabPanel("Dataset Selection",
             br(),
             div(class = "dataset-selection",
                 h3("Select Dataset"),
                 p("Choose a dataset from the available directories in the data folder:"),
                 fluidRow(
                   column(6,
                          selectInput("selected_dataset", 
                                      "Available Datasets:",
                                      choices = c("Loading..." = ""),
                                      width = "100%")
                   ),
                   column(6,
                          br(),
                          actionButton("load_dataset", "Load Dataset", 
                                       class = "btn-primary btn-lg",
                                       style = "margin-top: 5px;")
                   )
                 ),
                 # Loading progress section
                 div(id = "loading_section", style = "display: none; margin-top: 20px;",
                     div(class = "progress-container",
                         h4("Loading Dataset...", style = "color: #495057; margin-bottom: 15px;"),
                         div(class = "progress", style = "height: 25px; background-color: #e9ecef; border-radius: 12px; overflow: hidden;",
                             div(id = "progress_bar", class = "progress-bar progress-bar-striped progress-bar-animated", 
                                 style = "width: 0%; background-color: #007bff; transition: width 0.3s ease;")
                         ),
                         div(id = "progress_text", style = "text-align: center; margin-top: 10px; font-weight: 500; color: #495057;",
                             "Initializing...")
                     )
                 ),
                 div(id = "dataset_status", class = "status-message")
             ),
             br(),
             h4("Dataset Information"),
             p("Once you select and load a dataset, you can navigate to the other tabs to view:"),
             tags$ul(
               tags$li(strong("BGEE Summary Statistics:"), " Combined statistics from the Barcode Gene Extractor & Evaluator (BGEE) pipeline. "),
               tags$li(strong("Fastp Metrics:"), " Before and after read trimming statistics from Fastp (parsed from json files). "),
               tags$li(strong("FASTA Compare Results:"), " Barcode consensus sequence comparison and quality metrics. ")
             )
    ),
    
    # Barcoding Outcome Tab (3rd)
    tabPanel("Barcoding Outcome",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a dataset from the 'Dataset Selection' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               # Static description text area
               div(
                 h3("🚦 Barcoding Outcome Overview"),
                 p("Overview of BIN-compliant barcode extraction results for all Process IDs in the loaded dataset."),
                 style = "background-color: #f8f9fa; padding: 20px; border-radius: 8px; margin-bottom: 20px;"
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
    
    # Fastp Metrics Tab (4th)
    tabPanel("Fastp Metrics",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a dataset from the 'Dataset Selection' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               # Static description text area
               div(
                 h4("About Fastp Metrics"),
                 HTML("Fastp metrics are generated from JSON output files of the <a href='https://github.com/OpenGene/fastp' target='_blank'>fastp preprocessing tool</a>. These files contain detailed statistics about read quality before and after filtering (from the 'concat' mode of the BGEE pipeline)."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               # Field descriptions
               div(
                 class = "mb-3",
                 h4("Fastp JSON Field Descriptions"),
                 tags$ul(
                   tags$li(strong("Process ID:"), " Unique sample identifier."),
                   tags$li(strong("before_total_reads / after_total_reads:"), " Total reads before and after filtering."),
                   tags$li(strong("before_total_bases / after_total_bases:"), " Total bases before and after filtering."),
                   tags$li(strong("before_q20/q30_bases, after_q20/q30_bases:"), " Bases with quality ≥ Q20/Q30 before and after filtering."),
                   tags$li(strong("q20_rate / q30_rate:"), " Quality scores as percentages."),
                   tags$li(strong("read1_mean_length / read2_mean_length:"), " Mean PE read lengths before/after filtering."),
                   tags$li(strong("gc_content:"), " GC% content of reads before/after filtering."),
                   tags$li(strong("passed_filter_reads:"), " Reads that passed the quality filter."),
                   tags$li(strong("low_quality_reads / too_many_N_reads / too_short_reads / too_long_reads:"), " Now of reads failing filtering reasons."),
                   tags$li(strong("duplication_rate:"), " Estimated duplication rate in the reads.")
                 ),
                 style = "background-color: #fff; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;"
               ),
               
               downloadButton("downloadJsonData", "Download Filtered JSON Data", class = "download-btn btn-primary"),
               DTOutput("jsonDataTable"),
               
               # Bar chart with selectable y-axis
               hr(),
               h4("Process ID Overview"),
               fluidRow(
                 column(6, selectInput("json_bar_ycol", "Y Axis:", choices = NULL)),
                 column(6, div(style = "height: 34px;")) # Spacer for alignment
               ),
               plotlyOutput("jsonBarPlot", height = "400px"),
               
               # Interactive scatter plot
               hr(),
               h4("Interactive Scatter Plot"),
               fluidRow(
                 column(6, selectInput("json_xcol", "X Axis:", choices = NULL)),
                 column(6, selectInput("json_ycol", "Y Axis:", choices = NULL))
               ),
               plotlyOutput("jsonDataPlot", height = "500px")
             )
    ),
    
    # BGEE Summary Statistics Tab (5th)
    tabPanel("BGEE Summary Statistics",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a dataset from the 'Dataset Selection' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               # Static description text area
               div(
                 h4("About BGEE Summary Statistics"),
                 HTML("<a href='https://github.com/bge-barcoding/BGEE' target='_blank'>BGEE (Barcode Gene Extractor & Evaluator)</a> summary statistics are generated from the combined output of the MGE (MitoGeneExtractor) and fasta_cleaner steps of both 'concat' and 'merge' pre-processing modes of the workflow. The BGEE summary statistics include read processing metrics, barcode consensus sequence coverage information, and quality assessments. The BGEE workflow coordinates all elements of data pre-processing, barcode sequence recovery, and post-processing. Six MGE parameter combinations ('r' and 's' values) are utilised for each sample"),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               # Field descriptions
               div(
                 class = "mb-3",
                 h4("Summary Statistics Field Descriptions"),
                 tags$ul(
                   tags$li(strong("Filename:"), " Name of file/sequence header, including the Process (ID) and MGE parameters (mge_params)."),
                   tags$li(strong("process_id:"), " Unique Process identifier (ID) for the sample."),
                   tags$li(strong("mge_params:"), " Parameter string parsed from Filename. See MitoGeneExtractor GitHub repo for explanation or 'r' and 's'"),
                   tags$li(strong("n_reads_in:"), " Total number of reads input to, and processed by, MGE."),
                   tags$li(strong("n_reads_aligned:"), " Total number of reads aligned to the reference sequence by MGE."),
                   tags$li(strong("n_reads_skipped:"), " Total number of reads skipped by MGE, either due to low 'r', secondary alignments, frameshifts, etc."),
                   tags$li(strong("ref_length:"), " Length of reference sequence."),
                   tags$li(strong("cov_min:"), " Minimum coverage (depth) observed in the alignment."),
                   tags$li(strong("cov_max:"), " Maximum coverage (depth) observed in the alignment."),
                   tags$li(strong("cov_avg:"), " Mean coverage (depth) observed in the alignment."),
                   tags$li(strong("cov_med:"), " Median coverage (depth) observed in the alignment."),
                   tags$li(strong("cleaning_input_reads:"), " Total number of reads input into fasta_cleaner (from MGE alignment file)."),
                   tags$li(strong("cleaning_kept_reads:"), " Total number of reads kept by fasta_cleaner after filtering."),
                   tags$li(strong("cleaning_removed_human:"), " Total number of reads removed due to surpassing the sequence similarlity threshold to human COI."),
                   tags$li(strong("cleaning_removed_at:"), " Total number of reads removed due to surpassing the sequence AT% threshold."),
                   tags$li(strong("cleaning_ambig_bases:"), " Total number of ambiguous bases (i.e. not G, T, C, or G) in consensus sequence."),
                   tags$li(strong("cleaning_cov_percent:"), " Coverage percentage of reference (i.e. How much of the reference sequence had any coverage at all."),
                   tags$li(strong("cleaning_cov_avg:"), " Mean coverage (depth) observed in the alignment after cleaning. "),
                   tags$li(strong("cleaning_cov_max:"), " Maximum coverage (depth) observed in the alignment after cleaning. "),
                   tags$li(strong("cleaning_cov_min:"), " Minimum coverage (depth) observed in the alignment after cleaning. ")
                 ),
                 style = "background-color: #fff; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;"
               ),
               
               downloadButton("downloadSummaryStats", "Download Filtered Summary Stats", class = "download-btn btn-primary"),
               DTOutput("summaryStatsTable"),
               hr(),
               h4("Interactive Scatter Plot"),
               fluidRow(
                 column(6, selectInput("ycol", "Y Axis:", choices = NULL)),
                 column(6, selectInput("xcol", "X Axis:", choices = NULL))
               ),
               plotlyOutput("summaryStatsPlot", height = "500px"),
               
               # Bar chart for n_reads_in by Process ID with search functionality
               hr(),
               h4("Reads Processed by BGEE for 'concat' and 'merge' mode, for each Process ID"),
               
               # Process ID search section
               div(class = "process-search-section",
                   fluidRow(
                     column(6,
                            textInput("process_id_search", 
                                      "Search Process ID(s):",
                                      placeholder = "e.g., BSNHM001 (partial matching supported)",
                                      value = "")
                     ),
                     column(6,
                            div(style = "margin-top: 25px;",
                                span("💡 Tip: Leave empty to show all Process IDs, or enter partial text to filter", 
                                     style = "color: #6c757d; font-size: 0.9em;")
                            )
                     )
                   )
               ),
               
               plotlyOutput("summaryBarPlot", height = "400px")
             )
    ),
    
    # FASTA Compare Results Tab (6th)
    tabPanel("FASTA Compare Results",
             br(),
             conditionalPanel(
               condition = "output.dataset_loaded == false",
               div(p("Please select and load a dataset from the 'Dataset Selection' tab first."),
                   style = "color: #6c757d; font-style: italic;")
             ),
             conditionalPanel(
               condition = "output.dataset_loaded == true",
               
               # Static description text area
               div(
                 h4("About FASTA Compare Results"),
                 HTML("<a href='https://github.com/bge-barcoding/fasta_compare' target='_blank'>FASTA Compare</a> is the final post-processing step of the BGEE workflow. Results are generated by analysing and ranking consensus sequences from the 'concat' and 'merge' modes, as well as fasta_cleaner-generated barcode consensus sequences. The tool evaluates sequence quality based on factors like gaps, ambiguous bases, and barcode region contiguity to select the best barcode sequence for each process ID."),
                 style = "background-color: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; color: #495057;"
               ),
               
               # Field descriptions
               div(
                 class = "mb-3",
                 h4("FASTA Compare Field Descriptions"),
                 tags$ul(
                   tags$li(strong("file:"), " Source FASTA file path"),
                   tags$li(strong("process_id:"), " Extracted process identifier (e.g. BSNHM001-24)"),
                   tags$li(strong("parameters:"), " Extracted BGEE run parameters (e.g. r_1_s_100)"),
                   tags$li(strong("seq_id:"), " Full sequence identifier from FASTA header"),
                   tags$li(strong("length:"), " Full sequence length (including gaps and ambiguous bases)"),
                   tags$li(strong("leading_gaps:"), " Count of leading gap ('-') characters"),
                   tags$li(strong("trailing_gaps:"), " Count of trailing gap ('-') characters"),
                   tags$li(strong("internal_gaps:"), " Count of internal gap ('-') characters"),
                   tags$li(strong("ambiguous_bases:"), " Count of ambiguous bases ('N')"),
                   tags$li(strong("longest_stretch:"), " Longest continuous stretch without gaps or ambiguous bases"),
                   tags$li(strong("barcode_length:"), " Length of barcode region (fixed by --target)"),
                   tags$li(strong("barcode_ambiguous_bases:"), " Count of 'N' characters in barcode region"),
                   tags$li(strong("barcode_longest_stretch:"), " Longest continuous barcode stretch without gaps or 'N's"),
                   tags$li(strong("barcode_rank:"), " Barcode quality rank (1 = best, 6 = worst)"),
                   tags$li(strong("full_rank:"), " Full sequence quality rank (1 = best, 3 = worst)"),
                   tags$li(strong("best_sequence:"), " TRUE if this is the best sequence for the process_id"),
                   tags$li(strong("selected_full_fasta:"), " TRUE if this was selected for output to full FASTA"),
                   tags$li(strong("selected_barcode_fasta:"), " TRUE if this was selected for output to barcode FASTA")
                 ),
                 style = "background-color: #fff; padding: 15px; border-radius: 8px; border: 1px solid #dee2e6; margin-bottom: 20px;"
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
  
  # Reactive values to store loaded data
  values <- reactiveValues(
    summary_stats = data.frame(),
    fasta_compare = data.frame(),
    json_data = data.frame(),
    dataset_loaded = FALSE,
    current_dataset = "",
    bge_plates = data.frame(),
    search_results = data.frame()
  )
  
  # Load BGE plates data on startup
  observe({
    values$bge_plates <- load_bge_plates()
    if (nrow(values$bge_plates) == 0) {
      showNotification("Warning: BGE plates file not found or empty. Sample Finder may not work correctly.", 
                       type = "warning", duration = 5)
    }
  })
  
  # Real-time search functionality
  observe({
    # Trigger search when any search field changes
    plate_term <- input$plate_search %||% ""
    sample_term <- input$sample_search %||% ""
    process_term <- input$process_search %||% ""
    
    # Only search if at least one field has content
    if (plate_term != "" || sample_term != "" || process_term != "") {
      values$search_results <- search_bge_plates(values$bge_plates, plate_term, sample_term, process_term)
      render_search_results()
    } else {
      # Clear results when all fields are empty
      values$search_results <- data.frame()
      removeUI("#search_results_container > div")
    }
  })
  
  # Render search results as cards
  render_search_results <- function() {
    # Clear existing results
    removeUI("#search_results_container > div")
    
    if (nrow(values$search_results) == 0) {
      # Show no results message
      insertUI("#search_results_container", "beforeEnd",
               div(class = "no-results",
                   p("🔍 No results found"),
                   p("Try adjusting your search terms or check the spelling.")
               ))
      return()
    }
    
    # Group results by plate
    grouped_results <- values$search_results %>%
      group_by(Plate_Number) %>%
      summarise(
        samples = list(data.frame(Sample_ID = Sample_ID, Process_ID = Process_ID)),
        .groups = 'drop'
      )
    
    # Create header
    total_plates <- nrow(grouped_results)
    total_samples <- nrow(values$search_results)
    
    insertUI("#search_results_container", "beforeEnd",
             div(class = "search-results-header",
                 h4(paste("🔍 Search Results:", total_plates, "plate(s),", total_samples, "sample(s)")),
                 style = "margin-bottom: 20px;"
             ))
    
    # Create cards for each plate
    for (i in 1:nrow(grouped_results)) {
      plate_num <- grouped_results$Plate_Number[i]
      plate_samples <- grouped_results$samples[[i]]
      
      # Determine which items to highlight
      plate_term <- input$plate_search %||% ""
      sample_term <- input$sample_search %||% ""
      process_term <- input$process_search %||% ""
      
      # Create process items
      process_items <- list()
      for (j in 1:nrow(plate_samples)) {
        sample_id <- plate_samples$Sample_ID[j]
        process_id <- plate_samples$Process_ID[j]
        
        # Check if this item should be highlighted
        is_highlighted <- FALSE
        if (sample_term != "" && grepl(sample_term, sample_id, ignore.case = TRUE)) {
          is_highlighted <- TRUE
        }
        if (process_term != "" && grepl(process_term, process_id, ignore.case = TRUE)) {
          is_highlighted <- TRUE
        }
        
        item_class <- if (is_highlighted) "process-item highlighted" else "process-item"
        
        process_items[[j]] <- div(class = item_class,
                                  span(paste("•", process_id), style = "font-weight: 500;"),
                                  span(paste(" (", sample_id, ")", sep = ""), class = "sample-id")
        )
      }
      
      # Create the plate card
      card_content <- div(class = "plate-card",
                          div(class = "plate-header",
                              span("📋 Plate: ", style = "color: #007bff;"),
                              span(plate_num, style = "font-family: 'Courier New', monospace;")
                          ),
                          div(style = "margin-bottom: 10px;",
                              strong("Process IDs:")
                          ),
                          div(class = "process-list", process_items)
      )
      
      insertUI("#search_results_container", "beforeEnd", card_content)
    }
  }
  
  # Initialise available datasets
  observe({
    datasets <- get_available_datasets()
    if (length(datasets) > 0) {
      choices <- setNames(datasets, datasets)
      choices <- c("Select a dataset..." = "", choices)
    } else {
      choices <- c("No datasets found" = "")
    }
    
    updateSelectInput(session, "selected_dataset", choices = choices)
  })
  
  # Create output for conditionalPanel
  output$dataset_loaded <- reactive({
    values$dataset_loaded
  })
  outputOptions(output, "dataset_loaded", suspendWhenHidden = FALSE)
  
  # Load dataset when button is clicked
  observeEvent(input$load_dataset, {
    if (input$selected_dataset == "" || is.null(input$selected_dataset)) {
      showNotification("Please select a dataset first.", type = "warning")
      return()
    }
    
    # Disable the load button during loading
    updateActionButton(session, "load_dataset", label = "Loading...", icon = NULL)
    shinyjs::disable("load_dataset")
    
    # Hide status message and show loading section
    removeUI("#dataset_status div")
    shinyjs::show("loading_section")
    
    # Create a progress update function
    update_progress <- function(percent, message) {
      session$sendCustomMessage("updateProgress", list(percent = percent, message = message))
    }
    
    # Use future/promises for async loading (if available) or regular loading
    tryCatch({
      # Load the dataset with progress updates
      dataset_result <- load_dataset_with_progress(input$selected_dataset, update_progress)
      
      # Update reactive values
      values$summary_stats <- dataset_result$summary_stats
      values$fasta_compare <- dataset_result$fasta_compare
      values$json_data <- dataset_result$json_data
      values$dataset_loaded <- TRUE
      values$current_dataset <- input$selected_dataset
      
      # Hide loading section and show final status
      shinyjs::hide("loading_section")
      
      # Update status message
      status_class <- if (grepl("Successfully", dataset_result$status)) "status-success" else "status-error"
      insertUI("#dataset_status", "beforeEnd",
               div(class = paste("status-message", status_class),
                   p(dataset_result$status, style = "margin: 0;")))
      
      showNotification(paste("Dataset", input$selected_dataset, "loaded successfully!"), type = "message")
      
    }, error = function(e) {
      # Hide loading section on error
      shinyjs::hide("loading_section")
      
      # Show error message
      insertUI("#dataset_status", "beforeEnd",
               div(class = "status-message status-error",
                   p(paste("Error loading dataset:", e$message), style = "margin: 0;")))
      
      showNotification("Failed to load dataset. Please try again.", type = "error")
    })
    
    # Re-enable the load button
    updateActionButton(session, "load_dataset", label = "Load Dataset", icon = NULL)
    shinyjs::enable("load_dataset")
  })
  
  
  
  # ===== BARCODING OUTCOME TAB SERVER LOGIC =====
  # Reactive function to prepare outcome data
  # Reactive function to prepare outcome data
  prepare_outcome_data <- reactive({
    req(values$dataset_loaded, nrow(values$fasta_compare) > 0)
    
    if (!"process_id" %in% names(values$fasta_compare)) {
      return(data.frame())
    }
    
    # Group by process_id and get outcome information
    outcome_data <- values$fasta_compare %>%
      group_by(process_id) %>%
      summarise(
        # Check if any sequence was selected for barcode FASTA
        has_selected_barcode = any(tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
        
        # Get information from the selected sequence (where selected_barcode_fasta == yes)
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
        
        selected_seq_id = ifelse(
          any(tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          seq_id[tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t")][1],
          NA
        ),
        
        selected_parameters = ifelse(
          any(tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          parameters[tolower(as.character(selected_barcode_fasta)) %in% c("true", "yes", "1", "t")][1],
          NA
        ),
        
        # Fallback: Prioritize best_sequence = "yes", then min barcode_rank
        best_barcode_rank = ifelse(
          any(tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          barcode_rank[tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t")][1],
          min(barcode_rank, na.rm = TRUE)
        ),
        
        # Fallback: Get barcode_longest_stretch for the best sequence
        fallback_barcode_longest_stretch = ifelse(
          any(tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t"), na.rm = TRUE),
          barcode_longest_stretch[tolower(as.character(best_sequence)) %in% c("true", "yes", "1", "t")][1],
          barcode_longest_stretch[barcode_rank == min(barcode_rank, na.rm = TRUE)][1]
        ),
        
        .groups = 'drop'
      ) %>%
      # Handle infinite values
      mutate(
        best_barcode_rank = ifelse(is.infinite(best_barcode_rank), NA, best_barcode_rank),
        selected_barcode_longest_stretch = ifelse(is.na(selected_barcode_longest_stretch), 0, selected_barcode_longest_stretch),
        fallback_barcode_longest_stretch = ifelse(is.na(fallback_barcode_longest_stretch), 0, fallback_barcode_longest_stretch)
      )
    
    # Add coverage information
    if (nrow(values$summary_stats) > 0 && nrow(outcome_data) > 0) {
      # Create coverage lookup
      coverage_data <- data.frame()
      
      for (i in 1:nrow(outcome_data)) {
        process_id <- outcome_data$process_id[i]
        seq_id <- outcome_data$selected_seq_id[i]
        parameters <- outcome_data$selected_parameters[i]
        has_selected <- outcome_data$has_selected_barcode[i]
        
        coverage_avg <- NA
        
        # Only look up coverage if there's a selected barcode
        if (has_selected && !is.na(seq_id)) {
          # Clean seq_id to match Filename in summary_stats
          # Remove _fcleaner and _merge suffixes to get base filename
          target_filename <- seq_id
          target_filename <- gsub("_fcleaner_merge$", "", target_filename)
          target_filename <- gsub("_fcleaner$", "", target_filename)
          target_filename <- gsub("_merge$", "", target_filename)
          
          # Find matching row in summary stats using cleaned filename
          matching_rows <- values$summary_stats[values$summary_stats$Filename == target_filename, ]
          
          if (nrow(matching_rows) > 0) {
            matched_row <- matching_rows[1, ]
            
            # Check if parameters contain 'fcleaner'
            if (!is.na(parameters) && grepl("fcleaner", parameters, ignore.case = TRUE)) {
              # Use cleaning coverage stats
              coverage_avg <- ifelse("cleaning_cov_avg" %in% names(matched_row), matched_row$cleaning_cov_avg, NA)
            } else {
              # Use regular coverage stats
              coverage_avg <- ifelse("cov_avg" %in% names(matched_row), matched_row$cov_avg, NA)
            }
          }
        }
        
        coverage_data <- rbind(coverage_data, data.frame(
          process_id = process_id,
          coverage_avg = coverage_avg,
          stringsAsFactors = FALSE
        ))
      }
      
      # Merge coverage data
      outcome_data <- outcome_data %>%
        left_join(coverage_data, by = "process_id")
    } else {
      # Add empty coverage column if no summary stats
      outcome_data$coverage_avg <- NA
    }
    
    # Determine outcome category and create display columns
    outcome_data <- outcome_data %>%
      mutate(
        # CORRECTED LOGIC: Amber for ANY rank >= 4, Red only for no selected barcode AND all ranks < 4
        outcome_category = case_when(
          has_selected_barcode & (is.na(selected_barcode_rank) | selected_barcode_rank <= 3) ~ "success",
          # Check if ANY barcode rank >= 4 (selected or fallback)
          (!is.na(selected_barcode_rank) & selected_barcode_rank >= 4) | 
            (!is.na(best_barcode_rank) & best_barcode_rank >= 4) ~ "partial",
          # Only Red if no selected barcode AND all ranks < 4
          !has_selected_barcode & (!is.na(best_barcode_rank) & best_barcode_rank < 4) ~ "failed",
          # Default fallback
          TRUE ~ "failed"
        ),
        
        # Create display columns with fallback info for Red status
        outcome_status = case_when(
          outcome_category == "success" ~ "✅ Success",
          outcome_category == "partial" ~ "⚠️ Partial",
          TRUE ~ "❌ Failed"
        ),
        
        # Display barcode info: selected if available, fallback in parentheses if Red
        display_barcode_stretch = case_when(
          has_selected_barcode ~ as.character(selected_barcode_longest_stretch),
          !has_selected_barcode ~ paste0("(", fallback_barcode_longest_stretch, ")"),
          TRUE ~ NA_character_
        ),
        
        display_barcode_rank = case_when(
          has_selected_barcode ~ as.character(selected_barcode_rank),
          !has_selected_barcode ~ paste0("(", best_barcode_rank, ")"),
          TRUE ~ NA_character_
        )
      ) %>%
      # Select and rename columns for display (removed Coverage % and Coverage Type)
      select(
        `Process ID` = process_id,
        `Status` = outcome_status,
        `Barcode Longest Stretch` = display_barcode_stretch,
        `Best Barcode Rank` = display_barcode_rank,
        `Coverage Avg` = coverage_avg,
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
    
    # Remove the outcome_category column for display
    display_data <- outcome_data %>% select(-outcome_category)
    
    # Find numeric columns in the actual display data
    numeric_columns <- which(sapply(display_data, is.numeric))
    
    datatable(
      display_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-center', targets = c(1)), # Status column centered
          list(className = 'dt-right', targets = numeric_columns - 1) # Numeric columns right-aligned (0-indexed)
        ),
        rowCallback = JS(
          "function(row, data, index) {",
          "  var status = data[1];", # Status is in column 1 (0-indexed)
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
      # Only format numeric columns that actually exist
      {if(length(numeric_columns) > 0) formatRound(., columns = numeric_columns, digits = 2) else .}
  })
  
  # Render outcome summary
  output$outcome_summary <- renderUI({
    outcome_data <- prepare_outcome_data()
    
    if (nrow(outcome_data) == 0) {
      return(div())
    }
    
    # Calculate summary statistics
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
  
  # Download handler for outcome table
  output$downloadOutcomeTable <- downloadHandler(
    filename = function() {
      paste0("barcoding_outcome_", values$current_dataset, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      outcome_data <- prepare_outcome_data()
      if (nrow(outcome_data) > 0) {
        # Remove outcome_category column for export
        export_data <- outcome_data %>% select(-outcome_category)
        write.csv(export_data, file, row.names = FALSE)
      } else {
        write.csv(data.frame("No data available"), file, row.names = FALSE)
      }
    }
  )
  
  # ===== END BARCODING OUTCOME TAB SERVER LOGIC =====
  
  # Render Summary Stats table
  output$summaryStatsTable <- renderDT({
    req(values$dataset_loaded, nrow(values$summary_stats) > 0)
    
    datatable(
      values$summary_stats,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-right', targets = which(sapply(values$summary_stats, is.numeric)) - 1)
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    ) %>%
      formatRound(columns = which(sapply(values$summary_stats, is.numeric)), digits = 2)
  })
  
  # Update choices for summary stats column selectors
  observe({
    if (values$dataset_loaded && nrow(values$summary_stats) > 0) {
      numeric_cols <- names(values$summary_stats)[sapply(values$summary_stats, is.numeric)]
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
        updateSelectInput(session, "xcol", choices = numeric_cols, selected = numeric_cols[1])
      }
    }
  })
  
  # Function to generate 12 distinct colors for mge_params
  get_mge_colors <- function(mge_params_values) {
    # Create a palette of 12 distinct colors
    color_palette <- c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
      "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78"
    )
    
    # Create named vector mapping each unique mge_params to a color
    unique_params <- unique(mge_params_values)
    colors <- setNames(color_palette[1:length(unique_params)], unique_params)
    return(colors)
  }
  
  # Render summary stats scatter plot with mge_params color coding
  output$summaryStatsPlot <- renderPlotly({
    req(input$xcol, input$ycol, values$dataset_loaded, nrow(values$summary_stats) > 0)
    
    # Check if required columns exist
    if (!"mge_params" %in% names(values$summary_stats)) {
      # Fallback to original plot if mge_params doesn't exist
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
    } else {
      # Get color mapping for mge_params
      color_mapping <- get_mge_colors(values$summary_stats$mge_params)
      
      # Check if Filename column exists for hover
      if ("Filename" %in% names(values$summary_stats)) {
        # Create custom hover template with Filename and mge_params
        hover_template <- paste(
          "Filename: %{text}<br>",
          "MGE Params: %{customdata}<br>",
          input$xcol, ": %{x}<br>",
          input$ycol, ": %{y}",
          "<extra></extra>"
        )
        
        plot_ly(
          data = values$summary_stats,
          x = ~.data[[input$xcol]],
          y = ~.data[[input$ycol]],
          color = ~mge_params,
          colors = color_mapping,
          text = ~Filename,
          customdata = ~mge_params,
          type = 'scatter',
          mode = 'markers',
          marker = list(size = 10, line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
          hovertemplate = hover_template
        ) %>%
          layout(
            xaxis = list(title = input$xcol),
            yaxis = list(title = input$ycol),
            legend = list(
              orientation = "h",
              x = 0,
              y = -0.2,
              xanchor = 'left',
              yanchor = 'top'
            ),
            margin = list(t = 30, b = 120)
          )
      } else {
        # Create hover template without Filename
        hover_template <- paste(
          "MGE Params: %{customdata}<br>",
          input$xcol, ": %{x}<br>",
          input$ycol, ": %{y}",
          "<extra></extra>"
        )
        
        plot_ly(
          data = values$summary_stats,
          x = ~.data[[input$xcol]],
          y = ~.data[[input$ycol]],
          color = ~mge_params,
          colors = color_mapping,
          customdata = ~mge_params,
          type = 'scatter',
          mode = 'markers',
          marker = list(size = 10, line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
          hovertemplate = hover_template
        ) %>%
          layout(
            xaxis = list(title = input$xcol),
            yaxis = list(title = input$ycol),
            legend = list(
              orientation = "h",
              x = 0,
              y = -0.2,
              xanchor = 'left',
              yanchor = 'top'
            ),
            margin = list(t = 30, b = 120)
          )
      }
    }
  })
  
  # Render summary stats bar plot with Process ID search functionality
  output$summaryBarPlot <- renderPlotly({
    req(values$dataset_loaded, nrow(values$summary_stats) > 0)
    
    # Check if required columns exist
    if (!"process_id" %in% names(values$summary_stats) || !"n_reads_in" %in% names(values$summary_stats) || !"mge_params" %in% names(values$summary_stats)) {
      return(NULL)
    }
    
    # Prepare data with mode classification and get unique combinations
    plot_data <- values$summary_stats %>%
      mutate(mode = ifelse(grepl("_merge$", mge_params), "Merge", "Concat")) %>%
      select(process_id, n_reads_in, mode) %>%
      distinct() %>%
      filter(!is.na(n_reads_in))
    
    # Apply Process ID search filter
    search_term <- input$process_id_search %||% ""
    if (search_term != "") {
      plot_data <- plot_data %>%
        filter(grepl(search_term, process_id, ignore.case = TRUE))
    }
    
    # If no data after filtering, return empty plot
    if (nrow(plot_data) == 0) {
      return(
        plot_ly() %>%
          layout(
            title = list(text = "No Process IDs match your search criteria", x = 0.5),
            xaxis = list(title = "Process ID"),
            yaxis = list(title = "Number of Reads In"),
            margin = list(t = 50, b = 50)
          )
      )
    }
    
    # Sort and arrange data
    plot_data <- plot_data %>% arrange(process_id, mode)
    
    # Get unique process_ids for ordering
    unique_process_ids <- unique(plot_data$process_id)
    
    # Create separate data for each mode
    merge_data <- plot_data %>% filter(mode == "Merge")
    concat_data <- plot_data %>% filter(mode == "Concat")
    
    # Create the plot
    p <- plot_ly()
    
    # Add merge mode bars
    if (nrow(merge_data) > 0) {
      p <- p %>% add_bars(
        data = merge_data,
        x = ~factor(process_id, levels = unique_process_ids),
        y = ~n_reads_in,
        name = "Merge Mode",
        marker = list(color = 'rgba(51, 122, 183, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
        hovertemplate = ~paste("Process ID:", process_id, "<br>Mode: Merge<br>Reads In:", n_reads_in, "<extra></extra>")
      )
    }
    
    # Add concat mode bars
    if (nrow(concat_data) > 0) {
      p <- p %>% add_bars(
        data = concat_data,
        x = ~factor(process_id, levels = unique_process_ids),
        y = ~n_reads_in,
        name = "Concat Mode",
        marker = list(color = 'rgba(255, 99, 71, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
        hovertemplate = ~paste("Process ID:", process_id, "<br>Mode: Concat<br>Reads In:", n_reads_in, "<extra></extra>")
      )
    }
    
    # Layout with conditional x-axis labels
    show_labels <- length(unique_process_ids) <= 50  # Only show labels if 50 or fewer Process IDs
    
    p %>% layout(
      xaxis = list(title = "Process ID", showticklabels = show_labels),
      yaxis = list(title = "Number of Reads In"),
      barmode = 'group',
      margin = list(t = 30, b = 50),
      legend = list(x = 0.02, y = 0.98)
    )
  })
  
  # Render FASTA Compare table
  output$fastaCompareTable <- renderDT({
    req(values$dataset_loaded, nrow(values$fasta_compare) > 0)
    
    datatable(
      values$fasta_compare,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-right', targets = which(sapply(values$fasta_compare, is.numeric)) - 1)
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    )
  })
  
  # Update choices for FASTA compare column selectors
  observe({
    if (values$dataset_loaded && nrow(values$fasta_compare) > 0) {
      numeric_cols <- names(values$fasta_compare)[sapply(values$fasta_compare, is.numeric)]
      if (length(numeric_cols) > 0) {
        updateSelectInput(session, "fasta_xcol", choices = numeric_cols, selected = numeric_cols[1])
        updateSelectInput(session, "fasta_ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
      }
    }
  })
  
  # Render FASTA compare plot with seq_id hover
  output$fastaComparePlot <- renderPlotly({
    req(input$fasta_xcol, input$fasta_ycol, values$dataset_loaded, nrow(values$fasta_compare) > 0)
    
    # Check if seq_id column exists
    if ("seq_id" %in% names(values$fasta_compare)) {
      # Create custom hover template with seq_id
      hover_template <- paste(
        "Seq ID: %{text}<br>",
        input$fasta_xcol, ": %{x}<br>",
        input$fasta_ycol, ": %{y}",
        "<extra></extra>"
      )
      
      plot_ly(
        data = values$fasta_compare,
        x = ~.data[[input$fasta_xcol]],
        y = ~.data[[input$fasta_ycol]],
        text = ~seq_id,  # This provides the seq_id for the hover
        type = 'scatter',
        mode = 'markers',
        marker = list(size = 10, color = 'rgba(0, 180, 135, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
        hovertemplate = hover_template
      ) %>%
        layout(
          xaxis = list(title = input$fasta_xcol),
          yaxis = list(title = input$fasta_ycol),
          margin = list(t = 30)
        )
    } else {
      # Fallback to original plot if seq_id doesn't exist
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
    }
  })
  
  # Render JSON Data table
  output$jsonDataTable <- renderDT({
    req(values$dataset_loaded, nrow(values$json_data) > 0)
    
    datatable(
      values$json_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        searchHighlight = TRUE,
        dom = 'Bfrtip',
        buttons = list('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-right', targets = which(sapply(values$json_data, is.numeric)) - 1)
        )
      ),
      filter = 'top',
      rownames = FALSE,
      class = 'display compact',
      extensions = 'Buttons'
    ) %>%
      formatRound(columns = which(sapply(values$json_data, is.numeric)), digits = 4)
  })
  
  # Update choices for JSON column selectors (including bar chart)
  observe({
    if (values$dataset_loaded && nrow(values$json_data) > 0) {
      all_cols <- names(values$json_data)
      numeric_cols <- names(values$json_data)[sapply(values$json_data, is.numeric)]
      
      if (length(numeric_cols) > 0) {
        # For scatter plot
        updateSelectInput(session, "json_xcol", choices = numeric_cols, selected = numeric_cols[1])
        updateSelectInput(session, "json_ycol", choices = numeric_cols, selected = numeric_cols[min(2, length(numeric_cols))])
        
        # For bar chart - all columns except process_id
        bar_choices <- all_cols[all_cols != "process_id"]
        updateSelectInput(session, "json_bar_ycol", choices = bar_choices, selected = bar_choices[1])
      }
    }
  })
  
  # Render JSON bar plot (process_id on x-axis, selectable y-axis)
  output$jsonBarPlot <- renderPlotly({
    req(input$json_bar_ycol, values$dataset_loaded, nrow(values$json_data) > 0)
    
    # Check if required columns exist
    if (!"process_id" %in% names(values$json_data) || !input$json_bar_ycol %in% names(values$json_data)) {
      return(NULL)
    }
    
    # Get unique process_ids and their values, sort by process_id A-Z
    plot_data <- values$json_data %>%
      select(process_id, !!sym(input$json_bar_ycol)) %>%
      distinct() %>%
      arrange(process_id) %>%
      filter(!is.na(!!sym(input$json_bar_ycol)))
    
    # Get unique process_ids for proper factor levels
    unique_process_ids <- unique(plot_data$process_id)
    
    plot_ly(
      data = plot_data,
      x = ~factor(process_id, levels = unique_process_ids),
      y = ~.data[[input$json_bar_ycol]],
      type = 'bar',
      marker = list(color = 'rgba(255, 99, 71, 0.7)', line = list(width = 1, color = 'rgba(0,0,0,0.5)')),
      hovertemplate = ~paste("Process ID:", process_id, "<br>", input$json_bar_ycol, ":", .data[[input$json_bar_ycol]], "<extra></extra>")
    ) %>%
      layout(
        xaxis = list(title = "Process ID", showticklabels = FALSE),
        yaxis = list(title = input$json_bar_ycol),
        margin = list(t = 30, b = 50)
      )
  })
  
  # Render JSON scatter plot
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
  
  # Download handlers
  output$downloadSummaryStats <- downloadHandler(
    filename = function() {
      paste0("summary_stats_", values$current_dataset, "_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      filtered_data <- values$summary_stats
      if (!is.null(input$summaryStatsTable_rows_all)) {
        filtered_data <- filtered_data[input$summaryStatsTable_rows_all, ]
      }
      write.csv(filtered_data, file, row.names = FALSE)
    }
  )
  
  output$downloadFastaCompare <- downloadHandler(
    filename = function() {
      paste0("fasta_compare_", values$current_dataset, "_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      filtered_data <- values$fasta_compare
      if (!is.null(input$fastaCompareTable_rows_all)) {
        filtered_data <- filtered_data[input$fastaCompareTable_rows_all, ]
      }
      write.csv(filtered_data, file, row.names = FALSE)
    }
  )
  
  output$downloadJsonData <- downloadHandler(
    filename = function() {
      paste0("json_data_", values$current_dataset, "_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      filtered_data <- values$json_data
      if (!is.null(input$jsonDataTable_rows_all)) {
        filtered_data <- filtered_data[input$jsonDataTable_rows_all, ]
      }
      write.csv(filtered_data, file, row.names = FALSE)
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
