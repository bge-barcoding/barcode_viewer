# Barcoding Results Dashboard

A Shiny-based interactive dashboard for visualisation, filtering, and downloading statistics from the [Barcode Gene Extractor & Evaluator (BGEE)](https://github.com/bge-barcoding/BGEE) workflow.

## Features

### 🔍 **Plate-Based Search & Loading**
- **Plate Selection**: Select and browse results for specific plates, or one (or more) Process IDs.
- **Progress Tracking**: Real-time loading progress with detailed status updates.
- **Data Validation**: Comprehensive validation and error reporting for missing or corrupted data.

### 📊 **Barcoding Outcome Overview**
- **Success/Failure Classification**: Automated rank categorisation of barcoding results (see below for rankings).
  - ✅ **Success**: Barcode (rank 1-3) successfully recovered.
  - ⚠️ **Partial**: Rank 4+ (lower quality) recovered.
  - ❌ **Failed**: No barcode recovered.
- **BOLD Systems Integration**: Clickable Process ID links to BOLD database records.
- **Summary Statistics**: Success rates and outcome distributions.
- **Quality Metrics**: Barcode rank, longest stretch (without gaps or ambiguous bases), and sequence coverage across barcode regions calcualted.

### 📈 **Interactive Visualisations**
- **Scatter Plots**: Interactive plotly charts with customisable X/Y axes, including dynamic plot updates based on selected metrics.
- **Bar Charts**: Process ID overview with proportional display options.

### 📁 **Multi-Source Data Analysis**
- **Fastp Metrics**: Combined pre-processing quality statistics from JSON files.
  - Read counts before/after filtering.
  - Quality scores (Q20/Q30 rates).
  - Duplication rates and filtering results.
- **BGEE Summary Statistics**: Aggregated BGEE pipeline results.
- **FASTA Compare Results**: barcode consensus sequence ranking and selection based on BOLD BIN criteria.

### 💾 **Data Export & Management**
- **Filtered Downloads**: Export filtered results as CSV files or to the clipboard.
- **Interactive Tables**: Sortable, searchable data tables with pagination.


## Data Structure Requirements
The dashboard expects data in the following structure:
```
./data/
├── sample_metadata.csv          # Sample metadata with Process ID mapping
├── bgee_summary_stats/          # BGEE summary CSV files
├── fasta_compare_stats/         # FASTA compare CSV files
└── fastp_json/                  # Fastp JSON output files
```

### Sample Metadata Format
The metadata CSV should contain at minimum:
- `sample_id` (or variations): Sample identifiers in BGE_XXXXX_XXX format.
- `process_id` (or variations): Corresponding Process IDs for BOLD Systems.


## Quick Start
### Prerequisites and setup
1. [Download R and RStudio](https://posit.co/download/rstudio-desktop/)
2. Install required packages:
```r
install.packages(c("shiny", "DT", "jsonlite", "dplyr", "tidyr", "plotly", "shinyjs"))
```
### Running the Dashboard
1. Download or copy `launch_dash.R` to launch the dashboard app directly from this GitHub repository.
1. Open `launch_dash.R` in R/RStudio.
2. Click 'Run App' (top right of the text editor). The dashboard will open in your default web browser.


## Usage Workflow
1. **Load Plate**: Select a plate from the dropdown in the "Search Results" tab.
2. **View Outcomes**: Navigate to "Barcoding Outcome" to see an overview of barcode recovery success/failure.
3. **Explore Data**: Use other tabs to examine detailed metrics:
   - **Fastp Metrics**: Pre-processing quality statistics.
   - **BGEE Summary Statistics**: Pipeline summary results.  
   - **FASTA Compare Results**: Post-processing barcode analysis and selection.
4. **Export Results**: Download filtered data using the export buttons.


## Support
For issues or questions, please open an issue on this repository or check the [BGEE documentation](https://github.com/bge-barcoding/BGEE) for pipeline-specific questions.
