# Barcoding Results Dashboard
A Shiny-based interactive dashboard for visualisation, filtering, and downloading statistics from the [Barcode Gene Extractor & Evaluator (BGEE)](https://github.com/bge-barcoding/BGEE) workflow.

## Features

### 🔍 **Dual Search & Loading Options**
- **Plate Selection**: Select and browse results for Process IDs in specific plates (BGE_XXXXX format).
- **Project Selection**: Filter results by project codes to view all Process IDs containing specific project identifiers.
- **Progress Tracking**: Real-time loading progress and status updates.
- **Data Validation**: Comprehensive validation and error reporting for missing or corrupted data.

### 📊 **Barcoding Outcome Overview**
- **Three-Tier Classification**: Automated categorisation of barcoding results with traffic light system:
  - 🟢 **Green - Pass**: Both structural and taxonomic validation successful.
  - 🟡 **Amber - Partial**: Either structural OR taxonomic validation successful (still considered overall fail).
  - 🔴 **Red - Fail**: Neither structural nor taxonomic validation successful.
- **Validation Criteria**:
  - **Structural**:
	1. Filter to sequences with ambig_original == 0 (no ambiguous bases in original sequence)
	2. From those, pick ones with stop_codons == 0 (no stop codons)
	3. From those, pick ones with reading_frame in [1,2,3] (valid reading frames)
	4. From those, pick ones with min(ambig_basecount) (fewest Ns in HMM-extracted barcode region)
	5. From those, pick ones with max(nuc_basecount) (longest sequence)
	6. If multiple sequences remain, pick the one with highest cov_med (best median coverage as final tie-breaker)
  - **Taxonomic**: Exact word match between family rank in BOLD and taxonomic databases.
- **BOLD Systems Integration**: Clickable Process ID links to BOLD database records.
- **Summary Statistics**: Pass rates and outcome distributions with visual summaries.

### 🧬 **FASTA Sequences**
- **Barcode Sequences**: View and download consensus barcode sequences in FASTA format.
- **Smart Filtering**: Filter sequences by validation status (All/Pass/Partial/Fail).
- **Export Options**: Download FASTA files or copy sequences to clipboard.
- **Missing Data Handling**: Graceful handling of sequences with missing nucleotide data.

### 📈 **Interactive Visualisations**
- **Scatter Plots**: Interactive plotly charts with customisable X/Y axes across all data types.
- **Bar Charts**: Process ID overview with proportional display options for fastp metrics.
- **Dynamic Updates**: Real-time plot updates based on selected metrics and filters.

### 📁 **Multi-Source Data Analysis**
- **Fastp Metrics**: Pre-processing quality statistics from JSON files:
  - Read counts before/after filtering
  - Quality scores (Q20/Q30 rates)
  - Duplication rates and filtering results
- **BGEE Summary Statistics**: Aggregated pipeline results from MGE and fasta_cleaner steps.
- **Barcode Validation**: Structural and taxonomic validation results from barcode_validator tool.

### 💾 **Data Export & Management**
- **Filtered Downloads**: Export filtered results as CSV files with dynamic naming.
- **Interactive Tables**: Sortable, searchable data tables with pagination and filtering.
- **Multiple Formats**: Support for CSV downloads and clipboard operations.

## Supported Data
The dashboard can filter results by plate ID, with 451 plates supported:
```
e.g., BGE_00001, BGE_00002, BGE_00003, etc.
```
Or by the following (27) project codes:
```
BSNHM, NHMCG, BGETR, BSUIO, BGLIB, BSNTN, BGEGR, DTAUT, HMAUT, 
BGENL, DTKNU, BBIOP, BHNHM, UNIFI, DTULO, MEAMP, MUSBA, BGSNL, 
BGSNH, BGEPL, EBGEP, BSCRO, BIOSC, INVBG, BCEMI, ILECA, ALPFU
```

## Quick Start

### Prerequisites and Setup
1. [Download R and RStudio](https://posit.co/download/rstudio-desktop/)
2. Install required packages:
```r
install.packages(c("shiny", "DT", "jsonlite", "dplyr", "tidyr", "plotly", "shinyjs"))
```

### Running the Dashboard
1. Download or copy `launch_dash.R` to your local directory.
2. Open `launch_dash.R` in R/RStudio.
3. Click 'Run App' (top right of the text editor). The dashboard will open in your default web browser.
OR
1. Visit [the deployed app](https://schistodan.shinyapps.io/barcoding-dashboard/).


## Usage Workflow
1. **Select Data**: Choose either:
   - **Plate**: Select a specific plate from the dropdown (e.g., BGE_00001)
   - **Project**: Select a project code to view all Process IDs containing that identifier
2. **Load Results**: Click "Load Results" to import and process the data.
3. **View Outcomes**: Navigate to "Barcoding Outcome" to see:
   - Traffic light classification of all Process IDs
   - Summary statistics and pass rates
   - Clickable BOLD Systems links
4. **Explore Data**: Use additional tabs to examine detailed metrics:
   - **Barcodes**: View and download FASTA sequences with filtering options
   - **Fastp Metrics**: Pre-processing quality statistics and visualizations
   - **BGEE Summary Statistics**: Pipeline summary results with interactive plots
   - **Barcode Validation**: Detailed structural and taxonomic validation data
5. **Export Results**: Download filtered data using the export buttons on each tab.

## Technical Details

### Validation Logic
- **Structural Validation**: Requires ALL criteria to pass:
  - Nucleotide count ≥ 500bp
  - Ambiguous base count = 0
  - Stop codon count = 0
  - Reading frame ∈ {1, 2, 3}
- **Taxonomic Validation**: Requires exact word match (case-insensitive) between BOLD family name and any term in the taxonomic ID string (comma-separated)

### File Processing
- **Metadata**: Automatically detects and standardizes column names for Process ID and Sample ID
- **BGEE Data**: Processes CSV files with flexible column naming (ID, Process.ID, etc.)
- **Validation Data**: Handles TSV files with group_id or process_id columns
- **JSON Data**: Extracts Process IDs from filename patterns ([ProcessID]_fastp_report.json)

## Support

For issues or questions, please open an issue on this repository.
