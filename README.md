# Barcoding Results Dashboard
A Shiny-based interactive dashboard for visualisation, filtering, and downloading validated barcode sequences (and associated statistics) from the [Barcode Gene Extractor & Evaluator (BGEE)](https://github.com/bge-barcoding/BGEE) pipeline. Barcode validation is undertaken using the [barcode_validator](https://github.com/naturalis/barcode_validator/tree/main) tool. Both BGEE and barcode_validator were built for the [Biodiversity Genomics Europe (BGE) consortium](https://biodiversitygenomics.eu/).


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


## Features
### Search & loading options
- **Plate selection**: Select and browse results for Process IDs in specific plates (BGE_XXXXX format).
- **Project selection**: Filter results by project codes to view all Process IDs containing specific project identifiers (e.g. BSNHM, BSNTN, BIOSC, etc).
- **Progress tracking**: Data loading progress and status updates.
- **Data validation**: Logging and error reporting for missing data.

### Barcoding outcome overview
- 🟢 **Green - Pass**: Both structural and taxonomic validation successful.
- 🟡 **Amber - Partial**: Either structural OR taxonomic validation successful (still considered overall fail).
- 🔴 **Red - Fail**: Neither structural nor taxonomic validation successful.
- **barcode_validator validation criteria**:
  - **Structural**:
	1. Filter to sequences with ambig_original == 0 (no ambiguous bases in original sequence (evidence of chimeric sequences))
	2. From those, pick ones with stop_codons == 0 (no stop codons)
	3. From those, pick ones with reading_frame in [1,2,3] (valid reading frames)
	4. From those, pick ones with min(ambig_basecount) (fewest Ns in HMM-extracted barcode region (these Ns represent gaps in the barcode consensus sequence))
	5. From those, pick ones with max(nuc_basecount) (longest sequence)
	6. If multiple sequences remain, pick the one with highest cov_med (best median coverage as final tie-breaker)
  - **Taxonomic**: Exact word match between family rank in BOLD (expected taxonomy) and taxonomic databases (observed taxonomy).
- **BOLD systems integration**: Clickable Process ID links to corresponding BOLD database records.
- **Summary statistics**: Pass rates and outcome distributions with visual summaries.

### FASTA Sequences
- **Barcode sequences**: View and download consensus barcode sequences in FASTA format.
- **Export options**: Download FASTA files or copy sequences to clipboard.

### Interactive visualisation(s)
- **Scatter plots**: Interactive plotly charts with customisable X/Y axes across all data types.
- **Bar charts**: Process ID overview with proportional display options for fastp metrics.

### Supported Data
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

### Workflow usage 
1. **Select data**: Choose either:
   - **Plate**: Select a specific plate from the dropdown (e.g., BGE_00001)
   - **Project**: Select a project code to view all Process IDs containing that identifier
2. **Load results**: Click "Load Results" to import and process the data.
3. **View outcomes**: Navigate to "Barcoding Outcome" to see:
   - Traffic light classification of all Process IDs
   - Summary statistics and pass rates
   - Clickable BOLD Systems links
4. **Explore data**: Use additional tabs to examine detailed metrics:
   - **Barcodes**: View and download FASTA sequences with filtering options
   - **BGEE summary statistics**: Pipeline summary results with interactive plots
   - **Barcode validation**: Detailed structural and taxonomic validation data
5. **Export results**: Download filtered data using the export buttons on each tab.

## Support
For issues or questions, please open an issue on this repository.
