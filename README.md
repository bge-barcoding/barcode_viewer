# Barcode Viewer Dashboard - UPDATE README
A Shiny-based interactive dashboard for visualisation, filtering, and downloading statistics from the Barcode Gene Extractor & Evaluator (BGEE) workflow.


## Features
- **Sample Finder**: Search across datasets by plate number, sample ID, or process ID.
- **Dataset Management**: Load and explore available datasets of BGEE outputs.
- **Interactive Visualisations**: Interactive scatter plots and bar charts.
- **Data Export**: Subset, filter, copy and download filtered results as CSV or Excel files.
- **Multi-source Analysis**: 
  - Fastp pre-processing quality metrics (JSON parsing)
  - BGEE summary statistics (MGE + fasta_cleaner results for both pre-processing modes)
  - Fasta Compare results (Post-processing sequence ranking (based on BOLD BIN criteria) and selection)


## Quick Start
### Prerequisites
1. [Download R and RStudio](https://posit.co/download/rstudio-desktop/).
2. Install required packages (in R).
```r
install.packages(c("shiny", "DT", "jsonlite", "dplyr", "tidyr", "plotly", "shinyjs"))
```
### Running the Dashboard
1. Clone this repository
2. Open the launcher in R/RStudio and click 'Run' (top right of the text editor window). Or, open a new R/Rstudio script and run the launcher with (you might need to set your working directory to the project folder for this second method to work):
```r
source("Launch_dash.R")
```
The dashboard will open in your default web browser.


## Support
For issues or questions, please open an issue on this repository.
