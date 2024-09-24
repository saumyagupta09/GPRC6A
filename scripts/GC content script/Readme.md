# GC Content Analysis for Caviomorpha Species

This repository contains scripts for analyzing GC content in genomic sequences of Caviomorpha species. The analysis is performed in 100 bp windows, facilitating a detailed examination of GC content variation across the genome.

## Scripts

1. **`gc_content_windows.sh`**  
   This script generates GC content data in 100 bp windows from the provided genomic sequence. 

2. **`plot_gc_content.R`**  
   This R script visualizes the GC content data generated from the previous script.

## Usage Instructions

### Prerequisites

- Ensure you have the following software installed:
  - **Bash** for running the shell script
  - **R** and the **ggplot2** package for executing the R script

### Running the Scripts

1. **Execute the `gc_content_windows.sh` script**  
   To calculate GC content in 100 bp windows, run the following command in your terminal:

   ```bash
   ./gc_content_windows.sh <species_name> <output_file> <genome_file>
