# GC Content Analysis Toolkit

This toolkit generates publication-ready boxplots showing the distribution of mean GC content across mammalian orders. The analysis is based on genomic data from 206 mammalian species representing 18 different orders.

## Files Included

- `gc_content_boxplot.R` - Main R script for generating boxplots
- `mean_gccontent.txt` - Mean GC content data for each species
- `list_sps_mammals_with_order_sorted.txt` - Species-to-order mapping data
- `README.md` - This documentation file

## Requirements

### R Packages
- `ggplot2` - For creating plots
- `dplyr` - For data manipulation

Install required packages if not already installed:
```r
install.packages(c("ggplot2", "dplyr"))
```

## Usage

### Basic Usage
```bash
Rscript gc_content_boxplot.R
```

### What the Script Does

1. **Data Loading**: Reads GC content and species order data
2. **Data Cleaning**: Removes accession numbers from species names and merges datasets
3. **Visualization**: Creates a boxplot with species names shown for species with mean GC content > 0.47

### Output Files

The script generates 3 files in total:

- `Mammalian_GC_Content_Distribution.pdf`
- `Mammalian_GC_Content_Distribution.png`
- `Mammalian_GC_Content_Distribution.tiff`

## Plot Features

- **High Resolution**: All outputs are generated at 300 DPI
- **Color Coding**: Each mammalian order has a distinct color for both data points and x-axis labels
- **Reference Lines**: Subtle dotted lines at 0.4, 0.45, and 0.5 GC content for easy interpretation
- **Ordered by Median**: Orders are arranged by median GC content (lowest to highest)
- **Species Labeling**: Species with >0.47 mean GC content are labeled horizontally
- **Square Legend**: Legend uses square symbols for better visual consistency

## Data Summary

- **Total Species**: 206
- **Mammalian Orders**: 18
- **GC Content Range**: 0.355 - 0.510
- **High GC Species**: 6 species with >0.47 mean GC content

### High GC Content Species (>0.47)
- Acomys_russatus (Rodentia): 0.4843
- Capra_hircus (Artiodactyla): 0.4716
- Delphinapterus_leucas (Cetacea): 0.5101
- Eptesicus_fuscus (Chiroptera): 0.4795
- Octodon_degus (Caviomorpha): 0.4779
- Pipistrellus_kuhlii (Chiroptera): 0.4962

## Order Statistics

Orders are ranked by median GC content:

| Order | Species Count | Mean GC | Median GC | Min GC | Max GC |
|-------|---------------|---------|-----------|--------|--------|
| Primates | 31 | 0.412 | 0.408 | 0.405 | 0.450 |
| Eulipotyphla | 5 | 0.400 | 0.409 | 0.355 | 0.418 |
| Marsupialia | 8 | 0.417 | 0.417 | 0.405 | 0.433 |
| Perissodactyla | 5 | 0.426 | 0.425 | 0.422 | 0.438 |
| Cingulata | 1 | 0.425 | 0.425 | 0.425 | 0.425 |
| Scandentia | 1 | 0.426 | 0.426 | 0.426 | 0.426 |
| Afrotheria | 5 | 0.430 | 0.427 | 0.412 | 0.452 |
| Dermoptera | 2 | 0.428 | 0.428 | 0.417 | 0.438 |
| Eulipotphyla | 1 | 0.428 | 0.428 | 0.428 | 0.428 |
| Pholidota | 2 | 0.429 | 0.429 | 0.429 | 0.429 |
| Artiodactyla | 21 | 0.435 | 0.432 | 0.426 | 0.472 |
| Pilosa | 1 | 0.432 | 0.432 | 0.432 | 0.432 |
| Cetacea | 18 | 0.435 | 0.435 | 0.416 | 0.510 |
| Carnivora | 40 | 0.436 | 0.436 | 0.428 | 0.448 |
| Lagomorpha | 5 | 0.427 | 0.438 | 0.406 | 0.438 |
| Chiroptera | 23 | 0.443 | 0.444 | 0.417 | 0.496 |
| Rodentia | 34 | 0.445 | 0.447 | 0.417 | 0.484 |
| Caviomorpha | 3 | 0.467 | 0.466 | 0.458 | 0.478 |
