# GPRC6A: Evolutionary Analysis of G-Protein Coupled Receptor Family C Group 6 Member A

## Overview

This repository contains a comprehensive evolutionary analysis of the **GPRC6A gene** across diverse vertebrate taxa. The project investigates gene presence/absence patterns, evolutionary selection pressures, and genomic features across different taxonomic groups and ecological niches.

## Research Focus

The **GPRC6A gene** is a G-protein coupled receptor that plays important roles in:
- Calcium sensing and metabolism
- Bone development and homeostasis
- Energy metabolism
- Reproductive physiology

This project analyzes the evolutionary history and selective pressures acting on this gene across vertebrates, with particular attention to:
- Gene loss events in specific lineages
- Selection pressure variations (dN/dS ratios)
- Genomic context and synteny
- Phylogenetic distribution patterns

## Project Structure

### Core Analysis Directories

#### 1. **Absrel/** - Adaptive Branch-Site Random Effects Likelihood Analysis
- Contains results from HYPHY's aBSREL analysis for detecting episodic diversifying selection
- Organized by taxonomic groups (Primates, Rodents, Carnivores, etc.)
- Each group contains configuration files, results, and logs

#### 2. **CODEML/** - PAML Codon Substitution Analysis
- Results from PAML's codeml program for estimating dN/dS ratios
- Phylogenetic models and likelihood ratio tests
- Organized by taxonomic groups with control files and output data

#### 3. **HYPHY/** - Hypothesis Testing Using Phylogenies
- Results from HYPHY analyses including:
  - Branch-site models
  - Site-specific selection tests
  - Gene family evolution analyses
- Organized by taxonomic groups

#### 4. **DAMBE/** - Data Analysis in Molecular Biology and Evolution
- Codon alignment files in PHYLIP format
- Results from DAMBE analyses for sequence evolution

### Data and Alignment Directories

#### 5. **Alignment_and_CDS/**
- Coding sequence (CDS) alignments for each taxonomic group
- Protein alignments and codon alignments
- Phylogenetic trees in Newick format
- Species lists for each analysis

#### 6. **Alignments/**
- Multiple sequence alignments (MSA) generated with MAFFT
- Named alignments for different taxonomic groups

#### 7. **Annotations/**
- CDS sequences in FASTA format for each taxonomic group
- Longest isoform annotations
- CASR (Calcium Sensing Receptor) sequences for comparative analysis

### Taxonomic Group-Specific Analyses

#### 8. **Individual Taxa Directories**
- **Primates_Scandentia_Dermoptera/**
- **Rodentia_Lagomorpha/**
- **Carnivora_Pholidota/**
- **Chiroptera/**
- **Artiodactyla/**
- **Perissodactyla/**
- **Eulipotyphla/**
- **Afrotheria_Xenarthra/**
- **Cetacea/**
- **Caviomorpha/**
- **Ruminantia/**
- **Sirenia/**

Each contains:
- FASTA sequence files
- Genome BLAST results
- SRA BLAST analyses
- Alignment files
- PacBio data (where available)

### Specialized Analysis Tools

#### 9. **gbgc/** - GC-Biased Gene Conversion Analysis
- R scripts for analyzing GC content vs. selection pressure relationships
- PHAST analysis results for detecting selection
- Visualization scripts for GC* vs dN/dS plots

#### 10. **TOGA/** - Tool to infer Orthologs from Genome Alignments
- Orthology inference results for multiple species
- Gene presence/absence calls across genomes

#### 11. **Repeatmasker/**
- Repeat element analysis results
- Masked sequences for various taxa

#### 12. **Synteny_blasts/**
- Synteny conservation analysis
- BLAST results for syntenic regions across species

### Visualization and Results

#### 13. **Figure_7/**
- Phylogenetic ANOVA analysis scripts
- Terrestrial vs. aquatic species comparisons
- Tree collapsing and visualization tools
- Statistical analysis results

#### 14. **Independent_losses/**
- Analysis of independent gene loss events
- Species-specific loss investigations

#### 15. **No_loss/**
- Analysis of species that retained the GPRC6A gene
- Comparative genomics for retained genes

## Scripts and Tools

### Utility Scripts (`scripts/`)

#### **check_orf.sh**
- Validates open reading frames in FASTA sequences
- Identifies in-frame stop codons
- Reports ORF completeness

#### **exonwise.sh**
- Processes BLAST output to extract exon sequences
- Handles both forward and reverse strands
- Generates BED files for sequence extraction

#### **GC Content Analysis (`scripts/GC content script/`)**
- `gc_content_windows.sh`: Calculates GC content in 100bp windows
- Species-specific plotting scripts (R)
- Generates GC content visualizations

## Key Analyses Performed

### 1. **Selection Pressure Analysis**
- **dN/dS ratios** using PAML (CODEML) and HYPHY
- **Branch-site models** for detecting episodic selection
- **Site-specific selection** tests

### 2. **Gene Loss Detection**
- **TOGA analysis** for orthology inference
- **BLAST-based** gene presence/absence detection
- **Synteny analysis** for gene loss validation

### 3. **Phylogenetic Analysis**
- **Phylogenetic ANOVA** for trait associations
- **Tree collapsing** for visualization
- **Clade-specific** analyses

### 4. **Genomic Context Analysis**
- **GC content** variation analysis
- **Repeat element** identification
- **Synteny conservation** assessment

## Getting Started

### Prerequisites

- **Bioinformatics tools:**
  - BLAST+ suite
  - MAFFT for alignments
  - PAML (CODEML)
  - HYPHY
  - DAMBE
  - BEDTools
  - RepeatMasker

- **Programming languages:**
  - R (with packages: ggplot2, phytools, ape, dplyr)
  - Python
  - Bash

- **Data formats:**
  - FASTA sequences
  - BLAST output files
  - Phylogenetic trees (Newick format)
  - Multiple sequence alignments

### Basic Workflow

1. **Sequence Collection**: Gather CDS sequences for target species
2. **Alignment**: Generate multiple sequence alignments using MAFFT
3. **Phylogenetic Analysis**: Build trees and perform selection tests
4. **Gene Loss Detection**: Use TOGA and BLAST for orthology inference
5. **Statistical Analysis**: Perform phylogenetic ANOVA and other tests
6. **Visualization**: Generate plots and figures

## Results and Outputs

### Main Findings
- **Gene loss patterns** across vertebrate phylogeny
- **Selection pressure variations** between lineages
- **Ecological correlations** with gene presence/absence
- **Synteny conservation** patterns

### Key Output Files
- **Alignment files** (.aln, .fa, .phy)
- **Phylogenetic trees** (.nwk, .nexus)
- **Selection test results** (dN/dS ratios, p-values)
- **Statistical analysis results** (ANOVA tables, plots)
- **Gene presence/absence matrices**

## Research Applications

This project provides insights into:
- **Evolutionary genomics** of G-protein coupled receptors
- **Gene family evolution** and loss mechanisms
- **Adaptive evolution** in vertebrate lineages
- **Comparative genomics** methodologies
- **Phylogenetic comparative methods**

## Citation and References

If you use this project in your research, please cite the relevant publications and acknowledge the computational methods and analyses performed.

## Contributing

This is a research project repository. For questions about methodology or results, please contact the research team.

## License

This project is for research purposes. Please respect the intellectual property and cite appropriately when using the data or methods.

---

**Note**: This repository contains large genomic datasets and analysis results. Some files may be compressed or require specific bioinformatics tools to access and interpret.
