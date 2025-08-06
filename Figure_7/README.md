# Phylogenetic Analysis of Terrestrial vs Aquatic Species: Gene Status Comparison

This repository contains data, scripts, and output files used to investigate the presence/absence of the **GPRC6A gene** across different **trophic levels** (Aquatic, Carnivore, Herbivore, Omnivore) using **phylogenetic ANOVA**

## 📁 Folder Contents

| File | Description |
|------|-------------|
| `Terrestrial_Aquatic_species.nwk` | Phylogenetic tree in Newick format containing species with trophic and gene presence data. |
| `Terrestrial_Aquatic.tsv` | Metadata table including species, gene status (binary), and trophic level. |
| `collapsed_tree_with_traits.nexus` | Collapsed version of the phylogenetic tree, annotated with traits, ready for visualization in FigTree or similar tools. |
| `collapsed_clade_summary.tsv` | Summary of the clades that were collapsed, including species count, trait values, average branch lengths, and depths. |
| `data_analyses_script.r` | R script that performs the main statistical analysis (phylogenetic ANOVA) and generates plots. |
| `collapseTree.r` | R script that collapses monophyletic clades with uniform trait combinations and writes the annotated tree + summary file. |
| `bar_dot_plot.pdf` | Output plot showing gene presence across trophic levels (created by `data_analyses_script.r`). |

---

## 🧪 Analysis Overview

### 1. **Phylogenetic ANOVA**
Using the `phylANOVA()` function from the `phytools` package:
- Compares gene presence between **Aquatic** and **Terrestrial** species.
- Performs **pairwise comparisons** between **Carnivore**, **Herbivore**, and **Omnivore** groups.
- Simulates 1000 replicates to estimate significance.

### 2. **Plotting**
- A bar + dot plot of gene presence across trophic levels is created using `ggplot2`.
- Points represent individual species; bars show mean gene presence with standard error.

### 3. **Tree Collapsing**
The `collapseTree.r` script:
- Collapses clades with **uniform combinations** of gene status and trophic level.
- Produces a **cleaner tree** for visualization and downstream comparison.
- Outputs an **annotated NEXUS file** with traits and a **summary table**.

---

## ▶️ How to Run

To reproduce the full analysis:

1. **Ensure R is installed**, and install required packages:
   ```r
   install.packages(c("ape", "phytools", "ggplot2", "dplyr", "RColorBrewer"))


