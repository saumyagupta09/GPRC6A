#!/usr/bin/env Rscript

# GC Content Analysis: Boxplot Generation for Mammalian Orders
# This script generates publication-ready boxplots showing GC content distribution
# across mammalian orders with optional species labeling for high GC content species.

# Load required libraries
library(ggplot2)
library(dplyr)

# Read the data files
cat("Reading data files...\n")

# Read mean GC content data
gc_data <- read.table("mean_gccontent.txt", 
                     header = FALSE, 
                     sep = "\t", 
                     stringsAsFactors = FALSE, 
                     col.names = c("Species_ID", "Mean_GC_Content"))

# Read species order data
order_data <- read.table("list_sps_mammals_with_order_sorted.txt", 
                        header = FALSE, 
                        sep = "\t", 
                        stringsAsFactors = FALSE, 
                        col.names = c("Species", "Order"))

cat("Data loaded successfully!\n")
cat("GC data dimensions:", dim(gc_data), "\n")
cat("Order data dimensions:", dim(order_data), "\n")

# Clean species names in GC data by removing accession numbers
gc_data$Species_Clean <- gsub("_[XN]M_.*$", "", gc_data$Species_ID)

# Merge the datasets
merged_data <- merge(gc_data, order_data, by.x = "Species_Clean", by.y = "Species", all.x = TRUE)

# Check for any unmatched species
unmatched <- merged_data[is.na(merged_data$Order), ]
if(nrow(unmatched) > 0) {
  cat("Warning: Some species could not be matched to orders:\n")
  print(unmatched$Species_Clean)
}

# Remove unmatched species
merged_data <- merged_data[!is.na(merged_data$Order), ]

cat("Final dataset dimensions:", dim(merged_data), "\n")
cat("Number of orders:", length(unique(merged_data$Order)), "\n")

# Define a color palette for orders
order_colors <- c(
  "Artiodactyla" = "#1f77b4",   # Bright blue
  "Rodentia" = "#ff7f0e",       # Bright orange
  "Chiroptera" = "#2ca02c",     # Bright green
  "Lagomorpha" = "#d62728",     # Bright red
  "Cetacea" = "#9467bd",        # Bright purple
  "Pholidota" = "#8c564b",      # Brown
  "Perissodactyla" = "#e377c2", # Pink
  "Carnivora" = "#7f7f7f",      # Gray
  "Afrotheria" = "#a6cee3",     # Light blue
  "Dermoptera" = "#b2df8a",     # Light green
  "Eulipotyphla" = "#33a02c",   # Dark green
  "Marsupialia" = "#fb9a99",    # Light red
  "Pilosa" = "#fdbf6f",         # Light orange
  "Primates" = "#ff69b4",       # Hot pink
  "Scandentia" = "#f1c40f",     # Yellow
  "Caviomorpha" = "#6a3d9a",    # Dark purple
  "Cingulata" = "#ff7f0e",      # Orange
  "Eulipotphyla" = "#33a02c"    # Dark green (duplicate for consistency)
)

# Order the orders for better visualization (by median GC content)
order_medians <- merged_data %>%
  group_by(Order) %>%
  summarise(median_gc = median(Mean_GC_Content)) %>%
  arrange(median_gc)

merged_data$Order <- factor(merged_data$Order, levels = order_medians$Order)

# Create color mapping for x-axis labels
x_axis_colors <- order_colors[levels(merged_data$Order)]

# Create custom y-axis breaks for better readability (decimal format)
y_breaks <- seq(0.35, 0.52, by = 0.01)  # Every 0.01 from 0.35 to 0.52
y_labels <- format(y_breaks, nsmall = 2)  # Format as decimal with 2 decimal places

# Identify species with >0.47 mean GC content for labeling
high_gc_species <- merged_data[merged_data$Mean_GC_Content > 0.47, ]
cat("Species with >0.47 mean GC content:", nrow(high_gc_species), "\n")
if(nrow(high_gc_species) > 0) {
  cat("High GC species:", paste(high_gc_species$Species_Clean, collapse = ", "), "\n")
}

# Function to create the base plot
create_base_plot <- function(include_labels = TRUE) {
  p <- ggplot(merged_data, aes(x = Order, y = Mean_GC_Content, fill = Order)) +
    # Add boxplot
    geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    # Add individual data points
    geom_jitter(aes(color = Order), 
                width = 0.2, 
                height = 0, 
                size = 1.5, 
                alpha = 0.8) +
    # Add reference lines for better interpretation (subtle)
    geom_hline(yintercept = 0.4, linetype = "dotted", color = "gray70", alpha = 0.6) +
    geom_hline(yintercept = 0.45, linetype = "dotted", color = "gray70", alpha = 0.6) +
    geom_hline(yintercept = 0.5, linetype = "dotted", color = "gray70", alpha = 0.6) +
    # Set colors
    scale_fill_manual(values = order_colors, name = "Order") +
    scale_color_manual(values = order_colors, name = "Order") +
    # Custom y-axis with decimal format
    scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(0.35, 0.52),
      expand = c(0.01, 0.01)  # Small padding
    ) +
    # Labels and title
    labs(
      x = "Mammalian Order",
      y = "Mean GC Content",
      title = "Distribution of Mean GC Content"
    ) +
    # Theme and formatting
    theme_minimal(base_size = 12) +
    theme(
      # Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Grid
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.minor.y = element_line(color = "gray92", linewidth = 0.2),
      
      # Axes - Color code x-axis labels to match point colors
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, 
                                color = x_axis_colors),
      axis.text.y = element_text(size = 9, color = "gray30"),
      axis.title = element_text(size = 12, face = "bold"),
      
      # Y-axis ticks
      axis.ticks.y = element_line(color = "gray60", linewidth = 0.5),
      axis.ticks.length.y = unit(0.2, "cm"),
      
      # Title
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      
      # Legend - Position on the right with square symbols
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key = element_rect(color = NA, fill = "white"),
      legend.key.size = unit(0.4, "cm"),
      legend.box.background = element_rect(color = "gray90", fill = "white"),
      legend.box.margin = margin(5, 5, 5, 5),
      
      # Margins
      plot.margin = margin(20, 20, ifelse(include_labels, 40, 20), 20)
    ) +
    # Override legend symbols to be squares
    guides(fill = guide_legend(override.aes = list(shape = 15, size = 3)),
           color = guide_legend(override.aes = list(shape = 15, size = 3)))
  
  # Add species labels if requested
  if(include_labels && nrow(high_gc_species) > 0) {
    p <- p + geom_text(data = high_gc_species,
                       aes(x = Order, y = Mean_GC_Content, label = Species_Clean),
                       vjust = -1.2,
                       hjust = 0.5,
                       size = 3,
                       angle = 0,  # Horizontal text
                       color = "black",
                       fontface = "plain")  # Not italic
  }
  
  return(p)
}

# Create the plot with species names
cat("Creating GC content boxplot with species labels...\n")
p <- create_base_plot(include_labels = TRUE)

# Save the plot
cat("Saving plot...\n")

ggsave("Mammalian_GC_Content_Distribution.pdf", 
       plot = p, 
       width = 16, 
       height = 10, 
       units = "in", 
       device = "pdf",
       dpi = 300)

ggsave("Mammalian_GC_Content_Distribution.png", 
       plot = p, 
       width = 16, 
       height = 10, 
       units = "in", 
       device = "png",
       dpi = 300)

ggsave("Mammalian_GC_Content_Distribution.tiff", 
       plot = p, 
       width = 16, 
       height = 10, 
       units = "in", 
       device = "tiff",
       dpi = 300)

cat("Plot saved successfully!\n")
cat("Files created:\n")
cat("- Mammalian_GC_Content_Distribution.pdf/png/tiff\n")

# Print summary statistics
cat("\nSummary statistics by order:\n")
summary_stats <- merged_data %>%
  group_by(Order) %>%
  summarise(
    n_species = n(),
    mean_gc = round(mean(Mean_GC_Content), 4),
    median_gc = round(median(Mean_GC_Content), 4),
    sd_gc = round(sd(Mean_GC_Content), 4),
    min_gc = round(min(Mean_GC_Content), 4),
    max_gc = round(max(Mean_GC_Content), 4)
  ) %>%
  arrange(median_gc)

print(summary_stats)

cat("\nFeatures:\n")
cat("1. ✓ Species names for those with >0.47 mean GC content (horizontal text)\n")
cat("2. ✓ Simple title: 'Distribution of Mean GC Content'\n")
cat("3. ✓ Square symbols in legend\n")
cat("4. ✓ Order by median GC content\n")
cat("5. ✓ Color-coded x-axis labels and data points\n")

if(nrow(high_gc_species) > 0) {
  cat("\nSpecies with >0.47 mean GC content (labeled on plot):\n")
  for(i in 1:nrow(high_gc_species)) {
    cat(sprintf("- %s: %.4f (%s)\n", 
                high_gc_species$Species_Clean[i], 
                high_gc_species$Mean_GC_Content[i],
                high_gc_species$Order[i]))
  }
}

cat("\nScript completed successfully!\n")
