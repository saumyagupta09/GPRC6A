#!/usr/bin/env Rscript

library(ggplot2)

# Set input data files
input_file <- "Chinchilla_lanigera_gcpcnt.bed"  
rectangles_file <- "new_Chinchilla_lanigera_absolute_lengths.txt" 

# Read the input data
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                   col.names = c("ID", "Start", "End", "Value"))

# Calculate the midpoint
data$Midpoint <- (data$Start + data$End) / 2

# Read the additional rectangles info
rectangles <- read.table(rectangles_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                         col.names = c("ID", "Start", "End", "Value"))

# Filter rows where the ID column contains "exon"
rectangles_exons <- rectangles[grep("exon", rectangles$ID), ]

# Remove underscores from exon names
rectangles_exons$ID <- gsub("_", " ", rectangles_exons$ID)

# Extract unique start and end positions for labels
exon_positions <- unique(c(rectangles_exons$Start, rectangles_exons$End))

# Plot the data with rectangles and lines/points
p <- ggplot() +
  # Add rectangles for exon regions
  geom_rect(data = rectangles_exons, 
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.5, color = NA) +
  # Add text labels for exon names at the top of each rectangle
  geom_text(data = rectangles_exons, 
            aes(x = (Start + End) / 2, y = 0.8, label = ID),
            size = 3, vjust = -0.5) + 
  geom_line(data = data, aes(x = Midpoint, y = Value)) + 
  geom_point(data = data, aes(x = Midpoint, y = Value), size = 1.25) + 
  # Add a vertical line at x = -500
  geom_vline(xintercept = -500, color = "black") +
  # Set the X-axis limits and breaks at exon start and end positions
  scale_x_continuous(limits = c(0, 23672), 
                     breaks = exon_positions) +
  # Set the Y-axis breaks in increments of 0.05, limits set to 0.80
  scale_y_continuous(breaks = seq(0, 0.80, by = 0.05), limits = c(0, 0.80)) +
  labs(x = "Position (bp)", y = "GC content") +  # Updated Y-axis label
  ggtitle(expression(italic("Chinchilla lanigera"))) + # Italicized title
  theme_minimal() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 5, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8), # Borders on all sides
    plot.title = element_text(hjust = 0, size = 14, face = "bold") # Title positioning
  )
ggsave("Chinchilla_lanigera_gccontent_plot.pdf", plot = p, width = 10, height = 6)
