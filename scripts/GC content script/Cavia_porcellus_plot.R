#!/usr/bin/env Rscript
library(ggplot2)

# Set the correct path for your input data files
input_file <- "Cavia_porcellus_gcpcnt.bed"  # Update with actual path if needed
rectangles_file <- "new_Cavia_porcellus_absolute_lengths.txt"  # Update with actual path if needed

# Read the input data
data <- read.table(input_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                   col.names = c("ID", "Start", "End", "Value"))

# Calculate the midpoint
data$Midpoint <- (data$Start + data$End) / 2

# Read the additional file with rectangles info
rectangles <- read.table(rectangles_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
                         col.names = c("ID", "Start", "End", "Value"))

# Filter rows where the ID column contains "exon" but exclude "upstream," "intron," and "downstream"
rectangles_exons <- rectangles[grep("exon", rectangles$ID) & !grepl("upstream|intron|downstream", rectangles$ID), ]

# Remove underscores from exon names in the ID column
rectangles_exons$ID <- gsub("_", " ", rectangles_exons$ID)

# Extract unique start and end positions for labels
exon_positions <- unique(c(rectangles_exons$Start, rectangles_exons$End))

# Filter rows where the ID column contains "putative"
rectangles_putative <- rectangles[grep("putative", rectangles$ID, ignore.case = TRUE), ]

# Plot the data with rectangles and lines/points
p <- ggplot() +
  # Add rectangles for exon regions
  geom_rect(data = rectangles_exons, 
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),
            fill = "lightblue", alpha = 0.5, color = NA) +
  # Add rectangles for putative regions in light red
  geom_rect(data = rectangles_putative, 
            aes(xmin = Start, xmax = End, ymin = -Inf, ymax = Inf),
            fill = "lightcoral", alpha = 0.5, color = NA) +
  # Add text labels for exon names at the top of each rectangle
  geom_text(data = rectangles_exons, 
            aes(x = (Start + End) / 2, y = 0.90, label = ID),
            size = 2.5, vjust = -0.5) +
  geom_line(data = data, aes(x = Midpoint, y = Value)) + 
  geom_point(data = data, aes(x = Midpoint, y = Value), size = 1.25) +
  scale_x_continuous(limits = c(-500, 29000),  
                     breaks = exon_positions) +

  scale_y_continuous(breaks = seq(0, 0.90, by = 0.05), limits = c(0, 0.90)) +
  labs(x = "Position (bp)", y = "GC content") +
  ggtitle(expression(italic("Cavia porcellus"))) + 
  theme_minimal() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5, hjust = 1),  
    axis.text.y = element_text(size = 12),
    panel.grid.major = element_line(color = "gray80"),
    panel.grid.minor = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.title = element_text(hjust = 0, size = 14, face = "bold") 
  )
# Save the plot as a PDF
ggsave("Cavia_porcellus_gccontent_plot.pdf", plot = p, width = 10, height = 6)
