# ==============================================================================
# Publication Heatmap: Expanded Y-Axis with Sample ID Labels
# ==============================================================================
library(ggplot2)
library(dplyr)
library(viridis)

args = commandArgs(trailingOnly=TRUE)

# 1. Load Dataset
file_path <- args[1]
data <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)

# 2. Assign Abundance Column
if ("Mean Abundance" %in% colnames(data)) {
  data$Abundance <- as.numeric(data$`Mean Abundance`)
} else {
  data$Abundance <- as.numeric(data$`Median Abundance`)
}

# 3. Clean and Filter Data by Sample ID
plot_data <- data %>%
  filter(!is.na(`Query Sequence`) & !is.na(`Sample ID`)) %>%
  group_by(`Query Sequence`, `Sample ID`) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# 4. Group and Sort the Factors Sequentially
plot_data <- plot_data %>%
  arrange(`Sample ID`) %>%
  mutate(`Sample ID` = factor(`Sample ID`, levels = unique(`Sample ID`)))

# 5. Build Heatmap Plot
p <- ggplot(plot_data, aes(x = `Query Sequence`, y = `Sample ID`, fill = Abundance)) +
  geom_tile(color = "white", linewidth = 0.01, width=0.9) + # Fine grid lines to separate cells
  scale_fill_viridis_c(option = "magma", name = "Mean\nAbundance") +
  labs(
    x = "Query Sequence (Exon)",
    y = "Sample ID",
    title = "Exon Abundance Heatmap by Sample ID"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(size = 6, color = "black", family = "mono"), # Mono font keeps tags aligned
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

# Create output PDF name from input CSV name
output_pdf <- sub("\\.csv$", ".pdf", basename(file_path))

ggsave(output_pdf,
       plot = p,
       width = as.numeric(args[2]),
       height = as.numeric(args[3]),
       limitsize = FALSE)

# 6. Save with Explicitly Increased Height to Spread Out Rows
# Adjust height dynamically based on your row count if needed (e.g., height = 18 or 24)
#ggsave("heatmap.pdf", plot = p, width = 10, height = 30, limitsize = FALSE)
#ggsave("expanded_exon_heatmap.png", plot = p, width = 10, height = 20, dpi = 300, limitsize = FALSE)

cat(paste0("Successfully generated plot. Saved as '", output_pdf, "'.\n"))
