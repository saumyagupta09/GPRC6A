# ==============================================================================
# Publication Heatmap: Expanded Y-Axis with Combined Project/Sample Labels
# ==============================================================================
library(ggplot2)
library(dplyr)
library(viridis)

# 1. Load Dataset
file_path <- "metagraph_results_d6921258-15bd-4a4c-8c8f-ba8284de6a60_all.csv"
data <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)

# 2. Assign Abundance Column
if ("Mean Abundance" %in% colnames(data)) {
  data$Abundance <- as.numeric(data$`Mean Abundance`)
} else {
  data$Abundance <- as.numeric(data$`Median Abundance`)
}

# 3. Clean and Combine Project Title (Tissue) with Sample ID
plot_data <- data %>%
  filter(!is.na(`Query Sequence`) & !is.na(`Sample ID`) & !is.na(`Project Title`)) %>%
  # Create a combined label for the Y-axis
  mutate(Sample_Label = paste0("[", `Project Title`, "] ", `Sample ID`)) %>%
  group_by(`Query Sequence`, Sample_Label, `Project Title`) %>%
  summarize(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# 4. Group and Sort the Factors Sequentially
plot_data <- plot_data %>%
  arrange(`Project Title`, Sample_Label) %>%
  mutate(Sample_Label = factor(Sample_Label, levels = unique(Sample_Label)))

# 5. Build Heatmap Plot
p <- ggplot(plot_data, aes(x = `Query Sequence`, y = Sample_Label, fill = Abundance)) +
  geom_tile(color = "white", linewidth = 0.05) + # Fine grid lines to separate cells
  scale_fill_viridis_c(option = "magma", name = "Mean\nAbundance") +
  labs(
    x = "Query Sequence (Exon)",
    y = "Tissue Context & Sample ID",
    title = "Exon Abundance Heatmap by Project/Tissue"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", color = "black"),
    axis.text.y = element_text(size = 6, color = "black", family = "mono"), # Mono font keeps tags aligned
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13)
  )

# 6. Save with Explicitly Increased Height to Spread Out Rows
# Adjust height dynamically based on your row count if needed (e.g., height = 18 or 24)
ggsave("expanded_exon_heatmap.pdf", plot = p, width = 10, height = 20, limitsize = FALSE)
ggsave("expanded_exon_heatmap.png", plot = p, width = 10, height = 20, dpi = 300, limitsize = FALSE)

cat("Successfully generated plot. Saved as 'expanded_exon_heatmap.pdf' (Height: 20 inches).\n")
