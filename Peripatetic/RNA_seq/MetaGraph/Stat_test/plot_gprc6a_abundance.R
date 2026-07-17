#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
})

# ==================== 1. DATA LOADING AND MERGING ====================
file_bos   <- "Bos_taurus_all_exon_abundance.csv"
file_homo  <- "Homo_sapiens_all_exon_abundance.csv"
file_sus   <- "Sus_scrofa_all_exon_abundance.csv"

# Verify all mandatory input files exist before launching pipeline
missing_files <- c(file_bos, file_homo, file_sus)[!file.exists(c(file_bos, file_homo, file_sus))]
if (length(missing_files) > 0) {
  stop(paste("Error: Missing raw csv data files ->", paste(missing_files, collapse = ", ")))
}

cat("Loading raw SRA observation files...\n")
df_bos   <- read_csv(file_bos, show_col_types = FALSE)
df_homo  <- read_csv(file_homo, show_col_types = FALSE)
df_sus   <- read_csv(file_sus, show_col_types = FALSE)

# Combine the individual species files into a unified master observational dataframe
df_raw <- bind_rows(df_bos, df_homo, df_sus) %>%
  rename(Exon = `Query Sequence`, Abundance = `Median Abundance`) %>%
  mutate(
    Exon = factor(Exon, levels = c("exon_1", "exon_2", "exon_3", "exon_4", "exon_5", "exon_6")),
    Species = factor(Species, levels = c("Bos_taurus", "Sus_scrofa", "Homo_sapiens"))
  )

# ==================== 2. EXACT PAIRWISE WELCH'S T-TEST PIPELINE ====================
cat("\n===============================================================\n")
cat("   EXACT PAIRWISE WELCH'S T-TESTS (FDR ADJUSTED) PER EXON\n")
cat("===============================================================\n")

exons_list <- levels(df_raw$Exon)
species_list <- levels(df_raw$Species)

raw_tests_df <- data.frame()

# Loop through each exon to extract raw groups and apply pairwise testing
for (ex in exons_list) {
  df_ex <- df_raw %>% filter(Exon == ex)
  
  # Generate all unique species pairwise combinations
  pairs <- combn(species_list, 2, simplify = FALSE)
  
  for (pair in pairs) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    vec1 <- df_ex %>% filter(Species == s1) %>% pull(Abundance)
    vec2 <- df_ex %>% filter(Species == s2) %>% pull(Abundance)
    
    # Run test only if both populations have data points available
    if (length(vec1) > 1 & length(vec2) > 1) {
      t_result <- t.test(vec1, vec2, alternative = "two.sided", paired = FALSE, var.equal = FALSE)
      
      raw_tests_df <- rbind(raw_tests_df, data.frame(
        Exon = ex,
        Species1 = s1,
        Species2 = s2,
        Mean1 = mean(vec1),
        Mean2 = mean(vec2),
        Mean_Diff = mean(vec1) - mean(vec2),
        p_raw = t_result$p.value
      ))
    }
  }
}

# Apply Benjamini-Hochberg (FDR) correction globally across all comparisons
stats_final <- raw_tests_df %>%
  mutate(p_adj = p.adjust(p_raw, method = "BH")) %>%
  mutate(
    sig_flag = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ "ns"
    )
  )

# Output detailed summary text table directly to the console
for (ex in exons_list) {
  cat(paste("\n--- Exon Category:", ex, "---\n"))
  df_ex_print <- stats_final %>% filter(Exon == ex)
  for (i in 1:nrow(df_ex_print)) {
    row <- df_ex_print[i, ]
    cat(sprintf("  %-12s vs %-12s | Mean Diff: %8.2f | Adj p-value: %8.2e (%s)\n", 
                row$Species1, row$Species2, row$Mean_Diff, row$p_adj, row$sig_flag))
  }
}
cat("\n===============================================================\n\n")

# ==================== 3. COMPUTING SUMMARY STATISTICS FOR BOXPLOT LAYER ====================
# Calculate standard box limits from raw numbers for geom_boxplot(stat='identity')
df_summary <- df_raw %>%
  group_by(Exon, Species) %>%
  summarise(
    min = min(Abundance, na.rm = TRUE),
    Q1 = quantile(Abundance, 0.25, na.rm = TRUE),
    median = median(Abundance, na.rm = TRUE),
    Q3 = quantile(Abundance, 0.75, na.rm = TRUE),
    max = max(Abundance, na.rm = TRUE),
    mean = mean(Abundance, na.rm = TRUE),
    N = n(),
    .groups = "drop"
  )

# Establish width properties for the species box grouping layout spacing
dodge_width <- 0.85
dodge_setting <- position_dodge(width = dodge_width)
n_species <- length(species_list)

df_summary <- df_summary %>%
  mutate(
    species_idx = as.numeric(Species),
    exon_idx = as.numeric(Exon),
    box_x = exon_idx + (species_idx - (n_species + 1) / 2) * (dodge_width / n_species)
  )

# Define fixed target guide altitudes based on your categorical ceiling rules
df_summary <- df_summary %>%
  mutate(
    target_y = case_when(
      max == 88               ~ 40,
      max < 50                ~ 100,
      max >= 90 & max <= 200  ~ 50,
      max > 500               ~ 250,
      TRUE                    ~ max * 1.15
    )
  ) %>%
  mutate(
    text_line_y = target_y,
    text_val_y = target_y + (target_y * 0.12)
  )

# ==================== 4. GENERATING PLOT LABELS STACKING DATA ====================
stats_labels_df <- data.frame()

for (ex in exons_list) {
  df_ex_sum <- df_summary %>% filter(Exon == ex)
  df_ex_stat <- stats_final %>% filter(Exon == ex)
  
  # Set baseline height for stat labels directly above the data text line ceiling
  base_stat_y <- max(df_ex_sum$text_val_y) + 65
  step_y <- 85
  
  if (nrow(df_ex_stat) > 0) {
    for (i in 1:nrow(df_ex_stat)) {
      row <- df_ex_stat[i, ]
      
      # Clean species string formatting for plot layout space conservation
      s1_short <- gsub("_.*", "", row$Species1)
      s2_short <- gsub("_.*", "", row$Species2)
      
      lbl_text <- sprintf("%s vs %s\np = %8.1e\n(%s)", s1_short, s2_short, row$p_adj, row$sig_flag)
      current_label_y <- base_stat_y + ((i - 1) * step_y)
      
      stats_labels_df <- rbind(stats_labels_df, data.frame(
        Exon = ex,
        x_coord = as.numeric(factor(ex, levels = exons_list)),
        y_coord = current_label_y,
        label_str = lbl_text
      ))
    }
  }
}

# ==================== 5. GGPLOT DRAWING ENGINE ====================
p <- ggplot(df_summary, aes(x = Exon, fill = Species)) +
  # 1. Base Boxplots (Stat identity from summary layer)
  geom_boxplot(
    aes(ymin = min, lower = Q1, middle = median, upper = Q3, ymax = max),
    stat = "identity", 
    position = dodge_setting,
    color = "black", 
    width = 0.7,
    alpha = 0.7
  ) +
  # 2. Solid Mean Indicator Bar inside the box
  geom_errorbar(
    aes(x = box_x, ymin = mean, ymax = mean),
    inherit.aes = FALSE,
    color = "blue",
    linetype = "solid",
    width = 0.18
  ) +
  # 3. Solid Median Indicator Bar inside the box
  geom_errorbar(
    aes(x = box_x, ymin = median, ymax = median),
    inherit.aes = FALSE,
    color = "darkgreen",
    linetype = "solid",
    width = 0.18
  ) +
  
  # ==================== UNIFORM POINTER LINES ====================
  geom_segment(
    aes(x = box_x - 0.04, xend = box_x - 0.04, y = text_line_y, yend = median),
    color = "darkgreen", linetype = "dotted", linewidth = 0.45
  ) +
  geom_segment(
    aes(x = box_x + 0.04, xend = box_x + 0.04, y = text_line_y, yend = mean),
    color = "blue", linetype = "dotted", linewidth = 0.45
  ) +
  
  # ==================== TYPOGRAPHY METRIC DATA LABELS ====================
  geom_text(
    aes(x = box_x, y = max, label = max),
    vjust = -0.5, size = 3.3, color = "red", fontface = "bold"
  ) +
  geom_text(
    aes(x = box_x - 0.05, y = text_val_y, label = median),
    hjust = 1, vjust = 0.5, size = 3.3, color = "darkgreen", fontface = "bold"
  ) +
  geom_text(
    aes(x = box_x + 0.05, y = text_val_y, label = round(mean, 1)),
    hjust = 0, vjust = 0.5, size = 3.3, color = "blue", fontface = "bold"
  ) +
  geom_text(
    aes(x = box_x, y = min, label = min),
    vjust = 1.4, size = 3.3, color = "red", fontface = "bold"
  ) +
  
  # ==================== NON-OVERLAPPING STATISTICAL TIERS ====================
  geom_text(
    data = stats_labels_df,
    aes(x = x_coord, y = y_coord, label = label_str),
    inherit.aes = FALSE,
    vjust = 0,
    size = 2.6,
    lineheight = 0.85,
    color = "purple4",
    fontface = "bold.italic"
  ) +
  
  # ==================== SPACE ISOLATED LOWER N LABELS ====================
  geom_text(
    aes(x = box_x, y = -Inf, label = paste0("N=", N)),
    vjust = -1.5, size = 2.9, color = "black", fontface = "italic"
  ) +
  
  # ==================== CANVAS LAYOUT & STYLING ====================
  labs(
    x = "Exon Category",
    y = "RNA Mean Abundance Expression",
    fill = "Species"
  ) +
  scale_fill_brewer(palette = "Pastel1") + 
  theme_bw(base_size = 14) + 
  theme(
    panel.grid.major = element_line(color = "gray93"),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 12),
    axis.text.x = element_text(face = "bold", size = 12, margin = margin(t = 22)),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(55, 20, 25, 20)
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.45)))

output_name <- "gprc6a_pairwise_expression_plot.pdf"
ggsave(output_name, plot = p, device = cairo_pdf, width = 16, height = 9.5)
cat(paste("Successfully generated publication PDF chart:", output_name, "\n"))
