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

# Combine into a unified master observational dataframe with abbreviated names
df_raw <- bind_rows(df_bos, df_homo, df_sus) %>%
  rename(Exon = `Query Sequence`, Abundance = `Median Abundance`) %>%
  mutate(
    Species = case_when(
      Species == "Bos_taurus"   ~ "B. taurus",
      Species == "Sus_scrofa"   ~ "S. scrofa",
      Species == "Homo_sapiens" ~ "H. sapiens",
      TRUE                      ~ Species
    ),
    Exon = factor(Exon, levels = c("exon_1", "exon_2", "exon_3", "exon_4", "exon_5", "exon_6")),
    Species = factor(Species, levels = c("B. taurus", "S. scrofa", "H. sapiens"))
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
  pairs <- combn(species_list, 2, simplify = FALSE)
  
  for (pair in pairs) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    vec1 <- df_ex %>% filter(Species == s1) %>% pull(Abundance)
    vec2 <- df_ex %>% filter(Species == s2) %>% pull(Abundance)
    
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

# ==================== 4. GENERATING BRACKET PATHS & LABELS DATA ====================
brackets_df <- data.frame()
bracket_labels_df <- data.frame()

for (ex in exons_list) {
  df_ex_sum  <- df_summary %>% filter(Exon == ex)
  df_ex_stat <- stats_final %>% filter(Exon == ex)
  
  # CLEARANCE PROTECTION: Base height explicitly sits above the maximum value 
  # of BOTH the physical top whiskers (max) and the normal metrics labels (text_val_y)
  highest_element <- max(c(df_ex_sum$max, df_ex_sum$text_val_y), na.rm = TRUE)
  base_bracket_y <- highest_element + 50
  step_y <- 65 
  
  if (nrow(df_ex_stat) > 0) {
    for (i in 1:nrow(df_ex_stat)) {
      row <- df_ex_stat[i, ]
      
      x1 <- df_ex_sum %>% filter(Species == row$Species1) %>% pull(box_x)
      x2 <- df_ex_sum %>% filter(Species == row$Species2) %>% pull(box_x)
      
      y_bar <- base_bracket_y + ((i - 1) * step_y)
      y_tick <- y_bar - 12 
      
      # Determine color mapping based on test status: Significant (Dark Blue) vs Non-Significant (Slate Gray)
      bracket_color <- if (row$sig_flag == "ns") "#6B7280" else "#1E3A8A"
      
      # Build the classic bracket path sequence map
      brackets_df <- rbind(brackets_df, data.frame(
        Exon = ex,
        x = c(x1, x1, x2, x2),
        y = c(y_tick, y_bar, y_bar, y_tick),
        group = paste0(ex, "_", i),
        color_val = bracket_color
      ))
      
      lbl_text <- sprintf("p = %8.1e %s", row$p_adj, row$sig_flag)
      bracket_labels_df <- rbind(bracket_labels_df, data.frame(
        Exon = ex,
        x_mid = (x1 + x2) / 2,
        y_mid = y_bar + 5,
        label_str = lbl_text,
        color_val = bracket_color
      ))
    }
  }
}

# ==================== 5. GGPLOT DRAWING ENGINE ====================
p <- ggplot(df_summary, aes(x = Exon, fill = Species)) +
  # 1. Base Boxplots
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
  
  # ==================== NON-OVERLAPPING DYNAMIC BRACKETS LAYER ====================
  # Using identity scale via aesthetic parsing ensures colors display explicitly without an extra legend
  geom_path(
    data = brackets_df,
    aes(x = x, y = y, group = group, color = color_val),
    inherit.aes = FALSE,
    linewidth = 0.55
  ) +
  geom_text(
    data = bracket_labels_df,
    aes(x = x_mid, y = y_mid, label = label_str, color = color_val),
    inherit.aes = FALSE,
    vjust = 0,
    size = 2.6,
    fontface = "bold.italic"
  ) +
  scale_color_identity() + # Direct hex coloring interpreter injection
  
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
    legend.text = element_text(size = 12, face = "italic"),
    axis.text.x = element_text(face = "bold", size = 12, margin = margin(t = 22)),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(55, 20, 25, 20)
  ) +
  coord_cartesian(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.55))) # Expanded top ceiling to secure high tiers

output_name <- "gprc6a_pairwise_expression_plot.pdf"
ggsave(output_name, plot = p, device = cairo_pdf, width = 16, height = 9.5)
cat(paste("Successfully generated publication PDF chart:", output_name, "\n"))
