########################## Phylogenetic Analysis 
library(ape)
library(phytools)

tree <- read.tree("Terrestrial_Aquatic_species.nwk")
data <- read.table("Terrestrial_Aquatic.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

rownames(data) <- data$Species

data$TrophicGroup <- ifelse(data$TrophicLevel == "Aquatic", "Aquatic", "Terrestrial")
data$TrophicGroup <- as.factor(data$TrophicGroup)

grouped_trophic <- setNames(data$TrophicGroup, rownames(data))
gene_status <- setNames(data$gene_status, rownames(data))

#conduct Anova 
phylANOVA_Aq_Ter <- phylANOVA(tree, x = grouped_trophic, y = gene_status, nsim = 1000)  # this does overall ANOVA on Terrestrial and Aquatic group
print(phylANOVA_Aq_Ter)


run_pairwise_phylANOVA <- function(group1, group2) {


  subset_data <- data[data$TrophicLevel %in% c(group1, group2), ]
  subset_tree <- drop.tip(tree, setdiff(tree$tip.label, subset_data$Species))


  trophic <- setNames(as.factor(subset_data$TrophicLevel), subset_data$Species)
  gene_status_subset <- setNames(subset_data$gene_status, subset_data$Species)

  result <- phylANOVA(subset_tree, x = trophic, y = gene_status_subset, nsim = 1000)

  # Print output
  cat("\n", group1, "vs", group2, "\n")
  cat("F-value:", round(result$F, 3), "\n")
  cat("P-value (simulated):", result$Pf, "\n")

  # Only print posthoc if available and not all NA
  if (!is.null(result$posthoc) && any(!is.na(result$posthoc$p))) {
    cat("\nPairwise posthoc p-values (Holm corrected):\n")
    print(result$posthoc$p)
  }
}

run_pairwise_phylANOVA("Carnivore", "Herbivore")
run_pairwise_phylANOVA("Carnivore", "Omnivore")
run_pairwise_phylANOVA("Herbivore", "Omnivore")

subset_data <- data[data$TrophicLevel %in% c("Carnivore", "Herbivore", "Omnivore"), ]
subset_tree <- drop.tip(tree, setdiff(tree$tip.label, subset_data$Species))


trophic <- setNames(as.factor(subset_data$TrophicLevel), subset_data$Species)
gene_status_subset <- setNames(subset_data$gene_status, subset_data$Species)

phylANOVA(subset_tree, x = trophic, y = gene_status_subset, nsim = 1000)



################### Plotting

library(ggplot2)
library(dplyr)

df_summary <- data %>%
  group_by(TrophicLevel) %>%
  summarise(
    mean = mean(gene_status),
    se = sqrt((mean * (1 - mean)) / n())
  )

data$TrophicLevel <- factor(data$TrophicLevel, levels = c("Aquatic", "Carnivore", "Herbivore", "Omnivore"))
df_summary$TrophicLevel <- factor(df_summary$TrophicLevel, levels = levels(data$TrophicLevel))

mean_terrestrial <- mean(data$gene_status[data$TrophicGroup == "Terrestrial"])
mean_aquatic <- mean(data$gene_status[data$TrophicGroup == "Aquatic"])

dot_colors <- c(
  "Aquatic" = "slateblue4",
  "Carnivore" = "sienna4",  # dark orange/brown (sienna)
  "Herbivore" = "seagreen",
  "Omnivore" = "orchid4"
)
# Create plot
p <- ggplot() +
  annotate("rect", xmin = 1.5, xmax = 4.5, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "grey50") +

  geom_bar(data = df_summary, aes(x = TrophicLevel, y = mean, fill = TrophicLevel),
           stat = "identity", color = "black", width = 0.6, alpha = 0.8) +

  geom_errorbar(data = df_summary, aes(x = TrophicLevel, ymin = mean - se, ymax = mean + se),
                width = 0.15, color = "black") +

  geom_jitter(data = data, aes(x = TrophicLevel, y = gene_status, color = TrophicLevel),
              width = 0.15, height = 0.05, alpha = 0.6, size = 2) +

  geom_point(aes(x = 1, y = mean_aquatic), color = "red", size = 4, shape = 18) +
  geom_point(aes(x = 3, y = mean_terrestrial), color = "red", size = 4, shape = 18) +

scale_fill_manual(values = c(
  "Aquatic" = "slateblue1",
  "Carnivore" = "orange1",
  "Herbivore" = "springgreen3",
  "Omnivore" = "orchid1"
)) +
  scale_color_manual(values = dot_colors) +

  labs(
    title = "Gene Presence by Trophic Level",
    x = "Trophic Level",
    y = "Proportion of Species with GPRC6A"
  ) +

  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title.position = "plot",
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.major.y = element_line(color = "gray80", linewidth = 0.4),
    panel.grid.minor.y = element_line(color = "gray90", linewidth = 0.2)
  )

  
ggsave("bar_dot_plot.pdf", plot = p, width = 8, height = 6)
