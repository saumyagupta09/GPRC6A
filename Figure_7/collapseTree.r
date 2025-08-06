library(ape)
library(phytools)
library(RColorBrewer)

# === INPUT FILES ===
tree_file <- "Terrestrial_Aquatic_species.nwk"
metadata_file <- "Terrestrial_Aquatic.tsv"  # columns: species, gene_status, trophic_level

# === READ DATA ===
tree <- read.tree(tree_file)
metadata <- read.table(metadata_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Filter metadata for only species in tree
metadata <- metadata[metadata$species %in% tree$tip.label, ]
tree <- drop.tip(tree, setdiff(tree$tip.label, metadata$species))
metadata <- metadata[match(tree$tip.label, metadata$species), ]

# Combine traits
combined_trait <- paste(metadata$gene_status, metadata$trophic_level, sep="_")
names(combined_trait) <- metadata$species

# === Collapse clades with uniform traits and track summary ===
collapse_clades <- function(tree, traits) {
tip_traits <- traits[tree$tip.label]
collapse_info <- list()
trait_counter <- list()
  
  repeat {
    to_collapse <- list()
    node_list <- (length(tree$tip.label)+1):(length(tree$tip.label)+tree$Nnode)
    
    for (node in node_list) {
      desc <- getDescendants(tree, node)
      tip_desc <- desc[desc <= length(tree$tip.label)]
      if (length(tip_desc) < 2) next
      
      labels <- tree$tip.label[tip_desc]
      trait_vals <- traits[labels]
      
      if (length(unique(trait_vals)) == 1) {
        to_collapse[[length(to_collapse)+1]] <- list(
          tips=labels,
          trait=trait_vals[1],
          node=node
        )
      }
    }
    
    if (length(to_collapse) == 0) break
    
    collapse <- to_collapse[[1]]
    tips_to_remove <- collapse$tips
    trait_val <- collapse$trait
    mrca_node <- findMRCA(tree, tips_to_remove)
    
    # Count trait usage to ensure unique tip names
    if (is.null(trait_counter[[trait_val]])) {
      trait_counter[[trait_val]] <- 1
    } else {
      trait_counter[[trait_val]] <- trait_counter[[trait_val]] + 1
    }
    unique_id <- trait_counter[[trait_val]]
    
    n_collapsed <- length(tips_to_remove)
    gene_status <- strsplit(trait_val, "_")[[1]][1]
    trophic_level <- strsplit(trait_val, "_")[[1]][2]
    new_tip <- paste0("Collapsed_", trait_val, "_n", n_collapsed, "_id", unique_id)
    
    # === Branch length and depth ===
    subtree <- keep.tip(tree, tips_to_remove)
    avg_branch_length <- mean(subtree$edge.length)
    depths <- node.depth.edgelength(subtree)[1:length(tips_to_remove)]
    avg_depth <- mean(depths)
    
    # Record summary
    collapse_info[[length(collapse_info)+1]] <- data.frame(
      CollapsedTip = new_tip,
      N_Species = n_collapsed,
      GeneStatus = gene_status,
      TrophicLevel = trophic_level,
      AvgBranchLength = round(avg_branch_length, 5),
      AvgNodeDepth = round(avg_depth, 5),
      Species = paste(tips_to_remove, collapse = ",")
    )
    
    # Drop tips
    tree <- drop.tip(tree, tips_to_remove)
    
    # Find correct edge to bind tip
    edge_to_bind <- which.edge(tree, mrca_node)
    
    # Bind tip
    tree <- bind.tip(tree, new_tip, edge=edge_to_bind, position=0)
    
    # Add to trait list
    traits[new_tip] <- trait_val
  }
  
  collapsed_summary <- do.call(rbind, collapse_info)
  return(list(tree=tree, traits=traits, summary=collapsed_summary))
}

# Run collapsing
result <- collapse_clades(tree, combined_trait)
tree_collapsed <- result$tree
trait_map <- result$traits
collapsed_summary <- result$summary

# === Assign colors ===
unique_traits <- unique(trait_map)
colors <- setNames(brewer.pal(n = max(3, length(unique_traits)), "Set3"), unique_traits)
tip_colors <- colors[trait_map[tree_collapsed$tip.label]]

# === Plot ===
plot(tree_collapsed, show.tip.label=TRUE, cex=0.6)
tiplabels(pch=19, col=tip_colors, cex=0.5)

# === Prepare trait table for FigTree ===
trait_df <- data.frame(
  Taxon = tree_collapsed$tip.label,
  Trait = trait_map[tree_collapsed$tip.label],
  GeneStatus = sapply(strsplit(trait_map[tree_collapsed$tip.label], "_"), `[`, 1),
  TrophicLevel = sapply(strsplit(trait_map[tree_collapsed$tip.label], "_"), `[`, 2)
)

# === Write NEXUS with TRAITS block ===
write.nexus.with.traits <- function(tree, traits, file) {
  write.nexus(tree, file=file)
  lines <- readLines(file)
  insert_at <- grep("END;", lines)[1]
  
  trait_lines <- apply(traits, 1, function(row) paste(row, collapse="\t"))
  
  traits_block <- c(
    "",
    "BEGIN TRAITS;",
    paste0("    Dimensions NTax=", nrow(traits), ";"),
    "    Format labels=yes missing=? separator=Tab;",
    "    Matrix",
    trait_lines,
    ";",
    "END;"
  )
  
  lines <- append(lines, traits_block, after=insert_at)
  writeLines(lines, con=file)
}

# === OUTPUT FILES ===
nexus_file <- "collapsed_tree_with_traits.nexus"
summary_file <- "collapsed_clade_summary.tsv"

write.nexus.with.traits(tree_collapsed, trait_df, nexus_file)
write.table(collapsed_summary, file = summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

cat("✔ Collapsed tree written to:", nexus_file, "\n")
cat("✔ Clade summary written to:", summary_file, "\n")
