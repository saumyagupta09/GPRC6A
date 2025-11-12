#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compare a gene tree and species tree, compute topology differences (RF distance),
and plot readable, labeled trees with colored branches (no GUI dependencies).
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from ete3 import Tree
from io import StringIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from matplotlib.collections import LineCollection


def load_tree(tree_file):
    """Load and sanitize a Newick tree."""
    with open(tree_file, "r", encoding="utf-8", errors="ignore") as f:
        newick = f.read().strip()
    cleaned = newick.replace("//", "/").replace("[&", "").replace("]", "")
    try:
        return Tree(cleaned, format=1), cleaned
    except Exception:
        return Tree(cleaned, format=0), cleaned


def compare_splits(t1, t2):
    """Return internal splits for each tree."""
    def get_splits(tree):
        leaves = set(tree.get_leaf_names())
        splits = set()
        for node in tree.traverse():
            if not node.is_leaf():
                subset = frozenset(node.get_leaf_names())
                if 1 < len(subset) < len(leaves):
                    splits.add(subset)
        return splits

    return get_splits(t1), get_splits(t2)


def assign_branch_colors(tree, shared_splits, unique_color, shared_color):
    """Assign color per branch based on whether split is shared."""
    branch_colors = {}
    for node in tree.traverse():
        if node.is_leaf():
            branch_colors[node.name] = "black"
        else:
            subset = frozenset(node.get_leaf_names())
            branch_colors[node.name] = (
                shared_color if subset in shared_splits else unique_color
            )
    return branch_colors


def compute_coords(tree):
    """Compute coordinates for plotting a ladderized tree."""
    coords = {}
    leaves = tree.get_terminals()
    for i, leaf in enumerate(leaves):
        coords[leaf] = (tree.distance(leaf), i)

    for clade in tree.get_nonterminals(order="postorder"):
        child_coords = [coords[c][1] for c in clade.clades if c in coords]
        if not child_coords:
            continue
        y = sum(child_coords) / len(child_coords)
        x = tree.distance(clade)
        coords[clade] = (x, y)
    return coords


def draw_colored_tree(ax, newick_str, branch_colors, title):
    """Draw a simple, readable tree with labeled tips."""
    tree = Phylo.read(StringIO(newick_str), "newick")
    tree.ladderize()  # order branches for clarity
    coords = compute_coords(tree)

    lines = []
    colors = []

    for clade in tree.find_clades(order="level"):
        x, y = coords.get(clade, (0, 0))
        for child in clade.clades:
            if child not in coords:
                continue
            x2, y2 = coords[child]
            lines.append(((x, y), (x2, y2)))
            col = branch_colors.get(
                getattr(child, "name", ""),
                branch_colors.get(getattr(clade, "name", ""), "black"),
            )
            colors.append(col)

    # Draw branch segments
    lc = LineCollection(lines, colors=colors, linewidths=2)
    ax.add_collection(lc)

    # Draw tip labels
    for leaf in tree.get_terminals():
        x, y = coords[leaf]
        ax.text(x + 0.02, y, leaf.name, va="center", fontsize=9)

    ax.autoscale()
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.axis("off")


def compare_trees(gene_tree_file, species_tree_file, output_prefix):
    """Compare gene vs. species tree and plot results."""
    gene_tree, gene_newick = load_tree(gene_tree_file)
    species_tree, species_newick = load_tree(species_tree_file)

    # RF distance
    rf_results = gene_tree.robinson_foulds(species_tree, unrooted_trees=True)
    rf, max_rf, common_leaves = rf_results[:3]
    print(f"Robinson-Foulds distance: {rf} / {max_rf}")
    print(f"Common leaves ({len(common_leaves)}): {', '.join(sorted(common_leaves))}")

    # Shared splits
    splits_gene, splits_species = compare_splits(gene_tree, species_tree)
    shared_splits = splits_gene.intersection(splits_species)
    print(f"Shared splits: {len(shared_splits)}")

    # Assign colors
    gene_colors = assign_branch_colors(gene_tree, shared_splits, "red", "green")
    species_colors = assign_branch_colors(species_tree, shared_splits, "blue", "green")

    # Plot side-by-side readable trees
    fig, axes = plt.subplots(1, 2, figsize=(16, 9))
    draw_colored_tree(axes[0], gene_newick, gene_colors, "Gene Tree (Red = unique)")
    draw_colored_tree(axes[1], species_newick, species_colors, "Species Tree (Blue = unique)")

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_colored_trees.png", dpi=300)
    plt.close()

    # Output summary CSV
    summary = {
        "RF_distance": [rf],
        "Max_RF": [max_rf],
        "Shared_splits": [len(shared_splits)],
        "Common_leaves": [",".join(sorted(common_leaves))],
    }
    pd.DataFrame(summary).to_csv(f"{output_prefix}_summary.csv", index=False)

    print("\nOutputs written:")
    print(f"  {output_prefix}_summary.csv")
    print(f"  {output_prefix}_colored_trees.png")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare gene and species trees; highlight topological differences."
    )
    parser.add_argument("--gene_tree", required=True, help="Gene tree file (Newick)")
    parser.add_argument("--species_tree", required=True, help="Species tree file (Newick)")
    parser.add_argument("--output_prefix", required=True, help="Prefix for output files")
    args = parser.parse_args()
    compare_trees(args.gene_tree, args.species_tree, args.output_prefix)
