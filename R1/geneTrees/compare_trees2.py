#!/usr/bin/env python3
import argparse
from Bio import Phylo
from ete3 import Tree
import matplotlib.pyplot as plt
import dendropy
import csv
import os

def load_tree_ete(file):
    """Load Newick tree robustly into ETE3 Tree."""
    with open(file, "r", encoding="utf-8", errors="ignore") as f:
        newick = f.read().strip()
    try:
        return Tree(newick, format=1)
    except:
        return Tree(newick, format=0)

def load_tree_bio(file):
    """Load tree into Biopython Phylo Tree for visualization."""
    return Phylo.read(file, "newick")

def compare_trees(gene_tree_file, species_tree_file, output_prefix, outgroup=None):
    # --- Load trees in both frameworks
    gene_tree = load_tree_ete(gene_tree_file)
    species_tree = load_tree_ete(species_tree_file)
    gene_bio = load_tree_bio(gene_tree_file)
    species_bio = load_tree_bio(species_tree_file)

    # --- Optionally reroot both trees
    if outgroup:
        try:
            gene_tree.set_outgroup(outgroup)
            species_tree.set_outgroup(outgroup)
            gene_bio.root_with_outgroup(outgroup)
            species_bio.root_with_outgroup(outgroup)
        except Exception as e:
            print(f"Warning: Could not reroot trees with {outgroup}: {e}")

    # --- Compute Robinson-Foulds distance using ETE3
    rf_result = gene_tree.robinson_foulds(species_tree, unrooted_trees=True)
    rf_dist, max_rf, common_leaves, parts_t1, parts_t2, *_ = rf_result

    print(f"Robinson-Foulds distance: {rf_dist}/{max_rf}")
    print(f"Common leaves ({len(common_leaves)}): {common_leaves}")

    # --- Convert to DendroPy for quartet distance
    gene_den = dendropy.Tree.get(path=gene_tree_file, schema="newick")
    species_den = dendropy.Tree.get(path=species_tree_file, schema="newick")

    taxon_namespace = dendropy.TaxonNamespace()
    gene_den.migrate_taxon_namespace(taxon_namespace)
    species_den.migrate_taxon_namespace(taxon_namespace)

    # Compute normalized quartet distance (if supported)
    try:
        from dendropy.calculate.treecompare import symmetric_difference
        quartet_dist = symmetric_difference(gene_den, species_den)
    except Exception:
        quartet_dist = None

    if quartet_dist is not None:
        print(f"Quartet distance: {quartet_dist}")
    else:
        print("Quartet distance could not be computed (DendroPy version limitation).")

    # --- Identify unique and shared splits
    unique_gene = parts_t1 - parts_t2
    unique_species = parts_t2 - parts_t1
    shared = parts_t1 & parts_t2

    # --- Save summary CSV
    out_csv = f"{output_prefix}_split_summary.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Type", "Num_Taxa", "Taxa"])
        for split in shared:
            writer.writerow(["Shared", len(split), ";".join(sorted(map(str, split)))])
        for split in unique_gene:
            writer.writerow(["Gene_Unique", len(split), ";".join(sorted(map(str, split)))])
        for split in unique_species:
            writer.writerow(["Species_Unique", len(split), ";".join(sorted(map(str, split)))])

    print(f"Split summary written to {out_csv}")

    # --- Plot side-by-side phylogenies
    fig, axes = plt.subplots(1, 2, figsize=(32, 16), dpi=400)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, wspace=0.3)

    Phylo.draw(
        gene_bio,
        do_show=False,
        axes=axes[0],
        label_func=lambda x: x.name if x else "",
    )
    axes[0].set_title("Gene Tree")

    Phylo.draw(
        species_bio,
        do_show=False,
        axes=axes[1],
        label_func=lambda x: x.name if x else "",
    )
    axes[1].set_title("Species Tree")

    plt.tight_layout()
    out_png = f"{output_prefix}_trees.png"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"Side-by-side tree comparison saved to {out_png}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare gene and species trees (RF + visualization)")
    parser.add_argument("gene_tree", help="Path to gene tree in Newick format")
    parser.add_argument("species_tree", help="Path to species tree in Newick format")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--outgroup", help="Outgroup species name for rerooting", default=None)
    args = parser.parse_args()
    compare_trees(args.gene_tree, args.species_tree, args.output_prefix, args.outgroup)
