#!/usr/bin/env python3
"""
generate_lineage_mapping.py

Automatically groups taxa from a Newick tree into lineage clusters and
produces a mapping file usable by the contrastive learning script.

Usage:
    python generate_lineage_mapping.py --tree treefile.nwk --out mapping.tsv --n_lineages 3
"""

import argparse
from Bio import Phylo
import numpy as np
import pandas as pd
from itertools import count

def cluster_tree_by_depth(tree, n_lineages=3):
    """
    Splits a tree into roughly n_lineages major clades.
    Uses cumulative branch length to find split points.
    """
    # Convert to ultrametric (approximate) by scaling depths
    depths = tree.depths()
    tip_depths = np.array([depths[t] for t in tree.get_terminals()])
    split_thresholds = np.linspace(tip_depths.min(), tip_depths.max(), n_lineages+1)[1:-1]

    clusters = {}
    lineage_id = count(1)

    def assign_clade(clade, current_label=None):
        depth = depths[clade]
        if current_label is None:
            current_label = f"Lineage_{next(lineage_id)}"
        for sub in clade.clades:
            sub_depth = depths[sub]
            if sub.is_terminal():
                clusters[sub.name] = current_label
            else:
                # If subclade depth crosses threshold, start a new lineage
                if np.any(sub_depth > split_thresholds):
                    current_label = f"Lineage_{next(lineage_id)}"
                assign_clade(sub, current_label)

    assign_clade(tree.root)
    return clusters

def main():
    ap = argparse.ArgumentParser(description="Generate lineage mapping file from Newick tree.")
    ap.add_argument("--tree", required=True, help="Input Newick tree file")
    ap.add_argument("--out", required=True, help="Output mapping TSV file")
    ap.add_argument("--n_lineages", type=int, default=3, help="Approximate number of lineages")
    args = ap.parse_args()

    tree = Phylo.read(args.tree, "newick")
    mapping = cluster_tree_by_depth(tree, args.n_lineages)
    df = pd.DataFrame(list(mapping.items()), columns=["sequence_id", "lineage"])
    df.to_csv(args.out, sep="\t", index=False)
    print(f"[+] Saved {len(df)} entries to {args.out}")

if __name__ == "__main__":
    main()
