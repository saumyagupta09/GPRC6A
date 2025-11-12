#!/usr/bin/env python3
"""
msa_reliability.py

Input:
 - codon-aligned FASTA (length multiple of 3)
 - Newick tree (leaf names must match FASTA IDs)
Options:
 - sliding window (codons)
 - thresholds for gap_frac, entropy, parsimony steps
Outputs:
 - CSVs: per_site, per_sequence, per_clade, per_branch_site
 - BED-like codon regions of low reliability
 - Optionally a filtered FASTA with low-reliability codons masked (NNN) or removed while preserving frame
"""

import argparse
import math
from collections import Counter, defaultdict
from copy import deepcopy

import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

# -----------------------
# Utilities: entropy, codon grouping
# -----------------------
def shannon_entropy(counts):
    total = sum(counts)
    if total == 0:
        return 0.0
    ent = 0.0
    for c in counts:
        if c == 0:
            continue
        p = c / total
        ent -= p * math.log2(p)
    return ent

def codon_of_alignment(seq_str, codon_index):
    """Return 3-base codon string at codon_index (0-based)"""
    start = codon_index * 3
    return seq_str[start:start+3]

# -----------------------
# Fitch parsimony for codon states
# We'll treat each codon (3-mer) as a discrete state (strings), with '-' or gaps unified to a special state.
# Implement two-pass Fitch to get minimum changes and then assign specific states to nodes to map branch changes.
# -----------------------
def is_gap_codon(codon):
    return ('-' in codon) or (len(codon) < 3)

def fitch_up_pass(tree, leaf_state_map):
    """
    tree: Biopython Clade (rooted or unrooted)
    leaf_state_map: dict leaf_name -> state (string)
    returns: dict node -> set(states), total_changes (sum over internal merges)
    """
    node_sets = {}
    changes = 0

    def postorder(node):
        nonlocal changes
        if node.is_terminal():
            s = leaf_state_map.get(node.name, None)
            if s is None:
                node_sets[node] = set()
            else:
                node_sets[node] = set([s])
            return node_sets[node]
        child_sets = []
        for ch in node.clades:
            child_sets.append(postorder(ch))
        # Fitch: intersection if non-empty, else union and increment changes
        inter = set.intersection(*child_sets)
        if len(inter) > 0:
            node_sets[node] = inter
        else:
            node_sets[node] = set.union(*child_sets)
            changes += 1
        return node_sets[node]

    postorder(tree.root if hasattr(tree, 'root') else tree.clade)
    return node_sets, changes

def fitch_down_pass_assign(tree, node_sets, leaf_state_map):
    """
    Assign an actual state to each node (choose one arbitrarily from node_sets) consistent with Fitch.
    Then compute which branches have state changes.
    returns: dict node -> assigned_state, list of (parent_node, child_node) that have a change
    """
    assigned = {}
    changes_on_branches = []

    # choose root assignment: pick any element from its set (prefer majority from leaves if possible)
    root = tree.root if hasattr(tree, 'root') else tree.clade
    # pick a state for root: if leaf_state_map majority in root set, pick that
    root_choice = None
    root_candidates = node_sets[root]
    if len(root_candidates) == 0:
        root_choice = None
    else:
        # rank candidates by global frequency in leaves
        freq = Counter()
        for leaf, s in leaf_state_map.items():
            if s is not None:
                freq[s] += 1
        # pick highest freq that is in candidates, else arbitrary
        for k, _ in freq.most_common():
            if k in root_candidates:
                root_choice = k
                break
        if root_choice is None:
            root_choice = next(iter(root_candidates))
    assigned[root] = root_choice

    def preorder(node):
        for ch in node.clades:
            # choose child's state: if parent's assigned in child's set, pick it; else pick any from child's set
            child_set = node_sets[ch]
            if len(child_set) == 0:
                assigned[ch] = None
            else:
                if assigned[node] in child_set:
                    assigned[ch] = assigned[node]
                else:
                    assigned[ch] = next(iter(child_set))
            # mark change
            if assigned[node] != assigned[ch]:
                changes_on_branches.append((node, ch))
            preorder(ch)

    preorder(root)
    return assigned, changes_on_branches

# Map branches to IDs (string)
def branch_id(parent, child):
    # use node names when available, otherwise object ids
    p = getattr(parent, 'name', None) or f'NODE_{id(parent)}'
    c = getattr(child, 'name', None) or f'NODE_{id(child)}'
    return f'{p}--{c}'

# -----------------------
# Core analyzer
# -----------------------
class MSAReliability:
    def __init__(self, fasta_path, tree_path, window_codons=5):
        self.fasta_path = fasta_path
        self.tree_path = tree_path
        self.window_codons = window_codons

        self.records = list(SeqIO.parse(fasta_path, 'fasta'))
        if len(self.records) == 0:
            raise ValueError("No sequences found in FASTA.")
        # ensure all same length
        lengths = set(len(r.seq) for r in self.records)
        if len(lengths) != 1:
            raise ValueError("Sequences have unequal lengths.")
        self.L = lengths.pop()
        if self.L % 3 != 0:
            raise ValueError("Alignment length is not a multiple of 3 (not codon aligned).")
        self.n_codons = self.L // 3

        # map id -> seq_str
        self.seq_map = {r.id: str(r.seq).upper() for r in self.records}
        # load tree
        self.tree = Phylo.read(tree_path, 'newick')
        # make sure all leaf names in tree match sequences
        leaves = [term.name for term in self.tree.get_terminals()]
        missing = [l for l in leaves if l not in self.seq_map]
        if missing:
            raise ValueError(f"The following leaves are missing from the FASTA: {missing}")

        # containers for outputs
        self.per_site = []
        self.branch_names = []
        self.branch_site_changes = defaultdict(dict)  # branch_id -> codon_idx -> 1/0
        self.per_sequence = {}
        self.per_clade = {}  # computed later if clade mapping given

    def analyze(self):
        # precompute leaf list
        leaves = [term.name for term in self.tree.get_terminals()]
        # iterate codons
        for cidx in range(self.n_codons):
            # gather codon states for all leaves
            codon_states = {}
            for leaf in leaves:
                cod = codon_of_alignment(self.seq_map[leaf], cidx)
                if is_gap_codon(cod):
                    codon_states[leaf] = None  # treat as missing
                else:
                    codon_states[leaf] = cod

            # entropy (codon-level)
            counts = Counter([s for s in codon_states.values() if s is not None])
            ent = shannon_entropy(list(counts.values()))

            # gap fraction
            gap_frac = sum(1 for s in codon_states.values() if s is None) / len(leaves)

            # parsimony (Fitch)
            node_sets, fitch_changes = fitch_up_pass(self.tree, {k: (v if v is not None else None) for k, v in codon_states.items()})
            assigned, changes_on_branches = fitch_down_pass_assign(self.tree, node_sets, codon_states)

            # map branch changes
            branch_changes = set(branch_id(p, c) for (p, c) in changes_on_branches)
            for b in branch_changes:
                self.branch_site_changes[b][cidx] = 1
                if b not in self.branch_names:
                    self.branch_names.append(b)

            # site score: combine metrics (example simple scoring, you can tune)
            # higher entropy, higher gap frac, higher parsimony changes => less reliable
            score = (ent / (math.log2(61)+1e-9)) * 0.4 + gap_frac * 0.3 + (fitch_changes / (len(leaves)-1+1e-9)) * 0.3
            # normalized-ish score in [0..~1], lower is better reliability
            self.per_site.append({
                'codon_index': cidx,
                'start_pos': cidx*3,
                'end_pos': cidx*3+3,
                'gap_frac': gap_frac,
                'entropy': ent,
                'parsimony_changes': fitch_changes,
                'combined_unreliability_score': score
            })

        # per-sequence metrics
        seq_stats = {}
        threshold_bad = 0.5  # example threshold for codon considered 'bad' per-site; adjust via CLI
        bad_codons = set(idx for idx, s in enumerate(self.per_site) if s['combined_unreliability_score'] >= threshold_bad)
        for rec in self.records:
            seqid = rec.id
            seq_str = str(rec.seq)
            seq_gap_frac = sum(1 for i in range(self.L) if seq_str[i] == '-') / self.L
            # fraction of codons that are bad and where this sequence is non-gap
            carry_bad = 0
            total_non_gap_codons = 0
            for cidx in range(self.n_codons):
                cod = codon_of_alignment(seq_str, cidx)
                if is_gap_codon(cod):
                    continue
                total_non_gap_codons += 1
                if cidx in bad_codons:
                    carry_bad += 1
            frac_bad = (carry_bad / total_non_gap_codons) if total_non_gap_codons > 0 else np.nan
            seq_stats[seqid] = {'seq_gap_frac': seq_gap_frac, 'frac_bad_codons': frac_bad, 'total_non_gap_codons': total_non_gap_codons}

        self.per_sequence = seq_stats

    def export_csvs(self, out_prefix='msa_reliability'):
        df_site = pd.DataFrame(self.per_site)
        df_site.to_csv(f'{out_prefix}_per_site.csv', index=False)

        df_seq = pd.DataFrame.from_dict(self.per_sequence, orient='index').reset_index().rename(columns={'index':'seqid'})
        df_seq.to_csv(f'{out_prefix}_per_sequence.csv', index=False)

        # branch-site matrix
        # make large DataFrame: rows branches, columns codon_index, values 0/1
        branch_list = sorted(self.branch_names)
        rows = []
        for b in branch_list:
            row = {'branch': b}
            for c in range(self.n_codons):
                row[f'codon_{c}'] = int(self.branch_site_changes.get(b, {}).get(c, 0))
            rows.append(row)
        df_branch = pd.DataFrame(rows)
        df_branch.to_csv(f'{out_prefix}_per_branch_site.csv', index=False)
        print(f'Wrote {out_prefix}_per_site.csv, _per_sequence.csv, _per_branch_site.csv')

    def find_low_regions(self, score_threshold=0.5, min_codons=2):
        """Return list of (start_codon, end_codon, mean_score) for contiguous low score runs"""
        out = []
        in_run = False
        run_start = None
        for c in range(self.n_codons):
            s = self.per_site[c]['combined_unreliability_score']
            if s >= score_threshold:
                if not in_run:
                    in_run = True
                    run_start = c
            else:
                if in_run:
                    run_end = c - 1
                    if (run_end - run_start + 1) >= min_codons:
                        mean_score = np.mean([self.per_site[i]['combined_unreliability_score'] for i in range(run_start, run_end+1)])
                        out.append((run_start, run_end, mean_score))
                    in_run = False
        # finish
        if in_run:
            run_end = self.n_codons - 1
            if (run_end - run_start + 1) >= min_codons:
                mean_score = np.mean([self.per_site[i]['combined_unreliability_score'] for i in range(run_start, run_end+1)])
                out.append((run_start, run_end, mean_score))
        return out

    def mask_low_regions_in_fasta(self, low_regions, action='mask'):
        """
        action: 'mask' -> replace codon with NNN
                'remove' -> remove codon (keeps frame: but final sequences shorter; note downstream tools expect same length, so prefer mask)
        returns list of SeqRecord
        """
        new_records = []
        mask_indices = set()
        for (s, e, _) in low_regions:
            for c in range(s, e+1):
                mask_indices.add(c)
        for rec in self.records:
            s = str(rec.seq)
            out_chars = list(s)
            for c in mask_indices:
                start = c*3
                if action == 'mask':
                    out_chars[start:start+3] = list('NNN')
                elif action == 'remove':
                    out_chars[start:start+3] = ['','', '']  # not optimal
            new_seq = ''.join([ch for ch in out_chars if ch != ''])
            new_records.append(SeqRecord(Seq(new_seq), id=rec.id, description=''))
        return new_records


# -----------------------
# CLI
# -----------------------
def main():
    parser = argparse.ArgumentParser(description="Assess MSA reliability (codon-aware), export per-site, per-sequence, per-branch metrics, and mask/remove poor regions while preserving codon structure.")
    parser.add_argument('fasta', help='Codon-aligned FASTA')
    parser.add_argument('tree', help='Newick tree matching FASTA headers')
    parser.add_argument('--window', type=int, default=5, help='sliding window in codons')
    parser.add_argument('--score-threshold', type=float, default=0.5, help='combined unreliability threshold (0..1) for low regions')
    parser.add_argument('--min-codons', type=int, default=2, help='minimum contiguous codons to report as region')
    parser.add_argument('--mask-action', choices=['mask','remove'], default='mask', help="mask -> replace codon by NNN, remove -> delete codon (may break alignment length consistency)")
    parser.add_argument('--out-prefix', default='msa_reliability', help='output prefix')
    parser.add_argument('--plot', action='store_true', help='produce a quick plot of per-site combined score')
    args = parser.parse_args()

    R = MSAReliability(args.fasta, args.tree, window_codons=args.window)
    R.analyze()
    R.export_csvs(out_prefix=args.out_prefix)

    low_regions = R.find_low_regions(score_threshold=args.score_threshold, min_codons=args.min_codons)
    # write low regions to TSV (codon coordinates)
    with open(f'{args.out_prefix}_low_regions_codons.tsv', 'w') as fh:
        fh.write('start_codon\tend_codon\tmean_score\n')
        for s,e,m in low_regions:
            fh.write(f'{s}\t{e}\t{m:.5f}\n')
    print(f'Found {len(low_regions)} low-reliability regions (codons). Wrote {args.out_prefix}_low_regions_codons.tsv')

    # mask sequences if requested
    masked = R.mask_low_regions_in_fasta(low_regions, action=args.mask_action)
    SeqIO.write(masked, f'{args.out_prefix}_filtered.fasta', 'fasta')
    print(f'Wrote filtered FASTA: {args.out_prefix}_filtered.fasta (action={args.mask_action})')

    if args.plot:
        df_site = pd.DataFrame(R.per_site)
        plt.figure(figsize=(12,3))
        plt.plot(df_site['codon_index'], df_site['combined_unreliability_score'])
        plt.xlabel('Codon index')
        plt.ylabel('Combined unreliability score')
        plt.title('Per-codon combined unreliability')
        plt.axhline(args.score_threshold, color='red', linestyle='--')
        plt.tight_layout()
        plt.savefig(f'{args.out_prefix}_per_site_plot.png')
        print(f'Wrote plot {args.out_prefix}_per_site_plot.png')

if __name__ == '__main__':
    main()
