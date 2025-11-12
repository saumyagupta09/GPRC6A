#!/usr/bin/env python3
"""
phy_stats_sliding_window.py

Read a PHYLIP (sequential or interleaved) nucleotide alignment (codon alignment expected)
and a Newick tree, then compute sliding-window statistics including:
 - counts of synonymous changes per lineage (tip)
 - counts of nonsynonymous changes per lineage
 - counts of radical vs conservative amino-acid changes per lineage
 - Shannon entropy (amino-acid) per window
 - phylogenetically-weighted entropy per lineage (weights = exp(-dist/scale))
 - number of convergent/shared changes that violate the tree per window

Outputs a CSV table with one row per window.

Usage:
    python3 phy_stats_sliding_window.py -a alignment.phy -t tree.nwk -w 30 -s 3 --step 1 -o out.csv

Notes / assumptions:
 - Alignment must be nucleotide codons (length multiple of 3). Gaps ('-' or '.') or Ns in a codon cause that codon to be skipped.
 - Tree tip names must match sequence names in the PHYLIP file.
 - "Lineage" statistics are provided for every tip in the tree; a lineage's count is the sum of substitution events on branches along the path from root to that tip.
 - Radical vs conservative classification uses simple amino-acid group membership (polarity/charge groups). A change between groups = radical; within group = conservative.
 - Ancestral reconstruction uses Fitch parsimony on codon-level characters and a deterministic backtrace to obtain one plausible ancestral state per node.

Dependencies: only Python standard library.

"""

from __future__ import annotations
import argparse
import math
import csv
import sys
from collections import defaultdict, Counter, deque
from typing import List, Dict, Tuple, Optional, Set

# ----------------------------- Basic parsers ---------------------------------

def parse_phylip(path: str) -> Dict[str, str]:
    """Simple PHYLIP parser that handles sequential or interleaved PHYLIP.
    Returns dict name->sequence (no line breaks, uppercase)
    """
    with open(path) as f:
        header = f.readline().strip()
        if not header:
            raise ValueError("Empty PHYLIP file")
        parts = header.split()
        try:
            nseq = int(parts[0])
            seqlen = int(parts[1]) if len(parts) > 1 else None
        except Exception:
            # maybe missing header; try to parse as relaxed PHYLIP
            raise ValueError("PHYLIP header expected (NSeq NCols). Got: %r" % header)

        sequences = {}
        # Read blocks. We will read lines, collect names and sequences until we have nseq sequences with full length
        # Many PHYLIP variants: name and sequence on same line, or interleaved blocks.
        # Strategy: first read nseq lines — extract names if present then continue reading blocks until sequences reach seqlen.
        # Read first block
        first_block = []
        while len(first_block) < nseq:
            line = f.readline()
            if not line:
                break
            line = line.rstrip('\n')
            if not line.strip():
                continue
            first_block.append(line)
        # Parse first block lines: either "name seq" or just seq
        names = []
        seqs = []
        for L in first_block:
            if '\t' in L:
                name, seq = L.split(None, 1)
            else:
                parts = L.strip().split()
                if len(parts) == 1:
                    # no name, just sequence
                    name = None
                    seq = parts[0]
                else:
                    name = parts[0]
                    seq = ''.join(parts[1:])
            names.append(name)
            seqs.append(seq)
        # If names are None, create generic names
        for i, (n, s) in enumerate(zip(names, seqs)):
            if n is None:
                names[i] = f'seq{i+1}'
        # initialize sequences dict
        for n, s in zip(names, seqs):
            sequences[n] = s.replace(' ', '').replace('\t', '').replace('\r', '').upper()
        # read remaining blocks until lengths match seqlen
        if seqlen is not None:
            need_more = any(len(s) < seqlen for s in sequences.values())
            while need_more:
                # read next up to nseq non-empty lines
                block = []
                while len(block) < nseq:
                    line = f.readline()
                    if not line:
                        break
                    line = line.rstrip('\n')
                    if not line.strip():
                        continue
                    block.append(line)
                if not block:
                    break
                # append to sequences in order
                bi = 0
                for L in block:
                    parts = L.strip().split()
                    # block lines in interleaved PHYLIP often have only the sequence
                    seqpart = ''.join(parts) if len(parts) == 1 else ''.join(parts[1:])
                    # find next sequence to append to (order preserved)
                    key = list(sequences.keys())[bi]
                    sequences[key] += seqpart.replace(' ', '').replace('\t', '').upper()
                    bi += 1
                need_more = any(len(s) < seqlen for s in sequences.values())
        return sequences


# Simple Newick parser that supports branch lengths and internal node parentheses but does not require support for comments
class Node:
    def __init__(self, name: Optional[str] = None, length: float = 0.0):
        self.name = name
        self.length = length
        self.children: List[Node] = []
        self.parent: Optional[Node] = None

    def is_leaf(self):
        return len(self.children) == 0

    def traverse(self):
        yield self
        for c in self.children:
            yield from c.traverse()

    def __repr__(self):
        return f"Node({self.name},{self.length},children={len(self.children)})"


def parse_newick(s: str) -> Node:
    s = s.strip()
    if not s.endswith(';'):
        raise ValueError('Newick string must end with ;')
    s = s[:-1]
    i = 0

    def skip_whitespace():
        nonlocal i
        while i < len(s) and s[i].isspace():
            i += 1

    def parse_subtree(parent=None) -> Node:
        nonlocal i
        skip_whitespace()
        if s[i] == '(':
            # internal node
            i += 1
            node = Node()
            while True:
                child = parse_subtree(parent=node)
                child.parent = node
                node.children.append(child)
                skip_whitespace()
                if s[i] == ',':
                    i += 1
                    continue
                elif s[i] == ')':
                    i += 1
                    break
                else:
                    raise ValueError('Unexpected character in Newick: %r at %d' % (s[i], i))
            # now optional name and length
            name, length = parse_name_and_length()
            node.name = name
            node.length = length
            return node
        else:
            # leaf
            name, length = parse_name_and_length()
            node = Node(name=name, length=length)
            return node

    def parse_name_and_length() -> Tuple[Optional[str], float]:
        nonlocal i
        skip_whitespace()
        name = None
        length = 0.0
        # read name (until ':' or ',' or ')' or ';')
        start = i
        while i < len(s) and s[i] not in ':,()':
            i += 1
        token = s[start:i].strip()
        if token != '':
            name = token
        if i < len(s) and s[i] == ':':
            i += 1
            # read number
            start = i
            while i < len(s) and s[i] not in ',()':
                i += 1
            num = s[start:i].strip()
            try:
                length = float(num)
            except:
                length = 0.0
        return name, length

    root = parse_subtree()
    return root


# ---------------------------- Codon / AA helpers ------------------------------
# Standard genetic code mapping
CODON_TABLE = {
    # U/T used as T
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

AA_GROUPS = {
    'nonpolar': set(list('AVLIMFWP')),
    'polar_uncharged': set(list('STNQYC')),
    'positive': set(list('KRH')),
    'negative': set(list('DE'))
}

def aa_group(aa: str) -> str:
    for g, s in AA_GROUPS.items():
        if aa in s:
            return g
    return 'other'


def translate_codon(codon: str) -> str:
    c = codon.replace('U','T').upper()
    if len(c) != 3:
        return 'X'
    if any(ch not in 'ACGT' for ch in c):
        return 'X'
    return CODON_TABLE.get(c, 'X')


# ------------------------- Fitch parsimony (multi-state) ---------------------

def fitch_upward(tip_states: Dict[str, Set[str]], tree_root: Node) -> Dict[Node, Set[str]]:
    """Run Fitch's algorithm upward pass, returning set of possible states per node."""
    node_sets: Dict[Node, Set[str]] = {}

    # post-order traversal
    def postorder(node: Node):
        if node.is_leaf():
            st = tip_states.get(node.name, set(['X']))
            node_sets[node] = set(st)
            return
        for c in node.children:
            postorder(c)
        # intersection of children sets if non-empty else union
        sets = [node_sets[c] for c in node.children]
        inter = set.intersection(*sets)
        if inter:
            node_sets[node] = set(inter)
        else:
            node_sets[node] = set.union(*sets)

    postorder(tree_root)
    return node_sets


def fitch_backtrace(node_sets: Dict[Node, Set[str]], root: Node) -> Dict[Node, str]:
    """Deterministic backtrace: choose one state per node consistent with minimizing changes.
    Root picks arbitrary element from its set; children pick element that intersects with parent's choice if possible, otherwise arbitrary from its set.
    Returns dict node->state (single representative)
    """
    state_choice: Dict[Node, str] = {}
    # pick root
    root_choice = next(iter(node_sets[root]))
    state_choice[root] = root_choice

    # pre-order traversal
    stack = [root]
    while stack:
        nd = stack.pop()
        for c in nd.children:
            # if child's set contains parent's choice, pick that; else pick arbitrary
            if state_choice[nd] in node_sets[c]:
                state_choice[c] = state_choice[nd]
            else:
                state_choice[c] = next(iter(node_sets[c]))
            stack.append(c)
    return state_choice


def map_substitutions(root: Node, state_assign: Dict[Node, str]) -> List[Tuple[Node, Node, str, str]]:
    """Return list of substitutions as (parent_node, child_node, parent_state, child_state) for edges where they differ."""
    subs = []
    for node in root.traverse():
        for c in node.children:
            pstate = state_assign[node]
            cstate = state_assign[c]
            if pstate != cstate:
                subs.append((node, c, pstate, cstate))
    return subs


# --------------------------- Tree utilities ----------------------------------

def get_tip_nodes(root: Node) -> List[Node]:
    return [n for n in root.traverse() if n.is_leaf()]


def leaf_name_to_node_map(root: Node) -> Dict[str, Node]:
    return {n.name: n for n in root.traverse() if n.is_leaf()}


def node_to_root_path(node: Node) -> List[Node]:
    path = []
    while node is not None:
        path.append(node)
        node = node.parent
    return list(reversed(path))


def pairwise_distance(root: Node) -> Dict[Tuple[str,str], float]:
    tips = get_tip_nodes(root)
    names = [t.name for t in tips]
    # compute distances by finding LCA and summing branch lengths
    # precompute node->root cumulative lengths
    cumlen = {}
    def compute_cum(n: Node, acc: float):
        cumlen[n] = acc
        for c in n.children:
            compute_cum(c, acc + c.length)
    compute_cum(root, 0.0)
    # map name->node
    name2node = leaf_name_to_node_map(root)
    out = {}
    for i in range(len(names)):
        for j in range(i+1, len(names)):
            a = name2node[names[i]]
            b = name2node[names[j]]
            # find LCA by following parents sets
            ancestors_a = set(node_to_root_path(a))
            path_b = node_to_root_path(b)
            lca = None
            for nd in path_b:
                if nd in ancestors_a:
                    lca = nd
                    break
            if lca is None:
                dist = cumlen[a] + cumlen[b]
            else:
                dist = cumlen[a] + cumlen[b] - 2 * cumlen[lca]
            out[(names[i], names[j])] = dist
            out[(names[j], names[i])] = dist
    # self distances 0
    for n in names:
        out[(n,n)] = 0.0
    return out


def branchset_for_tip(root: Node, tip_name: str) -> List[Node]:
    name2node = leaf_name_to_node_map(root)
    node = name2node[tip_name]
    return node_to_root_path(node)[1:]  # exclude root

# ------------------------ Sliding window and main logic ----------------------

import argparse

def sliding_windows(aln_len_codons: int, window_codons: int, step_codons: int):
    i = 0
    while i + window_codons <= aln_len_codons:
        yield i, i + window_codons
        i += step_codons


def compute_stats(alignment: Dict[str,str], tree_root: Node, window_codons: int, step_codons: int, out_csv: str):
    taxa = list(alignment.keys())
    seqlen = len(next(iter(alignment.values())))
    if seqlen % 3 != 0:
        raise ValueError('Alignment length not divisible by 3; expected codon alignment')
    n_codons = seqlen // 3

    # precompute pairwise distances and scale
    pairdist = pairwise_distance(tree_root)
    # scale for weighting: mean pairwise distance
    dvals = [v for k,v in pairdist.items() if k[0] != k[1]]
    scale = sum(dvals)/len(dvals) if dvals else 1.0

    tip_nodes = get_tip_nodes(tree_root)
    tip_names = [t.name for t in tip_nodes]

    # prepare CSV header
    header = ['window_start_codon','window_end_codon','n_codons_considered','entropy_aa','convergence_events']
    # per-tip columns for counts
    for t in tip_names:
        header += [f'{t}_syn','{0}_nonsyn'.format(t), f'{t}_radical', f'{t}_conservative', f'{t}_phylo_weighted_entropy']

    with open(out_csv, 'w', newline='') as csvf:
        writer = csv.writer(csvf)
        writer.writerow(header)

        for ws, we in sliding_windows(n_codons, window_codons, step_codons):
            # per-tip event counters (events along branches aggregated to tip path)
            per_tip_syn = Counter()
            per_tip_nonsyn = Counter()
            per_tip_rad = Counter()
            per_tip_cons = Counter()
            # overall amino acid counts for entropy
            aa_counts = Counter()
            total_aa_obs = 0
            # convergence counter
            convergence_count = 0

            # for phylo-weighted entropy we'll compute per tip a weighted frequency across taxa in window
            phylo_weighted_entropy = {t: 0.0 for t in tip_names}

            # accumulate for each codon in window
            for codon_idx in range(ws, we):
                i0 = codon_idx*3
                codon_map = {}
                skip = False
                for tax in taxa:
                    cod = alignment[tax][i0:i0+3].upper()
                    if len(cod) < 3 or any(c in '-.' for c in cod) or 'N' in cod:
                        skip = True
                        break
                    codon_map[tax] = cod.replace('U','T')
                if skip:
                    continue
                # translate to AA
                aa_map = {tax: translate_codon(cd) for tax,cd in codon_map.items()}
                # count amino acids for entropy
                for a in aa_map.values():
                    if a == '*' or a == 'X':
                        continue
                    aa_counts[a] += 1
                    total_aa_obs += 1

                # Fitch parsimony on codon states
                tip_states = {tax: set([codon_map[tax]]) for tax in taxa}
                node_sets = fitch_upward(tip_states, tree_root)
                state_assign = fitch_backtrace(node_sets, tree_root)
                subs = map_substitutions(tree_root, state_assign)

                # For convergence detection: group by derived amino acid and find if same derived AA appears on >1 independent branch
                # We'll collect substitutions as (child_node, derived_aa)
                derived_by_branch = []
                for parent, child, pstate, cstate in subs:
                    # determine parent AA and child AA
                    paa = translate_codon(pstate)
                    caa = translate_codon(cstate)
                    # ignore stop or unknown
                    if caa in ('*','X'):
                        continue
                    # record event
                    derived_by_branch.append((child, caa, parent, pstate, cstate))
                    # classify syn/nonsyn
                    if paa == caa:
                        # synonymous
                        change_syn = True
                    else:
                        change_syn = False
                    # radical/conservative by aa groups
                    grp_p = aa_group(paa)
                    grp_c = aa_group(caa)
                    is_radical = (grp_p != grp_c)

                    # assign this event to descendant tips: we will add counts along each tip whose path includes this child branch
                    # compute descendant tips under child
                    desc = [leaf.name for leaf in child.traverse() if leaf.is_leaf()]
                    for tip in desc:
                        if change_syn:
                            per_tip_syn[tip] += 1
                        else:
                            per_tip_nonsyn[tip] += 1
                        if is_radical:
                            per_tip_rad[tip] += 1
                        else:
                            per_tip_cons[tip] += 1

                # convergence counting: if same derived aa appears in >=2 branches that are independent (neither is ancestor of the other), count as convergence event
                # a cheap test: if two branches' child nodes are not ancestor-descendant
                for i in range(len(derived_by_branch)):
                    for j in range(i+1, len(derived_by_branch)):
                        b1 = derived_by_branch[i][0]
                        b2 = derived_by_branch[j][0]
                        aa1 = derived_by_branch[i][1]
                        aa2 = derived_by_branch[j][1]
                        if aa1 != aa2:
                            continue
                        # check ancestor-descendant
                        path1 = set(node_to_root_path(b1))
                        path2 = set(node_to_root_path(b2))
                        if (b1 in path2) or (b2 in path1):
                            # not independent
                            continue
                        convergence_count += 1

            # compute entropy
            entropy = 0.0
            if total_aa_obs > 0:
                for aa, cnt in aa_counts.items():
                    p = cnt/total_aa_obs
                    entropy -= p * math.log2(p)

            # phylogenetically-weighted entropy per tip: for each tip, weight other taxa by exp(-dist/scale)
            # compute per tip weighted aa frequencies across taxa in window
            # For performance we recompute aa frequencies per taxon across window: we already counted per codon, but easier: aggregate amino-acids per taxon across window
            aa_by_taxon = {t: [] for t in taxa}
            for codon_idx in range(ws, we):
                i0 = codon_idx*3
                for tax in taxa:
                    cod = alignment[tax][i0:i0+3].upper()
                    if len(cod) < 3 or any(c in '-.' for c in cod) or 'N' in cod:
                        continue
                    aa = translate_codon(cod)
                    if aa in ('*','X'):
                        continue
                    aa_by_taxon[tax].append(aa)
            # for each tip compute weighted freq
            for tip in tip_names:
                # compute weights to other taxa
                weights = {}
                for other in taxa:
                    d = pairdist.get((tip, other), 0.0)
                    weights[other] = math.exp(-d/scale) if scale>0 else 1.0
                total_w = sum(weights.values())
                freq = Counter()
                tot = 0.0
                for other in taxa:
                    w = weights[other]
                    aas = aa_by_taxon[other]
                    for a in aas:
                        freq[a] += w
                        tot += w
                # normalize
                ent = 0.0
                if tot > 0:
                    for a,cnt in freq.items():
                        p = cnt / tot
                        if p>0:
                            ent -= p * math.log2(p)
                phylo_weighted_entropy[tip] = ent

            # prepare row
            n_considered = sum(1 for _ in aa_counts.values())
            row = [ws, we, total_aa_obs, entropy, convergence_count]
            for t in tip_names:
                row += [per_tip_syn[t], per_tip_nonsyn[t], per_tip_rad[t], per_tip_cons[t], phylo_weighted_entropy[t]]
            writer.writerow(row)


# ------------------------------- CLI ----------------------------------------

def main():
    parser = argparse.ArgumentParser(description='Compute sliding-window phylogenetic statistics from a codon PHYLIP alignment and a Newick tree')
    parser.add_argument('-a','--alignment', required=True, help='PHYLIP alignment file (codon alignment)')
    parser.add_argument('-t','--tree', required=True, help='Newick tree file with branch lengths')
    parser.add_argument('-w','--window', type=int, default=30, help='window size in codons (default 30)')
    parser.add_argument('-s','--step', type=int, default=3, help='step size in codons (default 3)')
    parser.add_argument('-o','--out', default='sliding_stats.csv', help='CSV output file')
    args = parser.parse_args()

    alignment = parse_phylip(args.alignment)
    with open(args.tree) as f:
        newick = f.read().strip()
    tree_root = parse_newick(newick)

    # simple check that tip names match
    names_in_aln = set(alignment.keys())
    names_in_tree = set(n.name for n in get_tip_nodes(tree_root))
    if not names_in_aln.issubset(names_in_tree) or not names_in_tree.issubset(names_in_aln):
        print('Warning: tip names in alignment and tree do not match exactly.', file=sys.stderr)
        print('In alignment but not in tree:', names_in_aln - names_in_tree, file=sys.stderr)
        print('In tree but not in alignment:', names_in_tree - names_in_aln, file=sys.stderr)

    compute_stats(alignment, tree_root, args.window, args.step, args.out)
    print('Wrote results to', args.out)

if __name__ == '__main__':
    main()
