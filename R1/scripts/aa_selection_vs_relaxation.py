#!/usr/bin/env python3
"""
aa_selection_vs_relaxation.py

Purpose:
    Takes a codon alignment in HyPhy PHYLIP format and a Newick tree file.
    Performs ancestral reconstruction at the amino-acid level using a
    Fitch parsimony approach (fast, alignment-driven). Infers amino-acid
    substitutions along branches, classifies substitutions by multiple
    physicochemical schemes (Grantham distance, polarity category changes,
    side-chain volume changes), computes site-wise entropy, and runs
    statistical tests comparing user-specified "test" branches vs the rest
    (reference branches) to infer signals consistent with relaxed selection
    vs positive selection.

Notes & caveats:
    - This implementation uses maximum-parsimony ancestral reconstruction
      (Fitch) on translated amino-acid sequences. It is fast and model-free
      but less powerful than likelihood-based ancestral reconstruction
      (e.g., PAML / HyPhy). Use when you want a portable, self-contained
      analysis. For the strongest inference combine with RELAX/aBSREL.

    - The script now implements several substitution classification schemes:
        1) Grantham distance (radical if > GRANTHAM_THRESHOLD)
        2) Polarity / charge category change (radical if category changes)
        3) Side-chain volume difference (radical if |delta| > VOLUME_THRESHOLD)

    - The script reports scheme-specific counts (per-branch and per-site),
      performs statistical tests per scheme (Fisher exact / Mann-Whitney U),
      and produces a combined heuristic interpretation that integrates
      entropy and scheme results.

Dependencies:
    biopython, numpy, pandas, scipy, matplotlib

Usage example:
    python aa_selection_vs_relaxation.py \
        --alignment alignment.phy --tree tree.nwk \
        --test-branches "(sp1,sp2),sp3" \
        --out-prefix analysis_out

Output:
    - CSV files with per-site metrics and per-branch substitution summaries
    - Simple plots (entropy distributions, radical fraction comparisons)
    - Interpretation text summarising scheme-specific and combined inferences

"""

import argparse
import math
import sys
from collections import defaultdict, Counter

import numpy as np
import pandas as pd
from Bio import AlignIO, Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from scipy import stats
import matplotlib.pyplot as plt

# ----------------------------- Constants & Property Tables -----------------------------
STANDARD_TABLE = CodonTable.unambiguous_dna_by_name["Standard"]

# Grantham distance matrix (complete symmetric mapping)
_GRANTHAM_MATRIX = {
    ('A','A'):0,('A','R'):112,('A','N'):111,('A','D'):126,('A','C'):195,('A','Q'):91,('A','E'):107,('A','G'):60,('A','H'):86,('A','I'):94,('A','L'):96,('A','K'):106,('A','M'):84,('A','F'):113,('A','P'):27,('A','S'):99,('A','T'):58,('A','W'):148,('A','Y'):112,('A','V'):64,
    ('R','A'):112,('R','R'):0,('R','N'):86,('R','D'):96,('R','C'):180,('R','Q'):43,('R','E'):54,('R','G'):125,('R','H'):29,('R','I'):97,('R','L'):102,('R','K'):26,('R','M'):91,('R','F'):97,('R','P'):103,('R','S'):110,('R','T'):71,('R','W'):101,('R','Y'):77,('R','V'):96,
    ('N','A'):111,('N','R'):86,('N','N'):0,('N','D'):23,('N','C'):139,('N','Q'):46,('N','E'):42,('N','G'):80,('N','H'):68,('N','I'):149,('N','L'):153,('N','K'):94,('N','M'):142,('N','F'):158,('N','P'):91,('N','S'):46,('N','T'):65,('N','W'):174,('N','Y'):143,('N','V'):133,
    ('D','A'):126,('D','R'):96,('D','N'):23,('D','D'):0,('D','C'):154,('D','Q'):61,('D','E'):45,('D','G'):94,('D','H'):81,('D','I'):168,('D','L'):172,('D','K'):101,('D','M'):160,('D','F'):177,('D','P'):108,('D','S'):65,('D','T'):85,('D','W'):181,('D','Y'):160,('D','V'):152,
    ('C','A'):195,('C','R'):180,('C','N'):139,('C','D'):154,('C','C'):0,('C','Q'):154,('C','E'):170,('C','G'):159,('C','H'):174,('C','I'):198,('C','L'):198,('C','K'):202,('C','M'):196,('C','F'):205,('C','P'):169,('C','S'):112,('C','T'):149,('C','W'):215,('C','Y'):194,('C','V'):192,
    ('Q','A'):91,('Q','R'):43,('Q','N'):46,('Q','D'):61,('Q','C'):154,('Q','Q'):0,('Q','E'):29,('Q','G'):87,('Q','H'):24,('Q','I'):109,('Q','L'):113,('Q','K'):53,('Q','M'):101,('Q','F'):116,('Q','P'):76,('Q','S'):68,('Q','T'):42,('Q','W'):130,('Q','Y'):99,('Q','V'):96,
    ('E','A'):107,('E','R'):54,('E','N'):42,('E','D'):45,('E','C'):170,('E','Q'):29,('E','E'):0,('E','G'):98,('E','H'):40,('E','I'):134,('E','L'):138,('E','K'):56,('E','M'):126,('E','F'):140,('E','P'):93,('E','S'):80,('E','T'):65,('E','W'):152,('E','Y'):122,('E','V'):121,
    ('G','A'):60,('G','R'):125,('G','N'):80,('G','D'):94,('G','C'):159,('G','Q'):87,('G','E'):98,('G','G'):0,('G','H'):98,('G','I'):135,('G','L'):138,('G','K'):127,('G','M'):127,('G','F'):153,('G','P'):42,('G','S'):56,('G','T'):59,('G','W'):184,('G','Y'):147,('G','V'):109,
    ('H','A'):86,('H','R'):29,('H','N'):68,('H','D'):81,('H','C'):174,('H','Q'):24,('H','E'):40,('H','G'):98,('H','H'):0,('H','I'):121,('H','L'):124,('H','K'):32,('H','M'):91,('H','F'):100,('H','P'):77,('H','S'):89,('H','T'):47,('H','W'):115,('H','Y'):83,('H','V'):84,
    ('I','A'):94,('I','R'):97,('I','N'):149,('I','D'):168,('I','C'):198,('I','Q'):109,('I','E'):134,('I','G'):135,('I','H'):121,('I','I'):0,('I','L'):5,('I','K'):102,('I','M'):10,('I','F'):21,('I','P'):95,('I','S'):142,('I','T'):89,('I','W'):61,('I','Y'):33,('I','V'):29,
    ('L','A'):96,('L','R'):102,('L','N'):153,('L','D'):172,('L','C'):198,('L','Q'):113,('L','E'):138,('L','G'):138,('L','H'):124,('L','I'):5,('L','L'):0,('L','K'):107,('L','M'):15,('L','F'):22,('L','P'):98,('L','S'):145,('L','T'):92,('L','W'):61,('L','Y'):36,('L','V'):32,
    ('K','A'):106,('K','R'):26,('K','N'):94,('K','D'):101,('K','C'):202,('K','Q'):53,('K','E'):56,('K','G'):127,('K','H'):32,('K','I'):102,('K','L'):107,('K','K'):0,('K','M'):95,('K','F'):102,('K','P'):103,('K','S'):121,('K','T'):78,('K','W'):110,('K','Y'):85,('K','V'):97,
    ('M','A'):84,('M','R'):91,('M','N'):142,('M','D'):160,('M','C'):196,('M','Q'):101,('M','E'):126,('M','G'):127,('M','H'):91,('M','I'):10,('M','L'):15,('M','K'):95,('M','M'):0,('M','F'):28,('M','P'):87,('M','S'):135,('M','T'):81,('M','W'):67,('M','Y'):36,('M','V'):21,
    ('F','A'):113,('F','R'):97,('F','N'):158,('F','D'):177,('F','C'):205,('F','Q'):116,('F','E'):140,('F','G'):153,('F','H'):100,('F','I'):21,('F','L'):22,('F','K'):102,('F','M'):28,('F','F'):0,('F','P'):114,('F','S'):155,('F','T'):103,('F','W'):40,('F','Y'):22,('F','V'):50,
    ('P','A'):27,('P','R'):103,('P','N'):91,('P','D'):108,('P','C'):169,('P','Q'):76,('P','E'):93,('P','G'):42,('P','H'):77,('P','I'):95,('P','L'):98,('P','K'):103,('P','M'):87,('P','F'):114,('P','P'):0,('P','S'):74,('P','T'):38,('P','W'):147,('P','Y'):110,('P','V'):68,
    ('S','A'):99,('S','R'):110,('S','N'):46,('S','D'):65,('S','C'):112,('S','Q'):68,('S','E'):80,('S','G'):56,('S','H'):89,('S','I'):142,('S','L'):145,('S','K'):121,('S','M'):135,('S','F'):155,('S','P'):74,('S','S'):0,('S','T'):58,('S','W'):177,('S','Y'):144,('S','V'):124,
    ('T','A'):58,('T','R'):71,('T','N'):65,('T','D'):85,('T','C'):149,('T','Q'):42,('T','E'):65,('T','G'):59,('T','H'):47,('T','I'):89,('T','L'):92,('T','K'):78,('T','M'):81,('T','F'):103,('T','P'):38,('T','S'):58,('T','T'):0,('T','W'):128,('T','Y'):92,('T','V'):69,
    ('W','A'):148,('W','R'):101,('W','N'):174,('W','D'):181,('W','C'):215,('W','Q'):130,('W','E'):152,('W','G'):184,('W','H'):115,('W','I'):61,('W','L'):61,('W','K'):110,('W','M'):67,('W','F'):40,('W','P'):147,('W','S'):177,('W','T'):128,('W','W'):0,('W','Y'):37,('W','V'):88,
    ('Y','A'):112,('Y','R'):77,('Y','N'):143,('Y','D'):160,('Y','C'):194,('Y','Q'):99,('Y','E'):122,('Y','G'):147,('Y','H'):83,('Y','I'):33,('Y','L'):36,('Y','K'):85,('Y','M'):36,('Y','F'):22,('Y','P'):110,('Y','S'):144,('Y','T'):92,('Y','W'):37,('Y','Y'):0,('Y','V'):64,
    ('V','A'):64,('V','R'):96,('V','N'):133,('V','D'):152,('V','C'):192,('V','Q'):96,('V','E'):121,('V','G'):109,('V','H'):84,('V','I'):29,('V','L'):32,('V','K'):97,('V','M'):21,('V','F'):50,('V','P'):68,('V','S'):124,('V','T'):69,('V','W'):88,('V','Y'):64,('V','V'):0
}

# Amino acid polarity / charge categories
# Categories: 'nonpolar', 'polar', 'positive', 'negative', 'aromatic'
AA_CATEGORY = {
    'A':'nonpolar','V':'nonpolar','L':'nonpolar','I':'nonpolar','M':'nonpolar','F':'aromatic','W':'aromatic','Y':'aromatic',
    'P':'nonpolar','G':'nonpolar',
    'S':'polar','T':'polar','C':'polar','N':'polar','Q':'polar',
    'K':'positive','R':'positive','H':'positive',
    'D':'negative','E':'negative'
}

# Side-chain volumes (approximate van der Waals volumes; units not important, only differences)
# Values are approximate but internally consistent for relative comparisons.
AA_VOLUME = {
    'A':  31.0, 'R': 124.0, 'N': 56.0, 'D': 54.0, 'C': 55.0,
    'Q': 85.0, 'E': 83.0, 'G': 3.0,  'H': 96.0, 'I': 111.0,
    'L': 111.0, 'K': 119.0, 'M': 105.0, 'F': 132.0, 'P': 32.0,
    'S': 32.0, 'T': 61.0, 'W': 170.0, 'Y': 136.0, 'V': 84.0
}

AA_LETTERS = set(AA_CATEGORY.keys())

# Thresholds (user adjustable later)
GRANTHAM_THRESHOLD = 100
VOLUME_THRESHOLD = 40.0

# ----------------------------- Utility functions -----------------------------

def grantham(a, b):
    if a == b:
        return 0
    if (a,b) in _GRANTHAM_MATRIX:
        return _GRANTHAM_MATRIX[(a,b)]
    if (b,a) in _GRANTHAM_MATRIX:
        return _GRANTHAM_MATRIX[(b,a)]
    return None

# Shannon entropy at a site given a list of residues
def shannon_entropy(residues):
    counts = Counter(residues)
    total = sum(counts.values())
    ent = 0.0
    for c in counts.values():
        p = c / total
        if p > 0:
            ent -= p * math.log2(p)
    return ent

# ----------------------------- Parsimony AR -----------------------------

def translate_alignment_codon_to_aa(aln_codon_records):
    """Translate codon sequences (aligned) to amino acid sequences.
    Expect sequences length multiple of 3; gaps ('-') are translated to '-'."""
    aa_records = []
    for rec in aln_codon_records:
        seq = str(rec.seq).upper()
        if len(seq) % 3 != 0:
            raise ValueError(f"Sequence length for {rec.id} not multiple of 3")
        aa_seq = []
        for i in range(0, len(seq), 3):
            codon = seq[i:i+3]
            if '-' in codon or 'N' in codon or 'n' in codon:
                aa_seq.append('-')
            else:
                try:
                    aa = Seq(codon).translate(table=STANDARD_TABLE, to_stop=False)
                    aa_seq.append(str(aa))
                except Exception:
                    aa_seq.append('X')
        aa_records.append(SeqRecord(Seq(''.join(aa_seq)), id=rec.id))
    return aa_records


def fitch_pass(tree, tip_aa_map, seq_length):
    """
    Perform Fitch parsimony to get sets at each node (post-order) and
    then assign a specific amino acid by backtracking (pre-order) to get
    one plausible reconstruction.

    tree: Bio.Phylo tree (Clade root)
    tip_aa_map: dict tip_name -> AA sequence (string)
    seq_length: number of amino acid sites

    Returns: dict node -> sequence (string). Node keys are clade objects.
    """
    # Ensure terminal names present
    for term in tree.get_terminals():
        if term.name not in tip_aa_map:
            raise ValueError(f"Tip {term.name} not found in alignment")

    # Initialize sets for each node
    sets = {}
    # Post-order traversal
    for node in tree.get_nonterminals(order='postorder') + tree.get_terminals():
        if node.is_terminal():
            seq = tip_aa_map[node.name]
            sets[node] = [set([aa]) if aa != '-' else set(['-']) for aa in seq]
        else:
            # children
            children = node.clades
            child_sets = [sets[ch] for ch in children]
            node_sets = []
            for i in range(seq_length):
                inter = set.intersection(*[cs[i] for cs in child_sets])
                if inter:
                    node_sets.append(set(inter))
                else:
                    union = set.union(*[cs[i] for cs in child_sets])
                    node_sets.append(set(union))
            sets[node] = node_sets

    # Backtracking to assign residues (choose arbitrary element from set, prefer consistency)
    assignments = {}
    root = tree.root
    # choose for root: pick for each site the element with highest frequency among tips if available
    root_seq = []
    for i in range(seq_length):
        s = sets[root][i]
        if '-' in s and len(s) > 1:
            s_no_gap = s - set(['-'])
            if s_no_gap:
                s = s_no_gap
        # pick most common among tips
        counts = Counter([tip_aa_map[t.name][i] for t in tree.get_terminals() if tip_aa_map[t.name][i] != '-'])
        chosen = None
        for c,_ in counts.most_common():
            if c in s:
                chosen = c
                break
        if not chosen:
            chosen = sorted(list(s))[0]
        root_seq.append(chosen)
    assignments[root] = ''.join(root_seq)

    # Pre-order traversal
    for node in tree.get_nonterminals(order='preorder'):
        parent_seq = assignments.get(node, None)
        for child in node.clades:
            if child in assignments:
                continue
            child_seq = []
            for i in range(seq_length):
                s = sets[child][i]
                p = parent_seq[i]
                if p in s:
                    child_seq.append(p)
                else:
                    child_seq.append(sorted(list(s))[0])
            assignments[child] = ''.join(child_seq)

    # terminals
    for t in tree.get_terminals():
        assignments[t] = tip_aa_map[t.name]

    return assignments

# ----------------------------- Substitution inference & classification -----------------------------

def infer_substitutions(tree, assignments):
    """Walk tree edges and list substitutions as tuples:
        (parent_node, child_node, site_index, parent_aa, child_aa)
    """
    subs = []
    # ensure we include root if present
    nodes = [tree.root] + tree.get_nonterminals(order='preorder') + tree.get_terminals()
    for parent in nodes:
        for child in getattr(parent, 'clades', []):
            parent_seq = assignments[parent]
            child_seq = assignments[child]
            for i, (pa, ca) in enumerate(zip(parent_seq, child_seq)):
                if pa == ca:
                    continue
                if pa == '-' or ca == '-':
                    # treat gaps as special; skip or record separately (skip here)
                    continue
                subs.append((parent, child, i, pa, ca))
    return subs


def classify_substitution(pa, ca):
    """Return classification dict for a substitution pa->ca under multiple schemes."""
    res = {}
    res['grantham'] = grantham(pa, ca)
    res['grantham_radical'] = (res['grantham'] is not None and res['grantham'] > GRANTHAM_THRESHOLD)

    # polarity/charge category change
    cat_pa = AA_CATEGORY.get(pa, 'other')
    cat_ca = AA_CATEGORY.get(ca, 'other')
    res['category_pa'] = cat_pa
    res['category_ca'] = cat_ca
    res['category_change'] = (cat_pa != cat_ca)

    # volume change
    vol_pa = AA_VOLUME.get(pa, None)
    vol_ca = AA_VOLUME.get(ca, None)
    if vol_pa is None or vol_ca is None:
        res['volume_diff'] = None
        res['volume_radical'] = False
    else:
        res['volume_diff'] = abs(vol_pa - vol_ca)
        res['volume_radical'] = (res['volume_diff'] > VOLUME_THRESHOLD)

    return res

# ----------------------------- CLI & Workflow -----------------------------

def parse_test_branches_arg(arg):
    parts = [p.strip() for p in arg.split(',') if p.strip()]
    return parts


def locate_nodes_for_tests(tree, parts):
    selected = set()
    for part in parts:
        if '(' in part or ')' in part:
            names = [n.strip() for n in part.replace('(', '').replace(')', '').split(',') if n.strip()]
            for node in tree.get_nonterminals(order='postorder'):
                tips = {t.name for t in node.get_terminals()}
                if set(names).issubset(tips):
                    selected.add(node)
                    break
        else:
            found = False
            for t in tree.get_terminals():
                if t.name == part:
                    selected.add(t)
                    found = True
                    break
            if not found:
                print(f"Warning: test branch spec '{part}' not found as a tip or clade")
    return selected


def run_analysis(aln_path, tree_path, test_branches_arg, out_prefix):
    # Read codon alignment
    aln = AlignIO.read(aln_path, 'phylip')
    seq_len_nt = aln.get_alignment_length()
    if seq_len_nt % 3 != 0:
        raise ValueError('Alignment length not multiple of 3 (codon alignment expected)')
    seq_len_aa = seq_len_nt // 3

    aa_records = translate_alignment_codon_to_aa(aln)
    tip_aa_map = {rec.id: str(rec.seq) for rec in aa_records}

    # Read tree
    tree = Phylo.read(tree_path, 'newick')
    try:
        _ = tree.root
    except Exception:
        tree.root_at_midpoint()

    assignments = fitch_pass(tree, tip_aa_map, seq_len_aa)
    subs = infer_substitutions(tree, assignments)

    # Collect scheme-specific counters
    branch_counters = defaultdict(lambda: defaultdict(int))
    site_subs = defaultdict(list)

    for parent, child, i, pa, ca in subs:
        cls = classify_substitution(pa, ca)
        key = f"{get_clade_label(parent)}->{get_clade_label(child)}"
        branch_counters[key]['total_subs'] += 1
        if cls['grantham_radical']:
            branch_counters[key]['grantham_radical'] += 1
        if cls['category_change']:
            branch_counters[key]['category_change'] += 1
        if cls['volume_radical']:
            branch_counters[key]['volume_radical'] += 1
        site_subs[i].append((key, pa, ca, cls))

    # Entropy per site
    site_entropy_extant = {}
    for i in range(seq_len_aa):
        residues = [seq[i] for seq in tip_aa_map.values() if seq[i] != '-']
        site_entropy_extant[i] = shannon_entropy(residues) if residues else 0.0

    # Identify test nodes
    parts = parse_test_branches_arg(test_branches_arg)
    test_nodes = locate_nodes_for_tests(tree, parts)

    def branch_is_test(child_node):
        child_tips = {t.name for t in child_node.get_terminals()}
        for tn in test_nodes:
            tn_tips = {t.name for t in tn.get_terminals()}
            if child_tips.issubset(tn_tips):
                return True
        return False

    # Helper to map branch key's child label to clade
    def child_clade_from_key(tree, key):
        # key format 'parent_label->child_label'
        child_label = key.split('->')[1]
        for cl in tree.find_clades():
            if cl.name == child_label:
                return cl
        for t in tree.get_terminals():
            if t.name == child_label:
                return t
        return None

    # Build per-site metrics
    site_records = []
    for i in range(seq_len_aa):
        subs_i = site_subs.get(i, [])
        test_subs = []
        ref_subs = []
        for s in subs_i:
            child_cl = child_clade_from_key(tree, s[0])
            if child_cl is None:
                ref_subs.append(s)
            elif branch_is_test(child_cl):
                test_subs.append(s)
            else:
                ref_subs.append(s)

        def summarize(sublist):
            n = len(sublist)
            gran = sum(1 for x in sublist if x[3]['grantham_radical'])
            cat = sum(1 for x in sublist if x[3]['category_change'])
            vol = sum(1 for x in sublist if x[3]['volume_radical'])
            return n, gran, cat, vol

        n_test, gran_test, cat_test, vol_test = summarize(test_subs)
        n_ref, gran_ref, cat_ref, vol_ref = summarize(ref_subs)
        ent = site_entropy_extant[i]
        site_records.append({'site': i+1, 'entropy': ent,
                              'n_test_subs': n_test, 'n_ref_subs': n_ref,
                              'grantham_test': gran_test, 'grantham_ref': gran_ref,
                              'category_test': cat_test, 'category_ref': cat_ref,
                              'volume_test': vol_test, 'volume_ref': vol_ref})

    site_df = pd.DataFrame(site_records)
    site_df.to_csv(f"{out_prefix}_site_metrics.csv", index=False)

    # Aggregate tests per scheme
    stats_out = []
    # Entropy comparison (sites with any test subs vs sites with any ref subs)
    sites_with_test = site_df[site_df['n_test_subs'] > 0]['entropy']
    sites_with_ref = site_df[site_df['n_ref_subs'] > 0]['entropy']
    if len(sites_with_test) > 0 and len(sites_with_ref) > 0:
        u_stat, p_entropy = stats.mannwhitneyu(sites_with_test, sites_with_ref, alternative='two-sided')
        stats_out.append({'scheme':'entropy', 'test':'sites_test_vs_ref', 'stat':u_stat, 'p_value':p_entropy})
    else:
        stats_out.append({'scheme':'entropy', 'test':'sites_test_vs_ref', 'stat':None, 'p_value':None})

    # For each substitution scheme, build contingency table radical/nonradical in test vs ref and run Fisher's exact
    def contingency_and_fisher(rad_test, nonrad_test, rad_ref, nonrad_ref):
        tot = rad_test + nonrad_test + rad_ref + nonrad_ref
        if tot == 0:
            return None, None
        contingency = np.array([[rad_test, nonrad_test], [rad_ref, nonrad_ref]])
        try:
            odds, p = stats.fisher_exact(contingency)
            return odds, p
        except Exception:
            return None, None

    # Grantham
    total_rad_test = site_df['grantham_test'].sum()
    total_sub_test = site_df['n_test_subs'].sum()
    total_nonrad_test = total_sub_test - total_rad_test
    total_rad_ref = site_df['grantham_ref'].sum()
    total_sub_ref = site_df['n_ref_subs'].sum()
    total_nonrad_ref = total_sub_ref - total_rad_ref
    odds_gr, p_gr = contingency_and_fisher(total_rad_test, total_nonrad_test, total_rad_ref, total_nonrad_ref)
    stats_out.append({'scheme':'grantham', 'test':'radical_fraction_test_vs_ref', 'oddsratio':odds_gr, 'p_value':p_gr})

    # Category (polarity/charge)
    total_cat_rad_test = site_df['category_test'].sum()
    total_cat_nonrad_test = total_sub_test - total_cat_rad_test
    total_cat_rad_ref = site_df['category_ref'].sum()
    total_cat_nonrad_ref = total_sub_ref - total_cat_rad_ref
    odds_cat, p_cat = contingency_and_fisher(total_cat_rad_test, total_cat_nonrad_test, total_cat_rad_ref, total_cat_nonrad_ref)
    stats_out.append({'scheme':'category_change', 'test':'radical_fraction_test_vs_ref', 'oddsratio':odds_cat, 'p_value':p_cat})

    # Volume
    total_vol_rad_test = site_df['volume_test'].sum()
    total_vol_nonrad_test = total_sub_test - total_vol_rad_test
    total_vol_rad_ref = site_df['volume_ref'].sum()
    total_vol_nonrad_ref = total_sub_ref - total_vol_rad_ref
    odds_vol, p_vol = contingency_and_fisher(total_vol_rad_test, total_vol_nonrad_test, total_vol_rad_ref, total_vol_nonrad_ref)
    stats_out.append({'scheme':'volume_change', 'test':'radical_fraction_test_vs_ref', 'oddsratio':odds_vol, 'p_value':p_vol})

    stats_df = pd.DataFrame(stats_out)
    stats_df.to_csv(f"{out_prefix}_stats_summary.csv", index=False)

    # Plots: entropy histogram and barplots of radical fractions per scheme
    plt.figure()
    plt.hist(site_df['entropy'].dropna(), bins=30)
    plt.title('Site entropy (extant sequences)')
    plt.xlabel('Shannon entropy (bits)')
    plt.ylabel('Count')
    plt.savefig(f"{out_prefix}_entropy_hist.png")
    plt.close()

    schemes = ['grantham','category_change','volume_change']
    rad_fracs_test = []
    rad_fracs_ref = []
    pvals = []
    for s in schemes:
        if s == 'grantham':
            rt = total_rad_test; st = total_sub_test; rf = total_rad_ref; sf = total_sub_ref
        elif s == 'category_change':
            rt = total_cat_rad_test; st = total_sub_test; rf = total_cat_rad_ref; sf = total_sub_ref
        else:
            rt = total_vol_rad_test; st = total_sub_test; rf = total_vol_rad_ref; sf = total_sub_ref
        frac_t = (rt / st) if st > 0 else 0
        frac_r = (rf / sf) if sf > 0 else 0
        rad_fracs_test.append(frac_t)
        rad_fracs_ref.append(frac_r)
        # p-value lookup
        pval = None
        if s == 'grantham': pval = p_gr
        if s == 'category_change': pval = p_cat
        if s == 'volume_change': pval = p_vol
        pvals.append(pval)

    x = np.arange(len(schemes))
    width = 0.35
    plt.figure(figsize=(8,4))
    plt.bar(x - width/2, rad_fracs_test, width, label='test')
    plt.bar(x + width/2, rad_fracs_ref, width, label='ref')
    plt.xticks(x, schemes)
    plt.ylabel('Fraction radical substitutions')
    plt.legend()
    plt.title('Radical substitution fractions by scheme')
    plt.savefig(f"{out_prefix}_radical_fractions_by_scheme.png")
    plt.close()

    # ----------------------------- Combined heuristic interpretation -----------------------------
    interpretation = []
    # Heuristic signals
    ent_p = stats_df.loc[stats_df['scheme']=='entropy','p_value'].values[0]
    ent_direction = None
    if ent_p is not None and ent_p < 0.05:
        mean_test_ent = sites_with_test.mean() if len(sites_with_test)>0 else 0
        mean_ref_ent = sites_with_ref.mean() if len(sites_with_ref)>0 else 0
        if mean_test_ent > mean_ref_ent:
            ent_direction = 'higher_in_test'
        else:
            ent_direction = 'higher_in_ref'

    # scheme significance
    sig_schemes = []
    for row in stats_out:
        if row.get('p_value') is not None and row.get('p_value') < 0.05:
            sig_schemes.append(row['scheme'])

    # Combined rules (heuristic):
    # - Positive selection signature: elevated fraction of radical substitutions in TEST across one or more schemes (significant), AND not accompanied by a broad increase in entropy alone. Also directional signal (excess radical substitutions) supports adaptation.
    # - Relaxation signature: increased entropy in TEST (significant) with no increase in radical substitution fractions (or even decrease), and overall increase in conservative substitutions.
    # - Mixed / uncertain: both entropy increase and radical increase, or neither significant.

    pos_evidence = any([s in ['grantham','category_change','volume_change'] for s in sig_schemes]) and not (ent_direction == 'higher_in_test' and len(sig_schemes)==0)
    relax_evidence = (ent_direction == 'higher_in_test') and not any([s in ['grantham','category_change','volume_change'] for s in sig_schemes])

    if pos_evidence:
        interpretation.append('Evidence consistent with positive selection: at least one substitution scheme shows a significantly higher fraction of radical substitutions on test branches compared to reference.')
    if relax_evidence:
        interpretation.append('Evidence consistent with relaxed selection: site entropy is significantly higher on test branches and radical substitution fractions are not elevated.')
    if not interpretation:
        interpretation.append('No clear signal from these heuristics. Results are ambiguous: consider combining with codon-model tests (RELAX, aBSREL) and likelihood-based ancestral reconstruction for stronger inference.')

    # Detailed notes to include
    interpretation.append('
Detailed results:')
    interpretation.append(f'Grantham radical fraction test p-value: {p_gr}')
    interpretation.append(f'Category-change (polarity/charge) p-value: {p_cat}')
    interpretation.append(f'Volume-change p-value: {p_vol}')
    interpretation.append(f'Entropy sites p-value (Mann-Whitney U): {ent_p}')

    with open(f"{out_prefix}_interpretation.txt", 'w') as fh:
        fh.write('
'.join(interpretation) + '
')

    # Save branch-level summary
    branch_rows = []
    for b, cnts in branch_counters.items():
        branch_rows.append({'branch': b, **cnts})
    branch_df = pd.DataFrame(branch_rows)
    branch_df.to_csv(f"{out_prefix}_branch_summary.csv", index=False)

    print(f"Wrote site metrics to {out_prefix}_site_metrics.csv")
    print(f"Wrote stats summary to {out_prefix}_stats_summary.csv")
    print(f"Wrote branch summary to {out_prefix}_branch_summary.csv")
    print(f"Wrote plots to {out_prefix}_entropy_hist.png and {out_prefix}_radical_fractions_by_scheme.png")
    print('Interpretation written to {}_interpretation.txt'.format(out_prefix))

# ----------------------------- Small helpers -----------------------------

def get_clade_label(clade):
    return clade.name if getattr(clade, 'name', None) else f"Node_{id(clade)}"

# ----------------------------- Main -----------------------------

def main():
    p = argparse.ArgumentParser(description='Infer relaxed vs positive selection signals from amino-acid substitutions and entropy.')
    p.add_argument('--alignment', required=True, help='Codon alignment in PHYLIP format (HyPhy phylip accepted)')
    p.add_argument('--tree', required=True, help='Tree in Newick format')
    p.add_argument('--test-branches', required=True, help='Comma-separated list of tip names or clade specs for test branches')
    p.add_argument('--out-prefix', default='aa_sel', help='Output file prefix')
    p.add_argument('--grantham-threshold', type=float, default=GRANTHAM_THRESHOLD, help='Grantham distance threshold to call radical (default 100)')
    p.add_argument('--volume-threshold', type=float, default=VOLUME_THRESHOLD, help='Volume difference threshold to call radical (default 40)')
    args = p.parse_args()

    global GRANTHAM_THRESHOLD, VOLUME_THRESHOLD
    GRANTHAM_THRESHOLD = args.grantham_threshold
    VOLUME_THRESHOLD = args.volume_threshold

    run_analysis(args.alignment, args.tree, args.test_branches, args.out_prefix)

if __name__ == '__main__':
    main()
