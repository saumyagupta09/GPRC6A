#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Branch-wise multi-hit-by-codon analysis via codon parsimony (no external deps).

What it does
------------
Given a codon alignment and a tree with branch lengths, this script:
  1) Reconstructs most-parsimonious codon states at every node (Sankoff DP).
  2) For each branch (parent -> child) and each codon site, computes the
     nucleotide Hamming distance between parent and child codons (0..3).
  3) Flags "multi-hit" codons where the distance >= 2.
  4) Summarizes per-branch counts and normalizes by branch length.

Normalization
-------------
Let L be branch length in substitutions per nucleotide site (as is standard for Newick).
Let N be the number of codon sites in the alignment.

We report:
  - MH_per_unit_length      = (# codons with dist >= 2) / L
  - MH_per_expected_nt_subs = (# codons with dist >= 2) / (3 * N * L)

The second one normalizes by the *expected number of nucleotide substitutions*
across the whole alignment on that branch (3 nucleotides per codon).

Inputs
------
- Alignment: FASTA or sequential PHYLIP containing codon sequences (length multiple of 3)
- Tree: Newick with branch lengths (required). Tip labels should match alignment names.
        Optional annotations in braces (e.g., Taxon{TEST}) are stripped by default.

Outputs (written next to the script / cwd)
------------------------------------------
- per_site_per_branch.tsv   : columns: branch, parent, child, child_is_tip, site, nt_distance (0..3), is_multi_hit
- per_branch_summary.tsv    : per-branch totals and length-normalized rates
- overall_summary.txt       : quick text summary

Usage
-----
  python3 mh_by_branch_codons.py --alignment codon_alignment.phy --tree Sus_scrofa_labelled.nwk

Optional flags:
  --no-strip-braces      (keep {...} in tip labels; default is to strip)
  --output-prefix PREFIX (default: 'mh_by_branch')

"""

import sys, re, math, argparse
from pathlib import Path
from typing import List, Dict, Tuple, Optional

# ----------------- IO helpers -----------------

def read_fasta(path: Path) -> Tuple[List[str], Dict[str,str]]:
    names, seqs = [], {}
    name, buf = None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                seqs[name] = "".join(buf).upper().replace(".", "-")
                names.append(name)
                buf = []
            name = line[1:].strip()
        else:
            buf.append(re.sub(r"\s+", "", line))
    if name is not None:
        seqs[name] = "".join(buf).upper().replace(".", "-")
        names.append(name)
    return names, seqs

def read_phylip_sequential(path: Path) -> Tuple[List[str], Dict[str,str]]:
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError("Empty PHYLIP file")
    m = re.match(r"\s*(\d+)\s+(\d+)", lines[0].strip())
    if not m:
        raise ValueError(f"Bad PHYLIP header: {lines[0]}")
    ntax, nchar = map(int, m.groups())
    i = 1
    names, seqs = [], {}
    seq_re = re.compile(r"^[ACGTURYKMSWBDHVNacgturykmswbdhvn\-\.\s]+$")
    while i < len(lines) and len(names) < ntax:
        while i < len(lines) and lines[i].strip() == "":
            i += 1
        if i >= len(lines): break
        name = lines[i].strip(); i += 1
        chunks = []
        while i < len(lines) and lines[i].strip() != "" and seq_re.match(lines[i]):
            chunks.append(re.sub(r"\s+", "", lines[i])); i += 1
        seq = "".join(chunks).upper().replace(".", "-")
        if len(seq) != nchar:
            raise ValueError(f"{name}: length {len(seq)} != nchar {nchar}")
        names.append(name); seqs[name] = seq
    return names, seqs

def detect_and_read_alignment(path: Path) -> Tuple[List[str], Dict[str,str]]:
    first = path.read_text(encoding="utf-8", errors="ignore")[:1]
    if first == ">":
        return read_fasta(path)
    else:
        return read_phylip_sequential(path)

# ----------------- Newick (with branch lengths) -----------------

class Node:
    __slots__=("name","children","parent","length")
    def __init__(self, name: Optional[str]=None, length: float=0.0):
        self.name = name
        self.children: List["Node"] = []
        self.parent: Optional["Node"] = None
        self.length = float(length)
    def is_leaf(self): return len(self.children)==0

def tokenize_newick(s: str):
    s = s.strip()
    if not s.endswith(";"): s += ";"
    i = 0
    WHITES = " \t\r\n"
    while i < len(s):
        c = s[i]
        if c in "(),:;":
            yield (c, c); i += 1
        elif c == "'":
            i += 1; buf=[]
            while i < len(s) and s[i] != "'":
                buf.append(s[i]); i += 1
            i += 1
            yield ("LABEL", "".join(buf))
        elif c in WHITES:
            i += 1
        else:
            buf=[]
            while i < len(s) and s[i] not in "(),:; \t\r\n":
                buf.append(s[i]); i += 1
            yield ("LABEL", "".join(buf))

def parse_newick(s: str) -> Node:
    toks = list(tokenize_newick(s))
    i = 0
    def parse_subtree() -> Node:
        nonlocal i
        if toks[i][0] == "(":
            i += 1
            n = Node()
            while True:
                child = parse_subtree()
                child.parent = n
                n.children.append(child)
                if toks[i][0] == ",":
                    i += 1; continue
                elif toks[i][0] == ")":
                    i += 1; break
                else:
                    raise ValueError(f"Unexpected token {toks[i]} inside group")
            # optional internal node label
            if toks[i][0] == "LABEL":
                n.name = toks[i][1]; i += 1
            # optional length
            if toks[i][0] == ":":
                i += 1
                if toks[i][0] == "LABEL":
                    try: n.length = float(toks[i][1])
                    except: n.length = 0.0
                    i += 1
            return n
        else:
            if toks[i][0] != "LABEL":
                raise ValueError("Expected leaf label")
            label = toks[i][1]; i += 1
            length = 0.0
            if toks[i][0] == ":":
                i += 1
                if toks[i][0] == "LABEL":
                    try: length = float(toks[i][1])
                    except: length = 0.0
                    i += 1
            return Node(name=label, length=length)
    root = parse_subtree()
    if toks[i][0] != ";":
        raise ValueError("Expected ';' at end of Newick")
    return root

def leaves_of(n: Node) -> List[Node]:
    out=[]
    def rec(x):
        if x.is_leaf(): out.append(x)
        else:
            for c in x.children: rec(c)
    rec(n); return out

def edges_with_lengths(root: Node) -> List[Tuple[Node,Node,float]]:
    edges=[]
    def rec(n):
        for c in n.children:
            edges.append((n,c,c.length))
            rec(c)
    rec(root); return edges

def postorder(n: Node) -> List[Node]:
    out=[]
    def rec(x):
        for c in x.children: rec(c)
        out.append(x)
    rec(n); return out

# ----------------- Codon helpers -----------------

def split_codons(seq: str) -> List[str]:
    return [seq[i:i+3] for i in range(0, len(seq), 3)]

def hamming3(a: str, b: str) -> int:
    return (a[0]!=b[0]) + (a[1]!=b[1]) + (a[2]!=b[2])

# ----------------- Sankoff (codon parsimony) -----------------

def parsimonious_edge_changes_for_site(root: Node,
                                       leaves: List[Node],
                                       codon_align: Dict[str,List[str]],
                                       site_idx: int
                                       ) -> Optional[Tuple[Dict[Tuple[Node,Node],int], Dict[Tuple[Node,Node],Tuple[str,str]]]]:
    """Return (changes, codon_on_edge) or None if all states unknown at this site."""
    # Candidate state set = observed leaf codons at this site (ignoring missing)
    states = set()
    for lf in leaves:
        cd = codon_align[lf.name][site_idx]
        if "-" in cd or "N" in cd or "?" in cd or len(cd)!=3:
            continue
        states.add(cd)
    if not states:
        return None
    states = list(states)

    # DP
    post_nodes = postorder(root)
    cost: Dict[Node, Dict[str,float]] = {}
    for n in post_nodes:
        if n.is_leaf():
            cd = codon_align[n.name][site_idx]
            if "-" in cd or "N" in cd or "?" in cd or len(cd)!=3:
                cost[n] = {s:0.0 for s in states}
            else:
                cost[n] = {s:(0.0 if s==cd else math.inf) for s in states}
        else:
            d = {}
            for s in states:
                total = 0.0
                for c in n.children:
                    cc = cost[c]
                    best = math.inf
                    for t in states:
                        v = cc[t] + hamming3(s,t)
                        if v < best: best = v
                    total += best
                d[s] = total
            cost[n] = d

    # Choose root state (min cost)
    root_state = min(cost[root], key=lambda s: cost[root][s])

    # Greedy pre-order to select child states and compute edge distances
    changes: Dict[Tuple[Node,Node],int] = {}
    codon_on_edge: Dict[Tuple[Node,Node],Tuple[str,str]] = {}

    def choose_child_state(parent_state: str, child: Node) -> str:
        cc = cost[child]
        best_t = None; best = math.inf
        for t in states:
            v = cc[t] + hamming3(parent_state, t)
            if v < best:
                best = v; best_t = t
        return best_t

    def pre(n: Node, parent_state: str):
        for c in n.children:
            child_state = choose_child_state(parent_state, c)
            changes[(n,c)] = hamming3(parent_state, child_state)
            codon_on_edge[(n,c)] = (parent_state, child_state)
            pre(c, child_state)

    pre(root, root_state)
    return changes, codon_on_edge

# ----------------- Main -----------------

def main():
    ap = argparse.ArgumentParser(description="Branch-wise multi-hit-by-codon via parsimony")
    ap.add_argument("--alignment", required=True, help="FASTA or sequential PHYLIP codon alignment")
    ap.add_argument("--tree", required=True, help="Newick tree with branch lengths")
    ap.add_argument("--no-strip-braces", action="store_true",
                    help="Do NOT strip {...} from tip labels (default: strip)")
    ap.add_argument("--output-prefix", default="mh_by_branch", help="Output file prefix")
    args = ap.parse_args()

    aln_path = Path(args.alignment)
    tree_path = Path(args.tree)
    strip_braces = not args.no_strip_braces
    prefix = args.output_prefix

    # read alignment
    names, seqs = detect_and_read_alignment(aln_path)
    # sanity: lengths multiple of 3
    L = len(next(iter(seqs.values())))
    if L % 3 != 0:
        raise SystemExit(f"Alignment length {L} is not a multiple of 3")
    n_codons = L // 3

    # read tree
    root = parse_newick(tree_path.read_text())
    # clean tip labels to match alignment names
    for lf in leaves_of(root):
        if lf.name and strip_braces:
            lf.name = re.sub(r"\{.*?\}", "", lf.name)

    # ensure taxa match
    leaf_names = {lf.name for lf in leaves_of(root)}
    missing = set(names) - leaf_names
    extra   = leaf_names - set(names)
    if missing:
        raise SystemExit(f"[Error] Tree is missing taxa present in alignment: {sorted(missing)}")
    if extra:
        # harmless but warn (e.g., tree has extra tips)
        sys.stderr.write(f"[Warning] Tree has extra tips not in alignment: {sorted(extra)}\n")

    # order leaves to match input order
    codon_align = {t: split_codons(seqs[t]) for t in names}
    leaves = leaves_of(root)

    # precompute edge list with lengths and names
    edges = edges_with_lengths(root)
    edge_names = []
    for p,c,blen in edges:
        pname = p.name if p.name else "INTERNAL"
        cname = c.name if c.name else "INTERNAL"
        edge_names.append((pname, cname, blen))

    # per-branch counters
    per_branch_mh = [0] * len(edges)       # # codon sites with dist >= 2
    per_branch_sh = [0] * len(edges)       # # codon sites with dist == 1
    per_branch_sumdist = [0] * len(edges)  # sum of distances (0..3) across sites

    # open per-site-per-branch writer
    site_file = Path(f"{prefix}.per_site_per_branch.tsv").open("w")
    site_file.write("branch\tparent\tchild\tchild_is_tip\tbranch_length\tsite\tnt_distance\tis_multi_hit\n")

    # main loop: sites
    for i in range(n_codons):
        res = parsimonious_edge_changes_for_site(root, leaves, codon_align, i)
        if res is None:
            # All missing at site; write zeros
            for idx,(p,c,blen) in enumerate(edges):
                child_is_tip = int(c.is_leaf())
                site_file.write(f"{idx}\t{p.name or 'INTERNAL'}\t{c.name or 'INTERNAL'}\t{child_is_tip}\t{blen:.6g}\t{i+1}\t-1\t0\n")
            continue
        changes, _ = res
        # edges are in the same (p,c) identity as 'edges' list
        for idx,(p,c,blen) in enumerate(edges):
            dist = changes.get((p,c), 0)
            is_mh = int(dist >= 2)
            is_sh = int(dist == 1)
            per_branch_mh[idx] += is_mh
            per_branch_sh[idx] += is_sh
            per_branch_sumdist[idx] += dist
            child_is_tip = int(c.is_leaf())
            site_file.write(f"{idx}\t{p.name or 'INTERNAL'}\t{c.name or 'INTERNAL'}\t{child_is_tip}\t{blen:.6g}\t{i+1}\t{dist}\t{is_mh}\n")

    site_file.close()

    # write per-branch summary with branch-length corrections
    N = n_codons
    br_path = Path(f"{prefix}.per_branch_summary.tsv")
    with br_path.open("w") as out:
        out.write("branch\tparent\tchild\tchild_is_tip\tbranch_length\tcodon_sites\t#MH(>=2)\t#SH(=1)\tSumDist(0..3)\tMH_per_unit_length\tMH_per_expected_nt_subs\n")
        for idx,(p,c,blen) in enumerate(edges):
            mh = per_branch_mh[idx]
            sh = per_branch_sh[idx]
            sd = per_branch_sumdist[idx]
            child_is_tip = int(c.is_leaf())
            # length-normalized rates
            if blen > 0.0:
                mh_per_len = mh / blen
                mh_per_exp_nt = mh / (3.0 * N * blen)
            else:
                mh_per_len = float('inf') if mh>0 else 0.0
                mh_per_exp_nt = float('inf') if mh>0 else 0.0
            out.write(f"{idx}\t{p.name or 'INTERNAL'}\t{c.name or 'INTERNAL'}\t{child_is_tip}\t{blen:.6g}\t{N}\t{mh}\t{sh}\t{sd}\t{mh_per_len:.6g}\t{mh_per_exp_nt:.6g}\n")

    # overall quick summary
    with Path(f"{prefix}.overall_summary.txt").open("w") as fh:
        total_edges = len(edges)
        total_mh = sum(per_branch_mh)
        fh.write(f"Taxa: {len(names)}\n")
        fh.write(f"Codon sites: {N}\n")
        fh.write(f"Branches: {total_edges}\n")
        fh.write(f"Multi-hit (>=2 nt) events counted as codons branches: {total_mh}\n")
        # Top 5 edges by MH
        top = sorted(
            [(i, per_branch_mh[i], edge_names[i]) for i in range(total_edges)],
            key=lambda x: -x[1]
        )[:5]
        fh.write("\nTop 5 branches by # multi-hit codons:\n")
        for i, mh, (pname, cname, blen) in top:
            fh.write(f"  idx={i}  {pname}->{cname}  blen={blen:.6g}  MH_codons={mh}\n")
        fh.write("\nSee per-branch details: " + str(br_path) + "\n")
        fh.write("Per-site branch table: " + f"{prefix}.per_site_per_branch.tsv" + "\n")

    print("[Done]")
    print(" - per-site table :", f"{prefix}.per_site_per_branch.tsv")
    print(" - per-branch sum :", f"{prefix}.per_branch_summary.tsv")
    print(" - summary        :", f"{prefix}.overall_summary.txt")

if __name__ == "__main__":
    main()
