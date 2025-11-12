#!/usr/bin/env python3
"""
Split a Newick tree + PHYLIP alignment into n phylogenetically contiguous groups
without introducing new branches or sequences. Preserves labels like "{fg}".

Features
- Strict PHYLIP reader (sequential first, interleaved fallback) honoring header ntaxa & seqlen.
- No new branches: groups are formed by cutting existing edges only (true subtrees).
- Taxa lists per group and PHYLIP-relaxed subset alignments.
- Robust ID matching: --id_match_mode exact|relaxed and optional --name_map CSV.

Requires: biopython
"""

import argparse, os, re, copy, csv
from collections import namedtuple
from typing import List, Tuple, Set, Dict, Optional

from Bio import Phylo, AlignIO
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Group = namedtuple("Group", ["tips"])

# ============ PHYLIP readers (strict; preserve labels) ============

_SEQ_CHARS_RE = re.compile(r"[A-Za-z\-\.\?\*]+")

def _clean_chunk(text: str) -> str:
    return "".join(_SEQ_CHARS_RE.findall(text or ""))

def _read_header(path: str) -> Tuple[int, int, List[str]]:
    with open(path, "r", encoding="utf-8") as fh:
        header = ""
        while True:
            header = fh.readline()
            if header == "": raise ValueError("Empty file or missing PHYLIP header.")
            if header.strip(): break
        try:
            ntaxa, seqlen = [int(x) for x in header.strip().split()[:2]]
        except Exception as e:
            raise ValueError(f"Failed to parse PHYLIP header: {header!r}") from e
        lines = [ln.rstrip("\n") for ln in fh]
    return ntaxa, seqlen, lines

def _split_id_and_chunk_preserve_labels(line: str) -> Tuple[str, str]:
    s = (line or "")
    if not s.strip(): return "", ""
    if s[0].isspace():  # chunk-only
        return "", _clean_chunk(s)
    m = re.search(r"\s([A-Za-z\-\.\?\*]+)", s)  # first seq block after whitespace
    if m:
        id_part = s[:m.start()].rstrip()
        seq_part = s[m.start():]
        return id_part, _clean_chunk(seq_part)
    return s.strip(), ""

def read_phylip_sequential_strict(path: str, trim_names: bool=False, allow_pad: bool=False) -> Tuple[MultipleSeqAlignment, int, int]:
    ntaxa, seqlen, lines = _read_header(path)
    records: List[SeqRecord] = []
    idx = 0
    def next_nonempty():
        nonlocal idx
        while idx < len(lines) and not lines[idx].strip(): idx += 1
        if idx >= len(lines): return None
        s = lines[idx]; idx += 1
        return s
    for _ in range(ntaxa):
        name_line = next_nonempty()
        if name_line is None: raise ValueError("Unexpected EOF while reading a taxon name.")
        taxon = name_line.strip() if not trim_names else name_line.strip()
        chunks, cur = [], 0
        while cur < seqlen:
            ln = next_nonempty()
            if ln is None: break
            chunk = _clean_chunk(ln)
            if not chunk:
                idx -= 1  # likely next name line
                break
            take = min(seqlen - cur, len(chunk))
            chunks.append(chunk[:take]); cur += take
        seq = "".join(chunks)
        if len(seq) != seqlen:
            if allow_pad and len(seq) < seqlen: seq += "-" * (seqlen - len(seq))
            else: raise ValueError(f"Strict PHYLIP: taxon '{taxon}' length {len(seq)} != header seqlen {seqlen}.")
        records.append(SeqRecord(Seq(seq), id=taxon, description=""))
    return MultipleSeqAlignment(records), ntaxa, seqlen

def read_phylip_interleaved_relaxed_strict(path: str, trim_names: bool=False, allow_pad: bool=False) -> Tuple[MultipleSeqAlignment, int, int]:
    ntaxa, seqlen, raw = _read_header(path)
    ids: List[str] = []; seqs: Dict[str, List[str]] = {}; idx = 0
    while idx < len(raw) and len(ids) < ntaxa:
        line = raw[idx]; idx += 1
        if not line.strip(): continue
        id_tok, chunk = _split_id_and_chunk_preserve_labels(line)
        if not id_tok: continue
        taxon = id_tok.strip() if not trim_names else id_tok.strip()
        ids.append(taxon); seqs[taxon] = [chunk]
    if len(ids) != ntaxa:
        raise ValueError(f"Interleaved parse error: expected {ntaxa} ids in first block, saw {len(ids)}.")
    def lens(): return [sum(len(x) for x in seqs[t]) for t in ids]
    while max(lens()) < seqlen and idx < len(raw):
        while idx < len(raw) and not raw[idx].strip(): idx += 1
        if idx >= len(raw): break
        peek = raw[idx]
        if peek and not peek[0].isspace():
            while idx < len(raw) and raw[idx].strip():
                line = raw[idx]; idx += 1
                id_tok, chunk = _split_id_and_chunk_preserve_labels(line)
                if not id_tok: continue
                t = id_tok.strip() if not trim_names else id_tok.strip()
                if t in seqs and chunk and sum(len(x) for x in seqs[t]) < seqlen:
                    seqs[t].append(chunk)
        else:
            block = []
            k = 0
            while idx < len(raw) and k < ntaxa:
                line = raw[idx]; idx += 1
                if not line.strip(): continue
                _, chunk = _split_id_and_chunk_preserve_labels(line)
                block.append(chunk); k += 1
            for t, ch in zip(ids, block):
                if sum(len(x) for x in seqs[t]) < seqlen and ch: seqs[t].append(ch)
    recs: List[SeqRecord] = []
    for t in ids:
        s = "".join(seqs[t])
        if len(s) != seqlen:
            if allow_pad and len(s) < seqlen: s += "-" * (seqlen - len(s))
            else: raise ValueError(f"Strict PHYLIP: taxon '{t}' length {len(s)} != header seqlen {seqlen}.")
        recs.append(SeqRecord(Seq(s), id=t, description=""))
    return MultipleSeqAlignment(recs), ntaxa, seqlen

def read_phylip_auto_strict(path: str, trim_names: bool=False, allow_pad: bool=False) -> Tuple[MultipleSeqAlignment, int, int]:
    try:
        return read_phylip_sequential_strict(path, trim_names=trim_names, allow_pad=allow_pad)
    except Exception:
        return read_phylip_interleaved_relaxed_strict(path, trim_names=trim_names, allow_pad=allow_pad)

# ============ ID matching (handles '{fg}' spacing etc.) ============

def _norm_variants(name: str) -> List[str]:
    variants = set()
    s = name
    variants.add(s)
    m = re.search(r"\s*\{[^}]+\}\s*$", s)
    if m:
        label = s[m.start():]
        base = s[:m.start()]
        variants.add(base.rstrip())                       # strip label
        variants.add(base + label.replace(" ", ""))      # tight
        variants.add(base.rstrip() + " " + label.strip())# spaced
    variants.add(s.replace("_", " "))
    variants.add(s.replace(" ", "_"))
    variants |= {v.lower() for v in list(variants)}
    variants |= {v.strip() for v in list(variants)}
    return [v for v in variants if v]

def load_name_map(path: Optional[str]) -> Dict[str, str]:
    if not path: return {}
    mapping = {}
    with open(path, "r", encoding="utf-8") as fh:
        rdr = csv.DictReader(fh)
        if "tree_id" not in rdr.fieldnames or "aln_id" not in rdr.fieldnames:
            raise SystemExit("ERROR: --name_map must be a CSV with headers: tree_id,aln_id")
        for row in rdr:
            t = row["tree_id"].strip(); a = row["aln_id"].strip()
            if t and a: mapping[t] = a
    return mapping

def build_id_mapping(tree_tips: Set[str], aln_ids: Set[str],
                     mode: str = "exact",
                     name_map: Optional[Dict[str,str]] = None) -> Dict[str,str]:
    name_map = name_map or {}
    mapping: Dict[str,str] = {}
    aln_ci = {a.lower(): a for a in aln_ids}
    for tip in tree_tips:
        if tip in name_map:
            a = name_map[tip]
            if a not in aln_ids:
                raise SystemExit(f"ERROR: name_map maps '{tip}' -> '{a}', but that aln_id is absent.")
            mapping[tip] = a; continue
        if tip in aln_ids:
            mapping[tip] = tip; continue
        if mode == "exact":
            continue
        # relaxed search
        candidates = set()
        for v in _norm_variants(tip):
            if v in aln_ids: candidates.add(v)
            if v.lower() in aln_ci: candidates.add(aln_ci[v.lower()])
        if len(candidates) == 1:
            mapping[tip] = next(iter(candidates))
        elif len(candidates) > 1:
            raise SystemExit(f"ERROR: Relaxed match for '{tip}' is ambiguous: {sorted(candidates)}. Use --name_map.")
    return mapping

# ============ Tree helpers (no new branches) ============

def get_all_tips(clade: Clade) -> List[str]:
    return [t.name for t in clade.get_terminals()]

def build_subtree_tipsets(tree: Tree) -> Dict[Clade, Set[str]]:
    tipset: Dict[Clade, Set[str]] = {}
    def dfs(node: Clade) -> Set[str]:
        if not node.clades:
            s = {node.name}; tipset[node] = s; return s
        s = set()
        for c in node.clades: s |= dfs(c)
        tipset[node] = s; return s
    dfs(tree.root); return tipset

def cut_edges_greedily(tree: Tree, n: int) -> List[Group]:
    """
    Partition tips into n components by cutting existing edges only.
    Greedy 'peel-off' heuristic with safe numeric tie-breaking.
    """
    all_tips = set(get_all_tips(tree.root))
    groups: List[Group] = [Group(all_tips)]
    if n <= 1 or len(all_tips) <= 1:
        return groups

    node_tipset = build_subtree_tipsets(tree)

    while len(groups) < n:
        # Pick largest peelable component
        groups.sort(key=lambda g: len(g.tips), reverse=True)
        G = next((g for g in groups if len(g.tips) >= 2), None)
        if G is None:
            break

        remaining_groups = n - (len(groups) - 1)
        target = max(1, round(len(G.tips) / 2)) if remaining_groups <= 2 else max(1, round(len(G.tips) / remaining_groups))

        best_key = None   # purely numeric key for comparison
        best_node = None  # the clade to peel

        for node, tips in node_tipset.items():
            if tips and tips.issubset(G.tips) and len(tips) < len(G.tips):
                size = len(tips)
                # Key: (imbalance from target, then smaller size to avoid huge peels when tied)
                cand_key = (abs(size - target), size)
                if (best_key is None) or (cand_key < best_key):
                    best_key = cand_key
                    best_node = node

        if best_node is None:
            break

        peeled = node_tipset[best_node]
        groups.remove(G)
        remainder = G.tips - peeled
        groups.append(Group(remainder))
        groups.append(Group(peeled))

    return groups[:n]


def write_group_subtree(original_tree: Tree, tipset: Set[str], out_path: str):
    t = copy.deepcopy(original_tree)
    keep = set(tipset)
    for tip in list(get_all_tips(t.root)):
        if tip not in keep:
            try: t.prune(t.find_any(name=tip))
            except Exception: pass
    Phylo.write(t, out_path, "newick")

def subset_alignment_by_mapping(aln: MultipleSeqAlignment, keep_tree_names: List[str], t2a: Dict[str,str]) -> MultipleSeqAlignment:
    keep_aln_ids = {t2a.get(t, t) for t in keep_tree_names}
    return MultipleSeqAlignment([rec for rec in aln if rec.id in keep_aln_ids])

# ============ CLI ============

def parse_args():
    p = argparse.ArgumentParser(description="Split tree+alignment into n groups (preserve labels like '{fg}', no new branches).")
    p.add_argument("--tree", required=True)
    p.add_argument("--aln", required=True)
    p.add_argument("--n_groups", type=int, required=True)
    p.add_argument("--outdir", default="splits")
    p.add_argument("--allow_extra_tree_tips", action="store_true")
    p.add_argument("--allow_extra_aln_seqs", action="store_true")
    p.add_argument("--trim_tip_names", action="store_true")
    p.add_argument("--allow_pad", action="store_true")
    p.add_argument("--id_match_mode", choices=["exact","relaxed"], default="exact",
                   help="Relaxed tries common variants (space before {…}, strip label, '_'↔' ', casefold).")
    p.add_argument("--name_map", help="CSV with headers tree_id,aln_id for explicit mappings.")
    return p.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Tree
    tree = Phylo.read(args.tree, "newick")
    if args.trim_tip_names:
        for t in tree.get_terminals():
            if t.name: t.name = t.name.strip()
    tree_tips = set(get_all_tips(tree.root))

    # Alignment
    aln, ntaxa_hdr, seqlen_hdr = read_phylip_auto_strict(args.aln, trim_names=args.trim_tip_names, allow_pad=args.allow_pad)
    aln_ids = set(rec.id for rec in aln)

    # Mapping
    explicit = load_name_map(args.name_map)
    t2a = build_id_mapping(tree_tips, aln_ids, mode=args.id_match_mode, name_map=explicit)

    present = {t for t in tree_tips if (t in aln_ids) or (t in t2a)}
    only_in_tree = tree_tips - present
    only_in_aln  = aln_ids - {t2a.get(t, t) for t in tree_tips}

    if only_in_tree and not args.allow_extra_tree_tips:
        raise SystemExit(
            f"ERROR: {len(only_in_tree)} tree tip(s) missing in alignment even after mapping: "
            f"{sorted(list(only_in_tree))[:5]}{' ...' if len(only_in_tree)>5 else ''}\n"
            "Try --id_match_mode relaxed or provide --name_map map.csv"
        )
    if only_in_aln and not args.allow_extra_aln_seqs:
        pass  # okay to ignore extras unless user forbids

    # Partition by cutting existing edges only
    n_target = max(1, min(args.n_groups, len(tree_tips)))
    groups = cut_edges_greedily(tree, n_target)

    # Write outputs
    for i, grp in enumerate(groups, start=1):
        tips = sorted(grp.tips & tree_tips)

        subtree_path = os.path.join(args.outdir, f"group_{i}.tree.nwk")
        write_group_subtree(tree, set(tips), subtree_path)

        subaln = subset_alignment_by_mapping(aln, tips, t2a)
        aln_path = os.path.join(args.outdir, f"group_{i}.phy")
        if len(subaln) > 0:
            AlignIO.write(subaln, aln_path, "phylip-relaxed")
        else:
            print(f"WARNING: group {i} has no matching alignment sequences; skipped {os.path.basename(aln_path)}.")

        taxa_path = os.path.join(args.outdir, f"group_{i}.taxa.txt")
        with open(taxa_path, "w", encoding="utf-8") as fh:
            for name in tips:
                fh.write(f"{name}\n")

        print(f"[group {i}] taxa={len(tips)} -> {os.path.basename(subtree_path)}, {os.path.basename(aln_path)}, {os.path.basename(taxa_path)}")

    sizes = [len(g.tips) for g in groups]
    print(f"\nDone. Created {len(groups)} groups from {len(tree_tips)} tree taxa.")
    print(f"PHYLIP header: ntaxa={ntaxa_hdr}, seqlen={seqlen_hdr}.")
    print("Group sizes:", sizes)

if __name__ == "__main__":
    main()
