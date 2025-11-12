#!/usr/bin/env python3
# coding: utf-8
"""
clans_add_groups.py - robust parser version (Python 3.6+)

Adds groups INSIDE <seqgroups>...</seqgroups>:
  * Protein-name groups (2nd whitespace token, e.g., 'CASR')
  * Per-sequence groups (one per node)

Indexing priority: sequence=<ID> -> nr=<ID> -> names= list order.
If <seqgroups> missing, creates one right after </seq>.
"""

import argparse
import re
import sys
from pathlib import Path
from collections import defaultdict
from colorsys import hsv_to_rgb
from typing import Dict, List, Tuple

# Block finders
RE_SEQ_BLOCK = re.compile(r"<\s*seq\s*>(?P<body>.*?)</\s*seq\s*>", re.IGNORECASE | re.DOTALL)
RE_SEQGROUPS_BLOCK = re.compile(r"(?s)<\s*seqgroups\s*>(.*?)(<\s*/\s*seqgroups\s*>)", re.IGNORECASE)
RE_END_SEQ = re.compile(r"</\s*seq\s*>", re.IGNORECASE)

# Very loose attribute matchers (work in any order, across lines)
RE_NAME_ANY = re.compile(r"\bname\s*=\s*([^\n\r><]*)", re.IGNORECASE)
RE_SEQ_ANY  = re.compile(r"\bsequence\s*=\s*(\d+)", re.IGNORECASE)
RE_NR_ANY   = re.compile(r"\bnr\s*=\s*(\d+)", re.IGNORECASE)

# Fallback: names= list (post-</seq>)
RE_NAMES_LIST = re.compile(r"\bnames\s*=\s*([^\n\r]+)")

def parse_args():
    p = argparse.ArgumentParser(description="Append protein-name and per-sequence groups inside <seqgroups>.")
    p.add_argument("input", help="Input CLANS file")
    p.add_argument("-o", "--output", help="Output file (default: <input>.tagged.clans)")
    p.add_argument("--alpha", type=int, default=255, help="RGBA alpha (0-255) for protein groups")
    p.add_argument("--min-size", type=int, default=1, help="Emit protein groups with at least this many members")
    p.add_argument("--no-individual", action="store_true", help="Skip one-per-sequence groups")
    p.add_argument("--prot-prefix", default="", help="Optional prefix for protein group names (e.g., 'prot:')")
    p.add_argument("--seq-prefix", default="", help="Optional prefix for per-sequence group names (e.g., 'seq:')")
    p.add_argument("--debug", action="store_true", help="Print parser diagnostics")
    return p.parse_args()

def read_text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace")

def write_text(p: Path, s: str):
    p.write_text(s, encoding="utf-8")

def extract_seq_block(txt: str) -> str:
    m = RE_SEQ_BLOCK.search(txt)
    return m.group("body") if m else ""

def choose_index_map(txt: str, debug: bool=False) -> Dict[int, str]:
    """
    Robust: parse the <seq>...</seq> block, split into chunks by '>',
    and for each chunk look for name=, and sequence= or nr= anywhere.
    If nothing found, fallback to names= list.
    """
    id2name = {}

    body = extract_seq_block(txt)
    if body:
        # split on '>' so each chunk looks like a node's attributes blob
        chunks = body.split('>')
        for ch in chunks:
            # ignore empty whitespace-only chunks
            if not ch.strip():
                continue
            # find name
            mname = RE_NAME_ANY.search(ch)
            if not mname:
                continue
            name = mname.group(1).strip()
            # try sequence=, else nr=
            midx = RE_SEQ_ANY.search(ch)
            if midx:
                idx = int(midx.group(1))
            else:
                midx = RE_NR_ANY.search(ch)
                if midx:
                    idx = int(midx.group(1))
                else:
                    # no explicit index in this chunk
                    continue
            id2name[idx] = name

    if id2name:
        if debug:
            sample = list(id2name.items())[:3]
            sys.stderr.write("Parsed {} nodes from <seq> block. Samples: {}\n".format(len(id2name), sample))
        return id2name

    # Fallback: names= list (semicolon-separated, ordered by index)
    m = RE_NAMES_LIST.search(txt)
    if m:
        raw = m.group(1).strip()
        items = [s for s in raw.split(";")]
        if items and items[-1] == "":
            items.pop()
        for i, nm in enumerate(items):
            id2name[i] = nm.strip()

    if debug:
        if id2name:
            sys.stderr.write("Used names= fallback with {} items.\n".format(len(id2name)))
        else:
            sys.stderr.write("Could not find <seq> block or names= list.\n")

    return id2name

def protein_from_name(full_name: str) -> str:
    """2nd whitespace token (e.g., 'CASR' in 'NP_000379.3 CASR [organism=...]')."""
    tokens = full_name.strip().split()
    if len(tokens) >= 2:
        return tokens[1]
    return re.sub(r"\[.*?\]", "", tokens[0]).strip() if tokens else "unnamed"

def make_colors(n: int, alpha: int) -> List[Tuple[int, int, int, int]]:
    colors = []
    if n <= 0:
        return colors
    for i in range(n):
        h = (i / float(n)) % 1.0
        r, g, b = hsv_to_rgb(h, 1.0, 1.0)
        colors.append((int(r * 255), int(g * 255), int(b * 255), alpha))
    if n == 1 and colors[0][0:3] == (0, 0, 0):
        colors[0] = (0, 0, 255, alpha)
    return colors

def serialize_group_block(name: str, color: Tuple[int,int,int,int], members: List[int]) -> str:
    r, g, b, a = color
    nums = ";".join(str(x) for x in sorted(members))
    return (
        "name={}\n"
        "type=7\n"
        "size=4\n"
        "hide=0\n"
        "color={};{};{};{}\n"
        "numbers={}\n".format(name, r, g, b, a, nums)
    )

def build_groups(id2name: Dict[int, str], alpha: int, min_size: int,
                 add_individual: bool, prot_prefix: str, seq_prefix: str) -> Tuple[str, int, int]:
    # Protein-name groups
    prot_groups = defaultdict(list)
    for idx, nm in id2name.items():
        prot = protein_from_name(nm)
        prot_groups[prot].append(idx)
    prot_groups = {k: v for k, v in prot_groups.items() if len(v) >= min_size}

    # Individual groups
    indiv_groups = {}
    if add_individual:
        for idx, nm in id2name.items():
            prot = protein_from_name(nm)
            label = "{}{}__seq{}".format(seq_prefix, prot, idx)
            indiv_groups[label] = [idx]

    # Colors
    prot_keys = sorted(prot_groups.keys())
    prot_colors = make_colors(len(prot_keys), alpha)
    prot_color_map = {k: c for k, c in zip(prot_keys, prot_colors)}

    indiv_keys = sorted(indiv_groups.keys())
    indiv_alpha = max(60, min(180, alpha // 2))
    indiv_colors = make_colors(len(indiv_keys), indiv_alpha)
    indiv_color_map = {k: c for k, c in zip(indiv_keys, indiv_colors)}

    # Serialize
    blocks = []
    for k in prot_keys:
        name = "{}{}".format(prot_prefix, k)
        blocks.append(serialize_group_block(name, prot_color_map[k], prot_groups[k]))
    for k in indiv_keys:
        blocks.append(serialize_group_block(k, indiv_color_map[k], indiv_groups[k]))

    return "\n".join(blocks), len(prot_keys), len(indiv_keys)

def inject_into_seqgroups(orig: str, blocks: str) -> str:
    """Insert blocks inside existing <seqgroups>...</seqgroups>, else create one after </seq>."""
    m = RE_SEQGROUPS_BLOCK.search(orig)
    if m:
        inner, end_tag = m.group(1), m.group(2)
        new_inner = (inner.rstrip() + "\n\n" +
                     "// ---- appended by clans_add_groups.py ----\n" +
                     blocks + "\n")
        start, stop = m.span(1)
        return orig[:start] + new_inner + orig[stop:]
    m2 = RE_END_SEQ.search(orig)
    if not m2:
        raise RuntimeError("Could not find </seq> nor <seqgroups> to place groups.")
    insert_pos = m2.end()
    new_block = ("\n<seqgroups>\n"
                 "// ---- created and appended by clans_add_groups.py ----\n" +
                 blocks + "\n"
                 "</seqgroups>\n")
    return orig[:insert_pos] + new_block + orig[insert_pos:]

def main():
    args = parse_args()
    in_path = Path(args.input)
    out_path = Path(args.output) if args.output else in_path.with_suffix(in_path.suffix + ".tagged.clans")

    txt = read_text(in_path)
    id2name = choose_index_map(txt, debug=args.debug)
    if not id2name:
        sys.stderr.write("ERROR: Could not extract (sequence/nr -> name) mappings from the CLANS file.\n")
        # Helpful hint:
        sys.stderr.write("HINT: Try running with --debug, or grep your file for 'sequence=' or 'nr=' or 'names='.\n")
        sys.exit(1)

    group_blocks, n_prot, n_indiv = build_groups(
        id2name,
        alpha=args.alpha,
        min_size=args.min_size,
        add_individual=(not args.no_individual),
        prot_prefix=args.prot_prefix,
        seq_prefix=args.seq_prefix,
    )

    try:
        new_txt = inject_into_seqgroups(txt, group_blocks)
    except RuntimeError as e:
        sys.stderr.write("ERROR: {}\n".format(e))
        sys.exit(2)

    write_text(out_path, new_txt)
    print("Wrote {} with {} protein groups and {} individual groups inside <seqgroups>."
          .format(out_path, n_prot, n_indiv))

if __name__ == "__main__":
    main()
