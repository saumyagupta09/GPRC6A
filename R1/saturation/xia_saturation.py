#!/usr/bin/env python3
"""
Xia-style substitution saturation analysis for codon alignments (PHYLIP sequential),
with optional filtering to drop high-contribution codons and/or taxa and write
cleaned alignments you can immediately re-use for inference.

Outputs:
  - out_prefix + xia_summary.tsv
  - out_prefix + codon_contrib.tsv
  - out_prefix + species_delta.tsv
  - [optional] out_prefix + mc_fss.tsv (if --mc_iters > 0)
  - [optional] PHYLIP alignments:
      * out_prefix + filtered_codons.phy
      * out_prefix + filtered_taxa.phy
      * out_prefix + filtered_combined.phy
  - [optional] manifests listing removals:
      * out_prefix + removed_codons.tsv
      * out_prefix + removed_taxa.tsv
"""

import argparse
import math
from pathlib import Path
from collections import Counter
import random
import csv
import sys

IUPAC_VALID = set("ACGTacgt?-Nn")  # ACGT as states; others => missing/ambiguous at that site

# ========================= PHYLIP (sequential) I/O =========================

import re
from pathlib import Path

# Accept IUPAC DNA + gap/missing (case-insensitive). We will strip any other chars.
_IUPAC_KEEP_RE = re.compile(r"[ACGTRYKMSWBDHVN\-\?]", re.IGNORECASE)

def read_sequential_phy(filepath):
    """
    Robust PHYLIP sequential reader:
      - Reads exactly ntaxa blocks (name + wrapped sequence lines).
      - For each taxon, keeps appending only valid DNA/IUPAC chars until nsites is reached.
      - If a line contains any stray chars, they are ignored.
      - If we overshoot nsites, extra chars are dropped (with a warning).
    Returns: taxa(list), seqs(dict name->str), nsites(int), wrap_len(int)
    """
    lines = Path(filepath).read_text().splitlines()
    if not lines:
        raise ValueError("Empty file.")

    # Header
    parts = lines[0].split()
    if len(parts) < 2:
        raise ValueError("Header must have at least two integers: 'ntaxa nsites'.")
    try:
        ntaxa_hdr, nsites_hdr = int(parts[0]), int(parts[1])
    except Exception as e:
        raise ValueError("Header parse failed for 'ntaxa nsites'.") from e

    taxa, seqs = [], {}
    i = 1  # pointer into lines
    wrap_len = None

    # read exactly ntaxa blocks
    for t in range(ntaxa_hdr):
        # skip blank lines between blocks
        while i < len(lines) and not lines[i].strip():
            i += 1
        if i >= len(lines):
            raise ValueError(f"Unexpected EOF while looking for taxon name #{t+1}.")

        # taxon name line (take as-is)
        name = lines[i].strip()
        i += 1
        if not name:
            raise ValueError(f"Empty taxon name at block #{t+1}.")
        if name in seqs:
            raise ValueError(f"Duplicate taxon name encountered: {name}")
        taxa.append(name)
        seq_accum = []

        # accumulate sequence chars until we reach nsites_hdr
        while len(seq_accum) < nsites_hdr:
            # skip blank lines
            while i < len(lines) and not lines[i].strip():
                i += 1
            if i >= len(lines):
                raise ValueError(
                    f"Unexpected EOF while reading sequence for '{name}': "
                    f"got {len(seq_accum)} of {nsites_hdr} required characters."
                )
            chunk_raw = lines[i].strip()
            i += 1

            # keep only valid DNA/IUPAC + gap/missing chars; ignore everything else silently
            kept = _IUPAC_KEEP_RE.findall(chunk_raw)
            if kept:
                if wrap_len is None:
                    wrap_len = len("".join(kept))  # first non-empty chunk length, as a hint
                # append but clamp if this chunk would overshoot nsites_hdr
                need = nsites_hdr - len(seq_accum)
                if len(kept) > need:
                    # warn once per taxon if we truncate
                    # (comment out the print if you prefer to stay silent)
                    print(f"Warning: extra characters beyond declared length for '{name}' "
                          f"were truncated.", file=sys.stderr)
                    kept = kept[:need]
                seq_accum.extend(kept)
            # if kept is empty (line had no valid bases), just continue

        seq = "".join(seq_accum).upper()
        seqs[name] = seq

    # optional: warn if there are trailing non-empty lines after reading ntaxa blocks
    trailing_nonempty = any(ln.strip() for ln in lines[i:])
    if trailing_nonempty:
        print("Warning: trailing lines after reading declared taxa were ignored.", file=sys.stderr)

    # Final sanity
    if len(taxa) != ntaxa_hdr:
        print(f"Warning: header ntaxa={ntaxa_hdr}, parsed {len(taxa)} taxa.", file=sys.stderr)

    lens = {len(seqs[t]) for t in taxa}
    if len(lens) != 1 or lens.pop() != nsites_hdr:
        raise ValueError("Failed to normalize sequences to declared length.")

    if nsites_hdr % 3 != 0:
        raise ValueError("Alignment length not divisible by 3 (required for codons).")

    if wrap_len is None or wrap_len <= 0:
        wrap_len = 60  # reasonable default for writing

    return taxa, seqs, nsites_hdr, wrap_len

def write_sequential_phy(taxa, seqs, wrap_len, out_path):
    """
    Writes PHYLIP sequential:
      header: 'ntaxa nsites'
      name line, then wrapped sequence lines (same wrap_len as input)
    """
    ntaxa = len(taxa)
    nsites = len(seqs[taxa[0]]) if taxa else 0
    with open(out_path, "w") as fh:
        fh.write(f"{ntaxa:4d}{nsites:8d}\n")
        for name in taxa:
            fh.write(f"{name}\n")
            s = seqs[name]
            for i in range(0, len(s), wrap_len):
                fh.write(s[i:i+wrap_len] + "\n")

# ========================= Entropy / Iss core ==============================

def shannon_entropy_from_counts(counts):
    total = sum(counts.values())
    if total == 0:
        return 0.0
    H = 0.0
    for k in ("A","C","G","T"):
        c = counts.get(k, 0)
        if c > 0:
            p = c / total
            H -= p * math.log(p)
    return H  # nats

def site_counts(column_chars):
    c = Counter([b for b in column_chars if b in ("A","C","G","T")])
    return c, sum(c.values())

def pool_base_freqs(matrix):
    c = Counter()
    for col in zip(*matrix):
        for b in col:
            if b in ("A","C","G","T"):
                c[b] += 1
    total = sum(c.values())
    if total == 0:
        return {"A":0.25,"C":0.25,"G":0.25,"T":0.25}
    return {b: c[b]/total for b in ("A","C","G","T")}

def expected_entropy_FSS(ntaxa, p):
    """
    Exact E[H] and Var[H] per-site under full substitution saturation (multinomial).
    Complexity ~ O(N^3). Fine for typical ntaxa (<= ~100).
    """
    N = ntaxa
    pA, pC, pG, pT = p["A"], p["C"], p["G"], p["T"]
    EH = 0.0
    EH2 = 0.0
    for na in range(N+1):
        for nc in range(N-na+1):
            for ng in range(N-na-nc+1):
                nt = N - na - nc - ng
                prob = (math.comb(N, na) *
                        math.comb(N-na, nc) *
                        math.comb(N-na-nc, ng) *
                        (pA**na) * (pC**nc) * (pG**ng) * (pT**nt))
                H = 0.0
                for cnt in (na, nc, ng, nt):
                    if cnt > 0:
                        frac = cnt / N
                        H -= frac * math.log(frac)
                EH += prob * H
                EH2 += prob * (H*H)
    varH = max(EH2 - EH*EH, 0.0)
    return EH, varH

def compute_Iss(taxa, seqs):
    names = taxa
    ntaxa = len(names)
    L = len(seqs[names[0]])
    matrix = [seqs[n] for n in names]

    H_sites = []
    pos_buckets = {0: [], 1: [], 2: []}
    valid_sites = 0

    for i in range(L):
        col = [s[i] for s in matrix]
        counts, obs = site_counts(col)
        if obs < 2:
            continue
        H = shannon_entropy_from_counts(counts)
        H_sites.append(H)
        valid_sites += 1
        pos_buckets[i % 3].append(H)

    if valid_sites == 0:
        raise ValueError("No usable sites (after removing all-missing/ambiguous).")

    H_bar = sum(H_sites) / valid_sites
    p = pool_base_freqs(matrix)
    HFSS, varHFSS = expected_entropy_FSS(ntaxa, p)
    Iss = H_bar / HFSS if HFSS > 0 else float('nan')
    pos_Iss = {pos: ( (sum(vals)/len(vals))/HFSS if vals and HFSS>0 else float('nan') )
               for pos, vals in pos_buckets.items()}

    return {
        "H_bar": H_bar,
        "HFSS": HFSS,
        "Var_HFSS": varHFSS,
        "Iss": Iss,
        "pos_Iss": pos_Iss,
        "valid_sites": valid_sites,
        "L": L,
        "pools": p
    }

# ========================= Contributions & filtering ========================

def per_codon_contrib(taxa, seqs, HFSS):
    names = taxa
    L = len(seqs[names[0]])
    n_codons = L // 3
    contrib = []
    for c in range(n_codons):
        Hs = []
        for offset in (0,1,2):
            i = 3*c + offset
            col = [seqs[n][i] for n in names]
            counts, obs = site_counts(col)
            if obs < 2:
                continue
            Hs.append(shannon_entropy_from_counts(counts))
        codon_Iss = (sum(Hs)/len(Hs))/HFSS if Hs and HFSS>0 else 0.0
        contrib.append({"codon_index": c+1, "Iss_contrib": codon_Iss})
    return contrib

def species_delta_Iss(taxa, seqs, base_Iss):
    deltas = []
    for t in taxa:
        taxa_sub = [x for x in taxa if x != t]
        seqs_sub = {x: seqs[x] for x in taxa_sub}
        stats = compute_Iss(taxa_sub, seqs_sub)
        deltas.append({"taxon": t, "delta_Iss": stats["Iss"] - base_Iss})
    return deltas

def build_filtered_by_codons(taxa, seqs, codons_to_drop):
    """Drop whole codons (1-based indices)."""
    L = len(next(iter(seqs.values())))
    n_codons = L // 3
    drop_set0 = {c-1 for c in codons_to_drop}
    keep_nt_cols = []
    for c in range(n_codons):
        if c not in drop_set0:
            keep_nt_cols.extend([3*c, 3*c+1, 3*c+2])
    if not keep_nt_cols:
        raise ValueError("All codons were dropped; cannot write an empty alignment.")
    out = {}
    for t in taxa:
        s = seqs[t]
        out[t] = "".join(s[i] for i in keep_nt_cols)
    return out

def build_filtered_by_taxa(taxa, seqs, taxa_to_drop):
    keep_taxa = [t for t in taxa if t not in set(taxa_to_drop)]
    if not keep_taxa:
        raise ValueError("All taxa were dropped; cannot write an empty alignment.")
    out = {t: seqs[t] for t in keep_taxa}
    return keep_taxa, out

# ========================= Monte Carlo FSS (optional) =======================

def mc_fss_pvalue(ntaxa, basefreqs, L_eff, H_bar_obs, iters=2000, seed=None):
    rng = random.Random(seed)
    pA, pC, pG, pT = basefreqs["A"], basefreqs["C"], basefreqs["G"], basefreqs["T"]
    cats = ("A","C","G","T")
    weights = (pA, pC, pG, pT)
    highs = 0
    lows = 0
    Hs = []
    for _ in range(iters):
        hsum = 0.0
        for _ in range(L_eff):
            draws = [rng.choices(cats, weights=weights, k=ntaxa).count(b) for b in cats]
            counts = {b:c for b,c in zip(cats, draws)}
            hsum += shannon_entropy_from_counts(counts)
        Hbar = hsum / L_eff
        Hs.append(Hbar)
        if Hbar >= H_bar_obs: highs += 1
        if Hbar <= H_bar_obs: lows += 1
    p_hi = (highs+1)/(iters+1)
    p_lo = (lows+1)/(iters+1)
    return p_hi, p_lo, (sum(Hs)/len(Hs))

# ========================= TSV helpers =====================================

def write_tsv(path, rows, header):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=header, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

# ========================= Main ============================================

def main():
    ap = argparse.ArgumentParser(description="Xia-style substitution saturation with optional codon/taxon filtering.")
    ap.add_argument("input", help="PHYLIP sequential codon alignment")
    ap.add_argument("--out_prefix", required=True, help="Prefix for output files")
    ap.add_argument("--codon_thresh", type=float, default=None,
                    help="Flag & drop codons with Iss_contrib >= this (e.g., 0.8)")
    ap.add_argument("--species_thresh", type=float, default=None,
                    help="Flag & drop taxa with ΔIss >= this (e.g., 0.02)")
    ap.add_argument("--mc_iters", type=int, default=0,
                    help="If >0, run Monte Carlo FSS test (e.g., 5000)")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for MC test")
    args = ap.parse_args()

    taxa, seqs, L, wrap_len = read_sequential_phy(args.input)
    base_stats = compute_Iss(taxa, seqs)

    # Summary
    summary_path = f"{args.out_prefix}xia_summary.tsv"
    with open(summary_path, "w") as f:
        f.write("metric\tvalue\n")
        f.write(f"ntaxa\t{len(taxa)}\n")
        f.write(f"L_total\t{base_stats['L']}\n")
        f.write(f"L_valid_sites\t{base_stats['valid_sites']}\n")
        f.write(f"H_bar\t{base_stats['H_bar']:.6f}\n")
        f.write(f"HFSS\t{base_stats['HFSS']:.6f}\n")
        f.write(f"Var_HFSS\t{base_stats['Var_HFSS']:.6f}\n")
        f.write(f"Iss\t{base_stats['Iss']:.6f}\n")
        f.write(f"Iss_pos1\t{base_stats['pos_Iss'][0]:.6f}\n")
        f.write(f"Iss_pos2\t{base_stats['pos_Iss'][1]:.6f}\n")
        f.write(f"Iss_pos3\t{base_stats['pos_Iss'][2]:.6f}\n")
        f.write(f"pools_A\t{base_stats['pools']['A']:.6f}\n")
        f.write(f"pools_C\t{base_stats['pools']['C']:.6f}\n")
        f.write(f"pools_G\t{base_stats['pools']['G']:.6f}\n")
        f.write(f"pools_T\t{base_stats['pools']['T']:.6f}\n")

    # Codon contributions
    codon_rows = per_codon_contrib(taxa, seqs, base_stats["HFSS"])
    if args.codon_thresh is not None:
        for r in codon_rows:
            r["flag"] = int(r["Iss_contrib"] >= args.codon_thresh)
    write_tsv(f"{args.out_prefix}codon_contrib.tsv", codon_rows,
              header=["codon_index","Iss_contrib"] + (["flag"] if args.codon_thresh is not None else []))

    # Species ΔIss
    species_rows = species_delta_Iss(taxa, seqs, base_stats["Iss"])
    if args.species_thresh is not None:
        for r in species_rows:
            r["flag"] = int(r["delta_Iss"] >= args.species_thresh)
    write_tsv(f"{args.out_prefix}species_delta.tsv", species_rows,
              header=["taxon","delta_Iss"] + (["flag"] if args.species_thresh is not None else []))

    # Optional Monte Carlo FSS test
    if args.mc_iters and args.mc_iters > 0:
        p_hi, p_lo, Hbar_sim_mean = mc_fss_pvalue(
            ntaxa=len(taxa),
            basefreqs=base_stats["pools"],
            L_eff=base_stats["valid_sites"],
            H_bar_obs=base_stats["H_bar"],
            iters=args.mc_iters,
            seed=args.seed
        )
        with open(f"{args.out_prefix}mc_fss.tsv","w") as f:
            f.write("stat\tvalue\n")
            f.write(f"Hbar_obs\t{base_stats['H_bar']:.6f}\n")
            f.write(f"Hbar_sim_mean\t{Hbar_sim_mean:.6f}\n")
            f.write(f"p_lo\t{p_lo:.6g}\n")
            f.write(f"p_hi\t{p_hi:.6g}\n")
        print(f"[MC FSS] p_lo={p_lo:.4g}, p_hi={p_hi:.4g}")

    # ----------------- Filtering & writing cleaned alignments -----------------

    wrote_any = False

    # Codon-filtered
    if args.codon_thresh is not None:
        codons_to_drop = [r["codon_index"] for r in codon_rows if r.get("flag", 0) == 1]
        if codons_to_drop:
            # Manifest of codons removed
            write_tsv(f"{args.out_prefix}removed_codons.tsv",
                      [{"codon_index": c} for c in codons_to_drop],
                      header=["codon_index"])
            seqs_c = build_filtered_by_codons(taxa, seqs, codons_to_drop)
            out_path_c = f"{args.out_prefix}filtered_codons.phy"
            write_sequential_phy(taxa, seqs_c, wrap_len, out_path_c)
            print(f"[FILTER] Wrote codon-filtered alignment: {out_path_c}  (dropped {len(codons_to_drop)} codons)")
            wrote_any = True
        else:
            print("[FILTER] No codons met the threshold; no codon-filtered file written.")

    # Taxon-filtered
    if args.species_thresh is not None:
        taxa_to_drop = [r["taxon"] for r in species_rows if r.get("flag", 0) == 1]
        if taxa_to_drop:
            write_tsv(f"{args.out_prefix}removed_taxa.tsv",
                      [{"taxon": t} for t in taxa_to_drop],
                      header=["taxon"])
            taxa_t, seqs_t = build_filtered_by_taxa(taxa, seqs, taxa_to_drop)
            out_path_t = f"{args.out_prefix}filtered_taxa.phy"
            write_sequential_phy(taxa_t, seqs_t, wrap_len, out_path_t)
            print(f"[FILTER] Wrote taxon-filtered alignment: {out_path_t}  (dropped {len(taxa_to_drop)} taxa)")
            wrote_any = True
        else:
            print("[FILTER] No taxa met the threshold; no taxon-filtered file written.")

    # Combined filter (apply codon + taxon)
    if args.codon_thresh is not None and args.species_thresh is not None:
        codons_to_drop = [r["codon_index"] for r in codon_rows if r.get("flag", 0) == 1]
        taxa_to_drop = [r["taxon"] for r in species_rows if r.get("flag", 0) == 1]
        if codons_to_drop or taxa_to_drop:
            taxa_ct, seqs_ct = build_filtered_by_taxa(taxa, seqs, taxa_to_drop) if taxa_to_drop else (taxa, seqs)
            seqs_ct2 = build_filtered_by_codons(taxa_ct, seqs_ct, codons_to_drop) if codons_to_drop else seqs_ct
            out_path_ct = f"{args.out_prefix}filtered_combined.phy"
            write_sequential_phy(taxa_ct, seqs_ct2, wrap_len, out_path_ct)
            print(f"[FILTER] Wrote combined-filtered alignment: {out_path_ct}  "
                  f"(dropped {len(taxa_to_drop)} taxa, {len(codons_to_drop)} codons)")
            wrote_any = True
        else:
            print("[FILTER] No taxa/codons met thresholds; no combined-filtered file written.")

    print(f"Summary: {summary_path}")
    print(f"Codons : {args.out_prefix}codon_contrib.tsv")
    print(f"Species: {args.out_prefix}species_delta.tsv")
    if args.mc_iters and args.mc_iters > 0:
        print(f"MC FSS : {args.out_prefix}mc_fss.tsv")
    if not wrote_any:
        print("[FILTER] No filtered alignments written (provide thresholds to enable).")

if __name__ == "__main__":
    main()
