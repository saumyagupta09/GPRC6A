#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, re, argparse, math
from collections import OrderedDict, Counter, defaultdict

# ---------- Config (tweak if desired) ----------
MIN_CHANGES_FOR_RELIABLE_OMEGA = 10   # require at least this many codon differences to trust ?^
ENTROPY_WINDOW = 1                    # per-site entropy (no windowing); keep for future tuning
LOW_ENTROPY_THRESHOLD = 0.2           # normalized Shannon entropy threshold for "conserved" (0..1)
RADICAL_GRANTHAM = 100                # >100 considered radical AA change
# Positive selection heuristic:
OMEGA_POS_CUTOFF = 1.2                # ?^ >= 1.2 and enough changes
# Relaxed selection heuristic:
OMEGA_RELAX_LO = 0.6                  # 0.6 <= ?^ <= 1.2 with other red flags
OMEGA_RELAX_HI = 1.2

# ---------- Genetic code and utilities ----------
STANDARD_CODE = {
    # U -> T mapping implicit; use DNA letters
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
IUPAC_DNA = set("ACGTRYKMSWBDHVN-?")

# Grantham distances (subset; fill all pairs by symmetry)
_GRANTHAM = {
    ('A','R'):112,('A','N'):111,('A','D'):126,('A','C'):195,('A','Q'):91,('A','E'):107,('A','G'):60,('A','H'):86,('A','I'):94,('A','L'):96,('A','K'):106,('A','M'):84,('A','F'):113,('A','P'):27,('A','S'):99,('A','T'):58,('A','W'):148,('A','Y'):112,('A','V'):64,
    ('R','N'):86,('R','D'):96,('R','C'):180,('R','Q'):43,('R','E'):54,('R','G'):125,('R','H'):29,('R','I'):97,('R','L'):102,('R','K'):26,('R','M'):91,('R','F'):97,('R','P'):103,('R','S'):110,('R','T'):71,('R','W'):101,('R','Y'):77,('R','V'):96,
    ('N','D'):23,('N','C'):139,('N','Q'):46,('N','E'):42,('N','G'):80,('N','H'):68,('N','I'):149,('N','L'):153,('N','K'):94,('N','M'):142,('N','F'):158,('N','P'):91,('N','S'):46,('N','T'):65,('N','W'):174,('N','Y'):143,('N','V'):133,
    ('D','C'):154,('D','Q'):61,('D','E'):45,('D','G'):94,('D','H'):81,('D','I'):168,('D','L'):172,('D','K'):101,('D','M'):160,('D','F'):177,('D','P'):108,('D','S'):65,('D','T'):85,('D','W'):181,('D','Y'):160,('D','V'):152,
    # ... for brevity we’ll auto-complete symmetry & self=0 below
}
def grantham(a, b):
    if a == b: return 0
    if (a,b) in _GRANTHAM: return _GRANTHAM[(a,b)]
    if (b,a) in _GRANTHAM: return _GRANTHAM[(b,a)]
    return 100  # neutral-ish default if missing

def is_pure_codon(c):
    return all(ch in "ACGT" for ch in c)

def translate(codon):
    if not is_pure_codon(codon): return None
    return STANDARD_CODE.get(codon, None)

# ---------- PHYLIP sequential reader (robust) ----------
SEQ_RE = re.compile(r"^[ACGTRYKMSWBDHVN\-\?]{10,}$", re.IGNORECASE)

def read_phylip_sequential(path):
    with open(path) as f:
        raw = [ln.rstrip("\r\n") for ln in f]
    if not raw: raise ValueError("Empty file.")
    p = raw[0].split()
    if len(p) < 2: raise ValueError("Header must be: '<ntaxa> <length>'")
    ntaxa, nsites = int(p[0]), int(p[1])

    taxa, seqs, i = [], OrderedDict(), 1
    for _ in range(ntaxa):
        while i < len(raw) and not raw[i].strip(): i += 1
        if i >= len(raw): raise ValueError("Unexpected EOF looking for taxon name.")
        name = raw[i].strip(); i += 1
        if not name: raise ValueError("Empty taxon name.")
        if name in seqs: raise ValueError(f"Duplicate taxon: {name}")
        taxa.append(name)
        acc = []
        while len(acc) < nsites:
            while i < len(raw) and not raw[i].strip(): i += 1
            if i >= len(raw):
                raise ValueError(f"EOF while reading {name}: have {len(acc)}/{nsites}.")
            chunk = raw[i].strip(); i += 1
            if SEQ_RE.match(chunk):
                need = nsites - len(acc)
                acc.extend(list(chunk.upper())[:need])
            else:
                # It's not a pure seq line; treat as noise OR a (badly wrapped) name leak.
                # We try not to misclassify: push pointer one step back and break (next loop will see it as a name).
                i -= 1
                break
        seqs[name] = "".join(acc)
    # basic sanity
    Ls = {len(s) for s in seqs.values()}
    if len(Ls) != 1:
        raise ValueError("Not all sequences have the same length.")
    return taxa, seqs, list(Ls)[0]

# ---------- Global codon-wise mask (drop codon if any taxon has '-' or '?') ----------
def codonwise_mask(seqs):
    names = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    if L % 3 != 0:
        L -= (L % 3)
        seqs = {k: v[:L] for k, v in seqs.items()}
    ncod = L // 3
    keep = [True]*ncod
    for j in range(ncod):
        for nm in names:
            codon = seqs[nm][3*j:3*j+3]
            if any(ch in "-?" for ch in codon) or len(codon) < 3:
                keep[j] = False; break
    if not any(keep): raise RuntimeError("All codons masked.")
    kept_idx = [j for j,k in enumerate(keep) if k]
    trimmed = {}
    for nm in names:
        s = []
        for j in kept_idx:
            s.append(seqs[nm][3*j:3*j+3])
        trimmed[nm] = "".join(s)
    return trimmed

# ---------- Consensus & entropy ----------
def majority_codon(seqs):
    names = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    ncod = L // 3
    cons = []
    for j in range(ncod):
        cnt = Counter(seqs[nm][3*j:3*j+3] for nm in names)
        codon, _ = cnt.most_common(1)[0]
        cons.append(codon)
    return "".join(cons)

def site_entropy(seqs):
    """Per-nucleotide entropy normalized in [0,1]."""
    names = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    ent = [0.0]*L
    for i in range(L):
        cnt = Counter(seqs[nm][i] for nm in names if seqs[nm][i] not in "-?")
        total = sum(cnt.values())
        if total == 0:
            ent[i] = 0.0
        else:
            H = 0.0
            for v in cnt.values():
                p = v/total
                H -= p*math.log(p, 2)
            # max entropy with A/C/G/T is 2 bits
            ent[i] = min(H/2.0, 1.0)
    return ent

# ---------- Fast NG-like per-taxon summary against consensus ----------
def per_taxon_stats(seqs):
    names = list(seqs.keys())
    L = len(next(iter(seqs.values())))
    assert L % 3 == 0
    ncod = L // 3
    cons = majority_codon(seqs)
    ent = site_entropy(seqs)
    # Precompute per-codon "is conserved" (all three positions low entropy)
    conserved_codon = []
    for j in range(ncod):
        e = ent[3*j:3*j+3]
        conserved_codon.append( all(x <= LOW_ENTROPY_THRESHOLD for x in e) )

    # Stats per taxon
    out = {}
    for nm in names:
        seq = seqs[nm]
        nonsyn = syn = 0
        stops = 0
        radical = 0
        conserved_hits = 0
        comparable = 0
        for j in range(ncod):
            c0 = cons[3*j:3*j+3]
            c1 = seq[3*j:3*j+3]
            aa0 = translate(c0)
            aa1 = translate(c1)
            if aa0 is None or aa1 is None:
                continue
            comparable += 1
            if aa1 == "*":
                stops += 1
            if c0 == c1:
                continue
            if aa0 == aa1:
                syn += 1
            else:
                nonsyn += 1
                if aa1 != "*" and aa0 != "*":
                    if grantham(aa0, aa1) > RADICAL_GRANTHAM:
                        radical += 1
                if conserved_codon[j]:
                    conserved_hits += 1
        # crude ?^ proxy:
        omega = (nonsyn + 1e-9) / (syn + 1e-9)
        out[nm] = {
            "codons_compared": comparable,
            "syn_changes": syn,
            "nonsyn_changes": nonsyn,
            "omega_hat": omega,
            "radical_changes": radical,
            "conserved_hits": conserved_hits,
            "stops": stops
        }
    return out

# ---------- Candidate ranking ----------
def shortlist(stats):
    pos_candidates = []    # likely positive selection
    relax_candidates = []  # likely relaxed selection
    for nm, st in stats.items():
        changes = st["syn_changes"] + st["nonsyn_changes"]
        if changes < MIN_CHANGES_FOR_RELIABLE_OMEGA:
            continue
        w = st["omega_hat"]
        # Positive selection: high ?^ with enough signal
        if w >= OMEGA_POS_CUTOFF and st["nonsyn_changes"] >= st["syn_changes"]+2:
            pos_candidates.append((nm, w, st))
        # Relaxed selection: ?^ near 1, but red flags (conserved hits, stops)
        if OMEGA_RELAX_LO <= w <= OMEGA_RELAX_HI and (st["conserved_hits"] >= 5 or st["stops"] > 0):
            relax_candidates.append((nm, w, st))
        # Very high stops also suggest relaxation/pseudogenization regardless of ?^
        if st["stops"] >= 1 and (nm, w, st) not in relax_candidates:
            relax_candidates.append((nm, w, st))
    # Sort
    pos_candidates.sort(key=lambda x: (-x[1], -(x[2]["nonsyn_changes"])))
    relax_candidates.sort(key=lambda x: (-x[2]["stops"], -x[2]["conserved_hits"], -x[1]))
    return pos_candidates, relax_candidates

# ---------- I/O ----------
def write_csv(stats, path):
    hdr = ["taxon","codons_compared","syn_changes","nonsyn_changes","omega_hat","radical_changes","conserved_hits","stops"]
    with open(path, "w") as out:
        out.write(",".join(hdr)+"\n")
        for nm, st in sorted(stats.items()):
            out.write(",".join([
                nm,
                str(st["codons_compared"]),
                str(st["syn_changes"]),
                str(st["nonsyn_changes"]),
                f"{st['omega_hat']:.4f}",
                str(st["radical_changes"]),
                str(st["conserved_hits"]),
                str(st["stops"]),
            ])+"\n")

def main():
    ap = argparse.ArgumentParser(description="Quick-screen candidates for relaxed/positive selection from a codon PHYLIP-seq MSA.")
    ap.add_argument("--input", required=True, help="PHYLIP sequential codon alignment")
    ap.add_argument("--out", default="selection_screen.csv", help="Output CSV")
    ap.add_argument("--no-mask", action="store_true", help="Disable global codon-wise gap masking")
    args = ap.parse_args()

    taxa, seqs, L = read_phylip_sequential(args.input)
    if L % 3 != 0:
        print(f"Note: input length {L} not divisible by 3; trimming trailing {L%3} nt.", file=sys.stderr)
        seqs = {k: v[:L - (L%3)] for k,v in seqs.items()}
        L -= (L%3)

    if not args.no_mask:
        seqs = codonwise_mask(seqs)

    stats = per_taxon_stats(seqs)
    write_csv(stats, args.out)

    # Shortlist
    pos, relax = shortlist(stats)

    print("\n=== Positive-selection candidates (screening) ===")
    if not pos:
        print("  (none at chosen thresholds)")
    else:
        for nm, w, st in pos[:10]:
            print(f"  {nm:25s}  ?{w:.2f}  dN={st['nonsyn_changes']} dS={st['syn_changes']}  conserved_hits={st['conserved_hits']}  radical={st['radical_changes']}")

    print("\n=== Relaxed-selection candidates (screening) ===")
    if not relax:
        print("  (none at chosen thresholds)")
    else:
        for nm, w, st in relax[:10]:
            flag = " +STOP" if st["stops"]>0 else ""
            print(f"  {nm:25s}  ?{w:.2f}  conserved_hits={st['conserved_hits']}  stops={st['stops']}{flag}  dN={st['nonsyn_changes']} dS={st['syn_changes']}")

    print(f"\nWrote: {args.out}")
    print("Notes:")
    print(" - ? here is a rough proxy vs. consensus (Nei Gojobori-style); use HyPhy/PAML for inference.")
    print(" - 'conserved_hits' counts nonsyn changes that land in low-entropy codons (more worrying).")
    print(" - Any internal stops are flagged; check for pseudogenization and consider RELAX.")

if __name__ == "__main__":
    # complete Grantham table symmetry
    keys = list(_GRANTHAM.keys())
    for a,b in keys:
        _GRANTHAM[(b,a)] = _GRANTHAM[(a,b)]
    for aa in set("ACDEFGHIKLMNPQRSTVWY"):
        _GRANTHAM[(aa,aa)] = 0
    main()
