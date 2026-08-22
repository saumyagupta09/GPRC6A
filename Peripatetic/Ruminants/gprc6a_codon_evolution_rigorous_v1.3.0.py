#!/usr/bin/env python3
"""
GPRC6A rigorous codon-evolution follow-up pipeline
=========================================

Purpose
-------
This is a downstream analysis stage for the GPRC6A comparative-genomics
pipeline and intactness refiner.  It reconstructs reference-anchored coding
remnants from raw first-pass BLAST HSPs plus exact nucleotide rescue HSPs from
the refiner when available, then aligns intact CDSs and disrupted
remnants with MACSE, and runs a panel of PAML/codeml models to distinguish:

  * continued purifying constraint,
  * relaxation toward neutrality after gene loss,
  * lineage-specific rate shifts, and
  * episodic positive selection.

The script is deliberately conservative about pseudogenes.  Missing sequence is
never filled with the reference allele: it is represented as N. Explicit target
deletions are retained as deletions, target insertions observed in BLAST HSPs
are retained, and disrupted/remnant sequences are passed to MACSE as
"less-reliable" sequences so frameshifts and stops can be aligned rather than
silently discarded.

A key design goal is alignment sensitivity analysis.  PAML is rerun on several
closely related PAML-safe versions of the MACSE codon alignment, a second
de-novo MACSE alignment with moderately altered pseudogene costs, and, unless
--skip-macse-refine is used, a refined MACSE alignment as well. The report
therefore states whether the qualitative evolutionary interpretation survives
small alignment/filtering changes.

Expected upstream layout
------------------------
OUTDIR/
  reference_exons.fa
  assemblies/<ACCESSION>/analysis.json
  assemblies/<ACCESSION>/blast.tsv.gz
  refined_intactness/
    tables/refined_assembly_summary.tsv
    tables/species_intactness.tsv
    tables/exon_rescue.tsv
    tables/event_validation.tsv
    tables/phylogenetic_gene_loss_origins.tsv
    species_tree.newick                  # if a tree was supplied to the refiner

Optional upstream layout
------------------------
If gprc6a_ml_gene_recovery.py has also been run, high-confidence intact CDSs
from ml_gene_recovery/recovered_species_cds.fasta are used as an optional
replacement only when the BLAST-remnant reconstruction is incomplete. Loss
lineages continue to use the reference-anchored remnant so disabling lesions
are not hidden by a predicted intact transcript model.

External programs
-----------------
  MACSE 2.x   : either a `macse` executable on PATH or --macse-jar FILE
  PAML/codeml : `codeml` on PATH or --codeml FILE
  HyPhy       : `hyphy` on PATH or --hyphy FILE, used for RELAX

Typical installation with Bioconda:
  conda install -c conda-forge -c bioconda macse paml hyphy pandas biopython

v1.3 rigorous additions
----------------------
  1. Foreground tests retain local ruminant phylogenetic context and explicit
     sister-lineage guardians so a terminal branch is not a collapsed ancestral path.
  2. Separate terminal, ancestral-Ruminantia, and inferred 0->1 loss-origin
     branches are tested.
  3. Independently intact CDSs are aligned and QCed by themselves first, and
     PAML-safe codon masks are generated only after each foreground taxon subset
     has been chosen.
  4. A codeml model with a nonzero return code, nonfinite likelihood, numerical
     failure, or convergence warning cannot enter any LRT.
  5. HyPhy RELAX independently tests relaxation (K<1) versus intensification
     (K>1) relative to intact reference branches.
  6. A simulation-calibration module runs constrained, relaxation,
     relaxation-with-lesions, and episodic-positive-selection scenarios through
     MACSE -> masking -> PAML. This estimates false-positive rates and asks whether
     relaxation alone can reproduce the observed branch-site statistic.

The default output subdirectory is codon_evolution_rigorous_v1_3. The previous
codon_evolution directory is never overwritten.

Example
-------
python gprc6a_codon_paml_followup.py \
    --outdir GPRC6A_Ruminantia \
    --tree ruminantia_species_tree.newick \
    --paml-workers 4 \
    --resume

If MACSE was installed as a jar instead of a wrapper executable:
python gprc6a_codon_paml_followup.py \
    --outdir GPRC6A_Ruminantia \
    --macse-jar /path/to/macse_v2.07.jar \
    --resume

Interpretive warning
--------------------
A branch-wide omega near 1 is consistent with relaxation but is not proof of
when gene inactivation occurred. A terminal branch can contain both a period of
functional constraint before the disabling lesion and a period of neutral drift
after it. Conversely, a branch-wide omega below 1 does not exclude an episode
of positive selection at a subset of sites. The pipeline therefore combines
branch, neutral-foreground, branch-site, free-ratio, site-model, and alignment
robustness evidence instead of using a single omega threshold.
"""

from __future__ import annotations

import argparse
import copy
import csv
import gzip
import html
import io
import json
import logging
import math
import os
import re
import random
import shlex
import shutil
import subprocess
import sys
import textwrap
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


VERSION = "1.3.0"
DNA = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}

BLASTN_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qseq", "sseq",
]

INTACT_SPECIES = {"INTACT_HIGH", "INTACT_LIKELY", "EXTERNAL_INTACT"}
LOSS_SPECIES_STRONG = {"DISRUPTED_HIGH", "DISRUPTED_LIKELY"}
LOSS_SPECIES_EXPLORATORY = {"DISRUPTED_POSSIBLE"}

INTACT_ASSEMBLY = {
    "INTACT_HIGH", "INTACT_LIKELY", "INTACT_POSSIBLE_ALTERNATIVE_MODEL",
    "INTACT_COPY_OR_ALT_LOCUS", "MIXED_OR_RESCUED",
}
LOSS_ASSEMBLY = {"DISRUPTED_STRONG", "DISRUPTED_SUPPORTED", "DISRUPTED_POSSIBLE"}

SUPPORTED_EVENT = {
    "SUPPORTED_BY_MINIPROT", "SUPPORTED_BY_MINIPROT_STOP", "SUPPORTED_MISSING_EXON",
    "SUPPORTED_BUT_ASSEMBLY_SENSITIVE",
}
RESCUED_EVENT = {
    "CONTRADICTED_BY_INTACT_MODEL", "RESCUED_BY_PROTEIN_MODEL",
    "RESCUED_BY_SPLICED_PROTEIN_MODEL", "ALTERNATIVE_INTACT_MODEL",
    "COMPENSATED_OR_ALIGNMENT_ARTIFACT",
}


# -----------------------------------------------------------------------------
# Generic utilities
# -----------------------------------------------------------------------------


def safe_name(s: Any) -> str:
    x = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s).strip()).strip("_")
    return x or "unnamed"


def norm_taxon(s: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower()).strip("_")


def species_key(name: str) -> str:
    p = str(name).strip().split()
    return " ".join(p[:2]) if len(p) >= 2 else str(name).strip()


def fnum(x: Any, default: float = float("nan")) -> float:
    try:
        if x is None or x == "":
            return default
        v = float(x)
        return v
    except Exception:
        return default


def inum(x: Any, default: int = 0) -> int:
    try:
        return int(float(x))
    except Exception:
        return default


def boolish(x: Any) -> bool:
    if isinstance(x, bool):
        return x
    return str(x).strip().lower() in {"1", "true", "yes", "y"}


def clamp(x: float, lo: float = 0.0, hi: float = 1.0) -> float:
    return max(lo, min(hi, x))


def wrap_fasta(seq: str, n: int = 80) -> str:
    return "\n".join(seq[i:i+n] for i in range(0, len(seq), n))


def read_json(path: Path) -> Dict[str, Any]:
    return json.loads(path.read_text())


def write_json(obj: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, default=str))


def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    g = pd.DataFrame() if df is None else df.copy()
    for c in g.columns:
        g[c] = g[c].map(lambda v: json.dumps(v, sort_keys=True) if isinstance(v, (dict, list, tuple, set)) else v)
    g.to_csv(path, sep="\t", index=False)


def run_cmd(cmd: Sequence[Any], cwd: Optional[Path] = None, check: bool = True,
            stdout_path: Optional[Path] = None, log_path: Optional[Path] = None) -> subprocess.CompletedProcess:
    cmd = [str(x) for x in cmd]
    logging.debug("RUN %s", shlex.join(cmd))
    stdout_handle = None
    try:
        if stdout_path:
            stdout_path.parent.mkdir(parents=True, exist_ok=True)
            stdout_handle = open(stdout_path, "wb")
            cp = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=stdout_handle,
                                stderr=subprocess.PIPE, check=False)
        else:
            cp = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, check=False)
        if log_path:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "ab") as fh:
                fh.write(("$ " + shlex.join(cmd) + "\n").encode())
                if cp.stdout and not stdout_path:
                    fh.write(cp.stdout)
                if cp.stderr:
                    fh.write(cp.stderr)
                fh.write(b"\n")
        if check and cp.returncode != 0:
            err = (cp.stderr or b"").decode(errors="replace")[-6000:]
            raise RuntimeError(f"Command failed ({cp.returncode}): {shlex.join(cmd)}\n{err}")
        return cp
    finally:
        if stdout_handle:
            stdout_handle.close()


def exe(path_or_name: Optional[str], default_name: str) -> Optional[str]:
    if path_or_name:
        p = Path(path_or_name).expanduser()
        if p.exists():
            return str(p.resolve())
        q = shutil.which(path_or_name)
        return q
    return shutil.which(default_name)


def load_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV, treating missing/blank/no-column files as an empty table.

    Some upstream summary tables are legitimately empty. In particular,
    pandas writes a DataFrame with zero rows *and zero columns* as a file
    containing only a newline. Such a file has non-zero size but raises
    pandas.errors.EmptyDataError when read back.
    """
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()

    # Fast-path whitespace-only files, including the one-newline files
    # produced by DataFrame.to_csv() for a 0-row, 0-column DataFrame.
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            prefix = fh.read(8192)
        if not prefix.strip():
            logging.info("Optional TSV is blank; treating as empty: %s", path)
            return pd.DataFrame()
    except OSError:
        pass

    try:
        return pd.read_csv(path, sep="\t", dtype=str).fillna("")
    except pd.errors.EmptyDataError:
        logging.info("Optional TSV has no parseable columns; treating as empty: %s", path)
        return pd.DataFrame()


def parse_jsonish_list(x: Any) -> List[str]:
    if isinstance(x, list):
        return [str(v) for v in x]
    s = str(x).strip()
    if not s:
        return []
    try:
        v = json.loads(s)
        if isinstance(v, list):
            return [str(z) for z in v]
    except Exception:
        pass
    if ";" in s:
        return [z.strip() for z in s.split(";") if z.strip()]
    return [s]


def bh_fdr(pvalues: Sequence[Optional[float]]) -> List[float]:
    vals = []
    for i, p in enumerate(pvalues):
        if p is None or not math.isfinite(float(p)):
            continue
        vals.append((i, min(1.0, max(0.0, float(p)))))
    out = [float("nan")] * len(pvalues)
    if not vals:
        return out
    vals.sort(key=lambda x: x[1])
    m = len(vals)
    raw = [0.0] * m
    for rank, (_, p) in enumerate(vals, 1):
        raw[rank-1] = p * m / rank
    for j in range(m-2, -1, -1):
        raw[j] = min(raw[j], raw[j+1])
    for (idx, _), q in zip(vals, raw):
        out[idx] = min(1.0, q)
    return out


def chi2_sf(x: float, df: int) -> float:
    if not math.isfinite(x) or x <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(x / 2.0))
    if df == 2:
        return math.exp(-x / 2.0)
    try:
        from scipy.stats import chi2  # optional
        return float(chi2.sf(x, df))
    except Exception:
        # Wilson-Hilferty approximation for uncommon fallback df values.
        z = ((x / df) ** (1 / 3) - (1 - 2 / (9 * df))) / math.sqrt(2 / (9 * df))
        return 0.5 * math.erfc(z / math.sqrt(2))


def aic(log_likelihood: float, np_: int) -> float:
    if not math.isfinite(log_likelihood) or np_ <= 0:
        return float("nan")
    return 2 * np_ - 2 * log_likelihood


# -----------------------------------------------------------------------------
# Reference model and exact first-pass remnant reconstruction
# -----------------------------------------------------------------------------


@dataclass
class RefExon:
    name: str
    index: int
    seq: str
    length: int
    cds_start: int
    cds_end: int


@dataclass
class AssemblySequence:
    accession: str
    organism_name: str
    species: str
    assembly_level: str
    assembly_status: str
    species_status: str
    sequence: str
    source: str
    callable_fraction: float
    reconstructed_fraction: float
    n_fraction: float
    explicit_deleted_nt: int
    insertion_nt: int
    sensitive_rescued_nt: int
    internal_stop_count_naive: int
    frame_remainder: int
    supported_events: int
    rescued_events: int
    quality_score: float
    exon_sequences: Dict[int, str]
    warnings: List[str]


def load_reference_exons(path: Path) -> Tuple[List[RefExon], str]:
    recs = list(SeqIO.parse(str(path), "fasta"))
    if not recs:
        raise ValueError(f"No reference exon FASTA records in {path}")

    def exon_num(rec: SeqRecord, fallback: int) -> int:
        for text in (rec.id, rec.description):
            m = re.search(r"exon[_\s-]*(\d+)", text, re.I)
            if m:
                return int(m.group(1))
        m = re.search(r"(\d+)$", rec.id)
        return int(m.group(1)) if m else fallback

    tmp = sorted([(exon_num(r, i+1), r) for i, r in enumerate(recs)], key=lambda x: x[0])
    out = []
    pos = 1
    for idx, rec in tmp:
        seq = str(rec.seq).upper().replace("U", "T")
        if re.search(r"[^ACGTN]", seq):
            raise ValueError(f"Unexpected character in reference exon {rec.id}")
        out.append(RefExon(rec.id, idx, seq, len(seq), pos, pos + len(seq) - 1))
        pos += len(seq)
    return out, "".join(e.seq for e in out)


def load_blast(path: Path) -> pd.DataFrame:
    """Load first-pass/refiner BLAST tabular output safely.

    BLAST can legitimately produce an output file with zero data rows. Some
    empty gzip/plain files have non-zero byte size, and pandas may return a
    zero-row DataFrame rather than raising EmptyDataError. On a zero-row frame,
    DataFrame.apply(axis=1) can return a DataFrame instead of a Series, causing
    assignment to one column to fail. Handle emptiness explicitly and avoid
    row-wise DataFrame.apply entirely.
    """
    base_cols = list(BLASTN_COLS)
    derived_cols = ["strand", "gmin", "gmax"]
    empty = lambda: pd.DataFrame(columns=base_cols + derived_cols)

    if not path.exists() or path.stat().st_size == 0:
        return empty()

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            names=base_cols,
            header=None,
            compression="gzip" if path.suffix == ".gz" else None,
        )
    except (pd.errors.EmptyDataError, EOFError):
        return empty()
    except Exception as e:
        logging.warning("Could not read BLAST %s: %s", path, e)
        return empty()

    if df.empty:
        return empty()

    for c in base_cols:
        if c not in df.columns:
            df[c] = pd.NA

    for c in ["length", "qstart", "qend", "sstart", "send", "qlen", "slen"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    for c in ["pident", "evalue", "bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # pandas-version-independent and safe for zero rows.
    df["strand"] = [
        "+" if int(sstart) <= int(send) else "-"
        for sstart, send in zip(df["sstart"].tolist(), df["send"].tolist())
    ]
    df["gmin"] = df[["sstart", "send"]].min(axis=1)
    df["gmax"] = df[["sstart", "send"]].max(axis=1)
    return df


def parse_hsp_alignment(row: pd.Series) -> Tuple[Dict[int, Dict[str, Any]], List[Dict[str, Any]]]:
    qseq = str(row.qseq).upper()
    sseq = str(row.sseq).upper()
    qpos = int(row.qstart) - 1
    scoord = int(row.sstart)
    sstep = 1 if str(row.strand) == "+" else -1
    calls: Dict[int, Dict[str, Any]] = {}
    insertions: List[Dict[str, Any]] = []
    ins_chars: List[str] = []
    ins_coords: List[int] = []
    ins_anchor: Optional[int] = None

    def flush() -> None:
        nonlocal ins_chars, ins_coords, ins_anchor
        if ins_chars:
            insertions.append({
                "anchor_qpos": int(ins_anchor or 0),
                "sequence": "".join(ins_chars),
                "length": len(ins_chars),
                "subject_start": min(ins_coords) if ins_coords else None,
                "subject_end": max(ins_coords) if ins_coords else None,
                "bitscore": fnum(row.bitscore, 0.0),
                "source": "within_hsp",
            })
        ins_chars, ins_coords, ins_anchor = [], [], None

    for qc, sc in zip(qseq, sseq):
        if qc != "-":
            flush()
            qpos += 1
            call = {
                "qbase": qc,
                "sbase": sc,
                "scoord": None,
                "bitscore": fnum(row.bitscore, 0.0),
                "pident": fnum(row.pident, 0.0),
            }
            if sc != "-":
                call["scoord"] = scoord
                scoord += sstep
            calls[qpos] = call
        elif sc != "-":
            if ins_anchor is None:
                ins_anchor = qpos
            ins_chars.append(sc)
            ins_coords.append(scoord)
            scoord += sstep
    flush()
    return calls, insertions


def exon_alias_mask(df: pd.DataFrame, exon: RefExon) -> pd.Series:
    if df.empty:
        return pd.Series([], dtype=bool)
    exact = df.qseqid.astype(str) == exon.name
    if exact.any():
        return exact
    patt = re.compile(rf"(?:^|[^0-9])(?:exon[_-]?)?{exon.index}(?:[^0-9]|$)", re.I)
    return df.qseqid.astype(str).map(lambda x: bool(patt.search(x)))


def reconstruct_exon(exon: RefExon, hsps: pd.DataFrame) -> Dict[str, Any]:
    best_calls: Dict[int, Dict[str, Any]] = {}
    insertion_map: Dict[Tuple[int, str], Dict[str, Any]] = {}
    if not hsps.empty:
        for _, row in hsps.iterrows():
            calls, ins = parse_hsp_alignment(row)
            for p, c in calls.items():
                if p < 1 or p > exon.length:
                    continue
                if p not in best_calls or c["bitscore"] > best_calls[p]["bitscore"]:
                    best_calls[p] = c
            for x in ins:
                key = (int(x["anchor_qpos"]), str(x["sequence"]))
                if key not in insertion_map or x["bitscore"] > insertion_map[key]["bitscore"]:
                    insertion_map[key] = x

    projected: List[Optional[str]] = [None] * exon.length
    kinds = ["missing"] * exon.length
    callable_nt = 0
    recovered_nt = 0
    deleted_nt = 0
    for p in range(1, exon.length+1):
        c = best_calls.get(p)
        if not c:
            continue
        s = str(c["sbase"]).upper()
        if s == "-":
            kinds[p-1] = "deletion"
            deleted_nt += 1
            continue
        projected[p-1] = s
        kinds[p-1] = "base"
        recovered_nt += 1
        if s in DNA:
            callable_nt += 1

    ins_by_anchor: Dict[int, List[Dict[str, Any]]] = defaultdict(list)
    for x in insertion_map.values():
        ins_by_anchor[int(x["anchor_qpos"])].append(x)

    parts: List[str] = []
    insertion_nt = 0
    for x in sorted(ins_by_anchor.get(0, []), key=lambda z: -fnum(z.get("bitscore"), 0))[:1]:
        parts.append(str(x["sequence"]).upper())
        insertion_nt += len(str(x["sequence"]))
    for p in range(1, exon.length+1):
        if kinds[p-1] == "deletion":
            pass
        elif kinds[p-1] == "base" and projected[p-1]:
            parts.append(str(projected[p-1]).upper())
        else:
            parts.append("N")
        for x in sorted(ins_by_anchor.get(p, []), key=lambda z: -fnum(z.get("bitscore"), 0))[:1]:
            parts.append(str(x["sequence"]).upper())
            insertion_nt += len(str(x["sequence"]))

    return {
        "sequence": "".join(parts),
        "callable_nt": callable_nt,
        "recovered_nt": recovered_nt,
        "deleted_nt": deleted_nt,
        "insertion_nt": insertion_nt,
        "reference_nt": exon.length,
    }


def naive_internal_stops(seq: str) -> int:
    s = seq.upper().replace("U", "T")
    n = (len(s) // 3) * 3
    stops = 0
    for i in range(0, max(0, n-3), 3):
        codon = s[i:i+3]
        if set(codon) <= DNA and codon in STOP_CODONS:
            stops += 1
    return stops


def reconstruct_assembly(analysis: Dict[str, Any], blast: pd.DataFrame, ref_exons: List[RefExon],
                         sensitive_blast: Optional[pd.DataFrame] = None) -> Tuple[str, Dict[int, str], Dict[str, Any]]:
    locus = analysis.get("locus") or {}
    contig = str(locus.get("contig", ""))
    strand = str(locus.get("strand", "+"))
    lstart = inum(locus.get("start"), 0)
    lend = inum(locus.get("end"), 0)

    g = blast.copy()
    if not g.empty and contig:
        g = g[g.sseqid.astype(str) == contig]
    if not g.empty:
        g = g[g.strand.astype(str) == strand]
    if not g.empty and lstart and lend:
        lo, hi = min(lstart, lend), max(lstart, lend)
        g = g[(g.gmax >= lo) & (g.gmin <= hi)]

    sg = sensitive_blast.copy() if sensitive_blast is not None else pd.DataFrame(columns=BLASTN_COLS)

    cds_parts: List[str] = []
    exon_seqs: Dict[int, str] = {}
    callable_nt = recovered_nt = deleted_nt = inserted_nt = ref_nt = 0
    sensitive_rescued_nt = 0
    for exon in ref_exons:
        h = g[exon_alias_mask(g, exon)] if not g.empty else pd.DataFrame(columns=g.columns)
        sh = sg[exon_alias_mask(sg, exon)] if not sg.empty else pd.DataFrame(columns=sg.columns)
        # The refiner's sensitive BLAST is against the retained local-locus FASTA,
        # so its subject coordinates are local rather than whole-genome coordinates.
        # We use it only as additional exact qseq/sseq evidence. Highest-bitscore
        # calls win at positions already covered by the first-pass alignment, while
        # previously missing positions can be rescued without reference imputation.
        combined = pd.concat([h, sh], ignore_index=True, sort=False) if not sh.empty else h
        r0 = reconstruct_exon(exon, h)
        r = reconstruct_exon(exon, combined)
        sensitive_rescued_nt += max(0, int(r["callable_nt"]) - int(r0["callable_nt"]))
        exon_seqs[exon.index] = r["sequence"]
        cds_parts.append(r["sequence"])
        callable_nt += r["callable_nt"]
        recovered_nt += r["recovered_nt"]
        deleted_nt += r["deleted_nt"]
        inserted_nt += r["insertion_nt"]
        ref_nt += r["reference_nt"]
    seq = "".join(cds_parts)
    metrics = {
        "callable_fraction": callable_nt / ref_nt if ref_nt else 0.0,
        "reconstructed_fraction": recovered_nt / ref_nt if ref_nt else 0.0,
        "n_fraction": seq.count("N") / len(seq) if seq else 1.0,
        "explicit_deleted_nt": deleted_nt,
        "insertion_nt": inserted_nt,
        "sensitive_rescued_nt": sensitive_rescued_nt,
        "frame_remainder": len(seq) % 3,
        "internal_stop_count_naive": naive_internal_stops(seq),
        "reference_nt": ref_nt,
        "sequence_nt": len(seq),
    }
    return seq, exon_seqs, metrics


# -----------------------------------------------------------------------------
# Upstream result discovery and representative selection
# -----------------------------------------------------------------------------


def discover_upstream(outdir: Path, refine_subdir: str) -> Dict[str, Any]:
    refine = outdir / refine_subdir
    tables = refine / "tables"
    species_df = load_tsv(tables / "species_intactness.tsv")
    assembly_df = load_tsv(tables / "refined_assembly_summary.tsv")
    exon_df = load_tsv(tables / "exon_rescue.tsv")
    event_df = load_tsv(tables / "event_validation.tsv")
    origins_df = load_tsv(tables / "phylogenetic_gene_loss_origins.tsv")
    if species_df.empty:
        raise RuntimeError(f"Missing or empty second-stage species table: {tables/'species_intactness.tsv'}")
    if assembly_df.empty:
        raise RuntimeError(f"Missing or empty second-stage assembly table: {tables/'refined_assembly_summary.tsv'}")

    analyses = []
    for p in sorted((outdir / "assemblies").glob("*/analysis.json")):
        try:
            d = read_json(p)
            d["_path"] = str(p)
            analyses.append(d)
        except Exception as e:
            logging.warning("Could not read %s: %s", p, e)
    if not analyses:
        raise RuntimeError(f"No first-pass assemblies/*/analysis.json under {outdir}")

    return {
        "species": species_df,
        "assemblies": assembly_df,
        "exons": exon_df,
        "events": event_df,
        "origins": origins_df,
        "analyses": analyses,
        "refine_root": refine,
    }


def load_ml_species_sequences(outdir: Path, ml_subdir: str) -> Dict[str, Dict[str, Any]]:
    path = outdir / ml_subdir / "recovered_species_cds.fasta"
    if not path.exists():
        return {}
    out: Dict[str, Dict[str, Any]] = {}
    for rec in SeqIO.parse(str(path), "fasta"):
        desc = rec.description
        msp = re.search(r"species=(.*?)\s+assembly=", desc)
        macc = re.search(r"assembly=([^\s]+)", desc)
        mqual = re.search(r"quality=([^\s]+)", desc)
        sp = msp.group(1).strip() if msp else rec.id.split("|")[0].replace("_", " ")
        seq = str(rec.seq).upper().replace("U", "T")
        out[norm_taxon(sp)] = {
            "species": sp,
            "accession": macc.group(1) if macc else "",
            "quality": mqual.group(1) if mqual else "",
            "sequence": seq,
            "path": str(path),
        }
    return out


def assembly_level_score(level: str) -> float:
    x = str(level).lower()
    if "complete" in x:
        return 4.0
    if "chromosome" in x:
        return 3.0
    if "scaffold" in x:
        return 2.0
    if "contig" in x:
        return 1.0
    return 0.0


def species_status_map(species_df: pd.DataFrame) -> Dict[str, str]:
    """Map species names to refiner species-level status.

    Store both the exact species_key and a normalized alias. This protects the
    follow-up against harmless whitespace/punctuation differences between the
    first-pass JSONs, refiner TSVs, and species-tree labels.
    """
    if species_df.empty:
        return {}
    out: Dict[str, str] = {}
    if "species_key" not in species_df.columns or "species_status" not in species_df.columns:
        logging.warning("species_intactness table lacks species_key/species_status columns")
        return out
    for _, r in species_df.iterrows():
        sp = str(r.get("species_key", "")).strip()
        st = str(r.get("species_status", "")).strip()
        if not sp:
            continue
        out[sp] = st
        out[norm_taxon(sp)] = st
    return out


def assembly_status_map(assembly_df: pd.DataFrame) -> Dict[str, str]:
    if assembly_df.empty:
        return {}
    return {str(r.assembly_accession): str(r.assembly_status) for _, r in assembly_df.iterrows()}


def event_counts(event_df: pd.DataFrame) -> Dict[str, Tuple[int, int]]:
    out: Dict[str, Tuple[int, int]] = {}
    if event_df.empty:
        return out
    for acc, g in event_df.groupby("assembly_accession"):
        st = set(g.validation_status.astype(str))
        sup = int(sum(x in SUPPORTED_EVENT for x in g.validation_status.astype(str)))
        resc = int(sum(x in RESCUED_EVENT for x in g.validation_status.astype(str)))
        out[str(acc)] = (sup, resc)
    return out


def exon_mean_quality(exon_df: pd.DataFrame) -> Dict[str, float]:
    out = {}
    if exon_df.empty:
        return out
    for acc, g in exon_df.groupby("assembly_accession"):
        vals = []
        for _, r in g.iterrows():
            vals.append(max(fnum(r.get("original_coverage"), 0),
                            fnum(r.get("protein_model_aa_coverage"), 0),
                            fnum(r.get("sensitive_blast_coverage"), 0)))
        out[str(acc)] = sum(vals)/len(vals) if vals else 0.0
    return out


def status_compatibility_score(species_status: str, assembly_status: str) -> float:
    if species_status in INTACT_SPECIES:
        order = {
            "EXTERNAL_REFSEQ_INTACT": 60,
            "INTACT_HIGH": 50, "INTACT_LIKELY": 45,
            "INTACT_POSSIBLE_ALTERNATIVE_MODEL": 38,
            "INTACT_COPY_OR_ALT_LOCUS": 30, "MIXED_OR_RESCUED": 25,
        }
        return float(order.get(assembly_status, 0))
    if species_status in LOSS_SPECIES_STRONG | LOSS_SPECIES_EXPLORATORY:
        order = {"DISRUPTED_STRONG": 50, "DISRUPTED_SUPPORTED": 45, "DISRUPTED_POSSIBLE": 35}
        return float(order.get(assembly_status, 0))
    return 10.0 if assembly_status not in {"UNRESOLVED_ERROR", "UNRESOLVED_NO_GENOME"} else 0.0


def build_assembly_sequences(outdir: Path, upstream: Dict[str, Any], ref_exons: List[RefExon],
                             ml_species: Dict[str, Dict[str, Any]]) -> List[AssemblySequence]:
    species_status = species_status_map(upstream["species"])
    assembly_status = assembly_status_map(upstream["assemblies"])
    ev_counts = event_counts(upstream["events"])
    exon_q = exon_mean_quality(upstream["exons"])
    rows: List[AssemblySequence] = []

    for i, a in enumerate(upstream["analyses"], 1):
        acc = str(a.get("assembly_accession", ""))
        org = str(a.get("organism_name", ""))
        sp = str(a.get("species_key") or species_key(org))
        sstat = species_status.get(sp, species_status.get(norm_taxon(sp), "UNRESOLVED"))
        astat = assembly_status.get(acc, "UNRESOLVED")
        blast_path = outdir / "assemblies" / safe_name(acc) / "blast.tsv.gz"
        if not blast_path.exists():
            # Some first-pass folders may not have been sanitized identically.
            ap = Path(str(a.get("_path", ""))).parent
            blast_path = ap / "blast.tsv.gz"
        blast = load_blast(blast_path)
        refine_root = Path(str(upstream.get("refine_root", outdir / "refined_intactness")))
        sensitive_path = refine_root / "assemblies" / safe_name(acc) / "sensitive_exon_blast.tsv"
        sensitive = load_blast(sensitive_path) if sensitive_path.exists() else pd.DataFrame(columns=BLASTN_COLS)
        warnings: List[str] = []
        if blast.empty:
            warnings.append("raw_blast_missing_or_empty")
        if not sensitive.empty:
            warnings.append("refiner_sensitive_exon_blast_incorporated")
        seq, exon_seqs, m = reconstruct_assembly(a, blast, ref_exons, sensitive)
        sup, resc = ev_counts.get(acc, (0, 0))
        qual = (
            100.0 * fnum(m.get("callable_fraction"), 0) +
            20.0 * exon_q.get(acc, 0.0) +
            3.0 * assembly_level_score(str(a.get("assembly_level", ""))) +
            status_compatibility_score(sstat, astat) -
            10.0 * fnum(m.get("n_fraction"), 1.0)
        )

        source = ("first_pass_plus_refiner_sensitive_remnant"
                  if inum(m.get("sensitive_rescued_nt"), 0) > 0 else
                  "first_pass_reference_anchored_remnant")
        ml = ml_species.get(norm_taxon(sp))
        # Use an ML-recovered exact CDS only for an intact species and only when
        # the reference-anchored reconstruction has meaningful missing data.
        # Disrupted species always retain the remnant sequence to avoid hiding lesions.
        if (sstat in INTACT_SPECIES and ml and
                str(ml.get("quality")) in {"RECOVERED_HIGH", "RECOVERED_LIKELY"} and
                fnum(m.get("callable_fraction"), 0) < 0.95):
            mseq = str(ml.get("sequence", "")).upper()
            if len(mseq) >= 300 and not re.search(r"[^ACGTN]", mseq):
                seq = mseq
                source = "ml_recovery_exact_cds"
                warnings.append("intact_sequence_replaced_by_high_confidence_ml_cds_due_to_missing_blast_bases")
                ref_len = sum(ex.length for ex in ref_exons)
                callable_ml = sum(base in DNA for base in seq)
                recovered_ml = len(seq) - seq.count("N")
                m["n_fraction"] = seq.count("N") / len(seq) if seq else 1.0
                m["callable_fraction"] = min(1.0, callable_ml / ref_len) if ref_len else 0.0
                m["reconstructed_fraction"] = min(1.0, recovered_ml / ref_len) if ref_len else 0.0
                m["frame_remainder"] = len(seq) % 3
                m["internal_stop_count_naive"] = naive_internal_stops(seq)

        rows.append(AssemblySequence(
            accession=acc,
            organism_name=org,
            species=sp,
            assembly_level=str(a.get("assembly_level", "")),
            assembly_status=astat,
            species_status=sstat,
            sequence=seq,
            source=source,
            callable_fraction=fnum(m.get("callable_fraction"), 0.0),
            reconstructed_fraction=fnum(m.get("reconstructed_fraction"), 0.0),
            n_fraction=fnum(m.get("n_fraction"), 1.0),
            explicit_deleted_nt=inum(m.get("explicit_deleted_nt"), 0),
            insertion_nt=inum(m.get("insertion_nt"), 0),
            sensitive_rescued_nt=inum(m.get("sensitive_rescued_nt"), 0),
            internal_stop_count_naive=inum(m.get("internal_stop_count_naive"), 0),
            frame_remainder=inum(m.get("frame_remainder"), 0),
            supported_events=sup,
            rescued_events=resc,
            quality_score=qual,
            exon_sequences=exon_seqs,
            warnings=warnings,
        ))
        if i % 25 == 0:
            logging.info("Reconstructed %d/%d assemblies", i, len(upstream["analyses"]))
    return rows


def load_external_intact_fasta(path: Optional[Path], ref_cds: str) -> Dict[str, AssemblySequence]:
    """Load independently curated intact outgroup GPRC6A CDS sequences.

    FASTA record IDs should be stable taxon IDs such as ``Homo_sapiens``.
    An optional ``|ACCESSION`` suffix is retained as provenance, e.g.
    ``Homo_sapiens|NM_001286312.2``.

    External sequences are accepted only when they are codon-complete, have no
    internal stop codon, have very little ambiguous sequence, and are broadly
    compatible with the reference CDS length. They are never used as loss
    foregrounds; they provide intact background branches for PAML.
    """
    if path is None:
        return {}
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"External intact FASTA not found: {path}")

    ref_len = len(ref_cds)
    out: Dict[str, AssemblySequence] = {}
    rejected = []
    for rec in SeqIO.parse(str(path), "fasta"):
        raw_id = str(rec.id)
        taxon_id, _, accession = raw_id.partition("|")
        species = taxon_id.replace("_", " ").strip()
        seq = re.sub(r"\s+", "", str(rec.seq).upper().replace("U", "T"))

        # Strip a conventional terminal stop codon. Internal stops remain fatal.
        if len(seq) >= 3 and seq[-3:] in STOP_CODONS:
            seq = seq[:-3]

        reasons = []
        if len(seq) < max(900, int(0.65 * ref_len)):
            reasons.append("too_short")
        if len(seq) > int(1.40 * ref_len):
            reasons.append("too_long")
        if len(seq) % 3 != 0:
            reasons.append("not_divisible_by_3")
        if re.search(r"[^ACGTN]", seq):
            reasons.append("non_IUPAC_ACGTN_character")
        nfrac = seq.count("N") / len(seq) if seq else 1.0
        if nfrac > 0.01:
            reasons.append("N_fraction_gt_0.01")
        stops = naive_internal_stops(seq)
        if stops:
            reasons.append(f"internal_stops_{stops}")

        if reasons:
            rejected.append((species, accession or raw_id, ";".join(reasons)))
            continue

        callable_fraction = min(1.0, sum(b in DNA for b in seq) / ref_len) if ref_len else 0.0
        reconstructed_fraction = min(1.0, (len(seq) - seq.count("N")) / ref_len) if ref_len else 0.0
        length_closeness = max(0.0, 1.0 - abs(len(seq) - ref_len) / max(ref_len, 1))
        qscore = 250.0 + 100.0 * callable_fraction + 30.0 * length_closeness - 20.0 * nfrac

        out[species] = AssemblySequence(
            accession=accession or raw_id,
            organism_name=species,
            species=species,
            assembly_level="curated_external_CDS",
            assembly_status="EXTERNAL_REFSEQ_INTACT",
            species_status="EXTERNAL_INTACT",
            sequence=seq,
            source=f"external_intact_fasta:{path.name}",
            callable_fraction=callable_fraction,
            reconstructed_fraction=reconstructed_fraction,
            n_fraction=nfrac,
            explicit_deleted_nt=0,
            insertion_nt=0,
            sensitive_rescued_nt=0,
            internal_stop_count_naive=0,
            frame_remainder=0,
            supported_events=0,
            rescued_events=0,
            quality_score=qscore,
            exon_sequences={},
            warnings=["independent_external_intact_background"],
        )

    if rejected:
        logging.warning(
            "Rejected %d external intact FASTA record(s): %s",
            len(rejected),
            " | ".join(f"{sp}:{acc}:{why}" for sp, acc, why in rejected[:12])
        )
    logging.info("Validated external intact background sequences: %d", len(out))
    return out

def choose_species_representatives(assemblies: List[AssemblySequence]) -> Dict[str, AssemblySequence]:
    groups: Dict[str, List[AssemblySequence]] = defaultdict(list)
    for a in assemblies:
        groups[a.species].append(a)
    reps = {}
    for sp, items in groups.items():
        reps[sp] = max(items, key=lambda x: (
            x.quality_score,
            x.callable_fraction,
            -x.n_fraction,
            assembly_level_score(x.assembly_level),
        ))
    return reps


# -----------------------------------------------------------------------------
# Sequence exports and background selection
# -----------------------------------------------------------------------------


def write_sequence_exports(root: Path, assemblies: List[AssemblySequence], reps: Dict[str, AssemblySequence],
                           ref_cds: str) -> None:
    sd = root / "sequence_recovery"
    sd.mkdir(parents=True, exist_ok=True)
    with open(sd / "all_assemblies_remnants.fasta", "w") as fh, \
         open(sd / "all_assemblies_reconstructed_exons.fasta", "w") as fe:
        for a in sorted(assemblies, key=lambda x: (x.species, x.accession)):
            aid = safe_name(a.species)
            fh.write(f">{aid}|{safe_name(a.accession)} species={a.species} assembly_status={a.assembly_status} source={a.source}\n{wrap_fasta(a.sequence)}\n")
            for exi, seq in sorted(a.exon_sequences.items()):
                fe.write(f">{aid}|{safe_name(a.accession)}|exon_{exi} species={a.species}\n{wrap_fasta(seq)}\n")
    with open(sd / "species_representatives_unaligned.fasta", "w") as fh:
        for sp, a in sorted(reps.items()):
            fh.write(f">{safe_name(sp)} species={sp} accession={a.accession} status={a.species_status} source={a.source}\n{wrap_fasta(a.sequence)}\n")
        fh.write(f">Sus_scrofa_REF species=Sus scrofa role=reference_outgroup\n{wrap_fasta(ref_cds)}\n")

    rows = []
    for a in assemblies:
        d = asdict(a)
        d.pop("sequence", None)
        d.pop("exon_sequences", None)
        rows.append(d)
    write_tsv(pd.DataFrame(rows), sd / "assembly_sequence_provenance.tsv")

    rep_rows = []
    for sp, a in reps.items():
        d = asdict(a)
        d.pop("sequence", None)
        d.pop("exon_sequences", None)
        d["selected_representative"] = True
        rep_rows.append(d)
    write_tsv(pd.DataFrame(rep_rows), sd / "species_representatives.tsv")


def tree_tip_norm_map(tree: Tree) -> Dict[str, List[Clade]]:
    m: Dict[str, List[Clade]] = defaultdict(list)
    for t in tree.get_terminals():
        m[norm_taxon(t.name or "")].append(t)
    return m


def standardized_tree(tree_path: Path, species_to_id: Dict[str, str]) -> Tuple[Tree, Dict[str, str]]:
    tree = Phylo.read(str(tree_path), "newick")
    by_norm = {norm_taxon(sp): sid for sp, sid in species_to_id.items()}
    matched: Dict[str, str] = {}
    for tip in tree.get_terminals():
        n = norm_taxon(tip.name or "")
        sid = by_norm.get(n)
        if sid:
            matched[sid] = str(tip.name)
            tip.name = sid
    return tree, matched


def make_unit_tree(tree: Tree) -> Tree:
    t = copy.deepcopy(tree)
    for c in t.find_clades(order="preorder"):
        if c is t.root:
            continue
        c.branch_length = 1.0
    return t


def topological_distance(tree: Tree, a: str, b: str) -> float:
    try:
        return float(make_unit_tree(tree).distance(a, b))
    except Exception:
        return float("inf")


def select_intact_background(reps: Dict[str, AssemblySequence], tree: Tree, focus_species: List[str],
                             min_callable: float, max_n: int, min_n: int,
                             fallback_min_callable: float = 0.70,
                             allow_conflicting_fallback: bool = True,
                             diagnostics_path: Optional[Path] = None) -> List[str]:
    """Select intact controls using explicit, auditable tiers.

    Tier 1 keeps the original conservative rule: species-level INTACT_HIGH or
    INTACT_LIKELY and high callable nucleotide reconstruction.

    If too few controls satisfy the requested threshold, the callable cutoff is
    stepped down only as far as ``fallback_min_callable`` while retaining the
    species-level intact requirement. This is justified because the refiner's
    intact call may come from a near-complete protein-to-genome model even when
    the reference-anchored BLAST nucleotide reconstruction is incomplete.

    Only if species-level intact controls are still insufficient may a
    CONFLICTING_ASSEMBLIES species enter as a fallback, and then only when the
    chosen representative assembly itself is strongly intact, has no supported
    disruptive event, no naive internal stop, and no net frame remainder.

    Every representative is written to a diagnostic table with the exact
    acceptance/rejection criteria. Strong/likely loss species are never used as
    backgrounds.
    """
    tree_tips = {t.name for t in tree.get_terminals()}
    fallback_min_callable = max(0.0, min(float(min_callable), float(fallback_min_callable)))

    diag_rows: List[Dict[str, Any]] = []
    for sp, r in sorted(reps.items()):
        sid = safe_name(sp)
        in_tree = sid in tree_tips
        strict_status = r.species_status in INTACT_SPECIES
        strong_rep_intact = r.assembly_status in {"INTACT_HIGH", "INTACT_LIKELY"}
        conflict_fallback_ok = (
            r.species_status == "CONFLICTING_ASSEMBLIES"
            and strong_rep_intact
            and r.supported_events == 0
            and r.internal_stop_count_naive == 0
            and r.frame_remainder == 0
        )
        diag_rows.append({
            "species": sp,
            "tree_id": sid,
            "in_tree": in_tree,
            "species_status": r.species_status,
            "assembly_status": r.assembly_status,
            "callable_fraction": r.callable_fraction,
            "reconstructed_fraction": r.reconstructed_fraction,
            "n_fraction": r.n_fraction,
            "internal_stop_count_naive": r.internal_stop_count_naive,
            "frame_remainder": r.frame_remainder,
            "supported_events": r.supported_events,
            "quality_score": r.quality_score,
            "strict_species_intact": strict_status,
            "conflicting_fallback_eligible": conflict_fallback_ok,
            "selected": False,
            "selection_tier": "",
            "rejection_reason": "",
        })

    # Log enough information to diagnose a genuinely sparse background set.
    status_counts = Counter(str(r.species_status) for r in reps.values())
    asm_counts = Counter(str(r.assembly_status) for r in reps.values())
    logging.info("Representative species-status counts: %s", dict(status_counts))
    logging.info("Representative assembly-status counts: %s", dict(asm_counts))

    strict_all = [sp for sp, r in reps.items()
                  if r.species_status in INTACT_SPECIES and safe_name(sp) in tree_tips]
    if strict_all:
        vals = sorted(reps[sp].callable_fraction for sp in strict_all)
        def q(frac: float) -> float:
            idx = min(len(vals)-1, max(0, int(round(frac*(len(vals)-1)))))
            return vals[idx]
        logging.info(
            "Species-level intact representatives in tree: %d; callable fraction min/q25/median/q75/max = %.3f/%.3f/%.3f/%.3f/%.3f",
            len(vals), vals[0], q(0.25), q(0.50), q(0.75), vals[-1]
        )
    else:
        logging.warning("No species-level INTACT_HIGH/INTACT_LIKELY representative matched the standardized tree")

    # Descending adaptive thresholds, never below fallback_min_callable.
    raw_thresholds = [float(min_callable), 0.80, 0.75, 0.70, float(fallback_min_callable)]
    thresholds: List[float] = []
    for th in raw_thresholds:
        if th > float(min_callable) + 1e-12 or th < fallback_min_callable - 1e-12:
            continue
        if not any(abs(th-x) < 1e-12 for x in thresholds):
            thresholds.append(th)
    thresholds.sort(reverse=True)

    pool: List[str] = []
    tier = ""
    used_threshold = float(min_callable)

    # Prefer species-level intact taxa. Lower only the sequence-completeness
    # threshold if required to reach the requested background count.
    best_strict: List[str] = []
    best_strict_threshold = float(min_callable)
    for th in thresholds:
        cand = [
            sp for sp, r in reps.items()
            if safe_name(sp) in tree_tips
            and r.species_status in INTACT_SPECIES
            and r.callable_fraction >= th
        ]
        cand.sort(key=lambda sp: reps[sp].quality_score, reverse=True)
        if len(cand) > len(best_strict):
            best_strict = cand
            best_strict_threshold = th
        if len(cand) >= min_n:
            pool = cand
            used_threshold = th
            tier = ("strict_species_intact"
                    if abs(th-float(min_callable)) < 1e-12
                    else f"species_intact_relaxed_callable_{th:.2f}")
            break

    if not pool:
        pool = list(best_strict)
        used_threshold = best_strict_threshold
        tier = (f"species_intact_best_available_{used_threshold:.2f}" if pool else "")

    # Only after exhausting species-level intact controls, admit carefully
    # constrained representatives from CONFLICTING_ASSEMBLIES.
    if len(pool) < min_n and allow_conflicting_fallback:
        conflict = [
            sp for sp, r in reps.items()
            if safe_name(sp) in tree_tips
            and sp not in pool
            and r.species_status == "CONFLICTING_ASSEMBLIES"
            and r.assembly_status in {"INTACT_HIGH", "INTACT_LIKELY"}
            and r.callable_fraction >= fallback_min_callable
            and r.supported_events == 0
            and r.internal_stop_count_naive == 0
            and r.frame_remainder == 0
        ]
        conflict.sort(key=lambda sp: reps[sp].quality_score, reverse=True)
        if conflict:
            pool.extend(conflict)
            tier = ((tier + "+") if tier else "") + "conflicting_species_strong_intact_assembly_fallback"
            logging.warning(
                "Using %d CONFLICTING_ASSEMBLIES taxon/taxa as fallback controls because their selected assembly is strongly intact; see background_candidate_diagnostics.tsv",
                len(conflict)
            )

    # Deduplicate while preserving quality ordering.
    pool = list(dict.fromkeys(pool))
    pool.sort(key=lambda sp: reps[sp].quality_score, reverse=True)

    if len(pool) < min_n:
        logging.warning(
            "Only %d acceptable intact background species remain after tiered selection; requested at least %d",
            len(pool), min_n
        )

    # Farthest-point sampling with the closest intact relative to every focal
    # lineage forced in first. This keeps local controls and broad phylogenetic spread.
    if len(pool) <= max_n:
        selected = list(pool)
    else:
        selected: List[str] = []
        for focal in focus_species:
            fid = safe_name(focal)
            if fid not in tree_tips:
                continue
            avail = [sp for sp in pool if safe_name(sp) in tree_tips]
            if not avail:
                continue
            nearest = min(avail, key=lambda sp: topological_distance(tree, fid, safe_name(sp)))
            if nearest not in selected:
                selected.append(nearest)
        if not selected and pool:
            selected.append(pool[0])

        while len(selected) < max_n:
            remaining = [sp for sp in pool if sp not in selected]
            if not remaining:
                break

            def diversity(sp: str) -> Tuple[float, float]:
                sid = safe_name(sp)
                ds = [topological_distance(tree, sid, safe_name(x))
                      for x in selected if safe_name(x) in tree_tips]
                d = min(ds) if ds else 0.0
                return d, reps[sp].quality_score

            selected.append(max(remaining, key=diversity))
        selected = selected[:max_n]

    selected_set = set(selected)
    for row in diag_rows:
        sp = row["species"]
        row["selected"] = sp in selected_set
        if sp in selected_set:
            if reps[sp].species_status in INTACT_SPECIES:
                row["selection_tier"] = (
                    "strict_species_intact"
                    if reps[sp].callable_fraction >= float(min_callable)
                    else f"species_intact_relaxed_callable_ge_{fallback_min_callable:.2f}"
                )
            else:
                row["selection_tier"] = "conflicting_species_strong_intact_assembly_fallback"
        else:
            reasons = []
            if not row["in_tree"]:
                reasons.append("not_in_standardized_tree")
            if row["species_status"] in LOSS_SPECIES_STRONG | LOSS_SPECIES_EXPLORATORY:
                reasons.append("loss_status_excluded")
            elif row["species_status"] not in INTACT_SPECIES and not row["conflicting_fallback_eligible"]:
                reasons.append("not_intact_background_status")
            if float(row["callable_fraction"]) < fallback_min_callable:
                reasons.append(f"callable_below_{fallback_min_callable:.2f}")
            if row["species_status"] == "CONFLICTING_ASSEMBLIES" and not row["conflicting_fallback_eligible"]:
                reasons.append("conflict_not_supported_by_clean_strong_intact_representative")
            row["rejection_reason"] = ";".join(reasons)

    if diagnostics_path is not None:
        diagnostics_path.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(pd.DataFrame(diag_rows), diagnostics_path)

    logging.info(
        "Selected %d intact background taxa using tier '%s' (requested callable %.2f; fallback floor %.2f)",
        len(selected), tier or "none", float(min_callable), fallback_min_callable
    )
    return selected


def write_macse_inputs(root: Path, reps: Dict[str, AssemblySequence], intact_background: List[str],
                       focus_species: List[str], ref_cds: str) -> Dict[str, Any]:
    ad = root / "alignment"
    ad.mkdir(parents=True, exist_ok=True)
    reliable = ad / "macse_reliable_intact.fasta"
    less = ad / "macse_less_reliable_loss.fasta"
    selected_meta = []
    ids: Set[str] = set()

    def unique_id(sp: str) -> str:
        base = safe_name(sp)
        x = base
        i = 2
        while x in ids:
            x = f"{base}_{i}"
            i += 1
        ids.add(x)
        return x

    species_to_id: Dict[str, str] = {}
    with open(reliable, "w") as fr:
        pid = "Sus_scrofa_REF"
        ids.add(pid)
        fr.write(f">{pid}\n{wrap_fasta(ref_cds)}\n")
        selected_meta.append({"taxon_id": pid, "species": "Sus scrofa", "role": "reference_outgroup",
                              "status": "INTACT_REFERENCE", "accession": "pig_reference", "source": "reference_exons.fa"})
        for sp in intact_background:
            r = reps[sp]
            sid = unique_id(sp)
            species_to_id[sp] = sid
            fr.write(f">{sid}\n{wrap_fasta(r.sequence)}\n")
            selected_meta.append({"taxon_id": sid, "species": sp, "role": "intact_background",
                                  "status": r.species_status, "accession": r.accession, "source": r.source})

    with open(less, "w") as fl:
        for sp in focus_species:
            if sp in species_to_id:  # e.g. cow is itself an intact focus and already reliable
                for row in selected_meta:
                    if row["species"] == sp:
                        row["role"] = "intact_background;focus"
                continue
            r = reps.get(sp)
            if not r:
                continue
            sid = unique_id(sp)
            species_to_id[sp] = sid
            fl.write(f">{sid}\n{wrap_fasta(r.sequence)}\n")
            selected_meta.append({"taxon_id": sid, "species": sp, "role": "loss_or_special_focus",
                                  "status": r.species_status, "accession": r.accession, "source": r.source})

    write_tsv(pd.DataFrame(selected_meta), ad / "selected_taxa.tsv")
    return {"reliable": reliable, "less": less, "species_to_id": species_to_id,
            "selected_meta": selected_meta, "pig_id": "Sus_scrofa_REF"}


# -----------------------------------------------------------------------------
# MACSE execution and codon-safe alignment perturbations
# -----------------------------------------------------------------------------


def resolve_macse(args) -> Tuple[List[str], str]:
    if args.macse_jar:
        jar = Path(args.macse_jar).expanduser().resolve()
        if not jar.exists():
            raise FileNotFoundError(f"MACSE jar not found: {jar}")
        java = exe(args.java, "java")
        if not java:
            raise RuntimeError("java not found; required for --macse-jar")
        return [java, f"-Xmx{args.java_memory}", "-jar", str(jar)], f"jar:{jar}"
    m = exe(args.macse, "macse")
    if not m:
        raise RuntimeError("MACSE not found. Install `macse` from Bioconda or provide --macse-jar.")
    return [m], m


def fasta_nonempty(path: Path) -> bool:
    return path.exists() and path.stat().st_size > 20 and any(True for _ in SeqIO.parse(str(path), "fasta"))


def run_macse(args, root: Path, inputs: Dict[str, Any]) -> Dict[str, Path]:
    base, macse_desc = resolve_macse(args)
    ad = root / "alignment"
    raw_nt = ad / "macse_raw_NT.fasta"
    raw_aa = ad / "macse_raw_AA.fasta"
    cmd = base + ["-prog", "alignSequences", "-seq", str(inputs["reliable"]),
                  "-out_NT", str(raw_nt), "-out_AA", str(raw_aa)]
    if fasta_nonempty(inputs["less"]):
        cmd += ["-seq_lr", str(inputs["less"]), "-fs_lr", str(args.fs_lr), "-stop_lr", str(args.stop_lr)]
    if not (args.resume and fasta_nonempty(raw_nt)):
        logging.info("Running MACSE alignSequences")
        run_cmd(cmd, check=True, log_path=ad / "macse_commands.log")
    else:
        logging.info("Reusing existing MACSE alignment: %s", raw_nt)

    out = {"baseline": raw_nt, "baseline_aa": raw_aa}

    # A second de-novo alignment with moderately stricter pseudogene penalties
    # directly tests whether small MACSE scoring changes alter the biological call.
    if not args.skip_macse_alt_cost and fasta_nonempty(inputs["less"]):
        alt_nt = ad / "macse_altcost_NT.fasta"
        alt_aa = ad / "macse_altcost_AA.fasta"
        cmd_alt = base + ["-prog", "alignSequences", "-seq", str(inputs["reliable"]),
                          "-out_NT", str(alt_nt), "-out_AA", str(alt_aa),
                          "-seq_lr", str(inputs["less"]),
                          "-fs_lr", str(args.alt_fs_lr), "-stop_lr", str(args.alt_stop_lr)]
        if not (args.resume and fasta_nonempty(alt_nt)):
            logging.info("Running MACSE alternate-cost alignment for sensitivity analysis")
            cp_alt = run_cmd(cmd_alt, check=False, log_path=ad / "macse_commands.log")
            if cp_alt.returncode != 0:
                logging.warning("MACSE alternate-cost alignment failed; continuing without it")
            elif fasta_nonempty(alt_nt):
                out["altcost"] = alt_nt
                out["altcost_aa"] = alt_aa
        else:
            out["altcost"] = alt_nt
            out["altcost_aa"] = alt_aa

    if not args.skip_macse_refine:
        ref_nt = ad / "macse_refined_NT.fasta"
        ref_aa = ad / "macse_refined_AA.fasta"
        cmd2 = base + ["-prog", "refineAlignment", "-align", str(raw_nt),
                       "-optim", "2", "-max_refine_iter", str(args.macse_refine_iter),
                       "-out_NT", str(ref_nt), "-out_AA", str(ref_aa),
                       "-seq", str(inputs["reliable"])]
        if fasta_nonempty(inputs["less"]):
            cmd2 += ["-seq_lr", str(inputs["less"]), "-fs_lr", str(args.fs_lr), "-stop_lr", str(args.stop_lr)]
        if not (args.resume and fasta_nonempty(ref_nt)):
            logging.info("Running MACSE refineAlignment for alignment sensitivity analysis")
            cp = run_cmd(cmd2, check=False, log_path=ad / "macse_commands.log")
            if cp.returncode != 0:
                logging.warning("MACSE refineAlignment failed; continuing with baseline alignment only")
            elif fasta_nonempty(ref_nt):
                out["refined"] = ref_nt
                out["refined_aa"] = ref_aa
        else:
            out["refined"] = ref_nt
            out["refined_aa"] = ref_aa
    write_json({"macse": macse_desc, "alignments": {k: str(v) for k, v in out.items()}}, ad / "macse_run.json")
    return out


def read_alignment(path: Path) -> Dict[str, str]:
    recs = {r.id: str(r.seq).upper().replace("U", "T").replace(".", "-") for r in SeqIO.parse(str(path), "fasta")}
    if not recs:
        raise ValueError(f"No sequences in alignment {path}")
    lengths = Counter(len(s) for s in recs.values())
    if len(lengths) != 1:
        mx = max(lengths)
        logging.warning("MACSE alignment has unequal sequence lengths %s; right-padding with gaps to %d", dict(lengths), mx)
        recs = {k: v + "-" * (mx-len(v)) for k, v in recs.items()}
    L = len(next(iter(recs.values())))
    if L % 3:
        pad = 3 - (L % 3)
        logging.warning("Alignment length %d not divisible by 3; padding %d columns", L, pad)
        recs = {k: v + "-" * pad for k, v in recs.items()}
    return recs


def is_valid_sense_codon(c: str) -> bool:
    return len(c) == 3 and set(c) <= DNA and c not in STOP_CODONS


def macse_alignment_qc(records: Dict[str, str], source_name: str) -> pd.DataFrame:
    rows = []
    L = len(next(iter(records.values())))
    ncod = L // 3
    for sid, seq in records.items():
        fs = seq.count("!")
        stops = 0
        terminal_stops = 0
        partial_gap = 0
        full_gap = 0
        ambiguous = 0
        sense = 0
        non_gap_indices = [i for i in range(ncod) if seq[3*i:3*i+3] != "---"]
        last_non_gap = non_gap_indices[-1] if non_gap_indices else -1
        for i in range(ncod):
            c = seq[3*i:3*i+3]
            if c == "---":
                full_gap += 1
            elif "!" in c:
                pass
            elif "-" in c:
                partial_gap += 1
            elif set(c) <= DNA:
                if c in STOP_CODONS:
                    if i == last_non_gap:
                        terminal_stops += 1
                    else:
                        stops += 1
                else:
                    sense += 1
            else:
                ambiguous += 1
        rows.append({
            "alignment_source": source_name, "taxon_id": sid, "aligned_nt": L,
            "aligned_codons": ncod, "frameshift_characters": fs,
            "internal_stop_codons": stops, "terminal_stop_codons": terminal_stops,
            "partial_gap_codons": partial_gap, "full_gap_codons": full_gap,
            "ambiguous_codons": ambiguous, "valid_sense_codons": sense,
            "valid_sense_fraction": sense / ncod if ncod else 0.0,
        })
    return pd.DataFrame(rows)


def alignment_problem_columns(records: Dict[str, str]) -> Tuple[Set[int], Dict[str, Set[int]]]:
    L = len(next(iter(records.values())))
    ncod = L // 3
    global_bad: Set[int] = set()
    per_seq: Dict[str, Set[int]] = defaultdict(set)
    for sid, seq in records.items():
        non_gap = [i for i in range(ncod) if seq[3*i:3*i+3] != "---"]
        last_non_gap = non_gap[-1] if non_gap else -1
        for i in range(ncod):
            c = seq[3*i:3*i+3]
            bad = False
            if "!" in c or ("-" in c and c != "---"):
                bad = True
            elif set(c) <= DNA and c in STOP_CODONS and i != last_non_gap:
                bad = True
            if bad:
                global_bad.add(i)
                per_seq[sid].add(i)
    return global_bad, per_seq


def build_paml_variant(records: Dict[str, str], occupancy: float, trim_flank: int,
                       strict_complete: bool = False) -> Tuple[Dict[str, str], Dict[str, Any]]:
    ids = list(records)
    L = len(next(iter(records.values())))
    ncod = L // 3
    bad, _ = alignment_problem_columns(records)
    trim: Set[int] = set()
    if trim_flank > 0:
        for i in bad:
            for j in range(max(0, i-trim_flank), min(ncod, i+trim_flank+1)):
                trim.add(j)

    keep: List[int] = []
    per_seq_valid = Counter()
    for i in range(ncod):
        if i in trim:
            continue
        valid = 0
        all_complete = True
        for sid in ids:
            c = records[sid][3*i:3*i+3]
            if is_valid_sense_codon(c):
                valid += 1
            else:
                all_complete = False
        frac = valid / len(ids) if ids else 0.0
        if strict_complete:
            if all_complete:
                keep.append(i)
        elif frac >= occupancy:
            keep.append(i)

    out: Dict[str, str] = {}
    for sid in ids:
        parts = []
        for i in keep:
            c = records[sid][3*i:3*i+3]
            if is_valid_sense_codon(c):
                parts.append(c)
                per_seq_valid[sid] += 1
            else:
                # PAML cannot accept stop codons and MACSE's ! is nonstandard.
                # Missing/frameshift/stop codons are represented as an unknown codon.
                parts.append("NNN")
        out[sid] = "".join(parts)

    stats = {
        "input_codons": ncod,
        "output_codons": len(keep),
        "occupancy_threshold": occupancy,
        "trim_flank_codons": trim_flank,
        "strict_complete": strict_complete,
        "problem_columns": len(bad),
        "trimmed_problem_neighborhood_columns": len(trim),
        "kept_column_indices_0based": keep,
        "valid_codons_per_taxon": dict(per_seq_valid),
    }
    return out, stats


def write_fasta_dict(records: Dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for sid, seq in records.items():
            fh.write(f">{sid}\n{wrap_fasta(seq)}\n")


def generate_alignment_variants(root: Path, macse_outputs: Dict[str, Path], min_codons: int) -> Tuple[Dict[str, Path], pd.DataFrame, pd.DataFrame]:
    vd = root / "alignment" / "variants"
    vd.mkdir(parents=True, exist_ok=True)
    variants: Dict[str, Path] = {}
    qc_frames = []
    summary = []

    sources = [("baseline", macse_outputs["baseline"])]
    if "refined" in macse_outputs:
        sources.append(("refined", macse_outputs["refined"]))
    if "altcost" in macse_outputs:
        sources.append(("altcost", macse_outputs["altcost"]))

    for source_name, path in sources:
        recs = read_alignment(path)
        qc_frames.append(macse_alignment_qc(recs, source_name))
        specs = [
            (f"{source_name}_mask50", 0.50, 0, False),
            (f"{source_name}_trim1_occ60", 0.60, 1, False),
        ]
        if source_name == "baseline":
            specs.extend([
                ("baseline_trim2_occ70", 0.70, 2, False),
                ("baseline_strict_complete", 1.00, 0, True),
            ])
        for name, occ, flank, strict in specs:
            v, st = build_paml_variant(recs, occ, flank, strict)
            st["variant"] = name
            st["source_alignment"] = source_name
            usable = st["output_codons"] >= min_codons
            st["usable_for_paml"] = usable
            path_out = vd / f"{name}.fasta"
            write_fasta_dict(v, path_out)
            if usable:
                variants[name] = path_out
            else:
                logging.warning("Skipping PAML variant %s: only %d codons remain (<%d)", name, st["output_codons"], min_codons)
            summary.append({k: v for k, v in st.items() if k not in {"kept_column_indices_0based", "valid_codons_per_taxon"}})
            write_json(st, vd / f"{name}.stats.json")

    qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()
    summ = pd.DataFrame(summary)
    write_tsv(qc, root / "alignment" / "macse_alignment_qc.tsv")
    write_tsv(summ, root / "alignment" / "alignment_variant_summary.tsv")
    return variants, qc, summ


# -----------------------------------------------------------------------------
# Tree pruning, outgroup attachment, and foreground branch labels
# -----------------------------------------------------------------------------


def prune_tree(tree: Tree, keep: Set[str]) -> Tree:
    t = copy.deepcopy(tree)
    for tip in list(t.get_terminals()):
        if tip.name not in keep:
            try:
                t.prune(tip)
            except Exception:
                pass
    return t


def add_pig_outgroup(tree: Tree, pig_id: str) -> Tree:
    t = copy.deepcopy(tree)
    if pig_id in {x.name for x in t.get_terminals()}:
        return t
    old = t.root
    pig = Clade(name=pig_id)
    # With clock=0 codeml treats the topology as unrooted.  A trifurcation at
    # the displayed root is therefore preferable to introducing an artificial
    # two-child root branch.  Attach pig as the third root branch when the
    # ruminant tree is bifurcating.
    if len(old.clades) >= 2:
        t.root = Clade(clades=[pig] + old.clades)
    else:
        old.branch_length = None
        t.root = Clade(clades=[pig, old])
    return t


def prune_with_pig_policy(base_tree: Tree, taxa: Set[str], pig_id: str) -> Tree:
    """Prune to taxa while preserving pig's supplied topology when available.

    Older pipeline versions assumed the input tree contained only Ruminantia and
    therefore attached the pig reference as a root outgroup. In the external
    background workflow the combined tree already contains ``Sus_scrofa_REF``
    at its correct Cetartiodactyla position. In that case we must preserve it.
    """
    base_tips = {x.name for x in base_tree.get_terminals()}
    if pig_id in base_tips:
        return prune_tree(base_tree, set(taxa))
    t = prune_tree(base_tree, set(taxa) - {pig_id})
    return add_pig_outgroup(t, pig_id)

def descendant_set(clade: Clade) -> Set[str]:
    return {t.name for t in clade.get_terminals() if t.name}


def find_exact_clade(tree: Tree, target_desc: Set[str]) -> Optional[Clade]:
    if not target_desc:
        return None
    matches = []
    for cl in tree.find_clades(order="preorder"):
        if descendant_set(cl) == target_desc:
            matches.append(cl)
    if not matches:
        return None
    return min(matches, key=lambda c: len(c.get_terminals()))


def paml_newick(tree: Tree, foreground_desc: Optional[Set[str]] = None) -> str:
    foreground_desc = foreground_desc or set()
    target_clade = find_exact_clade(tree, foreground_desc) if foreground_desc else None

    def rec(c: Clade) -> str:
        if c.clades:
            # Internal support/node labels from the species tree are omitted.
            # PAML branch labels are added explicitly below and are the only
            # internal labels needed for this analysis.
            s = "(" + ",".join(rec(x) for x in c.clades) + ")"
        else:
            s = safe_name(c.name or "taxon")
        if target_clade is c:
            s += " #1"
        return s
    return rec(tree.root) + ";\n"


def write_analysis_tree(base_tree: Tree, taxa: Set[str], pig_id: str, out: Path,
                        foreground_desc: Optional[Set[str]] = None) -> Tuple[Tree, Set[str]]:
    keep_requested = set(taxa)
    t = prune_with_pig_policy(base_tree, keep_requested, pig_id)
    present = {x.name for x in t.get_terminals()}
    if keep_requested - present:
        logging.debug("Tree missing taxa: %s", sorted(keep_requested-present))
    actual = {x.name for x in t.get_terminals()}
    fg = (foreground_desc or set()) & actual
    if foreground_desc and not find_exact_clade(t, fg):
        raise ValueError(f"Foreground descendants do not define one branch after pruning: {sorted(fg)}")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(paml_newick(t, fg if foreground_desc else None))
    return t, actual


# -----------------------------------------------------------------------------
# PAML control files, execution, and output parsing
# -----------------------------------------------------------------------------


def write_phylip(records: Dict[str, str], taxa: Sequence[str], path: Path) -> int:
    use = [(x, records[x]) for x in taxa if x in records]
    if len(use) < 3:
        raise ValueError("PAML needs at least 3 taxa")
    lens = {len(seq) for _, seq in use}
    if len(lens) != 1:
        raise ValueError(f"Unequal sequence lengths for PAML: {lens}")
    L = next(iter(lens))
    if L % 3:
        raise ValueError(f"PAML alignment length is not divisible by 3: {L}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f" {len(use)} {L}\n")
        for sid, seq in use:
            fh.write(f"{sid}  {seq}\n")
    return L // 3


def codeml_ctl_text(seqfile: str, treefile: str, outfile: str, model: int, nssites: int,
                    codon_freq: int, fix_omega: int = 0, omega: float = 0.4,
                    cleandata: int = 0, rate_ancestor: int = 0, ncatg: int = 8) -> str:
    return textwrap.dedent(f"""\
        seqfile = {seqfile}
        treefile = {treefile}
        outfile = {outfile}
        noisy = 3
        verbose = 1
        runmode = 0
        seqtype = 1
        CodonFreq = {codon_freq}
        clock = 0
        model = {model}
        NSsites = {nssites}
        icode = 0
        Mgene = 0
        fix_kappa = 0
        kappa = 2
        fix_omega = {fix_omega}
        omega = {omega}
        fix_alpha = 1
        alpha = 0
        Malpha = 0
        ncatG = {ncatg}
        getSE = 0
        RateAncestor = {rate_ancestor}
        Small_Diff = 5e-7
        cleandata = {cleandata}
        fix_blength = 0
        method = 0
        """)


def parse_lnl(text: str) -> Tuple[float, int]:
    matches = re.findall(r"lnL\([^\)]*np:\s*(\d+)\)\s*:\s*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", text)
    if not matches:
        matches = re.findall(r"lnL\([^\)]*\)\s*:\s*([+-]?\d+(?:\.\d+)?)", text)
        if matches:
            return float(matches[-1]), 0
        return float("nan"), 0
    np_, ll = matches[-1]
    return float(ll), int(np_)


def extract_newick_after_label(text: str, label: str) -> Optional[str]:
    lines = text.splitlines()
    for i, line in enumerate(lines):
        if label.lower() in line.lower():
            buff = []
            # Sometimes the Newick begins on the same line after the colon.
            after = line.split(":", 1)[1].strip() if ":" in line else ""
            if "(" in after:
                buff.append(after[after.find("("):])
            for j in range(i+1, min(len(lines), i+30)):
                if not buff and "(" not in lines[j]:
                    continue
                x = lines[j]
                if not buff and "(" in x:
                    x = x[x.find("("):]
                buff.append(x.strip())
                if ";" in x:
                    break
            if buff:
                s = "".join(buff)
                if ";" in s:
                    return s[:s.find(";")+1]
    return None


def parse_newick_string(s: Optional[str]) -> Optional[Tree]:
    if not s:
        return None
    try:
        return Phylo.read(io.StringIO(s), "newick")
    except Exception:
        return None


def branch_table_from_dn_ds(dn_tree: Optional[Tree], ds_tree: Optional[Tree]) -> pd.DataFrame:
    if not dn_tree or not ds_tree:
        return pd.DataFrame()
    dn = {}
    for c in dn_tree.find_clades(order="preorder"):
        if c is dn_tree.root:
            continue
        key = tuple(sorted(t.name for t in c.get_terminals() if t.name))
        dn[key] = fnum(c.branch_length, float("nan"))
    ds = {}
    for c in ds_tree.find_clades(order="preorder"):
        if c is ds_tree.root:
            continue
        key = tuple(sorted(t.name for t in c.get_terminals() if t.name))
        ds[key] = fnum(c.branch_length, float("nan"))
    rows = []
    for key in sorted(set(dn) | set(ds), key=lambda x: (len(x), x)):
        dnv = dn.get(key, float("nan")); dsv = ds.get(key, float("nan"))
        w = dnv / dsv if math.isfinite(dnv) and math.isfinite(dsv) and dsv > 0 else float("inf") if dsv == 0 and dnv > 0 else float("nan")
        label = key[0] if len(key) == 1 else f"MRCA[{len(key)}]:{key[0]}..{key[-1]}"
        rows.append({"descendant_key": "|".join(key), "n_descendants": len(key), "branch_label": label,
                     "dN": dnv, "dS": dsv, "omega": w})
    return pd.DataFrame(rows)


def parse_omega_vector(text: str) -> List[float]:
    pats = [
        r"w\s*\(dN/dS\)\s*for\s*branches\s*:\s*([^\n]+)",
        r"omega\s*\(dN/dS\)\s*=\s*([^\n]+)",
    ]
    for p in pats:
        m = re.search(p, text, re.I)
        if m:
            vals = []
            for x in re.findall(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", m.group(1)):
                try: vals.append(float(x))
                except Exception: pass
            if vals:
                return vals
    return []


def parse_beb(text: str) -> List[Dict[str, Any]]:
    idx = text.find("Bayes Empirical Bayes (BEB) analysis")
    if idx < 0:
        return []
    block = text[idx:idx+12000]
    rows = []
    for line in block.splitlines():
        m = re.match(r"\s*(\d+)\s+([A-Z*])\s+([01](?:\.\d+)?)\s*(\*\*|\*)?", line)
        if m:
            rows.append({"site": int(m.group(1)), "aa": m.group(2), "posterior": float(m.group(3)),
                         "stars": m.group(4) or ""})
    return rows


def parse_codeml_output(mlc: Path) -> Dict[str, Any]:
    if not mlc.exists():
        return {"status": "missing_output"}
    text = mlc.read_text(errors="replace")
    ll, np_ = parse_lnl(text)
    dn_newick = extract_newick_after_label(text, "dN tree")
    ds_newick = extract_newick_after_label(text, "dS tree")
    dn_tree = parse_newick_string(dn_newick)
    ds_tree = parse_newick_string(ds_newick)
    btab = branch_table_from_dn_ds(dn_tree, ds_tree)
    return {
        "status": "ok" if math.isfinite(ll) else "parse_incomplete",
        "lnL": ll,
        "np": np_,
        "AIC": aic(ll, np_),
        "omega_vector": parse_omega_vector(text),
        "dN_tree_newick": dn_newick or "",
        "dS_tree_newick": ds_newick or "",
        "branch_table": btab.to_dict(orient="records"),
        "beb_sites": parse_beb(text),
    }


def run_codeml_job(codeml: str, job_dir: Path, seq_records: Dict[str, str], taxa: Sequence[str],
                   tree: Tree, foreground_desc: Optional[Set[str]], model_name: str,
                   model: int, nssites: int, codon_freq: int, fix_omega: int = 0,
                   omega: float = 0.4, cleandata: int = 0, resume: bool = True,
                   rate_ancestor: int = 0, ncatg: int = 8) -> Dict[str, Any]:
    job_dir.mkdir(parents=True, exist_ok=True)
    mlc = job_dir / "mlc"
    if resume and mlc.exists() and mlc.stat().st_size > 100:
        parsed = parse_codeml_output(mlc)
        if math.isfinite(fnum(parsed.get("lnL"))):
            parsed.update({"model_name": model_name, "job_dir": str(job_dir), "codon_freq": codon_freq,
                           "returncode": 0, "resumed": True})
            return parsed

    seqfile = job_dir / "alignment.phy"
    treefile = job_dir / "tree.nwk"
    ctl = job_dir / "codeml.ctl"
    codons = write_phylip(seq_records, taxa, seqfile)
    treefile.write_text(paml_newick(tree, foreground_desc))
    ctl.write_text(codeml_ctl_text(seqfile.name, treefile.name, mlc.name, model, nssites,
                                   codon_freq, fix_omega, omega, cleandata, rate_ancestor, ncatg))
    cp = run_cmd([codeml, ctl.name], cwd=job_dir, check=False, log_path=job_dir / "run.log")
    parsed = parse_codeml_output(mlc)
    parsed.update({"model_name": model_name, "job_dir": str(job_dir), "codon_freq": codon_freq,
                   "returncode": cp.returncode, "resumed": False, "alignment_codons": codons})
    if cp.returncode != 0:
        parsed["status"] = "codeml_failed"
        parsed["stderr_tail"] = (cp.stderr or b"").decode(errors="replace")[-3000:]
    write_json({k: v for k, v in parsed.items() if k != "branch_table"}, job_dir / "parsed_summary.json")
    if parsed.get("branch_table"):
        write_tsv(pd.DataFrame(parsed["branch_table"]), job_dir / "branch_dnds.tsv")
    if parsed.get("beb_sites"):
        write_tsv(pd.DataFrame(parsed["beb_sites"]), job_dir / "beb_sites.tsv")
    return parsed


def branch_metric(parsed: Dict[str, Any], target_desc: Set[str]) -> Dict[str, float]:
    key = "|".join(sorted(target_desc))
    for r in parsed.get("branch_table", []) or []:
        if str(r.get("descendant_key")) == key:
            return {"dN": fnum(r.get("dN")), "dS": fnum(r.get("dS")), "omega": fnum(r.get("omega"))}
    return {"dN": float("nan"), "dS": float("nan"), "omega": float("nan")}


# -----------------------------------------------------------------------------
# Foreground definitions
# -----------------------------------------------------------------------------


@dataclass
class Foreground:
    foreground_id: str
    kind: str
    species_label: str
    descendant_species: List[str]
    descendant_ids: List[str]
    source: str
    species_status: str


def build_foregrounds(upstream: Dict[str, Any], reps: Dict[str, AssemblySequence], species_to_id: Dict[str, str],
                      include_possible: bool, always_test: Sequence[str]) -> List[Foreground]:
    allowed_loss = set(LOSS_SPECIES_STRONG)
    if include_possible:
        allowed_loss |= LOSS_SPECIES_EXPLORATORY
    out: List[Foreground] = []
    seen: Set[Tuple[str, Tuple[str, ...]]] = set()

    for sp, r in sorted(reps.items()):
        if r.species_status in allowed_loss:
            sid = species_to_id.get(sp)
            if not sid:
                continue
            key = ("terminal", (sid,))
            if key not in seen:
                out.append(Foreground(f"terminal_{safe_name(sp)}", "terminal", sp, [sp], [sid],
                                      "refined_species_loss_call", r.species_status))
                seen.add(key)

    # Always test named taxa such as Bos taurus even if their refined status is intact/unresolved.
    by_norm = {norm_taxon(sp): sp for sp in reps}
    for requested in always_test:
        sp = by_norm.get(norm_taxon(requested))
        if not sp or sp not in species_to_id:
            continue
        sid = species_to_id[sp]
        key = ("terminal", (sid,))
        if key not in seen:
            out.append(Foreground(f"terminal_{safe_name(sp)}", "terminal", sp, [sp], [sid],
                                  "explicit_focus", reps[sp].species_status))
            seen.add(key)

    # Add inferred shared loss-origin branches from the refiner where at least two
    # selected descendants remain available.
    origins = upstream.get("origins", pd.DataFrame())
    if origins is not None and not origins.empty:
        for i, row in origins.iterrows():
            desc_sp = parse_jsonish_list(row.get("descendants", ""))
            avail_sp = [sp for sp in desc_sp if sp in species_to_id]
            ids = sorted({species_to_id[sp] for sp in avail_sp})
            if len(ids) < 2:
                continue
            key = ("clade", tuple(ids))
            if key in seen:
                continue
            label = str(row.get("node", f"loss_origin_{i+1}"))
            out.append(Foreground(f"loss_origin_{i+1}_{safe_name(label)[:50]}", "clade", label,
                                  avail_sp, ids, "refiner_fitch_loss_origin", "INFERRED_LOSS_ORIGIN"))
            seen.add(key)
    return out


# -----------------------------------------------------------------------------
# PAML orchestration
# -----------------------------------------------------------------------------


def subset_records(records: Dict[str, str], taxa: Set[str]) -> Dict[str, str]:
    return {k: v for k, v in records.items() if k in taxa}


def read_variant_records(path: Path) -> Dict[str, str]:
    return {r.id: str(r.seq).upper() for r in SeqIO.parse(str(path), "fasta")}


def prepare_foreground_tree(base_tree: Tree, records: Dict[str, str], bg_ids: List[str], pig_id: str,
                            fg: Foreground) -> Tuple[Tree, List[str], Set[str]]:
    fg_ids = set(fg.descendant_ids) & set(records)
    taxa = set(bg_ids) | fg_ids | {pig_id}
    taxa &= set(records)
    if len(fg_ids) == 0:
        raise ValueError("Foreground absent from alignment")
    t = prune_with_pig_policy(base_tree, taxa, pig_id)
    actual = {x.name for x in t.get_terminals()} & set(records)
    fg_actual = fg_ids & actual
    if not find_exact_clade(t, fg_actual):
        raise ValueError(f"Foreground is not monophyletic/identifiable after pruning: {sorted(fg_actual)}")
    taxa_order = [x.name for x in t.get_terminals() if x.name in actual]
    return t, taxa_order, fg_actual


def model_job_specs() -> List[Tuple[str, int, int, int, float]]:
    # name, model, NSsites, fix_omega, omega
    return [
        ("branch_one_ratio", 0, 0, 0, 0.4),
        ("branch_two_ratio", 2, 0, 0, 0.8),
        ("branch_foreground_neutral", 2, 0, 1, 1.0),
        ("branch_site_alt", 2, 2, 0, 1.5),
        ("branch_site_null", 2, 2, 1, 1.0),
    ]


def select_run_matrix(variants: Sequence[str], codon_freqs: Sequence[int], depth: str) -> List[Tuple[str, int, str]]:
    """Return (variant, codon_freq, scope) where scope is full or reduced.

    full = all five foreground models.
    reduced = two-ratio + neutral + branch-site alt/null; one-ratio can be borrowed
              conceptually but is still cheap enough, so currently reduced omits only
              redundant site/global analyses.
    """
    primary = "baseline_mask50" if "baseline_mask50" in variants else list(variants)[0]
    out = []
    if depth == "exhaustive":
        for v in variants:
            for cf in codon_freqs:
                out.append((v, cf, "full"))
        return out
    # extensive: primary under all codon frequencies, all alignment variants under F3x4.
    for cf in codon_freqs:
        out.append((primary, cf, "full"))
    for v in variants:
        if v != primary:
            out.append((v, 2, "full"))
    # deduplicate
    seen = set(); uniq = []
    for x in out:
        if (x[0], x[1]) not in seen:
            uniq.append(x); seen.add((x[0], x[1]))
    return uniq


def run_paml_panel(args, root: Path, variants: Dict[str, Path], base_tree: Tree,
                   bg_species: List[str], species_to_id: Dict[str, str], pig_id: str,
                   foregrounds: List[Foreground]) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    codeml = exe(args.codeml, "codeml")
    if not codeml:
        raise RuntimeError("codeml not found. Install PAML or provide --codeml.")
    bg_ids = [species_to_id[x] for x in bg_species if x in species_to_id]
    codon_freqs = [int(x) for x in str(args.codon_freqs).split(",") if str(x).strip()]
    matrix = select_run_matrix(list(variants), codon_freqs, args.analysis_depth)
    paml_root = root / "paml"
    jobs = []

    # Global free-ratio models on all selected taxa, useful for descriptive dN/dS trees.
    global_specs = []
    primary_variant = "baseline_mask50" if "baseline_mask50" in variants else list(variants)[0]
    global_freqs = codon_freqs if args.analysis_depth in {"extensive", "exhaustive"} else [2]
    for cf in global_freqs:
        global_specs.append((primary_variant, cf))

    global_results = []
    for v, cf in global_specs:
        recs = read_variant_records(variants[v])
        taxa = set(recs)
        t = prune_with_pig_policy(base_tree, taxa, pig_id)
        actual = {x.name for x in t.get_terminals()} & set(recs)
        order = [x.name for x in t.get_terminals() if x.name in actual]
        for name, model, ns in [("global_M0", 0, 0), ("global_free_ratio", 1, 0)]:
            d = paml_root / v / f"F{cf}" / name
            parsed = run_codeml_job(codeml, d, recs, order, t, None, name, model, ns, cf,
                                    resume=args.resume, cleandata=0)
            global_results.append({
                "variant": v, "codon_freq": cf, "model": name,
                "lnL": parsed.get("lnL"), "np": parsed.get("np"), "AIC": parsed.get("AIC"),
                "status": parsed.get("status"), "job_dir": parsed.get("job_dir"),
            })
            if name == "global_free_ratio" and parsed.get("branch_table"):
                bdf = pd.DataFrame(parsed["branch_table"])
                bdf["variant"] = v; bdf["codon_freq"] = cf
                write_tsv(bdf, paml_root / v / f"F{cf}" / "free_ratio_branches.tsv")

    # Intact-only site models. These answer whether GPRC6A as a coding gene has
    # site classes with omega > 1, but are not lineage-specific tests.
    site_results = []
    for cf in global_freqs:
        recs_all = read_variant_records(variants[primary_variant])
        intact_taxa = set(bg_ids) | {pig_id}
        recs = subset_records(recs_all, intact_taxa)
        t = prune_with_pig_policy(base_tree, set(recs), pig_id)
        actual = {x.name for x in t.get_terminals()} & set(recs)
        order = [x.name for x in t.get_terminals() if x.name in actual]
        if len(order) >= 4:
            models = [
                ("site_M1a", 0, 1, 0, 0.4),
                ("site_M2a", 0, 2, 0, 1.5),
                ("site_M7", 0, 7, 0, 0.4),
                ("site_M8", 0, 8, 0, 1.5),
            ]
            parsed_by = {}
            for name, model, ns, fixw, w in models:
                d = paml_root / primary_variant / f"F{cf}" / "intact_site_models" / name
                p = run_codeml_job(codeml, d, recs, order, t, None, name, model, ns, cf,
                                   fixw, w, 0, args.resume, ncatg=10)
                parsed_by[name] = p
                site_results.append({"variant": primary_variant, "codon_freq": cf, "model": name,
                                     "lnL": p.get("lnL"), "np": p.get("np"), "AIC": p.get("AIC"),
                                     "BEB_sites_ge_0.95": sum(fnum(x.get("posterior"), 0) >= 0.95 for x in p.get("beb_sites", []) or []),
                                     "job_dir": p.get("job_dir"), "status": p.get("status")})
            for alt, null, label in [("site_M2a", "site_M1a", "M1a_vs_M2a"), ("site_M8", "site_M7", "M7_vs_M8")]:
                a, n = parsed_by.get(alt, {}), parsed_by.get(null, {})
                lr = max(0.0, 2*(fnum(a.get("lnL"), float("nan"))-fnum(n.get("lnL"), float("nan"))))
                p = chi2_sf(lr, 2) if math.isfinite(lr) else float("nan")
                site_results.append({"variant": primary_variant, "codon_freq": cf, "model": label,
                                     "LRT": lr, "p": p, "status": "comparison"})

    # Foreground jobs are embarrassingly parallel because codeml itself is single-threaded.
    def one_foreground_setting(fg: Foreground, v: str, cf: int) -> List[Dict[str, Any]]:
        recs = read_variant_records(variants[v])
        t, taxa_order, fg_ids = prepare_foreground_tree(base_tree, recs, bg_ids, pig_id, fg)
        if len(taxa_order) < 4:
            return [{"foreground_id": fg.foreground_id, "variant": v, "codon_freq": cf,
                     "status": "too_few_taxa"}]
        # Foreground codon information gate.
        informative = []
        for sid in fg_ids:
            seq = recs.get(sid, "")
            informative.append(sum(is_valid_sense_codon(seq[i:i+3]) for i in range(0, len(seq), 3)))
        min_info = min(informative) if informative else 0
        if min_info < args.min_foreground_codons:
            return [{"foreground_id": fg.foreground_id, "variant": v, "codon_freq": cf,
                     "status": "too_few_foreground_codons", "foreground_min_informative_codons": min_info}]

        parsed_by = {}
        rows = []
        for name, model, ns, fixw, w in model_job_specs():
            d = paml_root / v / f"F{cf}" / "foregrounds" / safe_name(fg.foreground_id) / name
            p = run_codeml_job(codeml, d, recs, taxa_order, t, fg_ids, name, model, ns, cf,
                               fixw, w, 0, args.resume, ncatg=8)
            parsed_by[name] = p
            metric = branch_metric(p, fg_ids) if name.startswith("branch_") and ns == 0 else {"dN": float("nan"), "dS": float("nan"), "omega": float("nan")}
            rows.append({
                "foreground_id": fg.foreground_id, "foreground_kind": fg.kind,
                "foreground_label": fg.species_label, "foreground_species_status": fg.species_status,
                "foreground_descendants": ";".join(fg.descendant_species),
                "foreground_ids": ";".join(sorted(fg_ids)), "variant": v, "codon_freq": cf,
                "model": name, "lnL": p.get("lnL"), "np": p.get("np"), "AIC": p.get("AIC"),
                "foreground_dN": metric["dN"], "foreground_dS": metric["dS"], "foreground_omega": metric["omega"],
                "BEB_sites_ge_0.95": sum(fnum(x.get("posterior"), 0) >= 0.95 for x in p.get("beb_sites", []) or []),
                "BEB_sites_ge_0.99": sum(fnum(x.get("posterior"), 0) >= 0.99 for x in p.get("beb_sites", []) or []),
                "foreground_min_informative_codons": min_info, "status": p.get("status"), "job_dir": p.get("job_dir"),
            })

        one = parsed_by.get("branch_one_ratio", {})
        two = parsed_by.get("branch_two_ratio", {})
        neu = parsed_by.get("branch_foreground_neutral", {})
        bsa = parsed_by.get("branch_site_alt", {})
        bsn = parsed_by.get("branch_site_null", {})
        metric = branch_metric(two, fg_ids)
        omega_vec = [fnum(x) for x in (two.get("omega_vector", []) or []) if math.isfinite(fnum(x))]
        background_omega = omega_vec[0] if len(omega_vec) >= 2 else float("nan")
        foreground_omega_vector = omega_vec[-1] if len(omega_vec) >= 2 else float("nan")
        lr_shift = max(0.0, 2*(fnum(two.get("lnL"), float("nan"))-fnum(one.get("lnL"), float("nan"))))
        lr_neutral = max(0.0, 2*(fnum(two.get("lnL"), float("nan"))-fnum(neu.get("lnL"), float("nan"))))
        lr_bs = max(0.0, 2*(fnum(bsa.get("lnL"), float("nan"))-fnum(bsn.get("lnL"), float("nan"))))
        p_bs_chi1 = chi2_sf(lr_bs, 1) if math.isfinite(lr_bs) else float("nan")
        p_bs_mix = 0.5*p_bs_chi1 if math.isfinite(p_bs_chi1) and lr_bs > 0 else 1.0
        rows.append({
            "foreground_id": fg.foreground_id, "foreground_kind": fg.kind,
            "foreground_label": fg.species_label, "foreground_species_status": fg.species_status,
            "foreground_descendants": ";".join(fg.descendant_species), "foreground_ids": ";".join(sorted(fg_ids)),
            "variant": v, "codon_freq": cf, "model": "COMPARISONS",
            "foreground_dN": metric["dN"], "foreground_dS": metric["dS"], "foreground_omega": metric["omega"],
            "background_omega": background_omega, "foreground_omega_vector": foreground_omega_vector,
            "omega_shift_from_background": metric["omega"] - background_omega if math.isfinite(metric["omega"]) and math.isfinite(background_omega) else float("nan"),
            "branch_shift_LRT": lr_shift, "branch_shift_p": chi2_sf(lr_shift, 1) if math.isfinite(lr_shift) else float("nan"),
            "foreground_vs_neutral_LRT": lr_neutral, "foreground_vs_neutral_p": chi2_sf(lr_neutral, 1) if math.isfinite(lr_neutral) else float("nan"),
            "branch_site_LRT": lr_bs, "branch_site_p_chi1": p_bs_chi1, "branch_site_p_mixture": p_bs_mix,
            "BEB_sites_ge_0.95": sum(fnum(x.get("posterior"), 0) >= 0.95 for x in bsa.get("beb_sites", []) or []),
            "BEB_sites_ge_0.99": sum(fnum(x.get("posterior"), 0) >= 0.99 for x in bsa.get("beb_sites", []) or []),
            "foreground_min_informative_codons": min_info,
            "status": "ok" if all(math.isfinite(fnum(x.get("lnL"))) for x in [one, two, neu, bsa, bsn]) else "incomplete_models",
        })
        return rows

    tasks = []
    with ThreadPoolExecutor(max_workers=max(1, args.paml_workers)) as pool:
        for fg in foregrounds:
            for v, cf, _scope in matrix:
                tasks.append(pool.submit(one_foreground_setting, fg, v, cf))
        foreground_rows = []
        done = 0
        for fut in as_completed(tasks):
            done += 1
            try:
                foreground_rows.extend(fut.result())
            except Exception as e:
                logging.exception("Foreground PAML job failed: %s", e)
            if done % 10 == 0 or done == len(tasks):
                logging.info("PAML foreground settings completed: %d/%d", done, len(tasks))

    fgdf = pd.DataFrame(foreground_rows)
    comp = fgdf[fgdf.model == "COMPARISONS"].copy() if not fgdf.empty else pd.DataFrame()
    if not comp.empty:
        # FDR correction within each alignment/codon-frequency sensitivity setting.
        # Report both conventional 50:50-mixture and conservative chi-square(1)
        # corrections.  The conservative q value is used for the final call.
        comp["branch_site_q"] = float("nan")
        comp["branch_site_q_mixture"] = float("nan")
        for (v, cf), idx in comp.groupby(["variant", "codon_freq"]).groups.items():
            q_cons = bh_fdr([fnum(comp.loc[i, "branch_site_p_chi1"]) for i in idx])
            q_mix = bh_fdr([fnum(comp.loc[i, "branch_site_p_mixture"]) for i in idx])
            for i, qq, qm in zip(idx, q_cons, q_mix):
                comp.loc[i, "branch_site_q"] = qq
                comp.loc[i, "branch_site_q_mixture"] = qm
        # write q values back into fgdf comparison rows
        for i, row in comp.iterrows():
            mask = ((fgdf.foreground_id == row.foreground_id) & (fgdf.variant == row.variant) &
                    (fgdf.codon_freq.astype(str) == str(row.codon_freq)) & (fgdf.model == "COMPARISONS"))
            fgdf.loc[mask, "branch_site_q"] = row.branch_site_q
            fgdf.loc[mask, "branch_site_q_mixture"] = row.branch_site_q_mixture

    return pd.DataFrame(global_results), pd.DataFrame(site_results), fgdf, comp


# -----------------------------------------------------------------------------
# Robustness-based evolutionary interpretation
# -----------------------------------------------------------------------------


def setting_call(row: pd.Series) -> str:
    if str(row.get("status")) != "ok":
        return "INSUFFICIENT"
    w = fnum(row.get("foreground_omega"))
    bg = fnum(row.get("background_omega"))
    ds = fnum(row.get("foreground_dS"))
    q = fnum(row.get("branch_site_q"), 1.0)
    p_neu = fnum(row.get("foreground_vs_neutral_p"), 1.0)
    p_shift = fnum(row.get("branch_shift_p"), 1.0)
    beb = inum(row.get("BEB_sites_ge_0.95"), 0)
    # A significant branch-site test is a distinct likelihood result and does
    # not require a stable branch-average dN/dS estimate.  Assess it before
    # rejecting settings with very small dS in the branch-wide model.
    if q < 0.05 and beb > 0:
        return "POSITIVE_SELECTION"
    if not math.isfinite(w):
        return "INSUFFICIENT"
    # Very small dS produces a division-by-near-zero omega; very high dS can
    # indicate synonymous saturation. Neither is used for a branch-wide call.
    if math.isfinite(ds) and ds < 1e-4:
        return "LOW_DS_UNSTABLE"
    if math.isfinite(ds) and ds > 3.0:
        return "HIGH_DS_SATURATED"
    if w > 1.0 and p_neu < 0.05:
        return "POSITIVE_SELECTION_BRANCH_AVERAGE_CANDIDATE"
    # Relaxation is stronger evidence than simply observing omega near one:
    # the intact background should be constrained, the foreground should move
    # toward one, and the two-ratio model should improve on one-ratio.
    if (math.isfinite(bg) and bg < 0.70 and 0.70 <= w <= 1.30 and
            w >= bg + 0.20 and p_neu >= 0.05 and p_shift < 0.05):
        return "RELAXED_TOWARD_NEUTRALITY"
    if 0.70 <= w <= 1.30 and p_neu >= 0.05:
        return "NEAR_NEUTRAL_WEAK_EVIDENCE"
    if w < 0.70 and p_neu < 0.05:
        return "PURIFYING_CONSTRAINT"
    if w < 1.0 and p_shift < 0.05:
        return "RATE_SHIFT_BELOW_ONE"
    return "AMBIGUOUS"


def summarize_foregrounds(comp: pd.DataFrame, foregrounds: List[Foreground]) -> pd.DataFrame:
    if comp.empty:
        return pd.DataFrame()
    comp = comp.copy()
    comp["setting_call"] = comp.apply(setting_call, axis=1)
    rows = []
    fg_meta = {x.foreground_id: x for x in foregrounds}
    for fid, g in comp.groupby("foreground_id"):
        valid = g[~g.setting_call.isin(["INSUFFICIENT", "LOW_DS_UNSTABLE", "HIGH_DS_SATURATED"])]
        calls = Counter(valid.setting_call.astype(str))
        top_call, top_n = calls.most_common(1)[0] if calls else ("INSUFFICIENT", 0)
        stability = top_n / len(valid) if len(valid) else 0.0
        omegas = [fnum(x) for x in valid.foreground_omega if math.isfinite(fnum(x))]
        bs_sig = sum(fnum(x, 1.0) < 0.05 for x in valid.branch_site_q) if len(valid) else 0
        near_neutral = sum(setting_call(r) == "RELAXED_TOWARD_NEUTRALITY" for _, r in valid.iterrows()) if len(valid) else 0
        pos = sum(setting_call(r).startswith("POSITIVE_SELECTION") for _, r in valid.iterrows()) if len(valid) else 0
        pur = sum(setting_call(r) == "PURIFYING_CONSTRAINT" for _, r in valid.iterrows()) if len(valid) else 0
        crosses_one = bool(omegas and min(omegas) < 1.0 < max(omegas))
        omega_dependent_calls = {"RELAXED_TOWARD_NEUTRALITY", "NEAR_NEUTRAL_WEAK_EVIDENCE",
                                 "PURIFYING_CONSTRAINT", "RATE_SHIFT_BELOW_ONE",
                                 "POSITIVE_SELECTION_BRANCH_AVERAGE_CANDIDATE"}
        if stability < 0.70 or (crosses_one and top_call in omega_dependent_calls):
            consensus = "ALIGNMENT_OR_MODEL_SENSITIVE"
        else:
            consensus = top_call
        meta = fg_meta.get(fid)
        rows.append({
            "foreground_id": fid,
            "foreground_kind": meta.kind if meta else g.foreground_kind.iloc[0],
            "foreground_label": meta.species_label if meta else g.foreground_label.iloc[0],
            "foreground_descendants": ";".join(meta.descendant_species) if meta else g.foreground_descendants.iloc[0],
            "species_status": meta.species_status if meta else g.foreground_species_status.iloc[0],
            "valid_sensitivity_settings": len(valid),
            "consensus_interpretation": consensus,
            "dominant_raw_call": top_call,
            "dominant_call_fraction": stability,
            "omega_median": float(pd.Series(omegas).median()) if omegas else float("nan"),
            "omega_min": min(omegas) if omegas else float("nan"),
            "omega_max": max(omegas) if omegas else float("nan"),
            "omega_crosses_1_across_settings": crosses_one,
            "branch_site_significant_settings": bs_sig,
            "positive_selection_call_settings": pos,
            "relaxed_neutral_call_settings": near_neutral,
            "purifying_call_settings": pur,
            "alignment_model_sensitive": consensus == "ALIGNMENT_OR_MODEL_SENSITIVE",
        })
    return pd.DataFrame(rows).sort_values(["foreground_kind", "foreground_label"])


# -----------------------------------------------------------------------------
# HTML report
# -----------------------------------------------------------------------------


def fmt(x: Any, n: int = 4) -> str:
    try:
        v = float(x)
        if math.isnan(v):
            return ""
        if math.isinf(v):
            return "∞"
        return f"{v:.{n}g}"
    except Exception:
        return html.escape(str(x)) if x is not None else ""


def html_table(df: pd.DataFrame, cols: Optional[List[str]] = None, max_rows: int = 500) -> str:
    if df is None or df.empty:
        return '<p class="muted">No records.</p>'
    g = df.copy()
    if cols:
        cols = [c for c in cols if c in g.columns]
        g = g[cols]
    g = g.head(max_rows)
    parts = ['<div class="tablewrap"><table><thead><tr>']
    parts.extend(f"<th>{html.escape(str(c))}</th>" for c in g.columns)
    parts.append("</tr></thead><tbody>")
    for _, r in g.iterrows():
        parts.append("<tr>")
        for c in g.columns:
            v = r[c]
            if isinstance(v, float):
                s = fmt(v)
            else:
                s = "" if v is None else str(v)
                if len(s) > 500:
                    s = s[:497] + "..."
                s = html.escape(s)
            parts.append(f"<td>{s}</td>")
        parts.append("</tr>")
    parts.append("</tbody></table></div>")
    return "".join(parts)


def badge(status: str) -> str:
    s = str(status)
    cls = "neutral"
    if "POSITIVE" in s:
        cls = "positive"
    elif "RELAX" in s or "NEUTRAL" in s:
        cls = "relaxed"
    elif "PURIFY" in s or "CONSTRAINT" in s:
        cls = "purifying"
    elif "SENSITIVE" in s or "AMBIG" in s or "INSUFFICIENT" in s:
        cls = "warn"
    return f'<span class="badge {cls}">{html.escape(s)}</span>'


def literature_interpretation_html() -> str:
    return """
    <p>The attached seminal-ribonuclease papers are used here as a methodological benchmark, not as evidence about GPRC6A itself. The 1996 analysis interpreted a high nonsilent:silent ratio on the ox seminal-RNase lineage as adaptive acceleration. The later 2007 codeml analysis found a best AIC branch grouping with approximately ω=1.6 on two recent bovine-lineage branches and approximately ω=0.3 over most of the remaining tree, and interpreted that recent bovine episode as positive selection.</p>
    <p>The same 2007 study is also a warning against a simplistic rule that a pseudogene must have ω≈1. In deer, known/shared disabling lesions were present even though estimated branch ω values were below one, and simulations showed large variance because the gene was short and the timing of inactivation within a branch was unknown. This report therefore calls relaxation only when branch-model, neutral-foreground, branch-site and robustness evidence are mutually consistent.</p>
    """


def write_report(root: Path, reps: Dict[str, AssemblySequence], bg_species: List[str], foregrounds: List[Foreground],
                 alignment_qc: pd.DataFrame, variant_summary: pd.DataFrame,
                 global_results: pd.DataFrame, site_results: pd.DataFrame,
                 fgdf: pd.DataFrame, comp: pd.DataFrame, fg_summary: pd.DataFrame,
                 tool_meta: Dict[str, Any], args) -> Path:
    report = root / "gprc6a_codon_evolution_report.html"
    rep_df = pd.DataFrame([{**{k:v for k,v in asdict(r).items() if k not in {"sequence","exon_sequences","warnings"}},
                            "warnings":";".join(r.warnings)} for r in reps.values()])
    bos = pd.DataFrame()
    if not fg_summary.empty:
        bos = fg_summary[fg_summary.foreground_label.astype(str).map(norm_taxon) == norm_taxon("Bos taurus")]

    css = """
    body{font-family:Arial,Helvetica,sans-serif;margin:0;background:#f4f6f8;color:#18212b;line-height:1.45}
    header{background:#172433;color:white;padding:25px 34px} main{max-width:1500px;margin:auto;padding:22px}
    section{background:white;margin:16px 0;padding:20px 24px;border-radius:10px;box-shadow:0 1px 5px #0002}
    h1,h2,h3{margin-top:.2em}.metrics{display:flex;gap:12px;flex-wrap:wrap}.metric{background:#eef2f6;padding:10px 14px;border-radius:8px;min-width:150px}
    .metric b{display:block;font-size:1.45em}.tablewrap{overflow:auto;max-height:680px;border:1px solid #d8dee6}
    table{border-collapse:collapse;width:100%;font-size:12px}th{position:sticky;top:0;background:#e9eef4;z-index:1}th,td{border:1px solid #d8dee6;padding:5px 7px;text-align:left;vertical-align:top}
    .badge{padding:3px 7px;border-radius:9px;font-size:11px;font-weight:bold;white-space:nowrap}.positive{background:#d9f5e5}.relaxed{background:#fff0c7}.purifying{background:#dcecff}.warn{background:#f7dce1}.neutral{background:#ececec}
    .alert{background:#fff7dd;border-left:5px solid #d5a20b;padding:12px}.good{background:#eaf7ee;border-left:5px solid #2d8b57;padding:12px}.muted{color:#66717d}
    code{background:#eef2f6;padding:1px 4px;border-radius:4px}details{margin:8px 0}
    """
    parts = [f"<!doctype html><html><head><meta charset='utf-8'><title>GPRC6A codon evolution</title><style>{css}</style></head><body>",
             f"<header><h1>GPRC6A codon-level evolutionary analysis</h1><div>Follow-up pipeline v{VERSION}</div></header><main>"]

    nloss = sum(r.species_status in LOSS_SPECIES_STRONG for r in reps.values())
    parts += ["<section><h2>Overview</h2><div class='metrics'>",
              f"<div class='metric'><b>{len(reps)}</b>species reconstructed</div>",
              f"<div class='metric'><b>{len(bg_species)}</b>intact background taxa</div>",
              f"<div class='metric'><b>{nloss}</b>strong/likely loss species</div>",
              f"<div class='metric'><b>{len(foregrounds)}</b>foreground tests</div>",
              f"<div class='metric'><b>{len(variant_summary)}</b>alignment variants</div>",
              "</div>",
              "<p class='alert'><b>Interpretation:</b> ω=dN/dS is a branch-average summary. A terminal loss branch can mix pre-loss constraint and post-loss drift. Positive selection is therefore not called from ω&gt;1 alone when synonymous information is sparse, and relaxation is not called from a single near-neutral estimate.</p>"]
    if not bos.empty:
        br = bos.iloc[0]
        parts.append("<h3>Bos taurus / cow</h3>")
        parts.append(f"<p>{badge(br.consensus_interpretation)} Median ω across valid sensitivity settings: <b>{fmt(br.omega_median)}</b>; range {fmt(br.omega_min)} to {fmt(br.omega_max)}. Positive-selection calls: {inum(br.positive_selection_call_settings)}/{inum(br.valid_sensitivity_settings)}; near-neutral/relaxation calls: {inum(br.relaxed_neutral_call_settings)}/{inum(br.valid_sensitivity_settings)}.</p>")
    else:
        parts.append("<p><b>Bos taurus:</b> no analyzable cow foreground was available in the selected sequence/tree set.</p>")
    parts.append("</section>")

    parts += ["<section><h2>What the attached RNase papers imply for interpretation</h2>", literature_interpretation_html(), "</section>"]

    parts += ["<section><h2>Sequence reconstruction and provenance</h2>",
              "<p>Loss/remnant sequences are reconstructed from the exact qseq/sseq evidence in the raw first-pass BLAST and, when present, the refiner's sensitive exon-rescue BLAST. Previously uncovered positions can therefore be filled by independent nucleotide rescue without borrowing a reference base. Positions still unrecovered become N, explicit target deletions remain absent, and target insertions are retained. A high-confidence ML-recovered CDS can replace an incomplete BLAST reconstruction only for species already classified as intact; loss lineages always retain their remnant sequence so disabling lesions are not hidden.</p>",
              html_table(rep_df, ["species","accession","species_status","assembly_status","assembly_level","source","callable_fraction","reconstructed_fraction","n_fraction","explicit_deleted_nt","insertion_nt","sensitive_rescued_nt","internal_stop_count_naive","supported_events","rescued_events","quality_score","warnings"], 1000),
              "<h3>Intact background selected for PAML</h3><p>" + html.escape("; ".join(bg_species)) + "</p></section>"]

    parts += ["<section><h2>MACSE alignment quality and sensitivity variants</h2>",
              "<p>Disrupted sequences are supplied to MACSE through <code>-seq_lr</code>. The pipeline also builds a second de-novo MACSE alignment with altered frameshift/stop penalties, unless disabled, and refines the baseline alignment. Frameshift-containing or stop-containing codons are converted to missing codons only in PAML-safe derivatives because codeml does not model stop codons or MACSE's <code>!</code> character. Every raw MACSE alignment is retained for inspection.</p>",
              html_table(variant_summary), "<h3>Per-taxon MACSE QC</h3>", html_table(alignment_qc, max_rows=2000), "</section>"]

    parts += ["<section><h2>Global PAML models</h2><h3>M0 and free-ratio branch models</h3>", html_table(global_results),
              "<h3>Intact-background site models</h3>", html_table(site_results), "</section>"]

    parts += ["<section><h2>Lineage-level interpretation</h2>",
              "<p><b>Positive selection</b> requires reproducible branch-site support and/or a branch-average foreground ω&gt;1 that significantly differs from the neutral foreground model, while remaining stable to alignment/model perturbation. Final branch-site calls use the conservative χ²₁ p value with within-setting FDR correction; the 50:50 mixture result is also reported. <b>Relaxed or near-neutral evolution</b> requires ω near one, no branch-site signal, and failure to reject a foreground ω of one across robust settings. <b>Purifying constraint</b> is assigned when foreground ω remains clearly below one and the neutral foreground model is rejected.</p>",
              html_table(fg_summary, ["foreground_id","foreground_kind","foreground_label","foreground_descendants","species_status","consensus_interpretation","valid_sensitivity_settings","dominant_call_fraction","omega_median","omega_min","omega_max","omega_crosses_1_across_settings","branch_site_significant_settings","positive_selection_call_settings","relaxed_neutral_call_settings","purifying_call_settings","alignment_model_sensitive"], 1000),
              "</section>"]

    if not comp.empty:
        cc = comp.copy(); cc["setting_call"] = cc.apply(setting_call, axis=1)
        parts += ["<section><h2>Alignment/model robustness matrix</h2>", html_table(cc, ["foreground_id","foreground_label","variant","codon_freq","foreground_min_informative_codons","foreground_dN","foreground_dS","foreground_omega","background_omega","omega_shift_from_background","branch_shift_p","foreground_vs_neutral_p","branch_site_LRT","branch_site_p_chi1","branch_site_p_mixture","branch_site_q","branch_site_q_mixture","BEB_sites_ge_0.95","setting_call","status"], 5000), "</section>"]

    parts += ["<section><h2>Model caveats</h2><ul>",
              "<li>PAML codon models assume an ancestral coding frame. MACSE preserves frameshift information for alignment, but PAML-safe analyses mask affected/stop codons. Deeply decayed pseudogene sequence can therefore become uninformative rather than genuinely neutral in the likelihood analysis.</li>",
              "<li>A terminal branch estimate describes the full history between its ancestral node and the sampled species. It cannot by itself date the disabling lesion.</li>",
              "<li>Very small dS makes ω unstable. Such settings are marked rather than interpreted as strong selection.</li>",
              "<li>Branch-site tests identify episodic selection at a subset of codons. Branch-wide ω can remain below one even when branch-site selection is real.</li>",
              "<li>Alignment sensitivity is explicitly measured. If the ω range crosses one or the dominant qualitative call is recovered in fewer than 70% of valid sensitivity settings, the final call is <b>ALIGNMENT_OR_MODEL_SENSITIVE</b>.</li>",
              "</ul></section>"]

    parts += ["<section><h2>Run metadata</h2><pre>" + html.escape(json.dumps(tool_meta, indent=2, default=str)) + "</pre>",
              "<p class='muted'>All raw codeml outputs, control files, branch dN/dS tables, alignments and machine-readable summaries are retained under this analysis directory.</p></section>",
              "</main></body></html>"]
    report.write_text("".join(parts), encoding="utf-8")
    return report


# -----------------------------------------------------------------------------
# CLI and main
# -----------------------------------------------------------------------------



# =============================================================================
# v1.3 rigorous analysis layer
# =============================================================================


def is_hybrid_species(sp: str) -> bool:
    s = str(sp)
    return " x " in s or " × " in s


def tree_parent_map(tree: Tree) -> Dict[int, Optional[Clade]]:
    parents: Dict[int, Optional[Clade]] = {id(tree.root): None}
    for p in tree.find_clades(order="preorder"):
        for c in p.clades:
            parents[id(c)] = p
    return parents


def branch_desc_key(c: Clade) -> Tuple[str, ...]:
    return tuple(sorted(t.name for t in c.get_terminals() if t.name))


def find_parent(tree: Tree, child: Clade) -> Optional[Clade]:
    return tree_parent_map(tree).get(id(child))


def choose_topology_guardians(tree: Tree, target_desc: Set[str], available_ids: Set[str],
                              id_quality: Dict[str, float]) -> Set[str]:
    """Choose sister-lineage tips needed to stop pruning from collapsing a target branch."""
    target = find_exact_clade(tree, target_desc)
    if target is None:
        return set()
    parent = find_parent(tree, target)
    if parent is None:
        return set()
    guardians: Set[str] = set()
    for sib in parent.clades:
        if sib is target:
            continue
        cand = [x.name for x in sib.get_terminals() if x.name in available_ids]
        if cand:
            guardians.add(max(cand, key=lambda x: id_quality.get(x, 0.0)))
    return guardians


def write_macse_inputs_v13(root: Path, reps: Dict[str, AssemblySequence], intact_background: List[str],
                            ruminant_context: List[str], focus_species: List[str], ref_cds: str) -> Dict[str, Any]:
    """Write a full alignment set containing intact controls plus phylogenetic context.

    Unlike v1.2, this does not align only one foreground against distant intact taxa.
    All sufficiently callable ruminants are retained by default, so terminal and
    ancestral foreground edges can remain distinct after tree pruning.
    """
    ad = root / "alignment"
    ad.mkdir(parents=True, exist_ok=True)
    reliable = ad / "macse_reliable_intact.fasta"
    less = ad / "macse_less_reliable_ruminants.fasta"
    ids: Set[str] = set()
    meta = []
    species_to_id: Dict[str, str] = {}

    def unique_id(sp: str) -> str:
        base = safe_name(sp); x = base; n = 2
        while x in ids:
            x = f"{base}_{n}"; n += 1
        ids.add(x); return x

    with open(reliable, "w") as fh:
        pig_id = "Sus_scrofa_REF"; ids.add(pig_id)
        fh.write(f">{pig_id}\n{wrap_fasta(ref_cds)}\n")
        meta.append({"taxon_id":pig_id,"species":"Sus scrofa","role":"reference_intact","status":"INTACT_REFERENCE"})
        for sp in intact_background:
            if sp not in reps: continue
            sid = unique_id(sp); species_to_id[sp] = sid; r = reps[sp]
            fh.write(f">{sid}\n{wrap_fasta(r.sequence)}\n")
            meta.append({"taxon_id":sid,"species":sp,"role":"external_intact_background","status":r.species_status,
                         "callable_fraction":r.callable_fraction,"quality_score":r.quality_score})

    selected_ruminants=[]
    for sp in list(ruminant_context) + list(focus_species):
        if sp in reps and sp not in selected_ruminants and sp not in species_to_id:
            selected_ruminants.append(sp)
    with open(less, "w") as fh:
        for sp in selected_ruminants:
            r=reps[sp]; sid=unique_id(sp); species_to_id[sp]=sid
            fh.write(f">{sid}\n{wrap_fasta(r.sequence)}\n")
            meta.append({"taxon_id":sid,"species":sp,"role":"ruminant_context_or_foreground","status":r.species_status,
                         "callable_fraction":r.callable_fraction,"quality_score":r.quality_score})
    write_tsv(pd.DataFrame(meta), ad/"selected_taxa_v13.tsv")
    return {"reliable":reliable,"less":less,"species_to_id":species_to_id,"selected_meta":meta,
            "pig_id":"Sus_scrofa_REF","ruminant_context":selected_ruminants}


def select_ruminant_context_v13(reps: Dict[str, AssemblySequence], base_tree: Tree, focus_species: List[str],
                                min_callable: float, include_hybrids: bool, mode: str,
                                anchor_count: int) -> List[str]:
    tips={t.name for t in base_tree.get_terminals()}
    external={sp for sp,r in reps.items() if r.species_status=="EXTERNAL_INTACT"}
    candidates=[]
    for sp,r in reps.items():
        if sp in external: continue
        if not include_hybrids and is_hybrid_species(sp): continue
        if safe_name(sp) not in tips: continue
        if r.callable_fraction >= min_callable:
            candidates.append(sp)
    for sp in focus_species:
        if sp in reps and (include_hybrids or not is_hybrid_species(sp)) and safe_name(sp) in tips and sp not in candidates:
            candidates.append(sp)
    if mode=="all_callable":
        return sorted(candidates)
    # Faster mode: quality-ranked anchors plus all focal species. Exact sister
    # guardians are added per foreground later.
    focus=set(focus_species)
    others=[sp for sp in candidates if sp not in focus]
    others.sort(key=lambda sp: reps[sp].quality_score, reverse=True)
    return sorted(focus | set(others[:max(0,anchor_count)]))


def run_intact_only_macse_qc_v13(args, root: Path, reliable_fasta: Path) -> Tuple[Path, pd.DataFrame, bool]:
    """Align only independently intact/reference CDSs before pseudogenes are introduced."""
    base,_=resolve_macse(args)
    d=root/"alignment"/"intact_only_qc"; d.mkdir(parents=True,exist_ok=True)
    nt=d/"intact_only_macse_NT.fasta"; aa=d/"intact_only_macse_AA.fasta"
    cmd=base+["-prog","alignSequences","-seq",str(reliable_fasta),"-out_NT",str(nt),"-out_AA",str(aa)]
    if not (args.resume and fasta_nonempty(nt)):
        run_cmd(cmd,check=True,log_path=d/"macse.log")
    recs=read_alignment(nt)
    qc=macse_alignment_qc(recs,"intact_only")
    write_tsv(qc,d/"intact_only_macse_qc.tsv")
    bad=qc[(pd.to_numeric(qc.frameshift_characters,errors="coerce").fillna(999)>args.intact_qc_max_frameshift) |
           (pd.to_numeric(qc.internal_stop_codons,errors="coerce").fillna(999)>args.intact_qc_max_internal_stop)]
    passed=bad.empty
    summary={"passed":passed,"n_taxa":len(qc),"max_frameshift":int(pd.to_numeric(qc.frameshift_characters,errors="coerce").max()),
             "max_internal_stop":int(pd.to_numeric(qc.internal_stop_codons,errors="coerce").max()),
             "failed_taxa":bad.taxon_id.astype(str).tolist() if len(bad) else []}
    write_json(summary,d/"intact_only_qc_summary.json")
    if not passed:
        msg=("Intact-only MACSE QC failed. At least one independently intact/reference CDS acquired more "
             "frameshift/stop characters than allowed. Inspect alignment/intact_only_qc before interpreting selection.")
        if args.allow_intact_qc_warning: logging.warning(msg)
        else: raise RuntimeError(msg+" Use --allow-intact-qc-warning only for exploratory continuation.")
    return nt,qc,passed


def fitch_loss_origins_v13(tree: Tree, reps: Dict[str, AssemblySequence], species_to_id: Dict[str,str],
                           pig_id: str) -> List[Tuple[Set[str], str]]:
    """Infer minimum-change 0->1 loss transitions with unresolved tips treated as unknown.

    Root state is preferentially functional (0) when compatible with Fitch sets,
    reflecting the independently intact external mammalian GPRC6A sequences.
    """
    id_to_species={sid:sp for sp,sid in species_to_id.items()}
    tip_sets={}
    for tip in tree.get_terminals():
        sid=tip.name
        if sid==pig_id: tip_sets[id(tip)]={0}; continue
        sp=id_to_species.get(sid)
        if sp and reps[sp].species_status in LOSS_SPECIES_STRONG:
            tip_sets[id(tip)]={1}
        elif sp and reps[sp].species_status in INTACT_SPECIES:
            tip_sets[id(tip)]={0}
        else:
            tip_sets[id(tip)]={0,1}
    state_sets={}
    def up(c):
        if not c.clades:
            s=set(tip_sets.get(id(c),{0,1})); state_sets[id(c)]=s; return s
        cs=[up(x) for x in c.clades]
        inter=set.intersection(*cs) if cs else {0,1}
        s=inter if inter else set.union(*cs)
        state_sets[id(c)]=s; return s
    up(tree.root)
    assigned={id(tree.root):(0 if 0 in state_sets[id(tree.root)] else min(state_sets[id(tree.root)]))}
    origins=[]
    def down(c):
        ps=assigned[id(c)]
        for ch in c.clades:
            opts=state_sets[id(ch)]
            cs=ps if ps in opts else (0 if 0 in opts else min(opts))
            assigned[id(ch)]=cs
            if ps==0 and cs==1:
                desc={t.name for t in ch.get_terminals() if t.name}
                origins.append((desc,"fitch_0_to_1"))
            down(ch)
    down(tree.root)
    return origins


def build_foregrounds_v13(upstream: Dict[str,Any], reps: Dict[str,AssemblySequence], species_to_id: Dict[str,str],
                           base_tree: Tree, pig_id: str, include_possible: bool,
                           always_test: Sequence[str], min_focus_callable: float) -> List[Foreground]:
    """Define terminal, ancestral-ruminant, and inferred-loss-origin foregrounds."""
    allowed=set(LOSS_SPECIES_STRONG)
    if include_possible: allowed |= LOSS_SPECIES_EXPLORATORY
    out=[]; seen=set(); tree_tips={t.name for t in base_tree.get_terminals()}
    by_norm={norm_taxon(sp):sp for sp in reps}
    # C: terminal branches.
    for sp,r in sorted(reps.items()):
        if r.species_status not in allowed or r.callable_fraction < min_focus_callable: continue
        if is_hybrid_species(sp): continue
        sid=species_to_id.get(sp)
        if sid and sid in tree_tips:
            out.append(Foreground(f"terminal_{safe_name(sp)}","terminal",sp,[sp],[sid],"refined_species_loss_call",r.species_status)); seen.add((sid,))
    for req in always_test:
        sp=by_norm.get(norm_taxon(req)); sid=species_to_id.get(sp) if sp else None
        if sp and sid in tree_tips and reps[sp].callable_fraction>=min_focus_callable and (sid,) not in seen:
            out.append(Foreground(f"terminal_{safe_name(sp)}","terminal",sp,[sp],[sid],"explicit_focus",reps[sp].species_status)); seen.add((sid,))
    # A: branch leading into Ruminantia. External taxa are excluded, pig is separate.
    rum_ids={sid for sp,sid in species_to_id.items() if reps.get(sp) and reps[sp].species_status!="EXTERNAL_INTACT" and sid in tree_tips}
    if len(rum_ids)>=2:
        try:
            mrca=base_tree.common_ancestor(list(rum_ids))
            desc={t.name for t in mrca.get_terminals() if t.name and t.name in rum_ids}
            # Require the MRCA not to contain external background tips.
            all_desc={t.name for t in mrca.get_terminals() if t.name}
            if desc and not (all_desc-desc):
                spdesc=[sp for sp,sid in species_to_id.items() if sid in desc]
                out.append(Foreground("ancestral_Ruminantia","ancestral_ruminant","Ancestral Ruminantia",spdesc,sorted(desc),
                                      "explicit_pre_ruminant_branch","ANCESTRAL_BRANCH"))
        except Exception as e:
            logging.warning("Could not define ancestral Ruminantia branch: %s",e)
    # B: infer 0->1 transitions directly on the combined tree, even when upstream origin table is empty.
    for i,(desc,method) in enumerate(fitch_loss_origins_v13(base_tree,reps,species_to_id,pig_id),1):
        ids=sorted(desc & set(species_to_id.values()))
        if not ids: continue
        species=[sp for sp,sid in species_to_id.items() if sid in ids]
        disrupted=sum(reps[sp].species_status in allowed for sp in species if sp in reps)
        if disrupted<1: continue
        key=tuple(ids)
        out.append(Foreground(f"inferred_loss_origin_{i}","inferred_loss_origin",f"Inferred loss origin {i}",species,ids,
                              method,"INFERRED_LOSS_ORIGIN"))
    return out


def foreground_context_v13(base_tree: Tree, all_records: Dict[str,str], fg: Foreground,
                           base_context_ids: Set[str], id_quality: Dict[str,float],
                           bg_ids: Set[str], pig_id: str) -> Tuple[Tree,List[str],Set[str],Set[str]]:
    """Build a topology-preserving foreground subset.

    The target edge is located in the full aligned tree before any pruning. For
    an internal foreground we retain at least one informative descendant from
    each immediate child of the target clade, which preserves the target node
    without forcing every poorly recovered pseudogene in that clade into PAML.
    At least one sister-lineage guardian outside the target is also retained so
    the edge leading into the target clade cannot collapse into its parent.
    """
    tree_tips={t.name for t in base_tree.get_terminals() if t.name}
    fg_full=set(fg.descendant_ids) & tree_tips
    if not fg_full: raise ValueError("Foreground absent from full standardized tree")
    full_target=find_exact_clade(base_tree,fg_full)
    if full_target is None:
        raise ValueError(f"Foreground does not identify one branch in full tree: {sorted(fg_full)}")
    avail=set(all_records)

    if not full_target.clades:
        fg_retained=fg_full & avail
    else:
        # Start with good context taxa already selected within the target.
        fg_retained=(fg_full & set(base_context_ids) & avail)
        # Preserve the internal target node by keeping at least one descendant
        # from each immediate child subtree of the original target clade.
        for child in full_target.clades:
            cand=[x.name for x in child.get_terminals() if x.name in avail and x.name in fg_full]
            if cand and not (set(cand)&fg_retained):
                fg_retained.add(max(cand,key=lambda x:id_quality.get(x,0.0)))
        if len(fg_retained)<2:
            extra=sorted((fg_full & avail)-fg_retained,key=lambda x:id_quality.get(x,0.0),reverse=True)
            fg_retained.update(extra[:max(0,2-len(fg_retained))])

    if not fg_retained:
        raise ValueError("No aligned foreground descendants remain")
    keep=(set(base_context_ids)|set(bg_ids)|{pig_id}|fg_retained) & avail
    keep |= choose_topology_guardians(base_tree,fg_full,avail,id_quality)
    t=prune_with_pig_policy(base_tree,keep,pig_id)
    actual={x.name for x in t.get_terminals()} & avail
    fg_actual=fg_retained & actual
    if not fg_actual or find_exact_clade(t,fg_actual) is None:
        raise ValueError(f"Target branch collapsed despite topology guardians: {fg.foreground_id}")
    # Internal targets must remain internal after pruning.
    target_after=find_exact_clade(t,fg_actual)
    if full_target.clades and (target_after is None or not target_after.clades):
        raise ValueError(f"Internal target collapsed to a terminal edge: {fg.foreground_id}")
    order=[x.name for x in t.get_terminals() if x.name in actual]
    return t,order,fg_actual,actual


def subset_alignment_variants_v13(root: Path, fg: Foreground, macse_outputs: Dict[str,Path], taxa: Set[str],
                                   min_codons: int) -> Tuple[Dict[str,Path],pd.DataFrame]:
    """Create masks after taxon subsetting, not before it."""
    d=root/"foregrounds"/safe_name(fg.foreground_id)/"alignment_variants"; d.mkdir(parents=True,exist_ok=True)
    variants={}; rows=[]
    sources=[("baseline",macse_outputs["baseline"])]
    if "refined" in macse_outputs: sources.append(("refined",macse_outputs["refined"]))
    if "altcost" in macse_outputs: sources.append(("altcost",macse_outputs["altcost"]))
    for src,path in sources:
        full=read_alignment(path); recs={k:v for k,v in full.items() if k in taxa}
        if len(recs)<4: continue
        specs=[(f"{src}_mask50",0.50,0,False),(f"{src}_trim1_occ60",0.60,1,False)]
        if src=="baseline": specs += [("baseline_trim2_occ70",0.70,2,False),("baseline_strict_complete",1.0,0,True)]
        for name,occ,flank,strict in specs:
            vv,st=build_paml_variant(recs,occ,flank,strict)
            st.update({"foreground_id":fg.foreground_id,"variant":name,"source_alignment":src,"subset_taxa":len(recs)})
            st["usable_for_paml"]=st["output_codons"]>=min_codons
            out=d/f"{name}.fasta"; write_fasta_dict(vv,out)
            write_json(st,d/f"{name}.stats.json")
            if st["usable_for_paml"]: variants[name]=out
            rows.append({k:v for k,v in st.items() if k not in {"kept_column_indices_0based","valid_codons_per_taxon"}})
    return variants,pd.DataFrame(rows)


def parse_codeml_validation_v13(job_dir: Path, returncode: Optional[int]) -> Dict[str,Any]:
    mlc=job_dir/"mlc"; runlog=job_dir/"run.log"
    parsed=parse_codeml_output(mlc)
    ll=fnum(parsed.get("lnL"))
    txt=""
    for p in (mlc,runlog):
        if p.exists():
            try: txt += "\n"+p.read_text(errors="replace")[-20000:]
            except Exception: pass
    fatal_markers=[]
    low=txt.lower()
    for marker in ["segmentation fault","floating point exception","nan;","nan ","fatal error","error: cannot","error: end","check convergence","not converged"]:
        if marker in low: fatal_markers.append(marker)
    finite=math.isfinite(ll)
    rc_ok=(returncode==0)
    status="ok" if rc_ok and finite and not fatal_markers else (
        "returncode_nonzero" if returncode not in (0,None) else "parse_or_numerical_failure")
    parsed.update({"status":status,"returncode":returncode,"finite_lnl":finite,
                   "fatal_markers":";".join(fatal_markers)})
    return parsed


def run_codeml_job_v13(codeml: str, job_dir: Path, seq_records: Dict[str,str], taxa: Sequence[str],
                        tree: Tree, foreground_desc: Optional[Set[str]], model_name: str,
                        model: int,nssites:int,codon_freq:int,fix_omega:int=0,omega:float=0.4,
                        cleandata:int=0,resume:bool=True,ncatg:int=8) -> Dict[str,Any]:
    """Strict codeml runner. Failed processes can never enter an LRT."""
    job_dir.mkdir(parents=True,exist_ok=True); summary=job_dir/"parsed_summary_v13.json"; mlc=job_dir/"mlc"
    if resume and summary.exists():
        old=read_json(summary)
        if old.get("status")=="ok" and mlc.exists():
            parsed=parse_codeml_validation_v13(job_dir,0)
            parsed.update({"model_name":model_name,"job_dir":str(job_dir),"codon_freq":codon_freq,"resumed":True})
            if parsed.get("status")=="ok": return parsed
    seqfile=job_dir/"alignment.phy"; treefile=job_dir/"tree.nwk"; ctl=job_dir/"codeml.ctl"
    codons=write_phylip(seq_records,taxa,seqfile)
    treefile.write_text(paml_newick(tree,foreground_desc))
    ctl.write_text(codeml_ctl_text(seqfile.name,treefile.name,mlc.name,model,nssites,codon_freq,fix_omega,omega,cleandata,0,ncatg))
    cp=run_cmd([codeml,ctl.name],cwd=job_dir,check=False,log_path=job_dir/"run.log")
    parsed=parse_codeml_validation_v13(job_dir,cp.returncode)
    parsed.update({"model_name":model_name,"job_dir":str(job_dir),"codon_freq":codon_freq,"resumed":False,"alignment_codons":codons})
    write_json({k:v for k,v in parsed.items() if k!="branch_table"},summary)
    if parsed.get("branch_table"): write_tsv(pd.DataFrame(parsed["branch_table"]),job_dir/"branch_dnds.tsv")
    if parsed.get("beb_sites"): write_tsv(pd.DataFrame(parsed["beb_sites"]),job_dir/"beb_sites.tsv")
    return parsed


def valid_model_set_v13(models: Sequence[Dict[str,Any]]) -> bool:
    return bool(models) and all(str(x.get("status"))=="ok" and math.isfinite(fnum(x.get("lnL"))) and x.get("returncode") in (0,None) for x in models)


def run_one_paml_setting_v13(args, root: Path, codeml: str, fg: Foreground, variant_name: str,
                              variant_path: Path, tree: Tree, taxa_order: List[str], fg_ids: Set[str],
                              codon_freq: int) -> Dict[str,Any]:
    recs=read_variant_records(variant_path)
    recs={k:v for k,v in recs.items() if k in set(taxa_order)}
    info=[sum(is_valid_sense_codon(recs[s][i:i+3]) for i in range(0,len(recs[s]),3)) for s in fg_ids if s in recs]
    min_info=min(info) if info else 0
    base={"foreground_id":fg.foreground_id,"foreground_kind":fg.kind,"foreground_label":fg.species_label,
          "foreground_descendants":";".join(fg.descendant_species),"foreground_species_status":fg.species_status,
          "variant":variant_name,"codon_freq":codon_freq,"foreground_min_informative_codons":min_info,
          "taxa_in_analysis":len(taxa_order)}
    if min_info<args.min_foreground_codons:
        return {**base,"status":"too_few_foreground_codons"}
    parsed={}
    for name,model,ns,fixw,w in model_job_specs():
        d=root/"paml"/safe_name(fg.foreground_id)/variant_name/f"F{codon_freq}"/name
        parsed[name]=run_codeml_job_v13(codeml,d,recs,taxa_order,tree,fg_ids,name,model,ns,codon_freq,fixw,w,0,args.resume)
    one=parsed["branch_one_ratio"]; two=parsed["branch_two_ratio"]; neu=parsed["branch_foreground_neutral"]
    bsa=parsed["branch_site_alt"]; bsn=parsed["branch_site_null"]
    status="ok" if valid_model_set_v13([one,two,neu,bsa,bsn]) else "MODEL_FAILED_NO_LRT"
    row={**base,"status":status,
         "model_statuses":";".join(f"{k}:{v.get('status')}" for k,v in parsed.items()),
         "model_returncodes":";".join(f"{k}:{v.get('returncode')}" for k,v in parsed.items())}
    if status!="ok": return row
    metric=branch_metric(two,fg_ids)
    om=[fnum(x) for x in (two.get("omega_vector") or []) if math.isfinite(fnum(x))]
    bg=om[0] if len(om)>=2 else float("nan")
    lr_shift=max(0.0,2*(fnum(two["lnL"])-fnum(one["lnL"])))
    lr_neu=max(0.0,2*(fnum(two["lnL"])-fnum(neu["lnL"])))
    lr_bs=max(0.0,2*(fnum(bsa["lnL"])-fnum(bsn["lnL"])))
    pbs=chi2_sf(lr_bs,1); pmix=0.5*pbs if lr_bs>0 else 1.0
    row.update({"foreground_dN":metric["dN"],"foreground_dS":metric["dS"],"foreground_omega":metric["omega"],
                "background_omega":bg,"omega_shift_from_background":metric["omega"]-bg if math.isfinite(metric["omega"]) and math.isfinite(bg) else float("nan"),
                "branch_shift_LRT":lr_shift,"branch_shift_p":chi2_sf(lr_shift,1),
                "foreground_vs_neutral_LRT":lr_neu,"foreground_vs_neutral_p":chi2_sf(lr_neu,1),
                "branch_site_LRT":lr_bs,"branch_site_p_chi1":pbs,"branch_site_p_mixture":pmix,
                "BEB_sites_ge_0.95":sum(fnum(x.get("posterior"),0)>=0.95 for x in bsa.get("beb_sites",[]) or []),
                "BEB_sites_ge_0.99":sum(fnum(x.get("posterior"),0)>=0.99 for x in bsa.get("beb_sites",[]) or [])})
    return row


def hyphy_newick_v13(tree: Tree, test_desc: Set[str], reference_tips: Set[str]) -> str:
    target=find_exact_clade(tree,test_desc)
    if target is None: raise ValueError("HyPhy test branch is not exact in tree")
    def rec(c):
        if c.clades: s="("+",".join(rec(x) for x in c.clades)+")"
        else: s=safe_name(c.name or "taxon")
        if c is target: s+="{test}"
        elif not c.clades and c.name in reference_tips: s+="{reference}"
        return s
    return rec(tree.root)+";\n"


def recursive_find_key_v13(obj: Any, key: str) -> Any:
    if isinstance(obj,dict):
        if key in obj: return obj[key]
        for v in obj.values():
            z=recursive_find_key_v13(v,key)
            if z is not None: return z
    elif isinstance(obj,list):
        for v in obj:
            z=recursive_find_key_v13(v,key)
            if z is not None: return z
    return None


def parse_relax_json_v13(path: Path) -> Dict[str,Any]:
    if not path.exists() or path.stat().st_size<10: return {"status":"missing_json"}
    try: data=json.loads(path.read_text())
    except Exception as e: return {"status":"json_parse_failed","error":str(e)}
    tr=recursive_find_key_v13(data,"test results")
    if not isinstance(tr,dict): tr={}
    K=fnum(tr.get("relaxation or intensification parameter"))
    p=fnum(tr.get("p-value")); lrt=fnum(tr.get("LRT"))
    status="ok" if all(math.isfinite(x) for x in [K,p,lrt]) else "parse_incomplete"
    return {"status":status,"K":K,"p":p,"LRT":lrt}


def run_relax_job_v13(args, root: Path, fg: Foreground, variant_name: str, variant_path: Path,
                       tree: Tree, fg_ids: Set[str], reference_ids: Set[str]) -> Dict[str,Any]:
    hy=exe(args.hyphy,"hyphy")
    if not hy: return {"foreground_id":fg.foreground_id,"variant":variant_name,"status":"hyphy_not_found"}
    d=root/"hyphy_relax"/safe_name(fg.foreground_id)/variant_name; d.mkdir(parents=True,exist_ok=True)
    recs=read_variant_records(variant_path); taxa=set(recs)
    refs=reference_ids & taxa
    if len(refs)<3: return {"foreground_id":fg.foreground_id,"variant":variant_name,"status":"too_few_reference_branches"}
    treefile=d/"tree_relax.nwk"; treefile.write_text(hyphy_newick_v13(tree,fg_ids,refs))
    out=d/"RELAX.json"; log=d/"RELAX.log"
    if args.resume and out.exists():
        parsed=parse_relax_json_v13(out)
        if parsed.get("status")=="ok":
            return {"foreground_id":fg.foreground_id,"foreground_label":fg.species_label,"foreground_kind":fg.kind,
                    "variant":variant_name,"returncode":0,"resumed":True,**parsed,"job_dir":str(d)}
    cmd=[hy,"relax","--alignment",str(variant_path),"--tree",str(treefile),"--test","test","--reference","reference",
         "--models",args.relax_models,"--output",str(out)]
    cp=run_cmd(cmd,check=False,log_path=log)
    parsed=parse_relax_json_v13(out)
    if cp.returncode!=0: parsed["status"]="hyphy_failed"
    return {"foreground_id":fg.foreground_id,"foreground_label":fg.species_label,"foreground_kind":fg.kind,
            "variant":variant_name,"returncode":cp.returncode,"resumed":False,**parsed,"job_dir":str(d)}


def summarize_relax_v13(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty: return pd.DataFrame()
    g=df.copy(); g["q"]=float("nan")
    for v,idx in g.groupby("variant").groups.items():
        vals=[fnum(g.loc[i,"p"],1.0) if str(g.loc[i,"status"])=="ok" else 1.0 for i in idx]
        for i,q in zip(idx,bh_fdr(vals)): g.loc[i,"q"]=q
    rows=[]
    for fid,x in g.groupby("foreground_id"):
        ok=x[x.status=="ok"].copy()
        calls=[]
        for _,r in ok.iterrows():
            if fnum(r.q,1)>=0.05: calls.append("NO_SIGNIFICANT_CHANGE")
            elif fnum(r.K)<1: calls.append("RELAXATION")
            else: calls.append("INTENSIFICATION")
        c=Counter(calls); top,n=c.most_common(1)[0] if c else ("INSUFFICIENT",0)
        frac=n/len(calls) if calls else 0
        rows.append({"foreground_id":fid,"foreground_label":x.foreground_label.iloc[0],"valid_relax_settings":len(ok),
                     "relax_consensus":top if frac>=0.67 else "RELAX_MODEL_SENSITIVE","relax_dominant_fraction":frac,
                     "K_median":float(pd.to_numeric(ok.K,errors="coerce").median()) if len(ok) else float("nan"),
                     "K_min":float(pd.to_numeric(ok.K,errors="coerce").min()) if len(ok) else float("nan"),
                     "K_max":float(pd.to_numeric(ok.K,errors="coerce").max()) if len(ok) else float("nan"),
                     "relaxation_settings":c.get("RELAXATION",0),"intensification_settings":c.get("INTENSIFICATION",0),
                     "no_change_settings":c.get("NO_SIGNIFICANT_CHANGE",0)})
    return pd.DataFrame(rows),g


def integrate_evidence_v13(paml_summary: pd.DataFrame, relax_summary: pd.DataFrame) -> pd.DataFrame:
    if paml_summary.empty: return pd.DataFrame()
    rmap={r.foreground_id:r for _,r in relax_summary.iterrows()} if relax_summary is not None and not relax_summary.empty else {}
    rows=[]
    for _,p in paml_summary.iterrows():
        rr=rmap.get(p.foreground_id); rc=str(rr.relax_consensus) if rr is not None else "NOT_RUN"
        pc=str(p.consensus_interpretation)
        if pc=="ALIGNMENT_OR_MODEL_SENSITIVE" or rc=="RELAX_MODEL_SENSITIVE": final="ALIGNMENT_OR_MODEL_SENSITIVE"
        elif pc.startswith("POSITIVE_SELECTION") and rc=="RELAXATION": final="MIXED_POSITIVE_SELECTION_AND_RELAXATION"
        elif pc.startswith("POSITIVE_SELECTION") and rc=="INTENSIFICATION": final="EPISODIC_POSITIVE_SELECTION_WITH_INTENSIFICATION"
        elif pc.startswith("POSITIVE_SELECTION"): final="EPISODIC_POSITIVE_SELECTION"
        elif rc=="RELAXATION": final="RELAXED_SELECTION"
        elif rc=="INTENSIFICATION": final="SELECTION_INTENSIFICATION"
        elif pc=="PURIFYING_CONSTRAINT": final="CONTINUED_PURIFYING_CONSTRAINT"
        elif pc in {"RELAXED_TOWARD_NEUTRALITY","NEAR_NEUTRAL_WEAK_EVIDENCE"}: final="RELAXATION_CANDIDATE"
        else: final="AMBIGUOUS"
        row=dict(p); row.update({"RELAX_consensus":rc,"RELAX_K_median":fnum(rr.K_median) if rr is not None else float("nan"),
                                 "integrated_interpretation":final})
        rows.append(row)
    return pd.DataFrame(rows)


def run_paml_and_relax_v13(args, root: Path, macse_outputs: Dict[str,Path], base_tree: Tree,
                           reps: Dict[str,AssemblySequence], species_to_id: Dict[str,str], pig_id: str,
                           bg_species: List[str], ruminant_context: List[str], foregrounds: List[Foreground]):
    codeml=exe(args.codeml,"codeml")
    if not codeml: raise RuntimeError("codeml not found")
    bg_ids={species_to_id[x] for x in bg_species if x in species_to_id}; ref_ids=set(bg_ids)|{pig_id}
    id_quality={sid:reps[sp].quality_score for sp,sid in species_to_id.items() if sp in reps}; id_quality[pig_id]=999
    base_context={species_to_id[sp] for sp in ruminant_context if sp in species_to_id}
    raw_records=read_alignment(macse_outputs["baseline"])
    codon_freqs=[int(x) for x in str(args.codon_freqs).split(",") if x.strip()]
    all_rows=[]; relax_jobs=[]; variant_qc=[]; fg_context_rows=[]
    per_fg={}
    for fg in foregrounds:
        try:
            t,order,fgids,actual=foreground_context_v13(base_tree,raw_records,fg,base_context,id_quality,bg_ids,pig_id)
            variants,qc=subset_alignment_variants_v13(root,fg,macse_outputs,actual,args.min_paml_codons)
            variant_qc.append(qc)
            per_fg[fg.foreground_id]=(t,order,fgids,actual,variants)
            fg_context_rows.append({"foreground_id":fg.foreground_id,"foreground_kind":fg.kind,"foreground_label":fg.species_label,
                                    "n_taxa":len(order),"n_foreground_descendants":len(fgids),"foreground_ids":";".join(sorted(fgids)),
                                    "taxa":";".join(order)})
        except Exception as e:
            logging.warning("Foreground context preparation failed %s: %s",fg.foreground_id,e)
    write_tsv(pd.DataFrame(fg_context_rows),root/"tables"/"foreground_contexts.tsv")
    if variant_qc: write_tsv(pd.concat(variant_qc,ignore_index=True),root/"tables"/"per_foreground_alignment_variants.tsv")

    tasks=[]
    with ThreadPoolExecutor(max_workers=max(1,args.paml_workers)) as pool:
        for fg in foregrounds:
            if fg.foreground_id not in per_fg: continue
            t,order,fgids,actual,variants=per_fg[fg.foreground_id]
            matrix=select_run_matrix(list(variants),codon_freqs,args.analysis_depth) if variants else []
            for v,cf,_ in matrix:
                tasks.append(pool.submit(run_one_paml_setting_v13,args,root,codeml,fg,v,variants[v],t,order,fgids,cf))
        for i,f in enumerate(as_completed(tasks),1):
            try: all_rows.append(f.result())
            except Exception as e: logging.exception("Strict PAML setting failed: %s",e)
            if i%20==0 or i==len(tasks): logging.info("Strict PAML settings: %d/%d",i,len(tasks))
    comp=pd.DataFrame(all_rows)
    if not comp.empty:
        comp["branch_site_q"]=float("nan"); comp["branch_site_q_mixture"]=float("nan")
        for (v,cf),idx in comp[comp.status=="ok"].groupby(["variant","codon_freq"]).groups.items():
            qs=bh_fdr([fnum(comp.loc[i,"branch_site_p_chi1"],1) for i in idx]); qm=bh_fdr([fnum(comp.loc[i,"branch_site_p_mixture"],1) for i in idx])
            for i,a,b in zip(idx,qs,qm): comp.loc[i,"branch_site_q"]=a; comp.loc[i,"branch_site_q_mixture"]=b
    paml_summary=summarize_foregrounds(comp,foregrounds) if not comp.empty else pd.DataFrame()

    relax_df=pd.DataFrame(); relax_summary=pd.DataFrame()
    if not args.skip_relax:
        if not exe(args.hyphy,"hyphy"):
            raise RuntimeError("HyPhy not found. Install `hyphy` or use --skip-relax for exploratory PAML-only execution.")
        rtasks=[]
        with ThreadPoolExecutor(max_workers=max(1,args.hyphy_workers)) as pool:
            for fg in foregrounds:
                if fg.foreground_id not in per_fg: continue
                t,order,fgids,actual,variants=per_fg[fg.foreground_id]
                names=[]
                if "baseline_mask50" in variants: names.append("baseline_mask50")
                if args.relax_variants=="all_mask50":
                    for z in ["refined_mask50","altcost_mask50"]:
                        if z in variants: names.append(z)
                for v in names:
                    rtasks.append(pool.submit(run_relax_job_v13,args,root,fg,v,variants[v],t,fgids,ref_ids))
            rrows=[]
            for i,f in enumerate(as_completed(rtasks),1):
                try:rrows.append(f.result())
                except Exception as e: logging.exception("RELAX job failed: %s",e)
                if i%10==0 or i==len(rtasks): logging.info("RELAX settings: %d/%d",i,len(rtasks))
        relax_df=pd.DataFrame(rrows)
        if not relax_df.empty:
            relax_summary,relax_df=summarize_relax_v13(relax_df)
    integrated=integrate_evidence_v13(paml_summary,relax_summary)
    return comp,paml_summary,relax_df,relax_summary,integrated,per_fg


# ----------------------------- Simulation calibration -------------------------

SENSE_CODONS=[a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG" if a+b+c not in STOP_CODONS]
AA_BY_CODON={c:str(Seq(c).translate()) for c in SENSE_CODONS}
TRANSITIONS={frozenset(("A","G")),frozenset(("C","T"))}
CODON_NEIGHBORS={}
for c in SENSE_CODONS:
    z=[]
    for j in range(3):
        for b in "ACGT":
            if b==c[j]: continue
            d=c[:j]+b+c[j+1:]
            if d not in AA_BY_CODON: continue
            mut=2.0 if frozenset((c[j],b)) in TRANSITIONS else 1.0
            nonsyn=AA_BY_CODON[c]!=AA_BY_CODON[d]
            z.append((d,mut,nonsyn))
    CODON_NEIGHBORS[c]=z
_CODON_NORM=sum(sum(mut*(1.0 if ns else 1.0) for _,mut,ns in CODON_NEIGHBORS[c]) for c in SENSE_CODONS)/len(SENSE_CODONS)


def simulate_codon_segment_v13(codon: str, length: float, omega: float, rng: random.Random) -> str:
    state=codon if codon in CODON_NEIGHBORS else rng.choice(SENSE_CODONS); t=0.0
    while True:
        rates=[mut*(omega if ns else 1.0)/_CODON_NORM for _,mut,ns in CODON_NEIGHBORS[state]]
        total=sum(rates)
        if total<=0: return state
        wait=rng.expovariate(total)
        if t+wait>=length: return state
        t+=wait; x=rng.random()*total; acc=0.0
        for (dest,_,_),rate in zip(CODON_NEIGHBORS[state],rates):
            acc+=rate
            if x<=acc:
                state=dest; break


def evolve_sequence_branch_v13(seq: str, length: float, omega: float, rng: random.Random,
                               selected: Optional[Set[int]]=None, selected_omega: Optional[float]=None) -> str:
    cod=[seq[i:i+3] for i in range(0,len(seq)-2,3)]; out=[]
    selected=selected or set()
    for i,c in enumerate(cod):
        w=selected_omega if i in selected and selected_omega is not None else omega
        out.append(simulate_codon_segment_v13(c,length,w,rng))
    return "".join(out)


def inject_terminal_lesions_v13(seq: str, n_stops: int, n_indels: int, rng: random.Random) -> str:
    s=list(seq); ncod=len(s)//3
    eligible=list(range(5,max(6,ncod-5))); rng.shuffle(eligible)
    for idx in eligible[:max(0,n_stops)]:
        stop=rng.choice(sorted(STOP_CODONS)); s[3*idx:3*idx+3]=list(stop)
    for _ in range(max(0,n_indels)):
        if len(s)>30:
            pos=rng.randrange(15,len(s)-15); del s[pos]
    return "".join(s)


def simulate_tree_alignment_v13(tree: Tree, target_desc: Set[str], root_seq: str, scenario: str,
                                 rng: random.Random, bg_omega: float, pos_omega: float,
                                 pos_fraction: float, loss_fraction: float,
                                 default_bl: float, lesion_stops: int, lesion_indels: int) -> Dict[str,str]:
    target=find_exact_clade(tree,target_desc)
    if target is None: raise ValueError("Simulation target branch absent")
    ncod=len(root_seq)//3; selected=set(rng.sample(range(ncod),max(1,int(round(ncod*pos_fraction))))) if scenario=="positive_selection" else set()
    tips={}
    def walk(c,parent_seq):
        for ch in c.clades:
            bl=fnum(ch.branch_length,default_bl)
            if not math.isfinite(bl) or bl<=0: bl=default_bl
            if ch is target and scenario.startswith("relaxation"):
                a=evolve_sequence_branch_v13(parent_seq,bl*loss_fraction,bg_omega,rng)
                child=evolve_sequence_branch_v13(a,bl*(1-loss_fraction),1.0,rng)
            elif ch is target and scenario=="positive_selection":
                child=evolve_sequence_branch_v13(parent_seq,bl,bg_omega,rng,selected,pos_omega)
            else:
                child=evolve_sequence_branch_v13(parent_seq,bl,bg_omega,rng)
            if not ch.clades:
                if ch is target and scenario=="relaxation_with_lesions":
                    child=inject_terminal_lesions_v13(child,lesion_stops,lesion_indels,rng)
                tips[ch.name]=child
            else: walk(ch,child)
    walk(tree.root,root_seq)
    if not tree.root.clades and tree.root.name: tips[tree.root.name]=root_seq
    return tips


def run_sim_macse_v13(args, d: Path, seqs: Dict[str,str], fg_ids: Set[str]) -> Path:
    base,_=resolve_macse(args); rel=d/"reliable.fa"; less=d/"less.fa"; nt=d/"macse_NT.fa"; aa=d/"macse_AA.fa"
    with open(rel,"w") as a,open(less,"w") as b:
        for sid,seq in seqs.items():
            fh=b if sid in fg_ids else a; fh.write(f">{sid}\n{wrap_fasta(seq)}\n")
    cmd=base+["-prog","alignSequences","-seq",str(rel),"-out_NT",str(nt),"-out_AA",str(aa)]
    if fasta_nonempty(less): cmd += ["-seq_lr",str(less),"-fs_lr",str(args.fs_lr),"-stop_lr",str(args.stop_lr)]
    run_cmd(cmd,check=True,log_path=d/"macse.log"); return nt


def simulation_scenarios_v13(args) -> List[Tuple[str,float]]:
    fr=[]
    for x in str(args.simulation_loss_fractions).split(","):
        try: fr.append(float(x))
        except Exception: pass
    out=[("constrained",1.0),("positive_selection",1.0)]
    for f in fr:
        if 0<=f<=1:
            out.append(("relaxation",f)); out.append(("relaxation_with_lesions",f))
    return out


def run_simulations_v13(args, root: Path, ref_cds: str, foregrounds: List[Foreground], per_fg: Dict[str,Any],
                         bg_species: List[str], species_to_id: Dict[str,str], pig_id: str) -> Tuple[pd.DataFrame,pd.DataFrame]:
    if args.skip_simulations or args.simulation_replicates<=0: return pd.DataFrame(),pd.DataFrame()
    codeml=exe(args.codeml,"codeml");
    if not codeml: return pd.DataFrame(),pd.DataFrame()
    by_label={norm_taxon(f.species_label):f for f in foregrounds if f.kind=="terminal"}
    fg=by_label.get(norm_taxon(args.simulation_focus))
    if fg is None or fg.foreground_id not in per_fg:
        logging.warning("Simulation focus %s is unavailable; skipping calibration",args.simulation_focus); return pd.DataFrame(),pd.DataFrame()
    tree,order,fgids,actual,_=per_fg[fg.foreground_id]
    # Use an intact reference-derived root sequence of requested size.
    codons=[ref_cds[i:i+3] for i in range(0,len(ref_cds)-2,3) if is_valid_sense_codon(ref_cds[i:i+3])]
    if not codons: raise RuntimeError("No valid reference codons for simulation")
    n=min(len(codons),max(120,args.simulation_codons)); root_seq="".join(codons[:n])
    simroot=root/"simulation_calibration"/safe_name(fg.foreground_id); simroot.mkdir(parents=True,exist_ok=True)
    rows=[]; rng_master=random.Random(args.simulation_seed)
    scenarios=simulation_scenarios_v13(args)
    for scenario,loss_fraction in scenarios:
        for rep in range(1,args.simulation_replicates+1):
            seed=rng_master.randrange(1,2**31-1); rng=random.Random(seed)
            d=simroot/f"{scenario}_f{loss_fraction:.2f}"/f"rep_{rep:04d}"; d.mkdir(parents=True,exist_ok=True)
            try:
                seqs=simulate_tree_alignment_v13(tree,fgids,root_seq,scenario,rng,args.simulation_background_omega,
                                                 args.simulation_positive_omega,args.simulation_positive_site_fraction,
                                                 loss_fraction,args.simulation_default_branch_length,
                                                 args.simulation_lesion_stops,args.simulation_lesion_indels)
                aln=run_sim_macse_v13(args,d,seqs,fgids)
                recs=read_alignment(aln); vv,st=build_paml_variant(recs,0.50,0,False); vpath=d/"baseline_mask50.fa"; write_fasta_dict(vv,vpath)
                if st["output_codons"]<args.min_paml_codons:
                    rows.append({"scenario":scenario,"loss_fraction":loss_fraction,"replicate":rep,"seed":seed,"status":"too_few_codons"}); continue
                r=run_one_paml_setting_v13(args,d,codeml,fg,"baseline_mask50",vpath,tree,order,fgids,2)
                r.update({"scenario":scenario,"loss_fraction":loss_fraction,"replicate":rep,"seed":seed,"simulated_codons":n,
                          "paml_mask_codons":st["output_codons"]})
                # For simulation calibration the raw single-replicate branch-site p is the relevant statistic.
                r["simulation_positive_call"]=bool(str(r.get("status"))=="ok" and fnum(r.get("branch_site_p_chi1"),1)<0.05 and inum(r.get("BEB_sites_ge_0.95"),0)>0)
                if args.simulation_run_relax and not args.skip_relax and exe(args.hyphy,"hyphy"):
                    refs={species_to_id[x] for x in bg_species if x in species_to_id}|{pig_id}
                    rr=run_relax_job_v13(args,d,fg,"sim_baseline_mask50",vpath,tree,fgids,refs)
                    r["RELAX_status"]=rr.get("status"); r["RELAX_K"]=rr.get("K"); r["RELAX_p"]=rr.get("p"); r["RELAX_LRT"]=rr.get("LRT")
                rows.append(r)
            except Exception as e:
                logging.exception("Simulation replicate failed %s %.2f rep %d",scenario,loss_fraction,rep)
                rows.append({"scenario":scenario,"loss_fraction":loss_fraction,"replicate":rep,"seed":seed,"status":"failed","error":str(e)})
    df=pd.DataFrame(rows); write_tsv(df,simroot/"simulation_replicates.tsv")
    sums=[]
    if not df.empty:
        for (sc,lf),g in df.groupby(["scenario","loss_fraction"]):
            ok=g[g.status=="ok"]
            lrt=pd.to_numeric(ok.get("branch_site_LRT"),errors="coerce") if len(ok) else pd.Series(dtype=float)
            w=pd.to_numeric(ok.get("foreground_omega"),errors="coerce") if len(ok) else pd.Series(dtype=float)
            pos=pd.Series(ok.get("simulation_positive_call",False)).astype(bool) if len(ok) else pd.Series(dtype=bool)
            sums.append({"scenario":sc,"loss_fraction":lf,"replicates_requested":len(g),"replicates_valid":len(ok),
                         "branch_site_positive_rate":float(pos.mean()) if len(pos) else float("nan"),
                         "branch_site_LRT_median":float(lrt.median()) if len(lrt) else float("nan"),
                         "branch_site_LRT_q95":float(lrt.quantile(.95)) if len(lrt) else float("nan"),
                         "omega_median":float(w.median()) if len(w) else float("nan")})
    sm=pd.DataFrame(sums); write_tsv(sm,simroot/"simulation_summary.tsv")
    return df,sm



def run_global_models_v13(args, root: Path, macse_outputs: Dict[str,Path], base_tree: Tree, pig_id: str) -> Tuple[pd.DataFrame,pd.DataFrame]:
    """Strict descriptive M0/free-ratio models on the full selected context.

    These are descriptive only. A failed free-ratio run is reported as failed and
    is never used to rescue or override a foreground inference.
    """
    codeml=exe(args.codeml,"codeml")
    if not codeml: return pd.DataFrame(),pd.DataFrame()
    raw=read_alignment(macse_outputs["baseline"]); vv,st=build_paml_variant(raw,0.50,0,False)
    if st["output_codons"]<args.min_paml_codons: return pd.DataFrame([{"status":"too_few_codons","codons":st["output_codons"]}]),pd.DataFrame()
    d=root/"paml"/"global_strict"; vf=d/"global_mask50.fa"; write_fasta_dict(vv,vf)
    t=prune_with_pig_policy(base_tree,set(vv),pig_id); order=[x.name for x in t.get_terminals() if x.name in vv]
    rows=[]; branches=[]
    for cf in [int(x) for x in str(args.codon_freqs).split(",") if x.strip()]:
        for name,model in [("M0",0),("free_ratio",1)]:
            pp=run_codeml_job_v13(codeml,d/f"F{cf}"/name,vv,order,t,None,name,model,0,cf,0,.4,0,args.resume)
            rows.append({"codon_freq":cf,"model":name,"status":pp.get("status"),"returncode":pp.get("returncode"),
                         "lnL":pp.get("lnL"),"np":pp.get("np"),"AIC":pp.get("AIC"),"alignment_codons":st["output_codons"]})
            if name=="free_ratio" and pp.get("status")=="ok" and pp.get("branch_table"):
                b=pd.DataFrame(pp["branch_table"]); b["codon_freq"]=cf; branches.append(b)
    return pd.DataFrame(rows),pd.concat(branches,ignore_index=True) if branches else pd.DataFrame()


def calibrate_simulations_v13(comp: pd.DataFrame, sim_df: pd.DataFrame, sim_summary: pd.DataFrame,
                               focus_label: str) -> pd.DataFrame:
    if sim_summary is None or sim_summary.empty: return sim_summary
    out=sim_summary.copy(); obs=float("nan")
    if comp is not None and not comp.empty:
        x=comp[(comp.foreground_label.astype(str).str.lower()==str(focus_label).lower()) &
               (comp.variant=="baseline_mask50") & (pd.to_numeric(comp.codon_freq,errors="coerce")==2) & (comp.status=="ok")]
        if len(x): obs=fnum(x.iloc[0].get("branch_site_LRT"))
    out["observed_focus_branch_site_LRT"]=obs; out["empirical_tail_P_LRT_ge_observed"]=float("nan")
    if math.isfinite(obs) and sim_df is not None and not sim_df.empty:
        for i,r in out.iterrows():
            g=sim_df[(sim_df.scenario==r.scenario) & (pd.to_numeric(sim_df.loss_fraction,errors="coerce")==fnum(r.loss_fraction)) & (sim_df.status=="ok")]
            vals=[fnum(x) for x in g.get("branch_site_LRT",[]) if math.isfinite(fnum(x))]
            if vals: out.loc[i,"empirical_tail_P_LRT_ge_observed"]=(1+sum(x>=obs for x in vals))/(len(vals)+1)
    return out

def write_report_v13(root: Path, integrated: pd.DataFrame, paml_comp: pd.DataFrame, relax_df: pd.DataFrame,
                      relax_summary: pd.DataFrame, intact_qc: pd.DataFrame, sim_summary: pd.DataFrame,
                      global_models: pd.DataFrame, global_branches: pd.DataFrame, intact_site: pd.DataFrame,
                      metadata: Dict[str,Any]) -> Path:
    def table(df,cols=None,n=500):
        if df is None or df.empty: return "<p><i>No rows.</i></p>"
        x=df.copy();
        if cols: x=x[[c for c in cols if c in x.columns]]
        return x.head(n).to_html(index=False,escape=True,border=1)
    cow=integrated[integrated.foreground_label.astype(str).str.lower()=="bos taurus"] if integrated is not None and not integrated.empty else pd.DataFrame()
    cowhtml=table(cow) if len(cow) else "<p>No cow row.</p>"
    parts=["<!doctype html><html><head><meta charset='utf-8'><title>GPRC6A rigorous evolution v1.3</title>",
           "<style>body{font-family:Arial;margin:30px;line-height:1.45}table{border-collapse:collapse;font-size:12px;display:block;overflow:auto;max-height:700px}th,td{padding:4px 6px;border:1px solid #ccc}section{margin:24px 0}code{background:#eee;padding:2px 4px}</style></head><body>",
           "<h1>GPRC6A rigorous codon-evolution analysis v1.3.0</h1>",
           "<p>This report preserves local ruminant topology for every foreground, creates codon masks after taxon subsetting, rejects failed codeml models before likelihood-ratio tests, adds HyPhy RELAX, and calibrates branch-site behavior with simulations.</p>",
           "<section><h2>Bos taurus</h2>",cowhtml,"</section>",
           "<section><h2>Integrated lineage interpretation</h2>",table(integrated,n=1000),"</section>",
           "<section><h2>HyPhy RELAX summary</h2>",table(relax_summary,n=1000),"</section>",
           "<section><h2>Intact-only MACSE QC</h2>",table(intact_qc,n=100),"</section>",
           "<section><h2>Simulation calibration</h2><p>Constrained simulations estimate false-positive behavior. Relaxation simulations test whether pseudogenization/relaxation alone can mimic a branch-site positive-selection signal. Positive-selection simulations estimate power. The empirical tail column compares the observed focus-lineage branch-site LRT with the simulated scenario.</p>",table(sim_summary,n=200),"</section>",
           "<section><h2>Global strict PAML models</h2>",table(global_models,n=100),"<h3>Free-ratio branch dN/dS</h3>",table(global_branches,n=1000),"</section>",
           "<section><h2>Intact-only PAML site models</h2>",table(intact_site,n=100),"</section>",
           "<section><h2>Strict PAML model comparisons</h2>",table(paml_comp,n=5000),"</section>",
           "<section><h2>RELAX per-setting results</h2>",table(relax_df,n=3000),"</section>",
           "<section><h2>Run metadata</h2><pre>",html.escape(json.dumps(metadata,indent=2,default=str)),"</pre></section></body></html>"]
    p=root/"gprc6a_codon_evolution_rigorous_report.html"; p.write_text("".join(parts)); return p


def run_intact_site_models_v13(args, root: Path, intact_alignment: Path, base_tree: Tree,
                                reference_ids: Set[str], pig_id: str) -> pd.DataFrame:
    codeml=exe(args.codeml,"codeml");
    if not codeml: return pd.DataFrame()
    raw=read_alignment(intact_alignment); vv,st=build_paml_variant(raw,1.0,0,True)
    if st["output_codons"]<args.min_paml_codons:
        vv,st=build_paml_variant(raw,0.80,0,False)
    d=root/"paml"/"intact_only_site_models"; v=d/"intact_paml_safe.fa"; write_fasta_dict(vv,v)
    t=prune_with_pig_policy(base_tree,set(vv),pig_id); order=[x.name for x in t.get_terminals() if x.name in vv]
    rows=[]
    for cf in [int(x) for x in str(args.codon_freqs).split(",") if x.strip()]:
        parsed={}
        for name,model,ns,fixw,w in [("M1a",0,1,0,.4),("M2a",0,2,0,1.5),("M7",0,7,0,.4),("M8",0,8,0,1.5)]:
            p=run_codeml_job_v13(codeml,d/f"F{cf}"/name,vv,order,t,None,name,model,ns,cf,fixw,w,0,args.resume,10)
            parsed[name]=p; rows.append({"codon_freq":cf,"model":name,"status":p.get("status"),"returncode":p.get("returncode"),"lnL":p.get("lnL"),"np":p.get("np"),"AIC":p.get("AIC")})
        for alt,null,label in [("M2a","M1a","M1a_vs_M2a"),("M8","M7","M7_vs_M8")]:
            a,n=parsed[alt],parsed[null]
            if valid_model_set_v13([a,n]):
                lr=max(0,2*(fnum(a["lnL"])-fnum(n["lnL"]))); rows.append({"codon_freq":cf,"model":label,"status":"comparison_ok","LRT":lr,"p":chi2_sf(lr,2)})
            else: rows.append({"codon_freq":cf,"model":label,"status":"comparison_not_run_model_failed"})
    return pd.DataFrame(rows)


def main_v13() -> int:
    args=parse_args(); args.outdir=args.outdir.expanduser().resolve(); root=args.outdir/args.evolution_subdir
    # Never silently write into the old v1.2 folder.
    if root.resolve()==(args.outdir/"codon_evolution").resolve():
        raise RuntimeError("v1.3 refuses to use the previous codon_evolution folder. Choose a new --evolution-subdir.")
    if root.exists() and any(root.iterdir()) and not args.resume:
        raise RuntimeError(f"New v1.3 output directory already exists and is non-empty: {root}. Use --resume or choose another --evolution-subdir.")
    root.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=getattr(logging,args.log_level),format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr),logging.FileHandler(root/"codon_evolution_v13.log")],force=True)
    logging.info("GPRC6A rigorous codon-evolution follow-up v%s",VERSION)
    ref_path=args.reference_exons.expanduser().resolve() if args.reference_exons else args.outdir/"reference_exons.fa"
    ref_exons,ref_cds=load_reference_exons(ref_path)
    if args.tree: tree_path=args.tree.expanduser().resolve()
    else:
        h=args.outdir/"external_outgroups"/"combined_species_tree.newick"; tree_path=h if h.exists() else args.outdir/args.refine_subdir/"species_tree.newick"
    if not tree_path.exists(): raise FileNotFoundError(tree_path)
    upstream=discover_upstream(args.outdir,args.refine_subdir); ml=load_ml_species_sequences(args.outdir,args.ml_subdir)
    assemblies=build_assembly_sequences(args.outdir,upstream,ref_exons,ml); reps=choose_species_representatives(assemblies)
    external_path=args.external_intact_fasta or (args.outdir/"external_outgroups"/"external_intact_gprc6a.fasta")
    external=load_external_intact_fasta(external_path if Path(external_path).exists() else None,ref_cds)
    overlap=set(reps)&set(external)
    if overlap: raise RuntimeError("External FASTA duplicates reconstructed species: "+", ".join(sorted(overlap)))
    reps.update(external); write_sequence_exports(root,assemblies,reps,ref_cds)

    # Preliminary standardized tree uses safe IDs.
    prelim_map={sp:safe_name(sp) for sp in reps}; prelim_tree,_=standardized_tree(tree_path,prelim_map)
    allowed=set(LOSS_SPECIES_STRONG)|(LOSS_SPECIES_EXPLORATORY if args.include_possible_loss else set())
    focus=[sp for sp,r in reps.items() if r.species_status in allowed and r.callable_fraction>=args.min_focus_callable and (args.include_hybrids or not is_hybrid_species(sp))]
    bynorm={norm_taxon(sp):sp for sp in reps}
    for x in args.always_test:
        sp=bynorm.get(norm_taxon(x))
        if sp and reps[sp].callable_fraction>=args.min_focus_callable and sp not in focus and (args.include_hybrids or not is_hybrid_species(sp)): focus.append(sp)
    bg_diag=root/"sequence_recovery"/"background_candidate_diagnostics.tsv"
    bg=select_intact_background(reps,prelim_tree,focus,args.min_intact_callable,args.max_intact_background,args.min_intact_background,
                               args.min_fallback_intact_callable,not args.no_conflicting_background_fallback,bg_diag)
    if len(bg)<3: raise RuntimeError("Need at least 3 independently intact background taxa")
    context=select_ruminant_context_v13(reps,prelim_tree,focus,args.context_min_callable,args.include_hybrids,args.context_mode,args.context_anchor_count)
    # Align an additional low-threshold guardian pool so exact sister lineages
    # remain available when a focal branch would otherwise collapse after
    # pruning. These taxa are not automatically admitted to the PAML context.
    external_species={sp for sp,r in reps.items() if r.species_status=="EXTERNAL_INTACT"}
    guardian_pool=[]
    prelim_tips={t.name for t in prelim_tree.get_terminals()}
    for sp,r in reps.items():
        if sp in external_species: continue
        if not args.include_hybrids and is_hybrid_species(sp): continue
        if safe_name(sp) not in prelim_tips: continue
        if r.callable_fraction>=args.guardian_min_callable: guardian_pool.append(sp)
    alignment_ruminants=sorted(set(context)|set(guardian_pool)|set(focus))
    inputs=write_macse_inputs_v13(root,reps,bg,alignment_ruminants,focus,ref_cds)
    species_to_id=inputs["species_to_id"]; pig_id=inputs["pig_id"]
    base_tree,_=standardized_tree(tree_path,species_to_id)
    # Remove unmatched tree tips before defining exact foreground clades. The
    # standardized_tree helper renames matched tips but intentionally leaves
    # unmatched tips in place; exact branch definitions require the analysis
    # tree to contain only sequences actually available to MACSE/PAML/RELAX.
    base_tree=prune_with_pig_policy(base_tree,set(species_to_id.values())|{pig_id},pig_id)

    # Suggestion 3a: intact-only alignment QC before adding pseudogenes.
    intact_aln,intact_qc,intact_ok=run_intact_only_macse_qc_v13(args,root,inputs["reliable"])
    # Main full-context MACSE alignment, then masks are built separately for each foreground.
    macse_outputs=run_macse(args,root,inputs)
    foregrounds=build_foregrounds_v13(upstream,reps,species_to_id,base_tree,pig_id,args.include_possible_loss,args.always_test,args.min_focus_callable)
    write_tsv(pd.DataFrame([asdict(x) for x in foregrounds]),root/"tables"/"foreground_definitions_v13.tsv")

    comp,paml_summary,relax_df,relax_summary,integrated,per_fg=run_paml_and_relax_v13(
        args,root,macse_outputs,base_tree,reps,species_to_id,pig_id,bg,context,foregrounds)
    intact_site=run_intact_site_models_v13(args,root,intact_aln,base_tree,{species_to_id[x] for x in bg if x in species_to_id}|{pig_id},pig_id)
    global_models,global_branches=run_global_models_v13(args,root,macse_outputs,base_tree,pig_id)
    sim_df,sim_summary=run_simulations_v13(args,root,ref_cds,foregrounds,per_fg,bg,species_to_id,pig_id)
    sim_summary=calibrate_simulations_v13(comp,sim_df,sim_summary,args.simulation_focus)
    if sim_summary is not None and not sim_summary.empty:
        write_tsv(sim_summary,root/"simulation_calibration"/"simulation_summary_calibrated.tsv")

    td=root/"tables"; td.mkdir(parents=True,exist_ok=True)
    write_tsv(comp,td/"foreground_strict_paml_comparisons.tsv"); write_tsv(paml_summary,td/"foreground_paml_summary.tsv")
    write_tsv(relax_df,td/"hyphy_relax_all_settings.tsv"); write_tsv(relax_summary,td/"hyphy_relax_summary.tsv")
    write_tsv(integrated,td/"integrated_evolutionary_interpretation.tsv"); write_tsv(intact_site,td/"intact_only_site_models.tsv")
    write_tsv(global_models,td/"global_strict_paml_models.tsv"); write_tsv(global_branches,td/"global_free_ratio_branches.tsv")
    metadata={"pipeline_version":VERSION,"output_directory":str(root),"previous_output_preserved":str(args.outdir/"codon_evolution"),
              "tree":str(tree_path),"intact_background":bg,"ruminant_context_count":len(context),"ruminant_context":context,
              "aligned_topology_guardian_pool_count":len(alignment_ruminants),"guardian_min_callable":args.guardian_min_callable,
              "foreground_count":len(foregrounds),"foreground_kinds":dict(Counter(f.kind for f in foregrounds)),
              "intact_only_macse_qc_passed":intact_ok,"codeml":exe(args.codeml,"codeml"),"hyphy":exe(args.hyphy,"hyphy"),
              "macse":exe(args.macse,"macse") if not args.macse_jar else args.macse_jar,
              "simulation_replicates_per_scenario":0 if args.skip_simulations else args.simulation_replicates,
              "six_rigorous_changes":["topology_preserving_ruminant_context","terminal_ancestral_and_inferred_loss_foregrounds",
                                       "per_foreground_post_subset_codon_masks_and_intact_only_MACSE_QC","strict_codeml_failure_gating",
                                       "HyPhy_RELAX","simulation_calibration_H0_relaxation_positive_selection"]}
    write_json(metadata,root/"run_metadata_v13.json")
    report=write_report_v13(root,integrated,comp,relax_df,relax_summary,intact_qc,sim_summary,global_models,global_branches,intact_site,metadata)
    summary={"pipeline_version":VERSION,"report":str(report),"output_directory":str(root),"integrated_counts":dict(Counter(integrated.integrated_interpretation)) if not integrated.empty else {},
             "simulation_summary":str(root/"simulation_calibration")}
    write_json(summary,root/"summary_v13.json")
    logging.info("Completed rigorous v1.3 analysis. Report: %s",report); return 0

def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Rigorous topology-preserving MACSE + PAML + HyPhy RELAX + simulation follow-up for GPRC6A.")
    p.add_argument("--outdir", required=True, type=Path, help="Root output directory from the first-pass/refiner pipeline")
    p.add_argument("--refine-subdir", default="refined_intactness")
    p.add_argument("--ml-subdir", default="ml_gene_recovery")
    p.add_argument("--evolution-subdir", default="codon_evolution_rigorous_v1_3", help="NEW output subdirectory; v1.3 never writes into the previous codon_evolution directory")
    p.add_argument("--reference-exons", type=Path, default=None, help="Pig GPRC6A exon FASTA; defaults to OUTDIR/reference_exons.fa")
    p.add_argument("--tree", type=Path, default=None, help="Species Newick. Defaults to refined_intactness/species_tree.newick if present")
    p.add_argument("--external-intact-fasta", type=Path, default=None,
                   help="Curated intact outgroup GPRC6A CDS FASTA. IDs should be Taxon_with_underscores|ACCESSION. Use with a --tree that contains those taxa.")

    p.add_argument("--macse", default=None, help="MACSE executable/wrapper; auto-detected as `macse` if omitted")
    p.add_argument("--macse-jar", default=None, help="MACSE jar, alternative to --macse")
    p.add_argument("--java", default=None)
    p.add_argument("--java-memory", default="4g")
    p.add_argument("--fs-lr", type=float, default=10.0, help="MACSE frameshift cost for less-reliable/pseudogene sequences")
    p.add_argument("--stop-lr", type=float, default=10.0, help="MACSE internal-stop cost for less-reliable/pseudogene sequences")
    p.add_argument("--skip-macse-refine", action="store_true", help="Skip the refined-alignment sensitivity replicate")
    p.add_argument("--macse-refine-iter", type=int, default=3)
    p.add_argument("--skip-macse-alt-cost", action="store_true", help="Skip the second de-novo MACSE alignment with altered pseudogene penalties")
    p.add_argument("--alt-fs-lr", type=float, default=15.0, help="Alternate MACSE less-reliable frameshift cost used only for robustness")
    p.add_argument("--alt-stop-lr", type=float, default=15.0, help="Alternate MACSE less-reliable internal-stop cost used only for robustness")

    p.add_argument("--codeml", default=None, help="PAML codeml executable")
    p.add_argument("--codon-freqs", default="1,2", help="Comma-separated PAML CodonFreq values; 1=F1x4, 2=F3x4, 3=Fcodon")
    p.add_argument("--analysis-depth", choices=["extensive", "exhaustive"], default="extensive",
                   help="Extensive runs all codon-frequency models on the primary alignment and F3x4 across perturbations; exhaustive crosses all")
    p.add_argument("--paml-workers", type=int, default=2, help="Concurrent single-threaded codeml jobs")
    p.add_argument("--min-paml-codons", type=int, default=120, help="Minimum alignment codons after filtering")
    p.add_argument("--min-foreground-codons", type=int, default=100, help="Minimum valid codons in every foreground taxon")

    p.add_argument("--min-intact-callable", type=float, default=0.85,
                   help="Preferred callable fraction for strict intact background taxa")
    p.add_argument("--min-fallback-intact-callable", type=float, default=0.70,
                   help="Lowest callable fraction allowed when strict intact controls are too sparse")
    p.add_argument("--no-conflicting-background-fallback", action="store_true",
                   help="Never use a CONFLICTING_ASSEMBLIES species as an intact control, even when its selected assembly is strongly intact")
    p.add_argument("--min-intact-background", type=int, default=6)
    p.add_argument("--max-intact-background", type=int, default=20)
    p.add_argument("--min-focus-callable", type=float, default=0.35)
    p.add_argument("--include-possible-loss", action="store_true", help="Also test DISRUPTED_POSSIBLE species as foregrounds")
    p.add_argument("--always-test", action="append", default=["Bos taurus"], help="Taxon to test even if not called lost; repeatable")
    p.add_argument("--context-min-callable", type=float, default=0.70,
                   help="Minimum callable fraction for ruminant phylogenetic context taxa; focal foregrounds and topology guardians can be retained below this threshold")
    p.add_argument("--context-mode", choices=["all_callable", "sister_anchors"], default="all_callable",
                   help="How much ruminant context to retain. all_callable is the rigorous default; sister_anchors is a faster fallback")
    p.add_argument("--context-anchor-count", type=int, default=12,
                   help="Maximum general ruminant anchors in sister_anchors mode, in addition to exact sister/topology guardians")
    p.add_argument("--guardian-min-callable", type=float, default=0.25,
                   help="Minimum callable fraction for extra ruminants aligned only as potential topology guardians; they are not automatically used as PAML context")
    p.add_argument("--include-hybrids", action="store_true",
                   help="Include hybrid taxa in phylogenetic selection analyses. Off by default because a bifurcating species tree is inappropriate for hybrids")
    p.add_argument("--intact-qc-max-frameshift", type=int, default=0,
                   help="Maximum MACSE frameshift characters allowed in any independently intact/reference sequence in the intact-only QC alignment")
    p.add_argument("--intact-qc-max-internal-stop", type=int, default=0,
                   help="Maximum internal stop codons allowed in any intact/reference sequence after intact-only MACSE alignment")
    p.add_argument("--allow-intact-qc-warning", action="store_true",
                   help="Continue despite intact-only MACSE QC failure. Not recommended for final inference")

    p.add_argument("--hyphy", default=None, help="HyPhy executable. Required for RELAX unless --skip-relax")
    p.add_argument("--skip-relax", action="store_true", help="Skip HyPhy RELAX analyses")
    p.add_argument("--hyphy-workers", type=int, default=2, help="Concurrent RELAX jobs")
    p.add_argument("--relax-variants", choices=["primary", "all_mask50"], default="all_mask50",
                   help="Run RELAX on baseline_mask50 only or on baseline/refined/altcost mask50 variants")
    p.add_argument("--relax-models", choices=["Minimal", "All"], default="Minimal")

    p.add_argument("--skip-simulations", action="store_true", help="Skip the simulation-calibration module")
    p.add_argument("--simulation-focus", default="Bos taurus", help="Foreground lineage used for calibration simulations")
    p.add_argument("--simulation-replicates", type=int, default=10,
                   help="Replicates per simulation scenario. Use 50-100+ for publication-grade calibration")
    p.add_argument("--simulation-codons", type=int, default=600)
    p.add_argument("--simulation-seed", type=int, default=1729)
    p.add_argument("--simulation-background-omega", type=float, default=0.30)
    p.add_argument("--simulation-positive-omega", type=float, default=2.50)
    p.add_argument("--simulation-positive-site-fraction", type=float, default=0.10)
    p.add_argument("--simulation-loss-fractions", default="0.25,0.50,0.75",
                   help="Fraction of the foreground branch spent under functional constraint before relaxation to omega=1")
    p.add_argument("--simulation-default-branch-length", type=float, default=0.05,
                   help="Branch length used when the supplied tree lacks a positive branch length")
    p.add_argument("--simulation-lesion-stops", type=int, default=2,
                   help="Internal stop codons injected in the relaxation_with_lesions terminal-foreground scenario")
    p.add_argument("--simulation-lesion-indels", type=int, default=2,
                   help="Single-nucleotide deletions injected in the relaxation_with_lesions terminal-foreground scenario")
    p.add_argument("--simulation-run-relax", action="store_true",
                   help="Also run HyPhy RELAX on every simulated replicate. Computationally expensive")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--log-level", choices=["DEBUG","INFO","WARNING","ERROR"], default="INFO")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.outdir = args.outdir.expanduser().resolve()
    root = args.outdir / args.evolution_subdir
    root.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=getattr(logging, args.log_level),
                        format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr), logging.FileHandler(root / "codon_evolution.log")])
    logging.info("GPRC6A codon-evolution follow-up v%s", VERSION)
    if args.stop_lr > 2 * args.fs_lr:
        logging.warning("Primary MACSE stop_lr (%s) exceeds 2*fs_lr (%s); MACSE may prefer two frameshifts over a stop", args.stop_lr, args.fs_lr)
    if args.alt_stop_lr > 2 * args.alt_fs_lr:
        logging.warning("Alternate MACSE stop_lr (%s) exceeds 2*fs_lr (%s); MACSE may prefer two frameshifts over a stop", args.alt_stop_lr, args.alt_fs_lr)

    ref_path = args.reference_exons.expanduser().resolve() if args.reference_exons else args.outdir / "reference_exons.fa"
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference exon FASTA not found: {ref_path}")
    ref_exons, ref_cds = load_reference_exons(ref_path)

    if args.tree:
        tree_path = args.tree.expanduser().resolve()
    else:
        helper_tree = args.outdir / "external_outgroups" / "combined_species_tree.newick"
        tree_path = helper_tree if helper_tree.exists() else args.outdir / args.refine_subdir / "species_tree.newick"
    if not tree_path.exists():
        raise FileNotFoundError(f"Species tree not found: {tree_path}. Supply --tree.")

    upstream = discover_upstream(args.outdir, args.refine_subdir)
    ml_species = load_ml_species_sequences(args.outdir, args.ml_subdir)
    logging.info("Optional ML-recovery species CDSs found: %d", len(ml_species))
    assemblies = build_assembly_sequences(args.outdir, upstream, ref_exons, ml_species)
    reps = choose_species_representatives(assemblies)

    external_path = args.external_intact_fasta
    if external_path is None:
        helper_fa = args.outdir / "external_outgroups" / "external_intact_gprc6a.fasta"
        if helper_fa.exists():
            external_path = helper_fa
    external_reps = load_external_intact_fasta(external_path, ref_cds)
    overlap = set(reps) & set(external_reps)
    if overlap:
        raise RuntimeError(
            "External intact FASTA duplicates ruminant species already reconstructed: "
            + ", ".join(sorted(overlap))
        )
    reps.update(external_reps)

    write_sequence_exports(root, assemblies, reps, ref_cds)
    logging.info("Species representatives: %d ruminant + %d external intact = %d total",
                 len(reps)-len(external_reps), len(external_reps), len(reps))

    # Determine focus lineages before background selection.
    allowed_loss = set(LOSS_SPECIES_STRONG)
    if args.include_possible_loss:
        allowed_loss |= LOSS_SPECIES_EXPLORATORY
    focus_species = [sp for sp, r in reps.items()
                     if r.species_status in allowed_loss and r.callable_fraction >= args.min_focus_callable]
    by_norm = {norm_taxon(sp): sp for sp in reps}
    for x in args.always_test:
        sp = by_norm.get(norm_taxon(x))
        if sp and reps[sp].callable_fraction >= args.min_focus_callable and sp not in focus_species:
            focus_species.append(sp)
    if not focus_species:
        raise RuntimeError("No analyzable loss/special-focus species passed sequence-quality filters")

    # Initial tree standardization uses IDs equal to safe species names.
    all_species_ids = {sp: safe_name(sp) for sp in reps}
    base_tree, matched = standardized_tree(tree_path, all_species_ids)
    if len(matched) < 3:
        raise RuntimeError("Too few reconstructed species names matched the supplied tree. Check tip naming.")

    bg_diag = root / "sequence_recovery" / "background_candidate_diagnostics.tsv"
    bg_species = select_intact_background(
        reps, base_tree, focus_species,
        args.min_intact_callable, args.max_intact_background,
        args.min_intact_background,
        fallback_min_callable=args.min_fallback_intact_callable,
        allow_conflicting_fallback=not args.no_conflicting_background_fallback,
        diagnostics_path=bg_diag,
    )
    if len(bg_species) < 3:
        raise RuntimeError(
            "Fewer than 3 defensible intact background species are available. "
            "This ruminant data set contains no species-level intact controls, so the analysis requires "
            "independent intact outgroup CDSs. Run gprc6a_fetch_intact_outgroups.py and then rerun this "
            "script with its external_intact_gprc6a.fasta and combined_species_tree.newick outputs. "
            f"Diagnostic table: {bg_diag}"
        )

    inputs = write_macse_inputs(root, reps, bg_species, focus_species, ref_cds)
    species_to_id = inputs["species_to_id"]
    # Re-standardize tree using exactly the final IDs written to FASTA.
    base_tree, matched2 = standardized_tree(tree_path, species_to_id)
    missing_focus = [sp for sp in focus_species if species_to_id.get(sp) not in {t.name for t in base_tree.get_terminals()}]
    if missing_focus:
        logging.warning("Focus species absent from supplied tree and cannot be tested: %s", "; ".join(missing_focus))

    macse_outputs = run_macse(args, root, inputs)
    variants, alignment_qc, variant_summary = generate_alignment_variants(root, macse_outputs, args.min_paml_codons)
    if not variants:
        raise RuntimeError("No PAML-safe alignment variant retained enough codons")

    foregrounds = build_foregrounds(upstream, reps, species_to_id, args.include_possible_loss, args.always_test)
    # Retain only foregrounds that exist as one branch in at least the standardized full tree.
    tree_tips = {t.name for t in base_tree.get_terminals()}
    valid_fg = []
    for fg in foregrounds:
        ids = set(fg.descendant_ids) & tree_tips
        if not ids:
            continue
        if fg.kind == "terminal" or find_exact_clade(base_tree, ids):
            fg.descendant_ids = sorted(ids)
            valid_fg.append(fg)
        else:
            logging.warning("Skipping inferred clade foreground not represented by one tree branch: %s", fg.foreground_id)
    foregrounds = valid_fg
    write_tsv(pd.DataFrame([asdict(x) for x in foregrounds]), root / "tables" / "foreground_definitions.tsv")

    global_results, site_results, fgdf, comp = run_paml_panel(
        args, root, variants, base_tree, bg_species, species_to_id, inputs["pig_id"], foregrounds)
    fg_summary = summarize_foregrounds(comp, foregrounds)

    td = root / "tables"; td.mkdir(parents=True, exist_ok=True)
    write_tsv(global_results, td / "global_paml_models.tsv")
    write_tsv(site_results, td / "intact_site_models.tsv")
    write_tsv(fgdf, td / "foreground_all_models.tsv")
    if not comp.empty:
        c = comp.copy(); c["setting_call"] = c.apply(setting_call, axis=1)
        write_tsv(c, td / "foreground_model_comparisons.tsv")
    write_tsv(fg_summary, td / "foreground_evolution_summary.tsv")

    # Collect primary free-ratio branch table into a global machine-readable table.
    free_tables = []
    for p in (root / "paml").glob("*/F*/free_ratio_branches.tsv"):
        g = load_tsv(p)
        if not g.empty:
            g["source_path"] = str(p)
            free_tables.append(g)
    free_df = pd.concat(free_tables, ignore_index=True) if free_tables else pd.DataFrame()
    write_tsv(free_df, td / "free_ratio_branches.tsv")

    tool_meta = {
        "pipeline_version": VERSION,
        "outdir": str(args.outdir),
        "reference_exons": str(ref_path),
        "tree": str(tree_path),
        "macse_executable": exe(args.macse, "macse") if not args.macse_jar else None,
        "macse_jar": args.macse_jar,
        "codeml": exe(args.codeml, "codeml"),
        "codon_freqs": args.codon_freqs,
        "analysis_depth": args.analysis_depth,
        "macse_primary_fs_lr": args.fs_lr,
        "macse_primary_stop_lr": args.stop_lr,
        "macse_alt_fs_lr": args.alt_fs_lr,
        "macse_alt_stop_lr": args.alt_stop_lr,
        "alignment_variants": sorted(variants),
        "intact_background": bg_species,
        "focus_species": focus_species,
        "foreground_count": len(foregrounds),
        "notes": [
            "Stop/frameshift codons are preserved in raw MACSE output but masked before codeml.",
            "Missing first-pass bases are N, never reference-imputed.",
            "Bos taurus is tested explicitly when present and sufficiently recovered.",
        ],
    }
    write_json(tool_meta, root / "run_metadata.json")
    report = write_report(root, reps, bg_species, foregrounds, alignment_qc, variant_summary,
                          global_results, site_results, fgdf, comp, fg_summary, tool_meta, args)
    summary = {
        "pipeline_version": VERSION,
        "species_reconstructed": len(reps),
        "intact_background": len(bg_species),
        "foregrounds": len(foregrounds),
        "alignment_variants": len(variants),
        "report": str(report),
        "consensus_counts": dict(Counter(fg_summary.consensus_interpretation)) if not fg_summary.empty else {},
    }
    write_json(summary, root / "summary.json")
    logging.info("Completed. Report: %s", report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main_v13())
