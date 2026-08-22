#!/usr/bin/env python3
"""
GPRC6A codon-evolution follow-up pipeline
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

Typical installation with Bioconda:
  conda install -c conda-forge -c bioconda macse paml pandas biopython

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


VERSION = "1.1.0"
DNA = set("ACGT")
STOP_CODONS = {"TAA", "TAG", "TGA"}

BLASTN_COLS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qseq", "sseq",
]

INTACT_SPECIES = {"INTACT_HIGH", "INTACT_LIKELY"}
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
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


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
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=BLASTN_COLS)
    try:
        df = pd.read_csv(path, sep="\t", names=BLASTN_COLS, header=None,
                         compression="gzip" if path.suffix == ".gz" else None)
    except Exception as e:
        logging.warning("Could not read BLAST %s: %s", path, e)
        return pd.DataFrame(columns=BLASTN_COLS)
    for c in ["length", "qstart", "qend", "sstart", "send", "qlen", "slen"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    for c in ["pident", "evalue", "bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["strand"] = df.apply(lambda r: "+" if int(r.sstart) <= int(r.send) else "-", axis=1)
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
    if species_df.empty:
        return {}
    return {str(r.species_key): str(r.species_status) for _, r in species_df.iterrows()}


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
        sstat = species_status.get(sp, "UNRESOLVED")
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
                m["n_fraction"] = seq.count("N") / len(seq) if seq else 1.0
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
                             min_callable: float, max_n: int, min_n: int) -> List[str]:
    candidates = [sp for sp, r in reps.items()
                  if r.species_status in INTACT_SPECIES and r.callable_fraction >= min_callable]
    candidates.sort(key=lambda sp: reps[sp].quality_score, reverse=True)
    if len(candidates) <= max_n:
        if len(candidates) < min_n:
            logging.warning("Only %d high-quality intact background species pass filters; requested at least %d", len(candidates), min_n)
        return candidates

    # Farthest-point sampling with the closest intact relative to every focal
    # lineage forced in first. This keeps local controls and broad phylogenetic spread.
    selected: List[str] = []
    tree_tips = {t.name for t in tree.get_terminals()}
    for focal in focus_species:
        fid = safe_name(focal)
        if fid not in tree_tips:
            continue
        avail = [sp for sp in candidates if safe_name(sp) in tree_tips]
        if not avail:
            continue
        nearest = min(avail, key=lambda sp: topological_distance(tree, fid, safe_name(sp)))
        if nearest not in selected:
            selected.append(nearest)
    if not selected and candidates:
        selected.append(candidates[0])

    while len(selected) < max_n:
        remaining = [sp for sp in candidates if sp not in selected]
        if not remaining:
            break
        def diversity(sp: str) -> Tuple[float, float]:
            sid = safe_name(sp)
            ds = [topological_distance(tree, sid, safe_name(x)) for x in selected if safe_name(x) in tree_tips]
            d = min(ds) if ds else 0.0
            return d, reps[sp].quality_score
        selected.append(max(remaining, key=diversity))
    return selected[:max_n]


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
    keep_rum = set(taxa) - {pig_id}
    t = prune_tree(base_tree, keep_rum)
    present = {x.name for x in t.get_terminals()}
    if keep_rum - present:
        logging.debug("Tree missing taxa: %s", sorted(keep_rum-present))
    t = add_pig_outgroup(t, pig_id)
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
    t = prune_tree(base_tree, taxa - {pig_id})
    t = add_pig_outgroup(t, pig_id)
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
        t = prune_tree(base_tree, taxa - {pig_id})
        t = add_pig_outgroup(t, pig_id)
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
        t = prune_tree(base_tree, set(recs)-{pig_id}); t = add_pig_outgroup(t, pig_id)
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


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Codon-level MACSE + PAML follow-up for GPRC6A intact/loss lineages.")
    p.add_argument("--outdir", required=True, type=Path, help="Root output directory from the first-pass/refiner pipeline")
    p.add_argument("--refine-subdir", default="refined_intactness")
    p.add_argument("--ml-subdir", default="ml_gene_recovery")
    p.add_argument("--evolution-subdir", default="codon_evolution")
    p.add_argument("--reference-exons", type=Path, default=None, help="Pig GPRC6A exon FASTA; defaults to OUTDIR/reference_exons.fa")
    p.add_argument("--tree", type=Path, default=None, help="Species Newick. Defaults to refined_intactness/species_tree.newick if present")

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

    p.add_argument("--min-intact-callable", type=float, default=0.85)
    p.add_argument("--min-intact-background", type=int, default=6)
    p.add_argument("--max-intact-background", type=int, default=20)
    p.add_argument("--min-focus-callable", type=float, default=0.35)
    p.add_argument("--include-possible-loss", action="store_true", help="Also test DISRUPTED_POSSIBLE species as foregrounds")
    p.add_argument("--always-test", action="append", default=["Bos taurus"], help="Taxon to test even if not called lost; repeatable")
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

    tree_path = args.tree.expanduser().resolve() if args.tree else args.outdir / args.refine_subdir / "species_tree.newick"
    if not tree_path.exists():
        raise FileNotFoundError(f"Species tree not found: {tree_path}. Supply --tree.")

    upstream = discover_upstream(args.outdir, args.refine_subdir)
    ml_species = load_ml_species_sequences(args.outdir, args.ml_subdir)
    logging.info("Optional ML-recovery species CDSs found: %d", len(ml_species))
    assemblies = build_assembly_sequences(args.outdir, upstream, ref_exons, ml_species)
    reps = choose_species_representatives(assemblies)
    write_sequence_exports(root, assemblies, reps, ref_cds)
    logging.info("Species representatives: %d", len(reps))

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

    bg_species = select_intact_background(reps, base_tree, focus_species,
                                          args.min_intact_callable, args.max_intact_background,
                                          args.min_intact_background)
    if len(bg_species) < 3:
        raise RuntimeError("Fewer than 3 intact background species are available; PAML comparison would be unreliable")

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
    raise SystemExit(main())
