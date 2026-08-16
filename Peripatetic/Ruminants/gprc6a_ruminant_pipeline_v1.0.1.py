#!/usr/bin/env python3
"""
GPRC6A comparative genomics pipeline for Ruminantia.

Purpose
-------
1. Download and analyze every genome assembly listed in a target TSV.
2. BLAST pig GPRC6A coding exons against each genome.
3. Select the most coherent orthologous GPRC6A locus.
4. Reconstruct each exon from BLAST HSPs and call:
   - substitutions and codon changes
   - premature stop codons
   - start-loss / stop-loss
   - frame-shifting and in-frame insertions/deletions
   - probable exon deletions versus unresolved missing sequence
   - potential noncanonical splice-site motifs
   - assembly-gap / N-content metrics
5. Detect candidate repeat insertions and optionally annotate them with RepeatMasker.
6. Aggregate all assemblies to species-level consensus.
7. Infer shared events and relative phylogenetic chronology using a supplied Newick tree.
8. Generate machine-readable TSV/JSON outputs and a self-contained HTML report.

Important interpretation note
-----------------------------
This pipeline distinguishes a *biological disruption* from a *failure to recover sequence*.
Missing BLAST sequence is NOT automatically called a deletion. A probable deletion requires
flanking evidence in a coherent locus and an assembled, low-N genomic interval.

The six supplied pig exons are expected to concatenate into the CDS. With the supplied
XM_021089187.1 exon FASTA, the concatenated CDS is 2787 nt (929 codons) and has the expected
terminal stop codon.

External executables
--------------------
Required for a full run:
  - blastn (NCBI BLAST+)
  - datasets (NCBI Datasets CLI), unless --genome-map supplies local FASTA paths
Optional:
  - RepeatMasker (for repeat family/class annotation)
  - dustmasker (fallback low-complexity annotation)

Python packages:
  pandas, biopython

Example
-------
python gprc6a_ruminant_pipeline.py \
  --targets-tsv ruminantia_assemblies.tsv \
  --exons Sus_scrofa_XM_021089187.1_exons.fa \
  --tree ruminantia_species_tree.newick \
  --outdir GPRC6A_ruminantia \
  --workers 2 --blast-threads 8 --resume

For a small test first:
python gprc6a_ruminant_pipeline.py ... --max-assemblies 3
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import html
import json
import logging
import math
import os
import re
import shutil
import statistics
import subprocess
import sys
import tempfile
import textwrap
import time
import zipfile
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Any

import pandas as pd
from Bio import Phylo, SeqIO
from Bio.Seq import Seq


VERSION = "1.0.1"
BLAST_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    "qlen", "slen", "qseq", "sseq",
]
DNA = set("ACGT")
DNA_N = set("ACGTN")


# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "unnamed"


def anchor_id(s: str) -> str:
    return safe_name(s).lower()


def run_cmd(cmd: Sequence[str], *, cwd: Optional[Path] = None,
            stdout=None, stderr=None, check: bool = True) -> subprocess.CompletedProcess:
    logging.debug("RUN: %s", " ".join(map(str, cmd)))
    return subprocess.run(
        list(map(str, cmd)), cwd=str(cwd) if cwd else None,
        stdout=stdout, stderr=stderr, check=check, text=False,
    )


def require_exe(name: str, required: bool = True) -> Optional[str]:
    p = shutil.which(name)
    if not p and required:
        raise RuntimeError(f"Required executable not found on PATH: {name}")
    return p


def md5_text(s: str) -> str:
    return hashlib.md5(s.encode("utf-8")).hexdigest()


def median_or_nan(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if x is not None and not pd.isna(x)]
    return statistics.median(vals) if vals else float("nan")


def mean_or_nan(xs: Iterable[float]) -> float:
    vals = [float(x) for x in xs if x is not None and not pd.isna(x)]
    return statistics.mean(vals) if vals else float("nan")


def frac(n: int, d: int) -> float:
    return n / d if d else float("nan")


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def species_key(name: str) -> str:
    """Collapse subspecies/infraspecific labels to binomial; retain hybrids as explicit labels."""
    name = str(name).strip()
    if " x " in name or " × " in name:
        return name.replace(" × ", " x ")
    p = name.split()
    return " ".join(p[:2]) if len(p) >= 2 else name


def confidence_grade(score: float) -> str:
    if score >= 0.85:
        return "HIGH"
    if score >= 0.65:
        return "MEDIUM"
    if score >= 0.40:
        return "LOW"
    return "UNRESOLVED"


def quote_tsv(v: Any) -> Any:
    if isinstance(v, (dict, list, tuple)):
        return json.dumps(v, sort_keys=True)
    return v


# -----------------------------------------------------------------------------
# Reference exons / CDS map
# -----------------------------------------------------------------------------

@dataclass
class RefExon:
    exon: str
    index: int
    seq: str
    length: int
    cds_start: int   # 1-based inclusive
    cds_end: int     # 1-based inclusive
    phase_at_start: int


def load_reference_exons(path: Path) -> Tuple[List[RefExon], str]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")

    def exon_num(rec_id: str, fallback: int) -> int:
        m = re.search(r"(?:exon[_-]?)?(\d+)$", rec_id, re.I)
        return int(m.group(1)) if m else fallback

    tmp = []
    for i, r in enumerate(records, 1):
        seq = str(r.seq).upper().replace("U", "T")
        bad = sorted(set(seq) - set("ACGTN"))
        if bad:
            raise ValueError(f"Unexpected characters in {r.id}: {bad}")
        tmp.append((exon_num(r.id, i), r.id, seq))
    tmp.sort(key=lambda x: x[0])

    exons: List[RefExon] = []
    offset = 0
    for idx, (_, rid, seq) in enumerate(tmp, 1):
        exons.append(RefExon(
            exon=rid,
            index=idx,
            seq=seq,
            length=len(seq),
            cds_start=offset + 1,
            cds_end=offset + len(seq),
            phase_at_start=offset % 3,
        ))
        offset += len(seq)

    cds = "".join(e.seq for e in exons)
    if len(cds) % 3 != 0:
        logging.warning("Concatenated reference exons are not divisible by 3: %d nt", len(cds))
    protein = str(Seq(cds).translate())
    internal_stops = [i + 1 for i, aa in enumerate(protein[:-1]) if aa == "*"]
    if internal_stops:
        logging.warning("Reference CDS contains internal stops at codons: %s", internal_stops)
    logging.info("Reference: %d exons, %d nt, %d codons, terminal AA=%s",
                 len(exons), len(cds), len(cds)//3, protein[-1] if protein else "NA")
    return exons, cds


def ref_cds_pos(exon: RefExon, exon_pos: int) -> int:
    return exon.cds_start + exon_pos - 1


def codon_of_cds_pos(cds_pos: int) -> Tuple[int, int]:
    return (cds_pos - 1)//3 + 1, (cds_pos - 1) % 3 + 1


# -----------------------------------------------------------------------------
# Target table and genome acquisition
# -----------------------------------------------------------------------------

def load_targets(path: Path, drop_paired_duplicates: bool = False) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"Assembly Accession", "Organism Name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Target TSV missing columns: {sorted(missing)}")
    if "Species Key" not in df.columns:
        df["Species Key"] = df["Organism Name"].map(species_key)
    if drop_paired_duplicates and "Assembly Paired Assembly Accession" in df.columns:
        # Prefer GCF over GCA for paired GenBank/RefSeq representations.
        accset = set(df["Assembly Accession"])
        drop = []
        for i, row in df.iterrows():
            acc = row["Assembly Accession"]
            pair = row.get("Assembly Paired Assembly Accession", "")
            if acc.startswith("GCA_") and pair.startswith("GCF_") and pair in accset:
                drop.append(i)
        if drop:
            logging.info("Dropping %d paired GCA records because matching GCF records are present", len(drop))
            df = df.drop(index=drop).reset_index(drop=True)
    return df


def load_genome_map(path: Optional[Path]) -> Dict[str, Path]:
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if not {"Assembly Accession", "Genome FASTA"}.issubset(df.columns):
        raise ValueError("--genome-map must contain 'Assembly Accession' and 'Genome FASTA' columns")
    return {r["Assembly Accession"]: Path(r["Genome FASTA"]).expanduser().resolve()
            for _, r in df.iterrows()}


def find_genomic_fna(root: Path) -> Optional[Path]:
    candidates = list(root.rglob("*.fna")) + list(root.rglob("*.fa")) + list(root.rglob("*.fasta"))
    if not candidates:
        return None
    # Prefer NCBI genomic naming, then largest FASTA.
    genomic = [p for p in candidates if "genomic" in p.name.lower()]
    cands = genomic or candidates
    return max(cands, key=lambda p: p.stat().st_size)


def _zip_is_valid(path: Path) -> bool:
    """Return True only for a non-trivial, readable ZIP archive."""
    try:
        return path.exists() and path.stat().st_size > 1000 and zipfile.is_zipfile(path)
    except OSError:
        return False


def _ncbi_ftp_base(accession: str) -> str:
    """
    Convert GCA_041121725.1 -> 
    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/041/121/725/
    """
    m = re.fullmatch(r"(GC[AF])_(\d{9})\.(\d+)", accession)
    if not m:
        raise ValueError(f"Unrecognized assembly accession format: {accession}")
    prefix, digits, _version = m.groups()
    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/"
        f"{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/"
    )


def _url_text(url: str, timeout: int = 90) -> str:
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "GPRC6A-ruminant-pipeline/1.0.1"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def _download_url(url: str, dest: Path, timeout: int = 180) -> None:
    """Stream an HTTP/HTTPS file to disk without loading it into memory."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    if tmp.exists():
        tmp.unlink()

    req = urllib.request.Request(
        url,
        headers={"User-Agent": "GPRC6A-ruminant-pipeline/1.0.1"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as out:
        shutil.copyfileobj(r, out, length=16 * 1024 * 1024)
    tmp.replace(dest)


def download_genome_from_ftp(accession: str, adir: Path,
                             retries: int = 3, backoff: float = 20.0) -> Path:
    """
    Fallback downloader using the official NCBI genomes/all HTTPS tree.

    This is intentionally independent of the Datasets API. It is useful when
    `datasets download` fails transiently even though the assembly is valid.
    """
    base = _ncbi_ftp_base(accession)
    last_error: Optional[Exception] = None

    for attempt in range(1, retries + 1):
        try:
            listing = _url_text(base)

            # NCBI assembly directory is normally ACCESSION_ASSEMBLYNAME/
            hrefs = re.findall(r'href=["\']([^"\']+/)["\']', listing, flags=re.I)
            candidates = [
                urllib.parse.unquote(h).rstrip("/")
                for h in hrefs
                if urllib.parse.unquote(h).rstrip("/").startswith(accession + "_")
            ]
            if not candidates:
                raise RuntimeError(
                    f"NCBI FTP directory found, but no assembly directory starts with {accession}_ at {base}"
                )

            assembly_dir = sorted(candidates)[0]
            assembly_url = urllib.parse.urljoin(base, assembly_dir + "/")
            files_html = _url_text(assembly_url)

            file_hrefs = re.findall(r'href=["\']([^"\']+)["\']', files_html, flags=re.I)
            genomic = [
                urllib.parse.unquote(h)
                for h in file_hrefs
                if urllib.parse.unquote(h).endswith("_genomic.fna.gz")
                and not urllib.parse.unquote(h).endswith("_rna_from_genomic.fna.gz")
            ]
            if not genomic:
                raise RuntimeError(
                    f"No *_genomic.fna.gz file listed for {accession} at {assembly_url}"
                )

            gz_name = sorted(genomic)[0]
            gz_url = urllib.parse.urljoin(assembly_url, gz_name)
            gz_path = adir / gz_name
            fna_path = adir / gz_name[:-3]

            logging.warning(
                "Falling back to NCBI genomes FTP/HTTPS for %s: %s",
                accession, gz_url
            )

            if not gz_path.exists() or gz_path.stat().st_size == 0:
                _download_url(gz_url, gz_path)

            tmp_fna = fna_path.with_suffix(fna_path.suffix + ".part")
            if tmp_fna.exists():
                tmp_fna.unlink()
            with gzip.open(gz_path, "rb") as inp, open(tmp_fna, "wb") as out:
                shutil.copyfileobj(inp, out, length=16 * 1024 * 1024)
            tmp_fna.replace(fna_path)

            if not fna_path.exists() or fna_path.stat().st_size == 0:
                raise RuntimeError(f"FTP FASTA decompression produced an empty file for {accession}")

            # Keep only the decompressed FASTA unless the user wants to inspect the cache.
            try:
                gz_path.unlink()
            except OSError:
                pass

            return fna_path

        except Exception as e:
            last_error = e
            if attempt < retries:
                delay = backoff * attempt
                logging.warning(
                    "FTP fallback attempt %d/%d failed for %s: %s. Retrying in %.0f s",
                    attempt, retries, accession, e, delay
                )
                time.sleep(delay)

    raise RuntimeError(
        f"NCBI FTP fallback failed for {accession} after {retries} attempts: {last_error}"
    )


def download_genome(accession: str, cache_dir: Path, datasets_exe: str,
                    retries: int = 4, backoff: float = 15.0,
                    ftp_fallback: bool = True) -> Path:
    """
    Download an assembly robustly.

    Strategy:
      1. Reuse a completed cached FASTA.
      2. Retry `datasets download` several times with exponential-ish backoff.
      3. Delete partial/corrupt ZIPs between attempts.
      4. Capture and log the real NCBI stderr message.
      5. Fall back to the official NCBI genomes FTP/HTTPS hierarchy.
    """
    adir = cache_dir / safe_name(accession)
    adir.mkdir(parents=True, exist_ok=True)
    marker = adir / ".download_complete"

    existing = find_genomic_fna(adir)
    if existing and existing.exists() and existing.stat().st_size > 0 and marker.exists():
        return existing

    # If a full FASTA exists from an earlier extraction but the marker was not
    # written due to a later interruption, it is safe to reuse it.
    if existing and existing.exists() and existing.stat().st_size > 1_000_000:
        marker.write_text(str(existing) + "\n")
        return existing

    zpath = adir / f"{accession}.zip"
    last_error = ""

    for attempt in range(1, retries + 1):
        # Never trust a ZIP left behind by a failed datasets invocation.
        if zpath.exists() and not _zip_is_valid(zpath):
            try:
                zpath.unlink()
            except OSError:
                pass

        cmd = [
            datasets_exe, "download", "genome", "accession", accession,
            "--include", "genome",
            "--filename", str(zpath),
            "--no-progressbar",
        ]

        logging.info(
            "Downloading %s with NCBI datasets (attempt %d/%d)",
            accession, attempt, retries
        )

        cp = subprocess.run(
            list(map(str, cmd)),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
            text=True,
        )

        stderr_text = (cp.stderr or "").strip()
        stdout_text = (cp.stdout or "").strip()

        if cp.returncode == 0 and _zip_is_valid(zpath):
            try:
                with zipfile.ZipFile(zpath) as zf:
                    zf.extractall(adir)
                fna = find_genomic_fna(adir)
                if not fna or not fna.exists() or fna.stat().st_size == 0:
                    raise RuntimeError(
                        f"No genomic FASTA found after extracting datasets package for {accession}"
                    )
                marker.write_text(str(fna) + "\n")
                return fna
            except Exception as e:
                last_error = f"ZIP extraction/validation failed: {type(e).__name__}: {e}"
        else:
            last_error = (
                f"datasets exit={cp.returncode}; "
                f"stderr={stderr_text[-3000:] or '<empty>'}; "
                f"stdout={stdout_text[-1000:] or '<empty>'}"
            )

        logging.warning(
            "NCBI datasets download failed for %s on attempt %d/%d: %s",
            accession, attempt, retries, last_error
        )

        # A failed invocation can leave a truncated ZIP. Remove it so the next
        # attempt does not immediately fail with BadZipFile.
        try:
            if zpath.exists():
                zpath.unlink()
        except OSError:
            pass

        if attempt < retries:
            delay = backoff * (2 ** (attempt - 1))
            logging.info("Retrying %s in %.0f s", accession, delay)
            time.sleep(delay)

    if ftp_fallback:
        try:
            fna = download_genome_from_ftp(
                accession, adir,
                retries=max(2, min(3, retries)),
                backoff=max(10.0, backoff),
            )
            marker.write_text(str(fna) + "\n")
            return fna
        except Exception as ftp_error:
            raise RuntimeError(
                f"Genome download failed for {accession}. "
                f"Datasets error: {last_error}. "
                f"FTP fallback error: {ftp_error}"
            ) from ftp_error

    raise RuntimeError(
        f"Genome download failed for {accession} after {retries} datasets attempts: {last_error}"
    )


# -----------------------------------------------------------------------------
# BLAST
# -----------------------------------------------------------------------------

def run_blast(exons_fasta: Path, genome_fasta: Path, out_gz: Path,
              blastn_exe: str, threads: int, force: bool = False) -> None:
    if out_gz.exists() and out_gz.stat().st_size > 20 and not force:
        return
    out_gz.parent.mkdir(parents=True, exist_ok=True)
    tmp = out_gz.with_suffix("")
    outfmt = "6 " + " ".join(BLAST_COLUMNS)
    cmd = [
        blastn_exe,
        "-task", "blastn",
        "-query", str(exons_fasta),
        "-subject", str(genome_fasta),
        "-out", str(tmp),
        "-outfmt", outfmt,
        "-evalue", "1e-3",
        "-word_size", "9",
        "-dust", "no",
        "-soft_masking", "false",
        "-max_target_seqs", "10000",
        "-max_hsps", "1000",
        "-num_threads", str(max(1, threads)),
    ]
    run_cmd(cmd)
    with open(tmp, "rb") as fi, gzip.open(out_gz, "wb", compresslevel=6) as fo:
        shutil.copyfileobj(fi, fo)
    tmp.unlink(missing_ok=True)


def load_blast(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLUMNS)
    try:
        df = pd.read_csv(path, sep="\t", names=BLAST_COLUMNS, header=None,
                         compression="gzip", dtype={"qseqid": str, "sseqid": str,
                                                    "qseq": str, "sseq": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=BLAST_COLUMNS)
    for c in ["pident", "evalue", "bitscore"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    for c in ["length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qlen", "slen"]:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df["strand"] = df.apply(lambda r: "+" if r.sstart <= r.send else "-", axis=1)
    df["gmin"] = df[["sstart", "send"]].min(axis=1)
    df["gmax"] = df[["sstart", "send"]].max(axis=1)
    return df


def union_coverage(intervals: Iterable[Tuple[int, int]], length: int) -> int:
    ints = sorted((max(1, min(a,b)), min(length, max(a,b))) for a,b in intervals)
    if not ints:
        return 0
    total = 0
    cs, ce = ints[0]
    for s,e in ints[1:]:
        if s <= ce + 1:
            ce = max(ce, e)
        else:
            total += ce - cs + 1
            cs, ce = s,e
    total += ce - cs + 1
    return total


def order_coherence(cluster: pd.DataFrame, exon_order: Dict[str, int]) -> float:
    med = []
    for q, g in cluster.groupby("qseqid"):
        if q not in exon_order:
            continue
        med.append((exon_order[q], float(statistics.median(((g.sstart + g.send)/2).tolist()))))
    if len(med) <= 1:
        return 0.5
    med.sort()
    strand = cluster.iloc[0]["strand"]
    coords = [x[1] for x in med]
    good = 0
    tot = 0
    for a,b in zip(coords, coords[1:]):
        tot += 1
        if (strand == "+" and b > a) or (strand == "-" and b < a):
            good += 1
    return good/tot if tot else 0.5


def locus_candidates(blast: pd.DataFrame, ref_exons: List[RefExon],
                     max_locus_gap: int = 2_000_000) -> List[Dict[str, Any]]:
    if blast.empty:
        return []
    exon_len = {e.exon: e.length for e in ref_exons}
    exon_order = {e.exon: e.index for e in ref_exons}
    candidates = []
    for (contig, strand), g in blast.groupby(["sseqid", "strand"], sort=False):
        g = g.sort_values("gmin").copy()
        clusters = []
        cur_idx = []
        cur_end = None
        for idx, r in g.iterrows():
            if cur_end is None or r.gmin - cur_end <= max_locus_gap:
                cur_idx.append(idx)
                cur_end = max(cur_end or r.gmax, r.gmax)
            else:
                clusters.append(cur_idx)
                cur_idx = [idx]
                cur_end = r.gmax
        if cur_idx:
            clusters.append(cur_idx)

        for ci, idxs in enumerate(clusters, 1):
            c = g.loc[idxs]
            cov_by_exon = {}
            bits_by_exon = {}
            for q, qg in c.groupby("qseqid"):
                if q not in exon_len:
                    continue
                cov_by_exon[q] = union_coverage(zip(qg.qstart, qg.qend), exon_len[q]) / exon_len[q]
                bits_by_exon[q] = float(qg.bitscore.max())
            exon_count = sum(v >= 0.20 for v in cov_by_exon.values())
            covsum = sum(cov_by_exon.values())
            coherence = order_coherence(c, exon_order)
            bits = sum(bits_by_exon.values())
            score = 1000*exon_count + 500*coherence + 300*covsum + bits/10
            candidates.append({
                "candidate_id": f"{contig}:{strand}:{ci}",
                "contig": contig,
                "strand": strand,
                "start": int(c.gmin.min()),
                "end": int(c.gmax.max()),
                "span": int(c.gmax.max() - c.gmin.min() + 1),
                "exon_count": exon_count,
                "coverage_sum": covsum,
                "order_coherence": coherence,
                "bitscore_sum": bits,
                "score": score,
                "indices": idxs,
            })
    candidates.sort(key=lambda x: x["score"], reverse=True)
    return candidates


# -----------------------------------------------------------------------------
# FASTA random access and sequence context
# -----------------------------------------------------------------------------

def get_seq(index, contig: str, start: int, end: int, strand: str = "+") -> str:
    if contig not in index:
        return ""
    rec = index[contig]
    start = max(1, int(start))
    end = min(len(rec.seq), int(end))
    if start > end:
        return ""
    seq = str(rec.seq[start-1:end]).upper()
    return reverse_complement(seq) if strand == "-" else seq


def n_stats(seq: str) -> Dict[str, Any]:
    seq = (seq or "").upper()
    n = len(seq)
    ncount = seq.count("N")
    runs = [len(m.group(0)) for m in re.finditer(r"N+", seq)]
    return {
        "length": n,
        "n_count": ncount,
        "n_fraction": frac(ncount, n) if n else float("nan"),
        "n_runs": len(runs),
        "longest_n_run": max(runs) if runs else 0,
    }


# -----------------------------------------------------------------------------
# HSP reconstruction
# -----------------------------------------------------------------------------

def parse_hsp_alignment(row: pd.Series) -> Tuple[Dict[int, Dict[str, Any]], List[Dict[str, Any]]]:
    """Return per-query-position calls and insertions within one HSP."""
    qseq = str(row.qseq).upper()
    sseq = str(row.sseq).upper()
    qpos = int(row.qstart) - 1
    scoord = int(row.sstart)
    sstep = 1 if row.strand == "+" else -1
    calls: Dict[int, Dict[str, Any]] = {}
    insertions: List[Dict[str, Any]] = []
    ins_chars = []
    ins_coords = []
    ins_anchor = None

    def flush_ins():
        nonlocal ins_chars, ins_coords, ins_anchor
        if ins_chars:
            seq = "".join(ins_chars)
            if row.strand == "-":
                # sseq from BLAST is already oriented to query, so do not reverse complement.
                pass
            insertions.append({
                "anchor_qpos": int(ins_anchor or 0),
                "sequence": seq,
                "length": len(seq),
                "subject_start": min(ins_coords) if ins_coords else None,
                "subject_end": max(ins_coords) if ins_coords else None,
                "source": "within_hsp",
                "bitscore": float(row.bitscore),
            })
        ins_chars, ins_coords, ins_anchor = [], [], None

    for qc, sc in zip(qseq, sseq):
        if qc != "-":
            flush_ins()
            qpos += 1
            call = {
                "qbase": qc,
                "sbase": sc,
                "scoord": None,
                "bitscore": float(row.bitscore),
                "pident": float(row.pident),
                "hsp_qstart": int(row.qstart),
                "hsp_qend": int(row.qend),
            }
            if sc != "-":
                call["scoord"] = scoord
                scoord += sstep
            calls[qpos] = call
        else:
            if sc != "-":
                if ins_anchor is None:
                    ins_anchor = qpos
                ins_chars.append(sc)
                ins_coords.append(scoord)
                scoord += sstep
    flush_ins()
    return calls, insertions


def choose_primary_hsps(hsps: pd.DataFrame, exon_len: int) -> pd.DataFrame:
    if hsps.empty:
        return hsps
    chosen = []
    covered = set()
    for idx, r in hsps.sort_values(["bitscore", "length"], ascending=False).iterrows():
        qset = set(range(max(1,int(r.qstart)), min(exon_len,int(r.qend))+1))
        novel = len(qset - covered)
        if not qset:
            continue
        if novel >= max(5, int(0.15*len(qset))) or not chosen:
            chosen.append(idx)
            covered |= qset
    return hsps.loc[chosen].sort_values("qstart") if chosen else hsps.iloc[0:0]


def oriented_subject_gap(prev: pd.Series, cur: pd.Series) -> Tuple[int, Optional[Tuple[int,int]]]:
    """Gap between two HSPs that are ordered by increasing query coordinate."""
    if prev.strand != cur.strand or prev.sseqid != cur.sseqid:
        return -1, None
    if prev.strand == "+":
        a = int(prev.send) + 1
        b = int(cur.sstart) - 1
    else:
        a = int(cur.sstart) + 1
        b = int(prev.send) - 1
    if b < a:
        return b-a+1, None
    return b-a+1, (a,b)


def merge_exon_calls(exon: RefExon, hsps: pd.DataFrame, fasta_index, contig: str,
                     strand: str, large_indel_threshold: int = 20) -> Dict[str, Any]:
    # Choose highest-bitscore observation at each query position.
    best_calls: Dict[int, Dict[str, Any]] = {}
    insertion_map: Dict[Tuple[int,str], Dict[str, Any]] = {}
    all_hsp_insertions = []
    for _, row in hsps.iterrows():
        calls, ins = parse_hsp_alignment(row)
        for qpos, c in calls.items():
            if qpos < 1 or qpos > exon.length:
                continue
            if qpos not in best_calls or c["bitscore"] > best_calls[qpos]["bitscore"]:
                best_calls[qpos] = c
        for x in ins:
            key = (x["anchor_qpos"], x["sequence"])
            if key not in insertion_map or x["bitscore"] > insertion_map[key]["bitscore"]:
                insertion_map[key] = x
            all_hsp_insertions.append(x)

    primary = choose_primary_hsps(hsps, exon.length)
    split_events = []
    if len(primary) >= 2:
        rows = [r for _,r in primary.iterrows()]
        for prev, cur in zip(rows, rows[1:]):
            qgap = int(cur.qstart) - int(prev.qend) - 1
            sgap, coords = oriented_subject_gap(prev, cur)
            if sgap < 0:
                continue
            delta = sgap - max(0, qgap)
            if delta >= large_indel_threshold and coords:
                seq = get_seq(fasta_index, contig, coords[0], coords[1], strand)
                # Remove the small query-uncovered component conservatively only from length interpretation,
                # not from sequence, because its exact alignment is unknown.
                split_events.append({
                    "kind": "candidate_large_insertion",
                    "anchor_qpos": int(prev.qend),
                    "query_gap": qgap,
                    "subject_gap": sgap,
                    "delta_length": delta,
                    "sequence": seq,
                    "subject_start": coords[0],
                    "subject_end": coords[1],
                    "source": "between_hsps",
                    "confidence": "MEDIUM",
                })
            elif -delta >= large_indel_threshold:
                split_events.append({
                    "kind": "candidate_large_deletion",
                    "anchor_qpos": int(prev.qend),
                    "query_gap": qgap,
                    "subject_gap": sgap,
                    "delta_length": delta,
                    "sequence": "",
                    "subject_start": None,
                    "subject_end": None,
                    "source": "between_hsps",
                    "confidence": "MEDIUM",
                })

    # Per-position categories.
    recovered = 0
    callable_bases = 0
    aligned_positions = len(best_calls)
    matches = 0
    mismatches = 0
    n_bases = 0
    deletion_positions = []
    substitutions = []
    genomic_coords = []
    projected = [None] * exon.length
    projected_kind = ["missing"] * exon.length

    for p in range(1, exon.length+1):
        c = best_calls.get(p)
        if not c:
            continue
        s = c["sbase"].upper()
        q = exon.seq[p-1]
        if s == "-":
            deletion_positions.append(p)
            projected_kind[p-1] = "deletion"
            continue
        recovered += 1
        projected[p-1] = s
        projected_kind[p-1] = "base"
        if c.get("scoord") is not None:
            genomic_coords.append(int(c["scoord"]))
        if s == "N":
            n_bases += 1
        if s in DNA:
            callable_bases += 1
            if s == q:
                matches += 1
            else:
                mismatches += 1
                cdspos = ref_cds_pos(exon, p)
                codon, codon_pos = codon_of_cds_pos(cdspos)
                substitutions.append({
                    "type": "substitution",
                    "exon": exon.exon,
                    "exon_index": exon.index,
                    "exon_pos": p,
                    "cds_pos": cdspos,
                    "codon": codon,
                    "codon_pos": codon_pos,
                    "ref": q,
                    "alt": s,
                    "subject_coord": c.get("scoord"),
                })

    # contiguous deletion blocks from explicit BLAST gaps
    del_blocks = []
    if deletion_positions:
        s = p0 = deletion_positions[0]
        for p in deletion_positions[1:]:
            if p == p0 + 1:
                p0 = p
            else:
                del_blocks.append((s,p0))
                s = p0 = p
        del_blocks.append((s,p0))

    # Deduplicate insertions and convert to event-like records.
    insertions = []
    for x in insertion_map.values():
        insertions.append({
            "kind": "insertion",
            **x,
            "confidence": "HIGH" if x["length"] >= 1 else "LOW",
        })
    insertions.extend(split_events)

    missing_positions = [p for p in range(1, exon.length+1) if p not in best_calls]
    missing_blocks = []
    if missing_positions:
        s = p0 = missing_positions[0]
        for p in missing_positions[1:]:
            if p == p0+1:
                p0=p
            else:
                missing_blocks.append((s,p0)); s=p0=p
        missing_blocks.append((s,p0))

    coverage = recovered / exon.length
    aligned_fraction = aligned_positions / exon.length
    callable_fraction = callable_bases / exon.length
    identity = matches / (matches+mismatches) if (matches+mismatches) else 0.0
    n_aligned_fraction = n_bases / recovered if recovered else float("nan")

    return {
        "exon": exon.exon,
        "exon_index": exon.index,
        "exon_length": exon.length,
        "coverage": coverage,
        "aligned_fraction": aligned_fraction,
        "callable_fraction": callable_fraction,
        "identity": identity,
        "n_aligned_fraction": n_aligned_fraction,
        "recovered_bases": recovered,
        "callable_bases": callable_bases,
        "aligned_positions": aligned_positions,
        "matches": matches,
        "mismatches": mismatches,
        "projected": projected,
        "projected_kind": projected_kind,
        "best_calls": best_calls,
        "substitutions": substitutions,
        "deletion_blocks": del_blocks,
        "insertions": insertions,
        "missing_blocks": missing_blocks,
        "genomic_min": min(genomic_coords) if genomic_coords else None,
        "genomic_max": max(genomic_coords) if genomic_coords else None,
        "primary_hsps": primary.to_dict(orient="records"),
    }


# -----------------------------------------------------------------------------
# Confidence, splice motifs, missing exon interpretation
# -----------------------------------------------------------------------------

def exon_confidence(ex: Dict[str, Any], locus: Dict[str, Any]) -> Tuple[float, str]:
    cov = ex["coverage"]
    ident = ex["identity"]
    nfrac = ex["n_aligned_fraction"]
    if pd.isna(nfrac):
        nfrac = 1.0 if cov == 0 else 0.0
    coherence = float(locus.get("order_coherence", 0.0))
    multi = min(1.0, float(locus.get("exon_count", 0))/6.0)
    score = 0.45*cov + 0.25*ident + 0.15*coherence + 0.10*multi + 0.05*(1-min(1.0,nfrac))
    return score, confidence_grade(score)


def splice_motifs(ex: Dict[str, Any], exon: RefExon, fasta_index, contig: str, strand: str,
                  n_exons: int) -> Dict[str, Any]:
    out = {"acceptor": "", "donor": "", "acceptor_status": "NA", "donor_status": "NA"}
    if ex["genomic_min"] is None or ex["coverage"] < 0.90:
        return out
    gmin, gmax = int(ex["genomic_min"]), int(ex["genomic_max"])
    # Oriented gene sequence: acceptor is 2 nt immediately 5' of exon; donor is 2 nt immediately 3'.
    if strand == "+":
        acceptor = get_seq(fasta_index, contig, gmin-2, gmin-1, "+")
        donor = get_seq(fasta_index, contig, gmax+1, gmax+2, "+")
    else:
        acceptor = get_seq(fasta_index, contig, gmax+1, gmax+2, "-")
        donor = get_seq(fasta_index, contig, gmin-2, gmin-1, "-")
    if exon.index > 1:
        out["acceptor"] = acceptor
        if "N" in acceptor or len(acceptor) != 2:
            out["acceptor_status"] = "UNRESOLVED"
        else:
            out["acceptor_status"] = "CANONICAL" if acceptor == "AG" else "NONCANONICAL"
    if exon.index < n_exons:
        out["donor"] = donor
        if "N" in donor or len(donor) != 2:
            out["donor_status"] = "UNRESOLVED"
        else:
            # GT is major; GC is a recognized minor donor. AT is reserved for U12 introns and is flagged separately.
            out["donor_status"] = "CANONICAL" if donor in {"GT", "GC"} else ("AT_MINOR" if donor == "AT" else "NONCANONICAL")
    return out


def expected_missing_interval(exon_i: int, exon_results: Dict[int, Dict[str, Any]],
                              strand: str) -> Optional[Tuple[int,int]]:
    prev = exon_results.get(exon_i-1)
    nxt = exon_results.get(exon_i+1)
    if not prev or not nxt:
        return None
    if prev.get("genomic_min") is None or nxt.get("genomic_min") is None:
        return None
    if strand == "+":
        a = int(prev["genomic_max"]) + 1
        b = int(nxt["genomic_min"]) - 1
    else:
        a = int(nxt["genomic_max"]) + 1
        b = int(prev["genomic_min"]) - 1
    if b < a:
        return None
    return a,b


# -----------------------------------------------------------------------------
# Variant/event calling
# -----------------------------------------------------------------------------

def make_event(event_type: str, exon: RefExon, exon_pos: Optional[int],
               details: Dict[str, Any], confidence: str, disruptive: bool = True) -> Dict[str, Any]:
    cdspos = ref_cds_pos(exon, exon_pos) if exon_pos else None
    codon, codon_pos = codon_of_cds_pos(cdspos) if cdspos else (None,None)
    d = {
        "event_type": event_type,
        "exon": exon.exon,
        "exon_index": exon.index,
        "exon_pos": exon_pos,
        "cds_pos": cdspos,
        "reference_codon": codon,
        "codon_pos": codon_pos,
        "confidence": confidence,
        "gene_disrupting": disruptive,
    }
    d.update(details)
    return d


def build_target_cds(ref_exons: List[RefExon], exon_results: Dict[int, Dict[str, Any]]) -> Tuple[str, List[Optional[int]]]:
    """Build target CDS in pig exon order, retaining insertions and omitting explicit deletions.

    Missing reference positions become N so an assembly gap is not mistaken for a deletion.
    Returns sequence and mapping from target base index to pig CDS position or insertion anchor.
    """
    seq_parts = []
    mapping: List[Optional[int]] = []
    for exon in ref_exons:
        ex = exon_results[exon.index]
        ins_by_anchor = defaultdict(list)
        for ins in ex["insertions"]:
            if ins.get("kind") == "insertion" and ins.get("source") == "within_hsp":
                ins_by_anchor[int(ins["anchor_qpos"])].append(ins)
        # Insertions before first ref base (anchor 0)
        for ins in sorted(ins_by_anchor.get(0, []), key=lambda x: -x.get("bitscore",0))[:1]:
            seq_parts.append(ins["sequence"])
            mapping.extend([exon.cds_start-1]*len(ins["sequence"]))
        for p in range(1, exon.length+1):
            kind = ex["projected_kind"][p-1]
            base = ex["projected"][p-1]
            cdspos = ref_cds_pos(exon,p)
            if kind == "deletion":
                pass
            elif kind == "base" and base:
                seq_parts.append(base)
                mapping.append(cdspos)
            else:
                seq_parts.append("N")
                mapping.append(cdspos)
            for ins in sorted(ins_by_anchor.get(p, []), key=lambda x: -x.get("bitscore",0))[:1]:
                seq_parts.append(ins["sequence"])
                mapping.extend([cdspos]*len(ins["sequence"]))
    return "".join(seq_parts), mapping


def call_events(ref_exons: List[RefExon], ref_cds: str,
                exon_results: Dict[int, Dict[str, Any]],
                fasta_index, contig: str, strand: str,
                locus_n: Dict[str, Any]) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    events: List[Dict[str, Any]] = []
    all_variants: List[Dict[str, Any]] = []

    # Small substitutions and explicit indels.
    for exon in ref_exons:
        ex = exon_results[exon.index]
        all_variants.extend(ex["substitutions"])
        grade = ex["confidence"]

        for s,e in ex["deletion_blocks"]:
            length = e-s+1
            typ = "frameshift_deletion" if length % 3 else "inframe_deletion"
            ev = make_event(typ, exon, s, {
                "length": length,
                "exon_end_pos": e,
                "sequence": exon.seq[s-1:e],
                "frame_delta": -length,
                "source": "aligned_hsp_gap",
            }, grade, disruptive=(length % 3 != 0))
            events.append(ev)
            all_variants.append({**ev, "type": typ})

        for ins in ex["insertions"]:
            if ins.get("kind") not in {"insertion", "candidate_large_insertion", "candidate_large_deletion"}:
                continue
            if ins.get("kind") == "candidate_large_deletion":
                length = abs(int(ins.get("delta_length",0)))
                typ = "candidate_frameshift_deletion" if length % 3 else "candidate_inframe_deletion"
                events.append(make_event(typ, exon, max(1,int(ins.get("anchor_qpos",1))), {
                    "length": length,
                    "frame_delta": -length,
                    "source": ins.get("source"),
                    "subject_gap": ins.get("subject_gap"),
                    "query_gap": ins.get("query_gap"),
                }, "MEDIUM", disruptive=(length % 3 != 0)))
                continue
            seq = str(ins.get("sequence", ""))
            length = int(ins.get("length", ins.get("delta_length", len(seq))))
            typ = "frameshift_insertion" if length % 3 else "inframe_insertion"
            if ins.get("kind") == "candidate_large_insertion":
                typ = "candidate_large_" + typ
            events.append(make_event(typ, exon, max(1,int(ins.get("anchor_qpos",1))), {
                "length": length,
                "sequence": seq,
                "sequence_md5": md5_text(seq) if seq else "",
                "frame_delta": length,
                "source": ins.get("source"),
                "subject_start": ins.get("subject_start"),
                "subject_end": ins.get("subject_end"),
            }, ins.get("confidence", grade), disruptive=(length % 3 != 0)))

        # Missing blocks are informational unless strong flanking evidence supports deletion.
        for s,e in ex["missing_blocks"]:
            length = e-s+1
            all_variants.append({
                "type": "unrecovered_sequence",
                "exon": exon.exon,
                "exon_index": exon.index,
                "exon_pos": s,
                "exon_end_pos": e,
                "length": length,
            })

    # Probable whole-exon deletion versus unresolved exon.
    for exon in ref_exons:
        ex = exon_results[exon.index]
        if ex["coverage"] >= 0.20:
            continue
        interval = expected_missing_interval(exon.index, exon_results, strand)
        if interval:
            seq = get_seq(fasta_index, contig, interval[0], interval[1], "+")
            ns = n_stats(seq)
            flanks_good = exon_results[exon.index-1]["coverage"] >= 0.70 and exon_results[exon.index+1]["coverage"] >= 0.70
            assembled = ns["length"] > 0 and ns["n_fraction"] <= 0.05 and ns["longest_n_run"] < 100
            if flanks_good and assembled:
                events.append(make_event("probable_exon_deletion", exon, 1, {
                    "length": exon.length,
                    "interval_start": interval[0],
                    "interval_end": interval[1],
                    "interval_n_fraction": ns["n_fraction"],
                    "reason": "flanking exons recovered on same locus and intervening genome interval is assembled",
                }, "HIGH" if ex["coverage"] == 0 else "MEDIUM", disruptive=True))
            else:
                events.append(make_event("exon_unresolved", exon, 1, {
                    "length": exon.length,
                    "interval_start": interval[0],
                    "interval_end": interval[1],
                    "interval_n_fraction": ns.get("n_fraction"),
                    "reason": "low exon recovery but deletion cannot be separated confidently from assembly/alignment failure",
                }, "LOW", disruptive=False))
        else:
            events.append(make_event("exon_unresolved", exon, 1, {
                "length": exon.length,
                "reason": "low recovery and insufficient flanking-exon evidence",
            }, "LOW", disruptive=False))

    # Potential splice motif disruptions.
    for exon in ref_exons:
        ex = exon_results[exon.index]
        sp = ex.get("splice", {})
        if sp.get("acceptor_status") == "NONCANONICAL":
            events.append(make_event("potential_splice_acceptor_disruption", exon, 1, {
                "motif": sp.get("acceptor"),
                "reason": "noncanonical 2-bp acceptor motif; requires manual/annotation validation",
            }, "MEDIUM", disruptive=True))
        if sp.get("donor_status") == "NONCANONICAL":
            events.append(make_event("potential_splice_donor_disruption", exon, exon.length, {
                "motif": sp.get("donor"),
                "reason": "noncanonical 2-bp donor motif; requires manual/annotation validation",
            }, "MEDIUM", disruptive=True))

    # Codon-level consequences from projected target bases in reference frame where no indel crosses the codon.
    ref_prot = str(Seq(ref_cds).translate())
    projected_ref = {}
    indel_ref_positions = set()
    for exon in ref_exons:
        ex = exon_results[exon.index]
        for p in range(1,exon.length+1):
            cdspos = ref_cds_pos(exon,p)
            if ex["projected_kind"][p-1] == "base":
                projected_ref[cdspos] = ex["projected"][p-1]
            if ex["projected_kind"][p-1] == "deletion":
                indel_ref_positions.add(cdspos)
        for ins in ex["insertions"]:
            if ins.get("kind") in {"insertion", "candidate_large_insertion"}:
                p = int(ins.get("anchor_qpos",1))
                if 1 <= p <= exon.length:
                    indel_ref_positions.add(ref_cds_pos(exon,p))

    for codon_idx in range(1, len(ref_cds)//3 + 1):
        ps = [3*codon_idx-2, 3*codon_idx-1, 3*codon_idx]
        if any(p in indel_ref_positions for p in ps):
            continue
        if not all(p in projected_ref and projected_ref[p] in DNA for p in ps):
            continue
        alt_codon = "".join(projected_ref[p] for p in ps)
        ref_codon = ref_cds[ps[0]-1:ps[-1]]
        alt_aa = str(Seq(alt_codon).translate())
        ref_aa = str(Seq(ref_codon).translate())
        if alt_aa == "*" and ref_aa != "*":
            # locate exon containing first base of codon
            exon = next(e for e in ref_exons if e.cds_start <= ps[0] <= e.cds_end)
            exonpos = ps[0] - exon.cds_start + 1
            events.append(make_event("premature_stop", exon, exonpos, {
                "target_codon": codon_idx,
                "ref_codon_seq": ref_codon,
                "alt_codon_seq": alt_codon,
                "ref_aa": ref_aa,
                "alt_aa": alt_aa,
            }, exon_results[exon.index]["confidence"], disruptive=True))

    # Start/stop loss.
    first3 = [projected_ref.get(i) for i in (1,2,3)]
    if all(x in DNA for x in first3):
        alt = "".join(first3)
        if alt != ref_cds[:3] and str(Seq(alt).translate()) != "M":
            events.append(make_event("start_lost", ref_exons[0], 1, {
                "ref_codon_seq": ref_cds[:3], "alt_codon_seq": alt,
            }, exon_results[1]["confidence"], disruptive=True))
    lastps = [len(ref_cds)-2, len(ref_cds)-1, len(ref_cds)]
    last3 = [projected_ref.get(i) for i in lastps]
    if all(x in DNA for x in last3):
        alt = "".join(last3)
        if str(Seq(ref_cds[-3:]).translate()) == "*" and str(Seq(alt).translate()) != "*":
            events.append(make_event("stop_lost", ref_exons[-1], ref_exons[-1].length-2, {
                "ref_codon_seq": ref_cds[-3:], "alt_codon_seq": alt,
            }, exon_results[ref_exons[-1].index]["confidence"], disruptive=True))

    # Translate reconstructed target to locate downstream stops, including consequences of frameshifts.
    target_cds, mapping = build_target_cds(ref_exons, exon_results)
    translatable_len = (len(target_cds)//3)*3
    target_prot = str(Seq(target_cds[:translatable_len]).translate()) if translatable_len >= 3 else ""
    frame_events = sorted(
        [e for e in events if "frameshift" in e["event_type"] and e.get("frame_delta")],
        key=lambda e: e.get("cds_pos") or 10**12,
    )
    seen_frameshift_blocks = set()
    for i, aa in enumerate(target_prot[:-1], 1):
        if aa != "*":
            continue
        tstart = (i-1)*3
        anchors = [p for p in mapping[tstart:tstart+3] if p]
        if not anchors:
            continue
        anchor = int(statistics.median(anchors))
        upstream = [e for e in frame_events if (e.get("cds_pos") or 0) <= anchor]
        cum = sum(int(e.get("frame_delta",0)) for e in upstream) % 3
        if cum:
            # Report the first stop reached in each distinct uncompensated frameshift block.
            block_key = tuple((e.get("cds_pos"), e.get("frame_delta")) for e in upstream)
            if block_key in seen_frameshift_blocks:
                continue
            seen_frameshift_blocks.add(block_key)
            exon = next((e for e in ref_exons if e.cds_start <= anchor <= e.cds_end), ref_exons[-1])
            ep = max(1, min(exon.length, anchor - exon.cds_start + 1))
            events.append(make_event("frameshift_induced_stop", exon, ep, {
                "target_codon": i,
                "reference_anchor_cds_pos": anchor,
                "cumulative_frame_offset": cum,
                "upstream_frameshift_events": len(upstream),
            }, "MEDIUM", disruptive=True))

    return events, all_variants


# -----------------------------------------------------------------------------
# One assembly
# -----------------------------------------------------------------------------

def analyze_assembly(row: Dict[str, str], args_dict: Dict[str, Any], ref_serialized: List[Dict[str,Any]], ref_cds: str,
                     genome_map: Dict[str,str]) -> Dict[str, Any]:
    args = argparse.Namespace(**args_dict)
    ref_exons = [RefExon(**x) for x in ref_serialized]
    accession = row["Assembly Accession"]
    organism = row["Organism Name"]
    species = row.get("Species Key") or species_key(organism)
    adir = Path(args.outdir) / "assemblies" / safe_name(accession)
    adir.mkdir(parents=True, exist_ok=True)
    analysis_json = adir / "analysis.json"
    if args.resume and analysis_json.exists() and analysis_json.stat().st_size > 100:
        try:
            with open(analysis_json) as f:
                previous = json.load(f)
            previous_status = previous.get("status", "")
            if previous_status in {"ok", "no_blast_hit"}:
                return previous
            logging.info(
                "Retrying %s because existing analysis.json has status=%s",
                accession, previous_status or "<missing>"
            )
        except Exception as e:
            logging.warning(
                "Could not reuse existing analysis.json for %s (%s); re-running assembly",
                accession, e
            )

    result: Dict[str, Any] = {
        "pipeline_version": VERSION,
        "assembly_accession": accession,
        "assembly_name": row.get("Assembly Name", ""),
        "organism_name": organism,
        "species_key": species,
        "assembly_level": row.get("Assembly Level", ""),
        "assembly_release_date": row.get("Assembly Release Date", ""),
        "status": "started",
        "error": "",
        "locus": {},
        "locus_candidates": [],
        "exons": [],
        "events": [],
        "all_variants": [],
    }
    try:
        # genome acquisition
        if accession in genome_map:
            genome = Path(genome_map[accession])
            if not genome.exists():
                raise FileNotFoundError(f"Genome map path does not exist: {genome}")
        else:
            genome = download_genome(
                accession,
                Path(args.outdir) / "cache" / "genomes",
                args.datasets_exe,
                retries=args.download_retries,
                backoff=args.download_backoff,
                ftp_fallback=not args.no_ftp_fallback,
            )
        result["genome_fasta"] = str(genome)

        blast_gz = adir / "blast.tsv.gz"
        run_blast(Path(args.exons), genome, blast_gz, args.blastn_exe, args.blast_threads, args.force_blast)
        blast = load_blast(blast_gz)
        candidates = locus_candidates(blast, ref_exons, args.max_locus_gap)
        result["locus_candidates"] = [{k:v for k,v in c.items() if k != "indices"} for c in candidates[:args.keep_locus_candidates]]
        if not candidates:
            result["status"] = "no_blast_hit"
            result["events"] = [{
                "event_type": "gene_unresolved_no_hit", "exon":"", "exon_index":None,
                "confidence":"LOW", "gene_disrupting":False,
                "reason":"No GPRC6A exon BLAST hit passed the search settings. This is not by itself evidence of gene loss."
            }]
            with open(analysis_json,"w") as f: json.dump(result,f,indent=2,default=str)
            return result

        best = candidates[0]
        best_df = blast.loc[best["indices"]].copy()
        locus = {k:v for k,v in best.items() if k != "indices"}
        result["locus"] = locus

        idx = SeqIO.index(str(genome), "fasta")
        try:
            locus_seq = get_seq(idx, best["contig"], max(1,best["start"]-args.locus_flank), best["end"]+args.locus_flank, "+")
            locus_ns = n_stats(locus_seq)
            locus.update({f"locus_{k}":v for k,v in locus_ns.items()})

            exon_results: Dict[int, Dict[str, Any]] = {}
            for exon in ref_exons:
                hsps = best_df[best_df.qseqid == exon.exon].copy()
                ex = merge_exon_calls(exon, hsps, idx, best["contig"], best["strand"], args.large_indel_threshold)
                score, grade = exon_confidence(ex, best)
                ex["confidence_score"] = score
                ex["confidence"] = grade
                sp = splice_motifs(ex, exon, idx, best["contig"], best["strand"], len(ref_exons))
                ex["splice"] = sp
                # local N context
                if ex["genomic_min"] is not None:
                    local = get_seq(idx, best["contig"], ex["genomic_min"]-args.exon_n_flank,
                                    ex["genomic_max"]+args.exon_n_flank, "+")
                    ex["genomic_n_context"] = n_stats(local)
                else:
                    ex["genomic_n_context"] = n_stats("")
                exon_results[exon.index] = ex

            events, allvars = call_events(ref_exons, ref_cds, exon_results, idx, best["contig"], best["strand"], locus_ns)

            # Strip bulky internals from JSON exon records while preserving reportable details.
            ex_out = []
            for exon in ref_exons:
                ex = exon_results[exon.index]
                keep = {k:v for k,v in ex.items() if k not in {"best_calls", "projected", "projected_kind", "primary_hsps"}}
                ex_out.append(keep)
            result["exons"] = ex_out
            result["events"] = events
            result["all_variants"] = allvars
            result["status"] = "ok"
        finally:
            idx.close()

        with open(analysis_json, "w") as f:
            json.dump(result, f, indent=2, default=str)

        if args.cleanup_genomes and accession not in genome_map:
            shutil.rmtree(Path(args.outdir)/"cache"/"genomes"/safe_name(accession), ignore_errors=True)
        return result
    except Exception as e:
        if logging.getLogger().isEnabledFor(logging.DEBUG):
            logging.exception("Assembly %s failed", accession)
        else:
            logging.error("Assembly %s failed: %s: %s", accession, type(e).__name__, e)
        result["status"] = "error"
        result["error"] = f"{type(e).__name__}: {e}"
        with open(analysis_json, "w") as f:
            json.dump(result, f, indent=2, default=str)
        return result


# -----------------------------------------------------------------------------
# Repeat annotation
# -----------------------------------------------------------------------------

def collect_insertions(results: List[Dict[str,Any]], min_len: int = 6) -> pd.DataFrame:
    rows=[]
    for r in results:
        for e in r.get("events",[]):
            if "insertion" not in e.get("event_type",""):
                continue
            seq = e.get("sequence","") or ""
            length = int(e.get("length",len(seq)) or 0)
            if length < min_len or not seq:
                continue
            rid = "INS_" + md5_text("|".join([r["assembly_accession"], str(e.get("exon_index")), str(e.get("exon_pos")), seq]))[:16]
            rows.append({
                "repeat_id": rid,
                "assembly_accession": r["assembly_accession"],
                "organism_name": r["organism_name"],
                "species_key": r["species_key"],
                "exon": e.get("exon"),
                "exon_index": e.get("exon_index"),
                "exon_pos": e.get("exon_pos"),
                "cds_pos": e.get("cds_pos"),
                "length": length,
                "sequence": seq,
                "event_type": e.get("event_type"),
                "confidence": e.get("confidence"),
            })
    return pd.DataFrame(rows)


def simple_repeat_features(seq: str) -> Dict[str,Any]:
    s = seq.upper()
    if not s:
        return {"gc_fraction":float("nan"),"entropy":float("nan"),"best_tandem_period":None,"best_tandem_match":0.0}
    gc=(s.count("G")+s.count("C"))/len(s)
    counts=Counter(c for c in s if c in DNA)
    total=sum(counts.values())
    entropy=0.0
    if total:
        for n in counts.values():
            p=n/total; entropy -= p*math.log2(p)
    best_p=None; best_match=0.0
    for p in range(1,min(21,max(2,len(s)//2+1))):
        comps=sum(1 for i in range(p,len(s)) if s[i] == s[i-p])
        denom=max(1,len(s)-p)
        m=comps/denom
        if m>best_match:
            best_match=m; best_p=p
    return {"gc_fraction":gc,"entropy":entropy,"best_tandem_period":best_p,"best_tandem_match":best_match}


def run_repeatmasker(insertions: pd.DataFrame, outdir: Path, threads: int,
                     repeatmasker_exe: Optional[str]) -> pd.DataFrame:
    if insertions.empty:
        return insertions.copy()
    d=outdir/"repeats"; d.mkdir(parents=True,exist_ok=True)
    fa=d/"insertions.fa"
    with open(fa,"w") as f:
        for _,r in insertions.iterrows():
            f.write(f">{r.repeat_id}\n{r.sequence}\n")

    ann = {rid:{"repeat_name":"","repeat_class_family":"","repeat_divergence":float("nan"),"repeat_coverage":0.0,"repeat_method":"simple"}
           for rid in insertions.repeat_id}

    if repeatmasker_exe:
        cmd=[repeatmasker_exe,"-species","mammal","-pa",str(max(1,threads)),"-dir",str(d),str(fa)]
        try:
            run_cmd(cmd)
            out=Path(str(fa)+".out")
            if not out.exists():
                out=d/(fa.name+".out")
            if out.exists():
                fragments=defaultdict(list)
                with open(out,errors="replace") as f:
                    for line in f:
                        if not line.strip() or line.lstrip().startswith(("SW","score")):
                            continue
                        parts=line.split()
                        if len(parts) < 11 or not re.match(r"^\d",parts[0]):
                            continue
                        try:
                            div=float(parts[1]); qname=parts[4]; qstart=int(parts[5]); qend=int(parts[6])
                            if parts[8] == "C":
                                rname=parts[9]; rclass=parts[10]
                            else:
                                rname=parts[8]; rclass=parts[9]
                        except Exception:
                            continue
                        fragments[qname].append((qstart,qend,rname,rclass,div))
                lens=dict(zip(insertions.repeat_id,insertions.length))
                for rid,fr in fragments.items():
                    if rid not in ann: continue
                    cov=union_coverage([(a,b) for a,b,*_ in fr],lens[rid])/lens[rid]
                    best=max(fr,key=lambda x:x[1]-x[0]+1)
                    ann[rid]={"repeat_name":best[2],"repeat_class_family":best[3],
                              "repeat_divergence":mean_or_nan(x[4] for x in fr),
                              "repeat_coverage":cov,"repeat_method":"RepeatMasker"}
        except Exception as e:
            logging.warning("RepeatMasker failed; retaining sequence-composition fallback: %s", e)

    rows=[]
    for _,r in insertions.iterrows():
        dct=r.to_dict()
        dct.update(simple_repeat_features(r.sequence))
        dct.update(ann[r.repeat_id])
        if not dct["repeat_name"]:
            # conservative fallback labels only sequence architecture, not TE family
            if dct["best_tandem_match"] >= 0.80 and dct["best_tandem_period"] and dct["best_tandem_period"] <= 20:
                dct["repeat_name"]="tandem-repeat-like"
                dct["repeat_class_family"]="simple/tandem"
            elif dct["entropy"] < 1.2:
                dct["repeat_name"]="low-complexity-like"
                dct["repeat_class_family"]="low_complexity"
        rows.append(dct)
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Aggregation
# -----------------------------------------------------------------------------

def flatten_results(results: List[Dict[str,Any]]) -> Tuple[pd.DataFrame,pd.DataFrame,pd.DataFrame,pd.DataFrame]:
    locus_rows=[]; exon_rows=[]; event_rows=[]; variant_rows=[]
    for r in results:
        base={
            "assembly_accession":r.get("assembly_accession"),
            "assembly_name":r.get("assembly_name"),
            "organism_name":r.get("organism_name"),
            "species_key":r.get("species_key"),
            "assembly_level":r.get("assembly_level"),
            "status":r.get("status"),
        }
        lr={**base, **{k:v for k,v in r.get("locus",{}).items() if not isinstance(v,(dict,list))}, "error":r.get("error","")}
        locus_rows.append(lr)
        for ex in r.get("exons",[]):
            row={**base}
            for k,v in ex.items():
                if k in {"substitutions","deletion_blocks","insertions","missing_blocks"}: continue
                if isinstance(v,dict):
                    for k2,v2 in v.items(): row[f"{k}_{k2}"]=v2
                elif not isinstance(v,list): row[k]=v
            exon_rows.append(row)
        for i,e in enumerate(r.get("events",[]),1):
            event_rows.append({**base,"assembly_event_id":f"{r.get('assembly_accession')}:E{i}",**e})
        for i,v in enumerate(r.get("all_variants",[]),1):
            variant_rows.append({**base,"assembly_variant_id":f"{r.get('assembly_accession')}:V{i}",**v})
    return pd.DataFrame(locus_rows),pd.DataFrame(exon_rows),pd.DataFrame(event_rows),pd.DataFrame(variant_rows)


def species_exon_summary(exon_df: pd.DataFrame) -> pd.DataFrame:
    if exon_df.empty: return pd.DataFrame()
    rows=[]
    for (sp,ei,exon),g in exon_df.groupby(["species_key","exon_index","exon"],dropna=False):
        good=g[g.status=="ok"]
        rows.append({
            "species_key":sp,"exon_index":ei,"exon":exon,
            "assemblies_total":g.assembly_accession.nunique(),
            "assemblies_ok":good.assembly_accession.nunique(),
            "coverage_median":median_or_nan(good.coverage),
            "coverage_max":max(good.coverage) if len(good) else float("nan"),
            "callable_fraction_median":median_or_nan(good.callable_fraction),
            "identity_median":median_or_nan(good.identity),
            "n_aligned_fraction_median":median_or_nan(good.n_aligned_fraction),
            "confidence_score_median":median_or_nan(good.confidence_score),
            "high_conf_assemblies":int((good.confidence=="HIGH").sum()) if len(good) else 0,
            "medium_or_high_assemblies":int(good.confidence.isin(["HIGH","MEDIUM"]).sum()) if len(good) else 0,
        })
    return pd.DataFrame(rows)


def event_signature(row: pd.Series) -> str:
    et=str(row.get("event_type", ""))
    exon=int(row.get("exon_index") or 0) if not pd.isna(row.get("exon_index")) else 0
    pos=row.get("exon_pos")
    pos=int(pos) if pos is not None and not pd.isna(pos) else 0
    if et in {"premature_stop","frameshift_induced_stop"}:
        cod=row.get("reference_codon")
        cod=int(cod) if cod is not None and not pd.isna(cod) else int(row.get("target_codon") or 0)
        return f"{et}|ex{exon}|codon{cod}"
    if "insertion" in et or "deletion" in et:
        ln=row.get("length"); ln=int(ln) if ln is not None and not pd.isna(ln) else 0
        # Pig exon coordinate anchors homologous indels across assemblies.
        return f"{et}|ex{exon}|pos{pos}|len{ln}"
    if et in {"probable_exon_deletion","exon_unresolved"}:
        return f"{et}|ex{exon}"
    if et.startswith("potential_splice_"):
        return f"{et}|ex{exon}|{row.get('motif','')}"
    if et in {"start_lost","stop_lost"}:
        return et
    return f"{et}|ex{exon}|pos{pos}"


def aggregate_species_events(event_df: pd.DataFrame, exon_df: pd.DataFrame) -> pd.DataFrame:
    if event_df.empty: return pd.DataFrame()
    e=event_df.copy()
    e["event_key"]=e.apply(event_signature,axis=1)
    rows=[]
    for (sp,key),g in e.groupby(["species_key","event_key"],dropna=False):
        # determine callable assemblies from exon confidence where possible
        exi=g.exon_index.dropna().iloc[0] if g.exon_index.notna().any() else None
        if exi is not None and not exon_df.empty:
            exi_int = int(float(exi))
            eg=exon_df[(exon_df.species_key==sp)&(pd.to_numeric(exon_df["exon_index"], errors="coerce")==exi_int)]
            callable_acc=set(eg[eg.confidence.isin(["HIGH","MEDIUM"])].assembly_accession)
        else:
            callable_acc=set(exon_df[exon_df.species_key==sp].assembly_accession) if not exon_df.empty else set()
        support_acc=set(g[g.confidence.isin(["HIGH","MEDIUM"])].assembly_accession)
        support_all=set(g.assembly_accession)
        callable_n=len(callable_acc)
        support_n=len(support_acc)
        support_fraction=support_n/callable_n if callable_n else float("nan")
        representative=g.sort_values("confidence",key=lambda s:s.map({"HIGH":0,"MEDIUM":1,"LOW":2,"UNRESOLVED":3}).fillna(9)).iloc[0]
        if support_n >= 2 and callable_n and support_fraction >= 0.6:
            conf="HIGH"
        elif support_n >= 1 and (not callable_n or support_fraction >= 0.5):
            conf="MEDIUM"
        else:
            conf="LOW"
        rows.append({
            "species_key":sp,"event_key":key,"event_type":representative.event_type,
            "exon":representative.get("exon"),"exon_index":representative.get("exon_index"),
            "exon_pos":representative.get("exon_pos"),"cds_pos":representative.get("cds_pos"),
            "reference_codon":representative.get("reference_codon"),
            "length":representative.get("length"),"motif":representative.get("motif"),
            "gene_disrupting":bool(representative.get("gene_disrupting",False)),
            "supporting_assemblies":len(support_all),"medium_high_supporting_assemblies":support_n,
            "callable_assemblies":callable_n,"support_fraction":support_fraction,
            "species_event_confidence":conf,
            "support_accessions":";".join(sorted(support_acc)),
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Phylogenetic inference
# -----------------------------------------------------------------------------

def clade_label(clade) -> str:
    if clade.name:
        return str(clade.name)
    tips=[t.name for t in clade.get_terminals()]
    if len(tips)==1:return tips[0]
    return f"MRCA[{len(tips)}]: {tips[0]} .. {tips[-1]}"


def assign_parents(tree) -> Tuple[Dict[int,Any],Dict[int,int]]:
    parents={id(tree.root):None}; depth={id(tree.root):0}
    for cl in tree.find_clades(order="preorder"):
        for child in cl.clades:
            parents[id(child)]=cl; depth[id(child)]=depth[id(cl)]+1
    return parents,depth


def fitch_event(tree, states: Dict[str,Optional[int]]) -> Dict[str,Any]:
    sets={}
    score=0
    for cl in tree.find_clades(order="postorder"):
        if cl.is_terminal():
            st=states.get(cl.name,None)
            sets[id(cl)]={0,1} if st is None else {int(st)}
        else:
            childsets=[sets[id(c)] for c in cl.clades]
            inter=set.intersection(*childsets) if childsets else {0,1}
            if inter:
                sets[id(cl)]=inter
            else:
                sets[id(cl)]=set.union(*childsets)
                score += 1
    # Root intact when equally parsimonious.
    assigned={id(tree.root): (0 if 0 in sets[id(tree.root)] else 1)}
    origins=[]; reversions=[]
    parents,depth=assign_parents(tree)
    for cl in tree.find_clades(order="preorder"):
        if cl is tree.root: continue
        p=parents[id(cl)]
        pst=assigned[id(p)]
        sset=sets[id(cl)]
        st=pst if pst in sset else (0 if 0 in sset else 1)
        assigned[id(cl)]=st
        if pst==0 and st==1:
            origins.append({"node":clade_label(cl),"depth":depth[id(cl)],"n_descendants":len(cl.get_terminals()),
                            "descendants":[t.name for t in cl.get_terminals()]})
        elif pst==1 and st==0:
            reversions.append({"node":clade_label(cl),"depth":depth[id(cl)],"n_descendants":len(cl.get_terminals()),
                               "descendants":[t.name for t in cl.get_terminals()]})
    return {"parsimony_score":score,"origins":origins,"reversions":reversions}


def phylogenetic_summary(species_events: pd.DataFrame, species_exons: pd.DataFrame,
                         tree_path: Path) -> pd.DataFrame:
    if species_events.empty or not tree_path.exists(): return pd.DataFrame()
    tree=Phylo.read(str(tree_path),"newick")
    tips={t.name for t in tree.get_terminals()}
    rows=[]
    for key,g in species_events.groupby("event_key"):
        rep=g.iloc[0]
        exi=rep.get("exon_index")
        states={sp:None for sp in tips}
        positives=set(g[g.species_event_confidence.isin(["HIGH","MEDIUM"])].species_key) & tips
        for sp in positives: states[sp]=1
        # High/medium exon recovery provides evidence of absence if the event is not called.
        if exi is not None and not pd.isna(exi) and not species_exons.empty:
            exi_int = int(float(exi))
            sg=species_exons[pd.to_numeric(species_exons["exon_index"], errors="coerce")==exi_int]
            for _,r in sg.iterrows():
                sp=r.species_key
                if sp not in tips or sp in positives: continue
                if (r.medium_or_high_assemblies or 0) >= 1 and (r.coverage_max or 0) >= 0.75:
                    states[sp]=0
        else:
            for sp in tips-positives: states[sp]=0
        inf=fitch_event(tree,states)
        # MRCA summary for all positive tips
        mrca=""; mrca_fraction=float("nan"); outside_pos=0
        if positives:
            m=tree.common_ancestor(list(positives)) if len(positives)>1 else next(t for t in tree.get_terminals() if t.name in positives)
            mrca=clade_label(m)
            desc={t.name for t in m.get_terminals()}
            callable_desc=[sp for sp in desc if states.get(sp) is not None]
            mrca_fraction=sum(states.get(sp)==1 for sp in callable_desc)/len(callable_desc) if callable_desc else float("nan")
            outside_pos=sum(states.get(sp)==1 for sp in tips-desc)
        rows.append({
            "event_key":key,"event_type":rep.event_type,"exon":rep.get("exon"),"exon_index":exi,
            "positive_species":len(positives),"positive_species_list":";".join(sorted(positives)),
            "mrca":mrca,"mrca_callable_support_fraction":mrca_fraction,"positive_species_outside_mrca":outside_pos,
            "parsimony_score":inf["parsimony_score"],
            "inferred_origins":" | ".join(x["node"] for x in inf["origins"]),
            "origin_details":json.dumps(inf["origins"],sort_keys=True),
            "inferred_reversions":" | ".join(x["node"] for x in inf["reversions"]),
            "chronology_note":"Relative phylogenetic chronology only; no absolute age inferred without a dated tree.",
        })
    return pd.DataFrame(rows).sort_values(["positive_species","event_key"],ascending=[False,True])


def annotate_repeat_event_keys(repeat_df: pd.DataFrame) -> pd.DataFrame:
    if repeat_df.empty: return repeat_df
    r=repeat_df.copy()
    def key(x):
        fam=x.get("repeat_class_family") or x.get("repeat_name") or "unclassified"
        # ±10 nt anchor tolerance approximated by 20-nt bins for split-HSP placement jitter.
        pos=int(x.get("exon_pos") or 0)
        bucket=int(round(pos/20)*20)
        return f"repeat_insertion|ex{int(x.exon_index)}|pos~{bucket}|{fam}"
    r["repeat_event_key"]=r.apply(key,axis=1)
    return r


def aggregate_species_repeats(repeat_df: pd.DataFrame, species_exons: pd.DataFrame) -> pd.DataFrame:
    if repeat_df.empty:
        return pd.DataFrame()
    rows=[]
    for (sp,key),g in repeat_df.groupby(["species_key","repeat_event_key"]):
        exi=int(float(g.exon_index.iloc[0]))
        eg=species_exons[(species_exons.species_key==sp)&(pd.to_numeric(species_exons["exon_index"],errors="coerce")==exi)] if not species_exons.empty else pd.DataFrame()
        callable_n=int(eg.medium_or_high_assemblies.iloc[0]) if len(eg) and not pd.isna(eg.medium_or_high_assemblies.iloc[0]) else 0
        support_acc=set(g.assembly_accession)
        support_n=len(support_acc)
        sf=support_n/callable_n if callable_n else float("nan")
        if support_n>=2 and callable_n and sf>=0.6:
            conf="HIGH"
        elif support_n>=1 and (not callable_n or sf>=0.5):
            conf="MEDIUM"
        else:
            conf="LOW"
        fam=g.repeat_class_family.replace("",pd.NA).dropna() if "repeat_class_family" in g.columns else pd.Series(dtype=object)
        rn=g.repeat_name.replace("",pd.NA).dropna() if "repeat_name" in g.columns else pd.Series(dtype=object)
        rows.append({
            "species_key":sp,"repeat_event_key":key,"exon_index":exi,
            "exon_pos_median":median_or_nan(g.exon_pos),"length_median":median_or_nan(g.length),
            "repeat_class_family":fam.mode().iloc[0] if len(fam) else "unclassified",
            "repeat_name":rn.mode().iloc[0] if len(rn) else "unclassified",
            "supporting_assemblies":support_n,"callable_assemblies":callable_n,
            "support_fraction":sf,"species_repeat_confidence":conf,
            "support_accessions":";".join(sorted(support_acc)),
        })
    return pd.DataFrame(rows)


def repeat_phylogeny(species_repeat_df: pd.DataFrame, species_exons: pd.DataFrame, tree_path: Path) -> pd.DataFrame:
    if species_repeat_df.empty or not tree_path.exists(): return pd.DataFrame()
    tree=Phylo.read(str(tree_path),"newick"); tips={t.name for t in tree.get_terminals()}
    rows=[]
    for key,g in species_repeat_df.groupby("repeat_event_key"):
        positives=set(g[g.species_repeat_confidence.isin(["HIGH","MEDIUM"])].species_key)&tips
        exi=int(float(g.exon_index.iloc[0]))
        states={sp:None for sp in tips}
        for sp in positives: states[sp]=1
        eg=species_exons[pd.to_numeric(species_exons["exon_index"], errors="coerce")==int(float(exi))] if not species_exons.empty else pd.DataFrame()
        for _,r in eg.iterrows():
            sp=r.species_key
            if sp in tips and sp not in positives and (r.medium_or_high_assemblies or 0)>=1 and (r.coverage_max or 0)>=0.75:
                states[sp]=0
        inf=fitch_event(tree,states)
        rows.append({
            "repeat_event_key":key,"exon_index":exi,
            "repeat_class_family":g.repeat_class_family.replace("",pd.NA).dropna().mode().iloc[0] if g.repeat_class_family.replace("",pd.NA).dropna().size else "unclassified",
            "repeat_name":g.repeat_name.replace("",pd.NA).dropna().mode().iloc[0] if g.repeat_name.replace("",pd.NA).dropna().size else "unclassified",
            "positive_species":len(positives),"positive_species_list":";".join(sorted(positives)),
            "parsimony_score":inf["parsimony_score"],
            "inferred_origins":" | ".join(x["node"] for x in inf["origins"]),
            "origin_details":json.dumps(inf["origins"],sort_keys=True),
            "chronology_note":"Relative phylogenetic chronology inferred from presence/absence; RepeatMasker divergence is reported separately when available.",
        })
    return pd.DataFrame(rows).sort_values(["positive_species","repeat_event_key"],ascending=[False,True])


# -----------------------------------------------------------------------------
# HTML report and SVG tree
# -----------------------------------------------------------------------------

def status_counts(species_events: pd.DataFrame) -> Dict[str,int]:
    if species_events.empty: return {}
    d=defaultdict(int)
    for sp,g in species_events.groupby("species_key"):
        d[sp]=int(((g.gene_disrupting==True)&g.species_event_confidence.isin(["HIGH","MEDIUM"])).sum())
    return dict(d)


def tree_svg(tree_path: Path, species_events: pd.DataFrame, species_exons: pd.DataFrame) -> str:
    if not tree_path.exists(): return "<p>Tree unavailable.</p>"
    tree=Phylo.read(str(tree_path),"newick")
    tips=tree.get_terminals(); y={id(t):30+i*22 for i,t in enumerate(tips)}
    depth={id(tree.root):0}
    for cl in tree.find_clades(order="preorder"):
        for c in cl.clades: depth[id(c)]=depth[id(cl)]+1
    def calc_y(cl):
        if cl.is_terminal(): return y[id(cl)]
        vals=[calc_y(c) for c in cl.clades]; y[id(cl)]=sum(vals)/len(vals); return y[id(cl)]
    calc_y(tree.root)
    maxd=max(depth.values()) if depth else 1
    xscale=150; x0=20
    counts=status_counts(species_events)
    cov={}
    if not species_exons.empty:
        for sp,g in species_exons.groupby("species_key"):
            cov[sp]=mean_or_nan(g.coverage_max)
    lines=[]; texts=[]
    for cl in tree.find_clades(order="preorder"):
        x=x0+depth[id(cl)]*xscale; yy=y[id(cl)]
        if cl.clades:
            cys=[y[id(c)] for c in cl.clades]
            lines.append(f'<line x1="{x}" y1="{min(cys):.1f}" x2="{x}" y2="{max(cys):.1f}" class="branch"/>')
            for c in cl.clades:
                cx=x0+depth[id(c)]*xscale; cy=y[id(c)]
                lines.append(f'<line x1="{x}" y1="{cy:.1f}" x2="{cx}" y2="{cy:.1f}" class="branch"/>')
            if cl.name and cl.name not in {"Ruminantia"}:
                texts.append(f'<text x="{x+3}" y="{yy-3:.1f}" class="internal">{html.escape(str(cl.name))}</text>')
        else:
            sp=cl.name; n=counts.get(sp,0); c=cov.get(sp,float("nan"))
            cls="disrupted" if n else ("uncertain" if pd.isna(c) or c<0.75 else "intact")
            label=html.escape(sp)
            texts.append(f'<a href="#species-{anchor_id(sp)}"><circle cx="{x+4}" cy="{yy:.1f}" r="5" class="{cls}"><title>{label}; disruptions={n}; mean max exon recovery={c:.3f}</title></circle><text x="{x+14}" y="{yy+4:.1f}" class="tip">{label}</text></a>')
    width=x0+(maxd+1)*xscale+340; height=50+len(tips)*22
    return f'''<div class="tree-scroll"><svg class="tree" width="{width}" height="{height}" viewBox="0 0 {width} {height}">
    {''.join(lines)}{''.join(texts)}</svg></div>'''


def fmt_pct(x) -> str:
    try:
        if pd.isna(x): return "NA"
        return f"{100*float(x):.1f}%"
    except Exception:return "NA"


def html_table(df: pd.DataFrame, cols: List[str], max_rows: Optional[int]=None, id_:str="") -> str:
    if df is None or df.empty: return "<p>No records.</p>"
    d=df[cols].copy() if set(cols).issubset(df.columns) else df.copy()
    if max_rows: d=d.head(max_rows)
    out=[f'<table id="{id_}"><thead><tr>']
    out += [f'<th>{html.escape(str(c))}</th>' for c in d.columns]
    out.append('</tr></thead><tbody>')
    for _,r in d.iterrows():
        out.append('<tr>')
        for c,v in r.items():
            if isinstance(v,float) and pd.isna(v): s=""
            elif isinstance(v,float): s=f"{v:.4g}"
            else:s=str(v)
            out.append(f'<td>{html.escape(s)}</td>')
        out.append('</tr>')
    out.append('</tbody></table>')
    return ''.join(out)


def exon_boxes(sp: str, g: pd.DataFrame, sev: pd.DataFrame) -> str:
    boxes=[]
    for i in range(1,7):
        x = g[pd.to_numeric(g["exon_index"], errors="coerce") == i] if (g is not None and not g.empty and "exon_index" in g.columns) else pd.DataFrame()
        if len(x):
            r=x.iloc[0]; cov=float(r.coverage_max) if not pd.isna(r.coverage_max) else 0.0
            conf=int(r.medium_or_high_assemblies or 0)
        else: cov=0.0; conf=0
        events = sev[pd.to_numeric(sev["exon_index"], errors="coerce") == i] if (sev is not None and not sev.empty and "exon_index" in sev.columns) else pd.DataFrame()
        ndis=int(((events.gene_disrupting==True)&events.species_event_confidence.isin(["HIGH","MEDIUM"])).sum()) if len(events) else 0
        cls="exon-bad" if ndis else ("exon-low" if cov<0.75 else "exon-good")
        boxes.append(f'<div class="exonbox {cls}"><b>Exon {i}</b><span>{fmt_pct(cov)} recovered</span><span>{conf} callable assemblies</span><span>{ndis} disruption(s)</span></div>')
    return '<div class="exonrow">'+''.join(boxes)+'</div>'


def write_html_report(outdir: Path, targets: pd.DataFrame, locus_df: pd.DataFrame, exon_df: pd.DataFrame,
                      event_df: pd.DataFrame, species_exons: pd.DataFrame, species_events: pd.DataFrame,
                      phylo_df: pd.DataFrame, repeats: pd.DataFrame, repeat_phylo_df: pd.DataFrame,
                      tree_path: Path, ref_exons: List[RefExon], ref_cds: str) -> Path:
    report=outdir/"report.html"
    species=sorted(set(targets["Species Key"]))
    failed=int((locus_df.status!="ok").sum()) if not locus_df.empty else len(targets)
    disruptive=species_events[(species_events.gene_disrupting==True)&species_events.species_event_confidence.isin(["HIGH","MEDIUM"])] if not species_events.empty else pd.DataFrame()
    shared=phylo_df[phylo_df.positive_species>=2] if not phylo_df.empty else pd.DataFrame()
    tree=tree_svg(tree_path,species_events,species_exons)
    generated=time.strftime("%Y-%m-%d %H:%M:%S")

    css='''
    body{font-family:Arial,Helvetica,sans-serif;margin:0;color:#18202a;background:#f5f7fa}header{background:#162332;color:white;padding:24px 34px;position:sticky;top:0;z-index:5}header h1{margin:0 0 5px 0;font-size:24px}main{max-width:1500px;margin:auto;padding:25px}.card{background:white;border:1px solid #dce2e8;border-radius:10px;padding:18px;margin:14px 0;box-shadow:0 1px 3px #0001}.metrics{display:flex;gap:12px;flex-wrap:wrap}.metric{background:#eef3f7;padding:12px 16px;border-radius:8px;min-width:150px}.metric b{display:block;font-size:22px}.toc a{margin-right:16px}table{border-collapse:collapse;width:100%;font-size:12px;margin:10px 0}th,td{border:1px solid #dfe4ea;padding:6px;vertical-align:top}th{background:#edf2f6;position:sticky;top:78px}.warn{background:#fff4d6;border-left:4px solid #d89400;padding:10px}.good{background:#e8f6ed}.bad{background:#fdeaea}.tree-scroll{overflow:auto;max-height:900px;border:1px solid #ddd;background:white}.branch{stroke:#657180;stroke-width:1}.tip{font-size:11px;fill:#111}.internal{font-size:9px;fill:#52606d}.intact{fill:#2f855a}.disrupted{fill:#c53030}.uncertain{fill:#b7791f}.exonrow{display:flex;gap:8px;flex-wrap:wrap}.exonbox{width:145px;padding:10px;border-radius:7px;border:1px solid #ccd3db}.exonbox span{display:block;font-size:11px;margin-top:3px}.exon-good{background:#e9f7ef}.exon-low{background:#fff5d9}.exon-bad{background:#fdecec}.species h3{margin-bottom:6px}details{margin:8px 0}code{background:#eef1f4;padding:2px 4px}.small{font-size:11px;color:#52606d}.pill{display:inline-block;padding:2px 7px;border-radius:12px;background:#e8edf2;margin:1px 3px;font-size:11px}a{color:#245b8a;text-decoration:none}a:hover{text-decoration:underline}
    '''
    parts=[f'''<!doctype html><html><head><meta charset="utf-8"><title>GPRC6A Ruminantia comparative genomics</title><style>{css}</style></head><body>
    <header><h1>GPRC6A comparative genomics across Ruminantia</h1><div>Reference: pig coding exons; generated {html.escape(generated)}</div></header><main>
    <div class="card toc"><b>Navigation:</b> <a href="#overview">Overview</a><a href="#tree">Phylogeny</a><a href="#shared">Shared disruptions</a><a href="#repeats">Repeat insertions</a><a href="#species">Per-species reports</a><a href="#methods">Methods / confidence</a></div>
    <section id="overview" class="card"><h2>Overview</h2>
    <div class="metrics"><div class="metric"><b>{len(targets)}</b>assembly targets</div><div class="metric"><b>{len(species)}</b>species/labels</div><div class="metric"><b>{len(ref_exons)}</b>pig exons</div><div class="metric"><b>{len(ref_cds)}</b>reference CDS nt</div><div class="metric"><b>{failed}</b>assemblies not fully analyzed</div><div class="metric"><b>{len(disruptive)}</b>species-event calls</div></div>
    <p class="warn"><b>Interpretation:</b> failure to recover an exon is not automatically a deletion. The pipeline calls a probable exon deletion only when flanking exons support the same locus and the intervening genomic interval is assembled with low N-content. Otherwise it reports the exon as unresolved.</p>
    <p>All assembly-level data are retained. Species-level calls summarize agreement across independent assemblies. Hybrids may be analyzed but are excluded from ordinary bifurcating phylogenetic inference unless represented in the supplied tree.</p>
    <p><b>Machine-readable data:</b> <a href="tables/assembly_locus_summary.tsv">loci</a> · <a href="tables/assembly_exon_metrics.tsv">exon metrics</a> · <a href="tables/assembly_events.tsv">assembly events</a> · <a href="tables/assembly_all_variants.tsv">all variants</a> · <a href="tables/species_events.tsv">species events</a> · <a href="tables/phylogenetic_event_summary.tsv">phylogenetic events</a> · <a href="tables/repeat_insertions.tsv">repeat insertions</a> · <a href="tables/species_repeat_summary.tsv">species repeat consensus</a></p></section>''']

    parts.append(f'<section id="tree" class="card"><h2>Phylogenetic overview</h2><p>Tip circles link to species reports. Red indicates at least one medium/high-confidence disruptive event; green indicates no such event with generally recoverable exons; amber indicates insufficient recovery.</p>{tree}</section>')

    parts.append('<section id="shared" class="card"><h2>Shared gene-disrupting events and inferred chronology</h2>')
    if not shared.empty:
        show=shared[shared.event_key.isin(set(disruptive.event_key))] if not disruptive.empty else shared
        parts.append(html_table(show,["event_key","event_type","exon","positive_species","mrca","mrca_callable_support_fraction","parsimony_score","inferred_origins","inferred_reversions","chronology_note"]))
    else: parts.append('<p>No event is currently supported in two or more species.</p>')
    parts.append('</section>')

    parts.append('<section id="repeats" class="card"><h2>Repeat insertions</h2>')
    if not repeat_phylo_df.empty:
        parts.append('<h3>Phylogenetic repeat-insertion groups</h3>')
        parts.append(html_table(repeat_phylo_df,["repeat_event_key","exon_index","repeat_class_family","repeat_name","positive_species","positive_species_list","parsimony_score","inferred_origins","chronology_note"]))
    if not repeats.empty:
        parts.append('<details><summary>All insertion/repeat annotations</summary>')
        parts.append(html_table(repeats,[c for c in ["repeat_id","assembly_accession","species_key","exon_index","exon_pos","length","repeat_name","repeat_class_family","repeat_divergence","repeat_coverage","repeat_method","best_tandem_period","best_tandem_match","entropy"] if c in repeats.columns]))
        parts.append('</details>')
    else: parts.append('<p>No insertion sequences at or above the configured minimum repeat length were recovered.</p>')
    parts.append('</section>')

    parts.append('<section id="species"><h2>Per-species reports</h2>')
    for sp in species:
        sg=species_exons[species_exons.species_key==sp] if not species_exons.empty else pd.DataFrame()
        se=species_events[species_events.species_key==sp] if not species_events.empty else pd.DataFrame()
        assemblies=locus_df[locus_df.species_key==sp] if not locus_df.empty else pd.DataFrame()
        raw_events=event_df[event_df.species_key==sp] if not event_df.empty else pd.DataFrame()
        parts.append(f'<article id="species-{anchor_id(sp)}" class="card species"><h3>{html.escape(sp)}</h3>')
        parts.append(exon_boxes(sp,sg,se))
        parts.append(f'<p><span class="pill">{len(assemblies)} assemblies</span><span class="pill">{len(se)} species-level events</span></p>')
        if not assemblies.empty and "assembly_accession" in assemblies.columns:
            blast_links = " · ".join(f'<a href="assemblies/{safe_name(str(a))}/blast.tsv.gz">{html.escape(str(a))}</a>' for a in assemblies.assembly_accession if str(a))
            if blast_links:
                parts.append(f'<p class="small"><b>Raw BLAST:</b> {blast_links}</p>')
        if not se.empty:
            parts.append('<h4>Species consensus events</h4>')
            parts.append(html_table(se,[c for c in ["event_key","event_type","exon","exon_pos","cds_pos","reference_codon","length","gene_disrupting","supporting_assemblies","callable_assemblies","support_fraction","species_event_confidence","support_accessions"] if c in se.columns]))
        if not sg.empty:
            parts.append('<details><summary>Exon recovery and confidence</summary>')
            parts.append(html_table(sg,[c for c in ["exon_index","exon","assemblies_total","assemblies_ok","coverage_median","coverage_max","callable_fraction_median","identity_median","n_aligned_fraction_median","confidence_score_median","high_conf_assemblies","medium_or_high_assemblies"] if c in sg.columns]))
            parts.append('</details>')
        if not assemblies.empty:
            parts.append('<details><summary>Assembly/locus details</summary>')
            parts.append(html_table(assemblies,[c for c in ["assembly_accession","assembly_name","organism_name","assembly_level","status","contig","strand","start","end","span","exon_count","order_coherence","score","locus_n_fraction","locus_longest_n_run","error"] if c in assemblies.columns]))
            parts.append('</details>')
        if not raw_events.empty:
            parts.append('<details><summary>All assembly-level event calls</summary>')
            parts.append(html_table(raw_events,[c for c in ["assembly_event_id","assembly_accession","organism_name","event_type","exon","exon_pos","cds_pos","reference_codon","length","confidence","gene_disrupting","motif","reason"] if c in raw_events.columns]))
            parts.append('</details>')
        parts.append('<p class="small"><a href="#tree">Back to tree</a> | <a href="#overview">Back to overview</a></p></article>')
    parts.append('</section>')

    parts.append('''<section id="methods" class="card"><h2>Methods and confidence framework</h2>
    <p><b>Ortholog locus:</b> BLAST HSPs are clustered by contig, strand and genomic proximity. Candidate loci are scored using the number of recovered exons, exon coverage, exon-order coherence and BLAST bit score. The best coherent locus is analyzed.</p>
    <p><b>Exon recovery:</b> sequence-recovered fraction is the fraction of pig exon positions with an aligned subject base. Callable fraction excludes N and gaps. Identity is computed over callable aligned A/C/G/T positions.</p>
    <p><b>Frameshifts:</b> explicit BLAST alignment gaps with length not divisible by three are called frameshifting. Large gaps between collinear HSPs are reported as candidate large indels and are deliberately assigned lower confidence unless independently corroborated.</p>
    <p><b>Premature stops:</b> codons are projected to the pig CDS coordinates. Stops created at callable codons before the reference terminal stop are reported. The reconstructed target CDS is additionally translated to identify stops downstream of frame-changing indels.</p>
    <p><b>Missing exons:</b> low recovery becomes a probable exon deletion only when adjacent exons identify the same locus and the genomic interval spanning the missing exon is assembled and low in N. Otherwise the result remains unresolved.</p>
    <p><b>Splice sites:</b> noncanonical two-base motifs at well-recovered exon boundaries are flagged as potential disruptions, not definitive loss, because transcript annotation and rare splice classes may require manual validation.</p>
    <p><b>Phylogenetic chronology:</b> event presence, confident absence and unknown states are mapped onto the supplied rooted species tree using Fitch parsimony. The report lists minimum-change origins/reversions and MRCA support. This is relative chronology. Absolute ages require a dated phylogeny.</p>
    <p><b>Repeats:</b> insertion sequences are passed to RepeatMasker when installed. Without RepeatMasker, the fallback can identify only tandem-repeat-like or low-complexity-like sequence architecture and does not claim a transposable-element family.</p>
    </section>''')
    parts.append('</main></body></html>')
    report.write_text(''.join(parts),encoding='utf-8')
    return report


# -----------------------------------------------------------------------------
# Output writer
# -----------------------------------------------------------------------------

def write_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    if df is None:
        df=pd.DataFrame()
    df2=df.copy()
    for c in df2.columns:
        df2[c]=df2[c].map(quote_tsv)
    df2.to_csv(path,sep='\t',index=False)


def load_existing_results(outdir: Path) -> List[Dict[str,Any]]:
    out=[]
    for p in sorted((outdir/"assemblies").glob("*/analysis.json")):
        try:
            out.append(json.loads(p.read_text()))
        except Exception as e:
            logging.warning("Cannot read %s: %s",p,e)
    return out


def aggregate_and_report(outdir: Path, targets: pd.DataFrame, ref_exons: List[RefExon], ref_cds: str,
                         tree_path: Path, repeatmasker_exe: Optional[str], repeat_threads:int) -> Dict[str,Path]:
    results=load_existing_results(outdir)
    locus,exon,events,variants=flatten_results(results)
    write_tsv(locus,outdir/"tables"/"assembly_locus_summary.tsv")
    write_tsv(exon,outdir/"tables"/"assembly_exon_metrics.tsv")
    write_tsv(events,outdir/"tables"/"assembly_events.tsv")
    write_tsv(variants,outdir/"tables"/"assembly_all_variants.tsv")

    sex=species_exon_summary(exon)
    sev=aggregate_species_events(events,exon)
    write_tsv(sex,outdir/"tables"/"species_exon_summary.tsv")
    write_tsv(sev,outdir/"tables"/"species_events.tsv")

    phy=phylogenetic_summary(sev,sex,tree_path)
    write_tsv(phy,outdir/"tables"/"phylogenetic_event_summary.tsv")

    ins=collect_insertions(results)
    rep=run_repeatmasker(ins,outdir,repeat_threads,repeatmasker_exe)
    rep=annotate_repeat_event_keys(rep) if not rep.empty else rep
    write_tsv(rep,outdir/"tables"/"repeat_insertions.tsv")
    sr=aggregate_species_repeats(rep,sex) if not rep.empty else pd.DataFrame()
    write_tsv(sr,outdir/"tables"/"species_repeat_summary.tsv")
    rphy=repeat_phylogeny(sr,sex,tree_path) if not sr.empty else pd.DataFrame()
    write_tsv(rphy,outdir/"tables"/"repeat_phylogenetic_summary.tsv")

    report=write_html_report(outdir,targets,locus,exon,events,sex,sev,phy,rep,rphy,tree_path,ref_exons,ref_cds)

    summary={
        "pipeline_version":VERSION,
        "assembly_targets":len(targets),
        "analysis_json_files":len(results),
        "status_counts":Counter(r.get("status","") for r in results),
        "species_keys":targets["Species Key"].nunique(),
        "species_event_calls":len(sev),
        "phylogenetic_event_groups":len(phy),
        "repeat_insertions":len(rep),
    }
    (outdir/"summary.json").write_text(json.dumps(summary,indent=2,default=str))
    return {"report":report,"tables":outdir/"tables","summary":outdir/"summary.json"}


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def parse_args():
    p=argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                              description="BLAST pig GPRC6A exons against many ruminant genomes and infer disruptive-event chronology.")
    p.add_argument("--targets-tsv",required=True,type=Path,help="Assembly table. Required columns: Assembly Accession, Organism Name. Ready-made ruminantia_assemblies.tsv is recommended.")
    p.add_argument("--exons",required=True,type=Path,help="Pig GPRC6A exon FASTA in CDS order")
    p.add_argument("--tree",required=True,type=Path,help="Rooted Newick species tree for phylogenetic inference")
    p.add_argument("--outdir",required=True,type=Path)
    p.add_argument("--genome-map",type=Path,default=None,help="Optional TSV with Assembly Accession and local Genome FASTA columns; bypasses NCBI download for mapped accessions")
    p.add_argument("--workers",type=int,default=1,help="Assemblies processed concurrently. Keep modest because mammalian genomes are IO-heavy.")
    p.add_argument("--blast-threads",type=int,default=4,help="Threads per blastn process")
    p.add_argument("--repeat-threads",type=int,default=4)
    p.add_argument("--max-assemblies",type=int,default=0,help="0=all; useful for a test run")
    p.add_argument("--start-at",type=int,default=0,help="0-based row offset after target filtering")
    p.add_argument("--num-shards",type=int,default=1,help="Optional deterministic sharding for HPC/job arrays")
    p.add_argument("--shard-index",type=int,default=0,help="0-based shard index")
    p.add_argument("--drop-paired-duplicates",action="store_true",help="Prefer GCF over a paired GCA representation. By default every accession is analyzed, including GCA/GCF pairs.")
    p.add_argument("--resume",action="store_true",help="Reuse per-assembly analysis.json and raw BLAST files")
    p.add_argument("--force-blast",action="store_true")
    p.add_argument("--cleanup-genomes",action="store_true",help="Delete downloaded genome package after successful assembly analysis; raw BLAST and analysis JSON are retained")
    p.add_argument("--download-retries",type=int,default=4,help="Number of NCBI datasets download attempts per assembly before FTP fallback")
    p.add_argument("--download-backoff",type=float,default=15.0,help="Initial retry delay in seconds; subsequent datasets attempts use exponential backoff")
    p.add_argument("--no-ftp-fallback",action="store_true",help="Disable direct fallback to the official NCBI genomes FTP/HTTPS tree")
    p.add_argument("--max-locus-gap",type=int,default=2_000_000,help="Maximum genomic distance joining BLAST hits into one candidate locus")
    p.add_argument("--locus-flank",type=int,default=10_000,help="Flank around selected locus for N-content metrics")
    p.add_argument("--exon-n-flank",type=int,default=100,help="Flank around recovered exon for local N metrics")
    p.add_argument("--large-indel-threshold",type=int,default=20)
    p.add_argument("--keep-locus-candidates",type=int,default=5)
    p.add_argument("--report-only",action="store_true",help="Skip genome analysis and rebuild tables/report from existing assembly JSON files")
    p.add_argument("--skip-report",action="store_true",help="Analyze assemblies but do not aggregate/report. Useful for concurrent HPC shards; run --report-only afterward.")
    p.add_argument("--log-level",choices=["DEBUG","INFO","WARNING","ERROR"],default="INFO")
    return p.parse_args()


def main():
    args=parse_args()
    args.outdir=args.outdir.resolve(); args.outdir.mkdir(parents=True,exist_ok=True)
    logging.basicConfig(level=getattr(logging,args.log_level),format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(sys.stderr),logging.FileHandler(args.outdir/"pipeline.log")])
    logging.info("GPRC6A pipeline v%s",VERSION)

    blastn_exe=require_exe("blastn", required=not args.report_only)
    datasets_exe=require_exe("datasets", required=False)
    repeatmasker_exe=require_exe("RepeatMasker", required=False)
    if not args.report_only and not datasets_exe and not args.genome_map:
        raise RuntimeError("NCBI datasets CLI is required unless --genome-map provides local genome FASTA files")
    if not repeatmasker_exe:
        logging.warning("RepeatMasker not found. Repeat family/class annotation will use conservative sequence-composition fallback only.")

    ref_exons,ref_cds=load_reference_exons(args.exons)
    targets=load_targets(args.targets_tsv,args.drop_paired_duplicates)
    if args.num_shards<1 or not (0<=args.shard_index<args.num_shards):
        raise ValueError("Invalid --num-shards/--shard-index")
    targets=targets.iloc[args.start_at:].reset_index(drop=True)
    if args.num_shards>1:
        targets=targets.iloc[[i for i in range(len(targets)) if i%args.num_shards==args.shard_index]].reset_index(drop=True)
    if args.max_assemblies>0:
        targets=targets.head(args.max_assemblies).copy()
    logging.info("Targets in this run: %d accessions, %d species keys",len(targets),targets["Species Key"].nunique())

    # Preserve exact run inputs.
    targets_file = args.outdir/(f"targets_used.shard{args.shard_index:03d}.tsv" if args.num_shards>1 and not args.report_only else "targets_used.tsv")
    targets.to_csv(targets_file,sep="\t",index=False)
    shutil.copy2(args.exons,args.outdir/"reference_exons.fa")
    shutil.copy2(args.tree,args.outdir/"species_tree.newick")

    if not args.report_only:
        gmap={k:str(v) for k,v in load_genome_map(args.genome_map).items()}
        args_dict=vars(args).copy()
        args_dict.update({"outdir":str(args.outdir),"exons":str(args.exons),"blastn_exe":blastn_exe,
                          "datasets_exe":datasets_exe or "datasets"})
        refser=[asdict(e) for e in ref_exons]
        rows=targets.to_dict(orient="records")
        if args.workers<=1:
            for i,row in enumerate(rows,1):
                logging.info("[%d/%d] %s | %s",i,len(rows),row["Assembly Accession"],row["Organism Name"])
                analyze_assembly(row,args_dict,refser,ref_cds,gmap)
        else:
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                fut={ex.submit(analyze_assembly,row,args_dict,refser,ref_cds,gmap):row for row in rows}
                done=0
                for f in as_completed(fut):
                    done+=1; row=fut[f]
                    try:
                        r=f.result(); logging.info("[%d/%d] completed %s: %s",done,len(rows),row["Assembly Accession"],r.get("status"))
                    except Exception:
                        logging.exception("Unhandled worker failure for %s",row["Assembly Accession"])

    if args.skip_report and not args.report_only:
        logging.info("Assembly analysis complete; aggregation/report skipped by --skip-report")
        print(args.outdir)
        return

    outputs=aggregate_and_report(args.outdir,targets,ref_exons,ref_cds,args.tree,repeatmasker_exe,args.repeat_threads)
    logging.info("HTML report: %s",outputs["report"])
    logging.info("Tables: %s",outputs["tables"])
    print(outputs["report"])


if __name__=="__main__":
    main()
