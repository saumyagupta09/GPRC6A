#!/usr/bin/env python3
"""
Fetch a conservative panel of independently annotated intact mammalian GPRC6A
CDSs from NCBI RefSeq and build a combined mammalian + ruminant species tree for
the codon/PAML follow-up.

The script does NOT infer intactness from gene name alone. Candidate RefSeq mRNA
records are inspected at the CDS level and rejected if they are partial, contain
internal stops, are off-frame, are highly ambiguous, or are grossly discordant
in length with the supplied GPRC6A reference CDS.

Default external panel:
  Homo sapiens
  Mus musculus
  Rattus norvegicus
  Oryctolagus cuniculus
  Canis lupus familiaris
  Felis catus
  Equus caballus

Sus scrofa is already present as the reference CDS in the main follow-up and is
therefore represented in the combined tree as Sus_scrofa_REF rather than
downloaded a second time.
"""

from __future__ import annotations
import argparse
import copy
import logging
import os
import re
import time
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import pandas as pd
from Bio import Entrez, Phylo, SeqIO
from Bio.Seq import Seq
from Bio.Phylo.BaseTree import Tree, Clade

VERSION = "1.0.0"
STOP = {"TAA","TAG","TGA"}

DEFAULT_SPECIES = [
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Oryctolagus cuniculus",
    "Canis lupus familiaris",
    "Felis catus",
    "Equus caballus",
]

def safe_name(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(s)).strip("_") or "unnamed"

def wrap(seq: str, n: int = 80) -> str:
    return "\n".join(seq[i:i+n] for i in range(0, len(seq), n))

def load_reference_cds(path: Path) -> str:
    recs = list(SeqIO.parse(str(path), "fasta"))
    if not recs:
        raise ValueError(f"No reference exon records in {path}")
    def key(r):
        m = re.search(r"(?:exon[_\s-]*)?(\d+)", r.id, re.I)
        return int(m.group(1)) if m else 10**9
    recs.sort(key=key)
    return "".join(str(r.seq).upper().replace("U","T") for r in recs)

def is_partial_location(feature) -> bool:
    s = str(feature.location)
    return "<" in s or ">" in s

def candidate_from_record(rec, species: str, ref_len: int) -> List[Dict[str,Any]]:
    out = []
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        q = feat.qualifiers
        gene = " ".join(q.get("gene", []))
        product = " ".join(q.get("product", []))
        text = (gene + " " + product).lower()
        if "gprc6a" not in text and "g protein-coupled receptor family c group 6 member a" not in text:
            continue
        seq = str(feat.extract(rec.seq)).upper().replace("U","T")
        terminal_stop = seq[-3:] in STOP if len(seq) >= 3 else False
        if terminal_stop:
            seq = seq[:-3]
        internal_stops = 0
        if len(seq) >= 3:
            n = len(seq) - (len(seq) % 3)
            aa = str(Seq(seq[:n]).translate())
            internal_stops = aa.count("*")
        nfrac = seq.count("N") / len(seq) if seq else 1.0
        length_ratio = len(seq) / ref_len if ref_len else 0.0
        accession = str(rec.id)
        complete = not is_partial_location(feat)
        valid = (
            complete
            and len(seq) % 3 == 0
            and internal_stops == 0
            and nfrac <= 0.005
            and 0.70 <= length_ratio <= 1.35
            and len(seq) >= 1500
        )
        prefix_score = 120 if accession.startswith("NM_") else (40 if accession.startswith("XM_") else 0)
        score = (
            prefix_score
            + (60 if complete else 0)
            + 50 * max(0.0, 1.0 - abs(1.0 - length_ratio))
            + 30 * (1.0 - min(1.0, nfrac))
            - 200 * internal_stops
        )
        out.append({
            "species":species, "accession":accession, "sequence":seq,
            "complete_cds":complete, "internal_stops":internal_stops,
            "n_fraction":nfrac, "length_nt":len(seq), "length_ratio":length_ratio,
            "score":score, "valid":valid, "product":product, "gene":gene,
        })
    return out

def esearch_ids(species: str, retmax: int = 30) -> List[str]:
    queries = [
        f'"{species}"[Organism] AND GPRC6A[Gene Name] AND refseq[filter] AND biomol_mrna[PROP]',
        f'"{species}"[Organism] AND GPRC6A[All Fields] AND refseq[filter] AND biomol_mrna[PROP]',
        f'"{species}"[Organism] AND "G protein-coupled receptor family C group 6 member A"[All Fields] AND refseq[filter]',
    ]
    for term in queries:
        for attempt in range(3):
            try:
                with Entrez.esearch(db="nuccore", term=term, retmax=retmax) as h:
                    ids = Entrez.read(h).get("IdList", [])
                if ids:
                    return list(ids)
                break
            except Exception as e:
                if attempt == 2:
                    logging.warning("ESearch failed for %s: %s", species, e)
                time.sleep(2*(attempt+1))
    return []

def fetch_records(ids: List[str]):
    if not ids:
        return []
    for attempt in range(3):
        try:
            with Entrez.efetch(db="nuccore", id=",".join(ids), rettype="gb", retmode="text") as h:
                return list(SeqIO.parse(h, "genbank"))
        except Exception as e:
            if attempt == 2:
                raise
            logging.warning("EFetch retry %d: %s", attempt+1, e)
            time.sleep(3*(attempt+1))
    return []

def fetch_best(species: str, ref_len: int) -> Dict[str,Any]:
    ids = esearch_ids(species)
    if not ids:
        raise RuntimeError(f"No RefSeq GPRC6A mRNA records found for {species}")
    recs = fetch_records(ids)
    candidates = []
    for rec in recs:
        candidates.extend(candidate_from_record(rec, species, ref_len))
    valid = [x for x in candidates if x["valid"]]
    if not valid:
        brief = "; ".join(
            f'{x["accession"]}:len={x["length_nt"]}:stops={x["internal_stops"]}:partial={not x["complete_cds"]}'
            for x in sorted(candidates, key=lambda z:z["score"], reverse=True)[:8]
        )
        raise RuntimeError(f"No high-quality intact CDS passed validation for {species}. Candidates: {brief}")
    return max(valid, key=lambda x:(x["score"], x["length_nt"]))

def combine_tree(ruminant_tree_path: Path) -> Tree:
    rt = Phylo.read(str(ruminant_tree_path), "newick")
    rum = copy.deepcopy(rt.root)
    rum.name = None

    def leaf(name):
        return Clade(name=name)

    rodents = Clade(clades=[leaf("Mus_musculus"), leaf("Rattus_norvegicus")])
    glires = Clade(clades=[leaf("Oryctolagus_cuniculus"), rodents])
    euarchontoglires = Clade(clades=[leaf("Homo_sapiens"), glires])

    carnivora = Clade(clades=[leaf("Canis_lupus_familiaris"), leaf("Felis_catus")])
    cetartiodactyla_selected = Clade(clades=[leaf("Sus_scrofa_REF"), rum])
    euungulata_selected = Clade(clades=[leaf("Equus_caballus"), cetartiodactyla_selected])
    laurasiatheria_selected = Clade(clades=[carnivora, euungulata_selected])

    root = Clade(clades=[euarchontoglires, laurasiatheria_selected])
    tree = Tree(root=root, rooted=False)
    for c in tree.find_clades():
        c.branch_length = None
        if c.clades:
            c.name = None
    return tree

def parse_args():
    p = argparse.ArgumentParser(description="Fetch validated intact GPRC6A mammalian outgroups for MACSE/PAML")
    p.add_argument("--outdir", type=Path, required=True, help="Main GPRC6A_Ruminantia output directory")
    p.add_argument("--email", required=True, help="Email sent to NCBI Entrez")
    p.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"), help="Optional NCBI API key")
    p.add_argument("--reference-exons", type=Path, default=None)
    p.add_argument("--ruminant-tree", type=Path, default=None)
    p.add_argument("--species", action="append", default=None,
                   help="Override default external panel; repeat for each species. Custom panels require your own combined tree.")
    p.add_argument("--sleep", type=float, default=0.35)
    p.add_argument("--force", action="store_true")
    p.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    return p.parse_args()

def main():
    args = parse_args()
    root = args.outdir.expanduser().resolve()
    out = root / "external_outgroups"
    out.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(level=getattr(logging,args.log_level),
                        format="%(asctime)s %(levelname)s %(message)s")

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    ref = args.reference_exons.expanduser().resolve() if args.reference_exons else root/"reference_exons.fa"
    rtree = args.ruminant_tree.expanduser().resolve() if args.ruminant_tree else root/"refined_intactness"/"species_tree.newick"
    if not ref.exists():
        raise FileNotFoundError(ref)
    if not rtree.exists():
        raise FileNotFoundError(rtree)

    ref_cds = load_reference_cds(ref)
    species_list = args.species if args.species else list(DEFAULT_SPECIES)

    fasta = out/"external_intact_gprc6a.fasta"
    meta = out/"external_intact_metadata.tsv"

    rows = []
    with open(fasta, "w") as fh:
        for i, sp in enumerate(species_list, 1):
            logging.info("[%d/%d] Fetching %s GPRC6A RefSeq CDS", i, len(species_list), sp)
            best = fetch_best(sp, len(ref_cds))
            sid = safe_name(sp)
            fh.write(f">{sid}|{best['accession']}\n{wrap(best['sequence'])}\n")
            rows.append({k:v for k,v in best.items() if k!="sequence"})
            logging.info("  selected %s | %d nt | ratio %.3f | Ns %.4f",
                         best["accession"], best["length_nt"], best["length_ratio"], best["n_fraction"])
            time.sleep(max(0.0,args.sleep))

    pd.DataFrame(rows).to_csv(meta, sep="\t", index=False)

    if args.species:
        logging.warning("Custom --species panel supplied. FASTA/metadata written, but automatic tree scaffold is disabled. Supply a biologically correct combined tree yourself.")
    else:
        tree = combine_tree(rtree)
        tree_path = out/"combined_species_tree.newick"
        Phylo.write(tree, str(tree_path), "newick")
        logging.info("Combined tree: %s", tree_path)

    logging.info("External intact FASTA: %s", fasta)
    logging.info("Metadata: %s", meta)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
