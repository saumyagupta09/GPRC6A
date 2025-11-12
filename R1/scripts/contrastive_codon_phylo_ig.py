#!/usr/bin/env python3
"""
contrastive_codon_phylo_ig.py

Contrastive (SimCLR) learning on codon alignments with:
 - PHYLIP (relaxed) alignment input (suitable for HyPhy)
 - Newick tree input -> patristic distances -> phylogenetic down-weighting
 - Integrated Gradients attribution in embedding space per sequence
 - Lineage-weighted per-site importance and lineage-level selection estimates (purifying / neutral / diversifying)

Usage:
    python contrastive_codon_phylo_ig.py --phy alignment.phy --map mapping.tsv --tree tree.nwk --outdir results --epochs 50

Required Python packages:
    numpy, pandas, biopython, torch, sklearn, matplotlib, scipy

Notes:
 - mapping.tsv: tab/space-separated file with two columns: sequence_name  lineage_label
 - Sequences in PHYLIP must be codon-aligned (length % 3 == 0) and mapping names must match PHYLIP names.
 - The script writes CSV and PNG outputs into --outdir.
"""
import os
import argparse
import math
import random
from collections import defaultdict

import numpy as np
import pandas as pd
from Bio import Phylo
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.manifold import TSNE
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import accuracy_score
from scipy.stats import chi2_contingency
import matplotlib.pyplot as plt

# ---------------------------
# Codon utilities & tables
# ---------------------------
CODON_LIST = [a+b+c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
CODON_TO_INDEX = {c:i for i,c in enumerate(CODON_LIST)}
INDEX_TO_CODON = {i:c for c,i in CODON_TO_INDEX.items()}

CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
AA_TO_CODONS = defaultdict(list)
for codon,aa in CODON_TABLE.items():
    AA_TO_CODONS[aa].append(codon)

# ---------------------------
# I/O helpers
# ---------------------------
def read_phylip_relaxed(path):
    """Read a relaxed PHYLIP (name seq on each line). Returns dict name->sequence."""
    seqs = {}
    with open(path) as f:
        header = f.readline()  # ignored beyond reading
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = ln.split()
            name = parts[0]
            seq = ''.join(parts[1:]).upper().replace('U','T')
            seqs[name] = seq
    return seqs

def read_mapping(path):
    """Simple mapping parser: seqname <tab/space> lineage"""
    mp = {}
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            parts = ln.split()
            if len(parts) >= 2:
                mp[parts[0]] = parts[1]
    return mp

def read_newick(path):
    if path is None:
        return None
    return Phylo.read(path, "newick")

def seq_to_codon_indices(seq):
    L = len(seq)
    if L % 3 != 0:
        raise ValueError("Sequence length not divisible by 3")
    out = []
    for i in range(0, L, 3):
        cod = seq[i:i+3]
        if any(ch not in "ATCG" for ch in cod):
            out.append(CODON_TO_INDEX['TTT'])
        else:
            out.append(CODON_TO_INDEX.get(cod, CODON_TO_INDEX['TTT']))
    return np.array(out, dtype=np.int64)

# ---------------------------
# Phylogenetic weighting
# ---------------------------
def patristic_distances(tree, names):
    """Return pairwise distance matrix between names using tree.distance(a,b)."""
    n = len(names)
    D = np.zeros((n,n), dtype=float)
    for i,a in enumerate(names):
        for j,b in enumerate(names):
            if i <= j:
                try:
                    d = tree.distance(a, b)
                except Exception:
                    # if names mismatch, try strip or fallback to large distance
                    d = float('inf')
                D[i,j] = d
                D[j,i] = d
    return D

def phylogenetic_weights(D, lam=1.0, eps=1e-12):
    """
    Given distance matrix D (n x n), compute weights w_i proportional to
        1 / (1 + sum_{j != i} exp(-lam * d_ij))
    Normalize to sum to 1.
    """
    n = D.shape[0]
    w = np.zeros(n, dtype=float)
    for i in range(n):
        # exclude self (d_ii = 0 typically)
        s = np.sum(np.exp(-lam * D[i,:])) - 1.0  # subtract exp(-lam*0)
        w[i] = 1.0 / (1.0 + s + eps)
    w = w / (w.sum() + eps)
    return w

# ---------------------------
# Augmentations
# ---------------------------
def random_mask(codon_indices, p=0.08):
    arr = codon_indices.copy()
    mask = np.random.rand(len(arr)) < p
    arr[mask] = CODON_TO_INDEX['TTT']  # use TTT as mask placeholder
    return arr

def random_synonymous(codon_indices, p=0.05):
    arr = codon_indices.copy()
    for i,cidx in enumerate(arr):
        if np.random.rand() < p:
            codon = INDEX_TO_CODON[int(cidx)]
            aa = CODON_TABLE.get(codon, 'X')
            choices = AA_TO_CODONS.get(aa, [codon])
            if len(choices) > 1:
                arr[i] = CODON_TO_INDEX[random.choice(choices)]
    return arr

def augment(codon_indices):
    if np.random.rand() < 0.5:
        arr = random_mask(codon_indices, p=0.08)
    else:
        arr = codon_indices.copy()
    if np.random.rand() < 0.6:
        arr = random_synonymous(arr, p=0.05)
    return arr

# ---------------------------
# Model & loss
# ---------------------------
class CodonEncoder(nn.Module):
    def __init__(self, codon_vocab=64, emb_dim=48, hidden=128, proj=64):
        super().__init__()
        self.emb = nn.Embedding(codon_vocab, emb_dim, padding_idx=CODON_TO_INDEX['TTT'])
        self.fc = nn.Sequential(nn.Linear(emb_dim, hidden), nn.ReLU(), nn.Linear(hidden, proj))
    def forward(self, x):
        # x: (batch, L)
        e = self.emb(x)               # (batch, L, emb_dim)
        m = e.mean(dim=1)             # mean pooling across positions
        out = self.fc(m)
        return nn.functional.normalize(out, dim=1)

def nt_xent_loss_weighted(z1, z2, weights, idxs_in_batch, temperature=0.5):
    """
    Compute per-sample NT-Xent loss and return weighted average according to 'weights' vector.
    - z1, z2: (batch, dim)
    - weights: full-sample weights (length N_total)
    - idxs_in_batch: indices (into full data) of samples in batch, length batch
    We'll compute per-sample loss for each of the 2*batch augmented samples and weight them using weights[idx].
    For simplicity, map each augmented sample back to its original sample weight.
    """
    device = z1.device
    batch = z1.size(0)
    z = torch.cat([z1, z2], dim=0)  # (2*batch, dim)
    sim = torch.matmul(z, z.T) / temperature  # (2b,2b)
    large_neg = -9e15
    diag_mask = torch.eye(2*batch, device=device).bool()
    sim = sim.masked_fill(diag_mask, large_neg)

    # positives: i-th in z1 matches i-th in z2 (offset)
    pos = torch.cat([torch.diag(sim, batch), torch.diag(sim, -batch)], dim=0)  # (2b,)
    exp_sim = torch.exp(sim)
    denom = torch.sum(exp_sim, dim=1)  # (2b,)

    losses = -torch.log(torch.exp(pos) / denom + 1e-12)  # (2b,)

    # compute weights per augmented view: for first batch positions use weights[idxs_in_batch],
    # for second batch the same mapping
    w_batch = torch.tensor(weights[np.array(idxs_in_batch)], dtype=torch.float, device=device)  # (batch,)
    w2 = torch.cat([w_batch, w_batch], dim=0)  # (2b,)
    if w2.sum().item() == 0:
        # fallback unweighted mean
        return losses.mean()
    weighted_loss = torch.sum(losses * w2) / (w2.sum() + 1e-12)
    return weighted_loss

# ---------------------------
# Integrated Gradients (embedding-space IG)
# ---------------------------
def integrated_gradients(model, classifier, seq_idx, baseline_idx=None, steps=25, device='cpu'):
    """
    seq_idx: 1-D numpy array of codon indices (length L)
    classifier: torch.nn.Linear mapping embedding->class logits (expects model.fc output)
    Returns: 1-D numpy array of length L with per-site attribution magnitudes
    """
    model.eval()
    if baseline_idx is None:
        baseline_idx = np.array([CODON_TO_INDEX['TTT']] * len(seq_idx), dtype=np.int64)
    x = torch.tensor(seq_idx[None,:], dtype=torch.long, device=device)
    b = torch.tensor(baseline_idx[None,:], dtype=torch.long, device=device)
    emb_b = model.emb(b).detach()   # (1, L, emb_dim)
    emb_x = model.emb(x).detach()   # (1, L, emb_dim)
    total_grad = torch.zeros_like(emb_x).to(device)
    for alpha in np.linspace(0.0, 1.0, steps+1)[1:]:
        inter = emb_b + (emb_x - emb_b) * float(alpha)
        inter = inter.detach().requires_grad_(True)
        m = inter.mean(dim=1)  # (1, emb_dim)
        out = model.fc(m)      # (1, proj_dim)
        logits = classifier(out)
        pred = logits.argmax(dim=1)
        loss = nn.CrossEntropyLoss()(logits, pred)
        loss.backward()
        total_grad += inter.grad.clone()
        model.zero_grad()
    avg_grad = total_grad / steps
    attr = (emb_x - emb_b) * avg_grad
    site_attr = torch.sum(torch.abs(attr), dim=2).cpu().numpy()[0]  # (L,)
    return site_attr

# ---------------------------
# High-level pipeline
# ---------------------------
def run_pipeline(phy_path, map_path, tree_path, outdir,
                 epochs=50, batch_size=32, device='cpu', lam=1.0, ig_steps=25):
    os.makedirs(outdir, exist_ok=True)
    # read inputs
    seqs = read_phylip_relaxed(phy_path)
    mapping = read_mapping(map_path)
    tree = None
    if tree_path:
        tree = read_newick(tree_path)

    # filter sequences to those in mapping
    seqs = {n:seqs[n] for n in seqs if n in mapping}
    if len(seqs) == 0:
        raise RuntimeError("No sequences from PHYLIP matched names in mapping file.")

    names = list(seqs.keys())
    L_nt = len(next(iter(seqs.values())))
    if any(len(s) != L_nt for s in seqs.values()):
        raise RuntimeError("Aligned sequences have differing lengths.")
    if L_nt % 3 != 0:
        raise RuntimeError("Sequences length not divisible by 3 (not codon-aligned).")
    L = L_nt // 3

    # codon matrix
    codon_matrix = {n: seq_to_codon_indices(seqs[n]) for n in names}
    # compute phylogenetic distances and weights
    if tree is not None:
        D = patristic_distances(tree, names)
        # convert any infinite distances to large number
        D[np.isinf(D)] = np.nanmax(D[np.isfinite(D)]) * 2.0 if np.any(np.isfinite(D)) else 1e6
        weights = phylogenetic_weights(D, lam=lam)
    else:
        # uniform weights
        weights = np.ones(len(names), dtype=float) / len(names)

    # map names to index
    name_to_idx = {n:i for i,n in enumerate(names)}

    # training data arrays
    X = np.stack([codon_matrix[n] for n in names])  # shape (N, L)
    Y = np.array([mapping[n] for n in names])
    le = LabelEncoder(); y_enc = le.fit_transform(Y)

    device = torch.device(device if torch.cuda.is_available() and device.startswith('cuda') else 'cpu')

    # model
    model = CodonEncoder()
    model.to(device)
    opt = optim.Adam(model.parameters(), lr=1e-3)

    N = X.shape[0]
    # precompute sample global weights in same order as names
    sample_weights = weights  # numpy array length N

    # training loop (weighted)
    for epoch in range(epochs):
        model.train()
        perm = np.random.permutation(N)
        epoch_losses = []
        for i in range(0, N, batch_size):
            batch_idx = perm[i:i+batch_size]
            batch = X[batch_idx]  # numpy shape (b, L)
            # two augmentations
            a1 = np.stack([augment(b) for b in batch])
            a2 = np.stack([augment(b) for b in batch])
            a1_t = torch.tensor(a1, dtype=torch.long, device=device)
            a2_t = torch.tensor(a2, dtype=torch.long, device=device)
            z1 = model(a1_t)
            z2 = model(a2_t)
            loss = nt_xent_loss_weighted(z1, z2, sample_weights, batch_idx, temperature=0.5)
            opt.zero_grad()
            loss.backward()
            opt.step()
            epoch_losses.append(loss.item())
        if (epoch+1) % max(1, epochs//5) == 0 or epoch == 0:
            print(f"[SimCLR] Epoch {epoch+1}/{epochs}  mean loss: {np.mean(epoch_losses):.4f}")

    # save model
    torch.save(model.state_dict(), os.path.join(outdir, "simclr_model.pt"))

    # embeddings (full set)
    model.eval()
    with torch.no_grad():
        X_t = torch.tensor(X, dtype=torch.long, device=device)
        embeds = model(X_t).cpu().numpy()

    emb_df = pd.DataFrame(embeds, index=names)
    emb_df['lineage'] = Y
    emb_csv = os.path.join(outdir, "embeddings.csv")
    emb_df.to_csv(emb_csv)

    # t-SNE (adjust perplexity)
    perp = min(30, max(2, max(2, N//3)))
    tsne = TSNE(n_components=2, perplexity=perp, random_state=0)
    z2d = tsne.fit_transform(embeds)
    plt.figure(figsize=(6,5))
    unique_labels = np.unique(Y)
    for lab in unique_labels:
        sel = (Y == lab)
        plt.scatter(z2d[sel,0], z2d[sel,1], label=str(lab), s=20)
    plt.legend()
    plt.title("t-SNE of SimCLR embeddings")
    plt.tight_layout()
    tsne_path = os.path.join(outdir, "tsne.png")
    plt.savefig(tsne_path)
    plt.close()

    # downstream classifier (logistic) to inspect separability
    clf = LogisticRegression(max_iter=2000)
    clf.fit(embeds, y_enc)
    preds = clf.predict(embeds)
    acc = accuracy_score(y_enc, preds)
    print("Downstream classifier accuracy (self):", acc)

    # convert sklearn classifier to torch Linear for attribution
    classifier = nn.Linear(embeds.shape[1], len(le.classes_)).to(device)
    classifier.weight.data = torch.tensor(clf.coef_, dtype=torch.float)
    classifier.bias.data = torch.tensor(clf.intercept_, dtype=torch.float)

    # per-sample Integrated Gradients and grad-saliency
    grad_sal_all = np.zeros((N, L), dtype=float)
    ig_all = np.zeros((N, L), dtype=float)

    print("Computing per-sequence attributions (Integrated Gradients)...")
    for i, name in enumerate(names):
        seq_idx = X[i]
        # grad-saliency (single step)
        seq_t = torch.tensor(seq_idx[None,:], dtype=torch.long, device=device)
        emb_layer = model.emb(seq_t).detach().requires_grad_(True)
        m = emb_layer.mean(dim=1)
        out = model.fc(m)
        logits = classifier(out)
        target = torch.tensor([y_enc[i]], dtype=torch.long, device=device)
        loss = nn.CrossEntropyLoss()(logits, target)
        loss.backward()
        g = emb_layer.grad.detach().cpu().numpy()[0]  # (L, emb)
        sal = np.sum(np.abs(g), axis=1)
        grad_sal_all[i,:] = sal
        model.zero_grad()

        # Integrated Gradients
        ig = integrated_gradients(model, classifier, seq_idx, baseline_idx=None, steps=ig_steps, device=device)
        ig_all[i,:] = ig

    # Lineage-weighted per-site importance
    # For each lineage, compute weighted mean IG across sequences in lineage (weights normalized within lineage)
    lineage_names = np.unique(Y)
    lineage_site_matrix = np.zeros((len(lineage_names), L), dtype=float)
    lineage_weights = {}
    for li, lineage in enumerate(lineage_names):
        inds = [i for i,name in enumerate(names) if mapping[names[i]] == lineage]
        if len(inds) == 0:
            continue
        # use global phylogenetic weights for sequences (sample_weights)
        w = sample_weights[np.array(inds)]
        if np.sum(w) == 0:
            w = np.ones_like(w) / len(w)
        else:
            w = w / np.sum(w)
        lineage_weights[lineage] = w
        lineage_site_matrix[li,:] = np.sum(ig_all[inds,:] * w[:,None], axis=0)  # weighted mean

    # also compute global weighted site importance (weighted across all sequences by sample_weights)
    global_site_importance = np.sum(ig_all * sample_weights[:,None], axis=0)

    # Save site importance (global and per-lineage)
    site_df = pd.DataFrame({
        'site': np.arange(1, L+1),
        'global_ig': global_site_importance,
        'global_grad_saliency': np.sum(grad_sal_all * sample_weights[:,None], axis=0)
    })
    # append lineage columns
    for li,lineage in enumerate(lineage_names):
        site_df[f'ig_{lineage}'] = lineage_site_matrix[li,:]
    site_csv = os.path.join(outdir, "site_importance_weighted.csv")
    site_df.to_csv(site_csv, index=False)

    # lineage-level selection estimates
    # For each lineage, summarize mean IG across sites and derive selection type.
    # We'll compute z-score of lineage mean vs global mean and classify:
    #   z > 1 -> diversifying; z < -1 -> purifying; else neutral
    lineage_summary = []
    global_mean = np.mean(global_site_importance)
    global_std = np.std(global_site_importance) + 1e-12
    for li,lineage in enumerate(lineage_names):
        vals = lineage_site_matrix[li,:]
        mean = float(np.mean(vals))
        median = float(np.median(vals))
        var = float(np.var(vals))
        # z-score
        z = (mean - global_mean) / global_std
        if z > 1.0:
            inferred = "diversifying"
        elif z < -1.0:
            inferred = "purifying"
        else:
            inferred = "neutral"
        lineage_summary.append({
            'lineage': lineage,
            'mean_ig': mean,
            'median_ig': median,
            'var_ig': var,
            'zscore_vs_global': z,
            'inferred_selection': inferred
        })
    lineage_df = pd.DataFrame(lineage_summary)
    lineage_csv = os.path.join(outdir, "lineage_selection_summary.csv")
    lineage_df.to_csv(lineage_csv, index=False)

    # For the top IG sites per lineage, do chi-square test of codon counts across lineages (for a quick sanity check)
    topk = max(5, int(0.05 * L))
    stats_rows = []
    for li,lineage in enumerate(lineage_names):
        top_sites = np.argsort(-lineage_site_matrix[li,:])[:topk]
        for s in top_sites:
            # contingency table: rows=lineages, cols=codon counts (reduced to non-zero columns)
            table = []
            for lab in lineage_names:
                inds_lab = [i for i,name in enumerate(names) if mapping[names[i]] == lab]
                cod_counts = [int(np.sum(X[inds_lab, s] == c)) for c in range(64)]
                table.append(cod_counts)
            table = np.array(table)
            nonzero = table.sum(axis=0) > 0
            if np.sum(nonzero) < 2:
                pval = 1.0
            else:
                try:
                    chi, pval, dof, exp = chi2_contingency(table[:, nonzero])
                except Exception:
                    pval = float('nan')
            stats_rows.append({'lineage': lineage, 'site': s+1, 'pval': pval})
    stats_df = pd.DataFrame(stats_rows)
    stats_csv = os.path.join(outdir, "site_stats.csv")
    stats_df.to_csv(stats_csv, index=False)

    # Heatmap: lineage x site matrix (IG)
    plt.figure(figsize=(max(6, L*0.15), max(3, len(lineage_names)*0.6)))
    im = plt.imshow(lineage_site_matrix, aspect='auto', interpolation='nearest')
    plt.colorbar(im, label='weighted IG attribution')
    plt.yticks(np.arange(len(lineage_names)), lineage_names)
    plt.xticks(np.arange(0, L, max(1, L//20)), np.arange(1, L+1)[::max(1, L//20)])
    plt.xlabel('codon site')
    plt.title('Lineage-weighted per-site Integrated Gradients')
    heatmap_path = os.path.join(outdir, "lineage_site_heatmap.png")
    plt.tight_layout()
    plt.savefig(heatmap_path)
    plt.close()

    # final summary
    summary = {
        'embeddings': emb_csv,
        'tsne': tsne_path,
        'site_importance': site_csv,
        'site_stats': stats_csv,
        'lineage_summary': lineage_csv,
        'heatmap': heatmap_path,
        'model': os.path.join(outdir, "simclr_model.pt"),
        'downstream_accuracy_self': float(acc)
    }
    print("Pipeline finished. Outputs written to:", outdir)
    return summary

# ---------------------------
# CLI
# ---------------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Contrastive codon SimCLR with phylogenetic weighting and IG attributions")
    p.add_argument('--phy', required=True, help='PHYLIP alignment (relaxed format: name seq per line)')
    p.add_argument('--map', required=True, help='Mapping file: seqname <tab> lineage')
    p.add_argument('--tree', required=False, help='Newick tree file (optional but recommended)')
    p.add_argument('--outdir', default='results', help='Output directory')
    p.add_argument('--epochs', type=int, default=50, help='Training epochs (SimCLR)')
    p.add_argument('--batch_size', type=int, default=32, help='Batch size')
    p.add_argument('--device', default='cpu', help='Device string (cpu or cuda:0)')
    p.add_argument('--lam', type=float, default=1.0, help='Phylogenetic decay lambda for weighting')
    p.add_argument('--ig_steps', type=int, default=25, help='Integrated Gradients steps')
    args = p.parse_args()

    summary = run_pipeline(
        phy_path=args.phy,
        map_path=args.map,
        tree_path=args.tree,
        outdir=args.outdir,
        epochs=args.epochs,
        batch_size=args.batch_size,
        device=args.device,
        lam=args.lam,
        ig_steps=args.ig_steps
    )
    print("Summary of outputs:")
    for k,v in summary.items():
        print(f" - {k}: {v}")
