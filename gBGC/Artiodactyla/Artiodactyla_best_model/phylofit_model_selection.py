#!/usr/bin/env python3
import os
import sys
import subprocess
import re

if len(sys.argv) != 3:
    print("Usage: ./phylofit_model_selection.py <alignment_file> <tree_file>")
    sys.exit(1)

alignment_file = sys.argv[1]
tree_file = sys.argv[2]

models = ["JC69", "F81", "HKY85", "REV", "UNREST", "SSREV"]
results = []

output_dir = "phylofit_model_test"
os.makedirs(output_dir, exist_ok=True)

for model in models:
    print(f"Running phyloFit with model {model}...")
    out_root = os.path.join(output_dir, model + "_fit")

    cmd = [
        "phyloFit",
        "--msa", alignment_file,
        "--tree", tree_file,
        "--subst-mod", model,
        "--out-root", out_root,
        "--msa-format", "FASTA"
    ]

    # Capture stdout
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    output = proc.stdout

    # Parse stdout for numpar and log-likelihood
    lnL = None
    numpar = None
    for line in output.splitlines():
        if "numpar" in line:
            m = re.search(r"numpar\s*=\s*(\d+)", line)
            if m:
                numpar = int(m.group(1))
        if "log(likelihood)" in line:
            m = re.search(r"log\(likelihood\)\s*=\s*([-+]?\d*\.\d+|\d+)", line)
            if m:
                lnL = float(m.group(1))

    if lnL is not None and numpar is not None:
        AIC = 2 * numpar - 2 * lnL
        results.append((model, lnL, numpar, AIC))
    else:
        print(f"Warning: Could not extract lnL or numpar for {model}")

if results:
    print("\n=== Model Selection Results ===")
    for model, lnL, numpar, AIC in results:
        print(f"{model:<7} | lnL: {lnL:<12} | numpar: {numpar:<4} | AIC: {AIC:<12}")

    best_model = min(results, key=lambda x: x[3])
    print(f"\nBest model based on AIC: {best_model[0]} with AIC={best_model[3]}")
else:
    print("No models were successfully parsed.")
