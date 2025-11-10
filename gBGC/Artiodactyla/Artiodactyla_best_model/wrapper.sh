#!/bin/bash

# Usage: ./run_best_phylofit.sh <alignment_file> <tree_file> <out_root>

ALIGNMENT=$1
TREE=$2
OUT_ROOT=$3

if [ -z "$ALIGNMENT" ] || [ -z "$TREE" ] || [ -z "$OUT_ROOT" ]; then
    echo "Usage: $0 <alignment_file> <tree_file> <out_root>"
    exit 1
fi

# Temporary file to capture model selection output
MODEL_SELECTION_OUTPUT=$(mktemp)

# Step 1: Run the Python script to select the best model
python3 phylofit_model_selection.py "$ALIGNMENT" "$TREE" > "$MODEL_SELECTION_OUTPUT"

# Step 2: Extract the best model from the Python script output
BEST_MODEL=$(grep "Best model" "$MODEL_SELECTION_OUTPUT" | awk -F": " '{print $2}' | awk '{print $1}')

if [ -z "$BEST_MODEL" ]; then
    echo "Error: Could not determine the best model from model selection output."
    exit 1
fi

echo "Best model selected: $BEST_MODEL"

# Step 3: Run phyloFit with the best model
phyloFit --msa "$ALIGNMENT" \
         --tree "$TREE" \
         --subst-mod "$BEST_MODEL" \
         --out-root "$OUT_ROOT" \
         --msa-format FASTA

"echo phyloFit Done"

# Step 4: Run phastBias for each species in the alignment
# Extract species names from FASTA headers (lines starting with ">")
for SPECIES in $(grep "^>" "$ALIGNMENT" | sed 's/>//'); do
    OUT_GFF="${SPECIES}_bias.gff"
    OUT_POST="${SPECIES}_bias.post"
    echo "Running phastBias for $SPECIES..."
    phastBias --bgc 3 \
              --output-tracts "$OUT_GFF" \
              --posteriors full \
              "$ALIGNMENT" \
              "${OUT_ROOT}.mod" \
              "$SPECIES" > "$OUT_POST"
done

echo "phastBias Done!"
