#!/bin/bash
# bil-atlas-repro - Step 3: Global Pangenomics
# Principle: Define the core and accessory gene pool of the species.

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
ANN_DIR="$BASE_DIR/../../results/annotations" # Reusing annotations from Stage 1 if available
OUT_DIR="$BASE_DIR/results/02.pangenome"

mkdir -p $OUT_DIR

echo ">>> Running Global Pangenome (Panaroo)..."
# Logic: Construct the graph for 4,098 genomes. 
# Requires GFF files from all genomes.
# panaroo -i $ANN_DIR/*/*.gff -o $OUT_DIR --clean-mode strict -a core --threads 32
