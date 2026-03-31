#!/bin/bash
# bil-atlas-repro - Step 2: High-Resolution Species Delineation
# Principle: Use ANI to confirm B. infantis vs B. longum separation.

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
GENOME_DIR="$BASE_DIR/data/raw_genomes"
OUT_DIR="$BASE_DIR/results/01.taxonomy"

mkdir -p $OUT_DIR

echo ">>> Calculating All-vs-All ANI Matrix (skani)..."
# Logic: Fast comparison of 4,098^2 pairs to define species boundaries (95% threshold).
conda run -n env-ecogenomics skani triangle $GENOME_DIR/*.fa -o $OUT_DIR/all_vs_all_ani.txt --threads 32

echo ">>> ANI calculation complete."
