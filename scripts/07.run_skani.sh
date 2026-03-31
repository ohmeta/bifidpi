#!/bin/bash
# blongpi - Step 7: ANI Clustering (skani)
# Requirements: skani (found in env-ecogenomics)

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
MAG_DIR="$BASE_DIR/data/mags"
OUT_DIR="$BASE_DIR/results/ani"

mkdir -p $OUT_DIR

echo ">>> Running skani for Average Nucleotide Identity (ANI) calculation..."

# Using conda run to execute from the correct environment
conda run -n env-ecogenomics skani triangle $MAG_DIR/*.fa* -o $OUT_DIR/skani_ani.txt --threads 32

echo ">>> skani complete. Results in $OUT_DIR/skani_ani.txt"
