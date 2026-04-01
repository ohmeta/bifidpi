#!/bin/bash
# blongpi - Step 7: ANI Screening (skani)
# Requirements: skani (found in env-ecogenomics)

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
MAG_DIR="$BASE_DIR/data/mags"
BIL_DIR="/mnt/store/omics/SDB/bilgenomes/BIL_genome"
OUT_DIR="$BASE_DIR/results/ani"

mkdir -p $OUT_DIR

echo ">>> Running skani to screen 28 MAGs against 4,098 global genomes..."

# Using 'dist' instead of 'triangle' because we are comparing two sets (Search mode)
# -q: query (your MAGs), -r: reference (the 4,098 genomes)
conda run -n env-ecogenomics skani dist -q "$MAG_DIR"/*.fa* -r "$BIL_DIR"/*.fa* -o "$OUT_DIR/mag_vs_bil_ani.txt" -t 32

echo ">>> Screening complete. Results in $OUT_DIR/mag_vs_bil_ani.txt"
