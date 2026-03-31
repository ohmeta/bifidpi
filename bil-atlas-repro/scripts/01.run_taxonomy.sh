#!/bin/bash
# bil-atlas-repro - Step 1: Quality Control & Taxonomy
# Principles from Shao et al. (2026) Cell

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
GENOME_DIR="$BASE_DIR/data/raw_genomes"
OUT_DIR="$BASE_DIR/results/01.taxonomy"

mkdir -p $OUT_DIR

echo ">>> Phase 1: High-throughput Taxonomy Assignment (GTDB-Tk)..."
# gtdbtk classify_wf --genome_dir $GENOME_DIR --out_dir $OUT_DIR/gtdbtk --cpus 32 -x fa

echo ">>> Phase 2: Quality Assessment (CheckM2)..."
# checkm2 predict --input $GENOME_DIR --output-directory $OUT_DIR/checkm2 --threads 32 -x fa
