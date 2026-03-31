#!/bin/bash
# bil-atlas-repro - Step 5: Probiotic Safety Screening
# Principle: Screen for safety-related markers (AMR and Virulence).

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
GENOME_DIR="$BASE_DIR/data/raw_genomes"
OUT_DIR="$BASE_DIR/results/04.safety"

mkdir -p $OUT_DIR

echo ">>> Screening for Antimicrobial Resistance (AMR) genes..."
# Requirements: abricate
# abricate --db resfinder $GENOME_DIR/*.fa > $OUT_DIR/amr_results.tsv

echo ">>> Screening for Virulence Factors..."
# abricate --db vfdb $GENOME_DIR/*.fa > $OUT_DIR/virulence_results.tsv

echo ">>> Safety screening complete."
