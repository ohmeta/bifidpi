#!/bin/bash
# blongpi - Step 8: CAZyme Profiling (dbCAN)
# Requirements: run_dbcan (found in env-blongpi)

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
ANN_DIR="$BASE_DIR/results/annotations"
OUT_DIR="$BASE_DIR/results/dbcan"

mkdir -p $OUT_DIR

echo ">>> Running dbCAN for Carbohydrate-Active enZymes (CAZyme) annotation..."

# Loop through all FAA (protein) files in annotation results
for faa in $ANN_DIR/*/*.faa; do
    id=$(basename $faa .faa)
    echo "Processing $id..."
    
    # Create sample-specific output directory
    sample_out="$OUT_DIR/$id"
    mkdir -p $sample_out
    
    # Run dbCAN
    # -t: tools (hmmer and diamond are standard)
    # --db_dir: Assuming dbCAN database is pre-configured or in standard location
    run_dbcan $faa protein --out_dir $sample_out --threads 32 --tools hmmer diamond
done

echo ">>> dbCAN complete. Results in $OUT_DIR"
