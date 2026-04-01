#!/bin/bash
# blongpi - Step 1: Foundation
# Requirements: bakta, panaroo, iqtree

# Get the absolute path of the blongpi directory
BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
MAG_DIR="$BASE_DIR/data/mags"
RES_DIR="$BASE_DIR/results"

mkdir -p $RES_DIR/annotations $RES_DIR/pangenome $RES_DIR/phylogeny

# 1. Annotation
echo ">>> Running Bakta for local MAGs..."
for f in $MAG_DIR/*.fa; do
    id=$(basename $f .fa)
    conda run -n env-blongpi prokka --outdir $RES_DIR/annotations/$id --prefix $id --genus Bifidobacterium --species longum --cpus 32 $f
done

# 1b. Optional: Annotate selected references from Cell Atlas
REF_LIST="$BASE_DIR/data/selected_references.txt"
if [ -f "$REF_LIST" ]; then
    echo ">>> Found selected references. Annotating..."
    while read -r ref_path; do
        id=$(basename "$ref_path" .fa)
        if [ ! -d "$RES_DIR/annotations/$id" ]; then
            echo "Annotating reference: $id"
            conda run -n env-blongpi prokka --outdir $RES_DIR/annotations/$id --prefix $id --genus Bifidobacterium --species longum --cpus 32 "$ref_path"
        fi
    done < "$REF_LIST"
fi

# 2. Pangenome (Panaroo)
echo ">>> Running Panaroo..."
# Collect all GFF files (including local and reference)
GFFS=$(ls $RES_DIR/annotations/*/*.gff)
# Suppress Biopython warnings and adjust core threshold for MAGs
export PYTHONWARNINGS="ignore::BiopythonDeprecationWarning"
conda run -n env-blongpi panaroo -i $GFFS -o $RES_DIR/pangenome --clean-mode strict -a core --core_threshold 0.5 --threads 32

# 3. Phylogeny
echo ">>> Running IQ-TREE..."
conda run -n env-blongpi iqtree -s $RES_DIR/pangenome/core_gene_alignment.aln -m GTR+G -nt AUTO -bb 1000 -redo -pre $RES_DIR/phylogeny/blongum_core
