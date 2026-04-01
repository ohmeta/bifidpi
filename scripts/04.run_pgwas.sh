#!/bin/bash
# bifidpi - Step 4: Pangenome GWAS (Scoary)
# Requirements: scoary

BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
PA_FILE="$BASE_DIR/results/pangenome/gene_presence_absence.csv"
TRAIT_FILE="$BASE_DIR/data/metadata_traits.csv"
OUT_DIR="$BASE_DIR/results/pgwas"
FILTERED_PA="$BASE_DIR/results/pgwas/gene_presence_absence_filtered.csv"

mkdir -p $OUT_DIR

echo ">>> Preparing metadata for Scoary..."
# Scoary needs a simple CSV with: Name,Trait1,Trait2
# We extract this from your detailed analysis table
METADATA="/mnt/store/users/zhujie/HOLA/assay/results/tables/b_longum_detailed_analysis.tsv"

# Create Scoary trait file: Column 1 is ID, Column 2 is Group (Binary 0/1)
echo "Name,is_exposure" > $TRAIT_FILE
tail -n +2 $METADATA | awk -F'\t' '{
    n=split($1,a,"/");
    id=a[n];
    gsub(".fa.gz$","",id);
    group=($11=="Exposure"?1:0);
    print id","group
}' | sort -u >> $TRAIT_FILE

# Filter PA file to only include genomes with metadata
echo ">>> Filtering PA file to match trait file..."
# Get trait names (first column of trait file, skip header)
tail -n +2 $TRAIT_FILE | cut -d',' -f1 > /tmp/valid_genomes.txt
# Use awk to filter: determine valid column indices from header, then keep only those in data rows
awk -F',' -v OFS=',' -v traitfile="/tmp/valid_genomes.txt" '
BEGIN {
    while ((getline line < traitfile) > 0) {
        valid[line] = 1
    }
}
NR==1 {
    # Always keep columns 1-3 (Gene, Non-unique Gene name, Annotation)
    keep[1]=1; keep[2]=1; keep[3]=1
    # Determine which sample columns to keep
    for (i=4; i<=NF; i++) {
        if ($i in valid) {
            keep[i] = 1
        }
    }
    # Print header
    first = 1
    for (i=1; i<=NF; i++) {
        if (i in keep) {
            if (!first) printf ","; printf "%s", $i; first=0
        }
    }
    printf "\n"
    next
}
{
    first = 1
    for (i=1; i<=NF; i++) {
        if (i in keep) {
            if (!first) printf ","; printf "%s", $i; first=0
        }
    }
    printf "\n"
}' $PA_FILE > $FILTERED_PA

echo ">>> Running Scoary P-GWAS..."
# -g: pangenome matrix, -t: traits, -c: correction method
source /home/zhujie/.conda/envs/env-base/etc/profile.d/conda.sh && conda activate env-blongpi && scoary -g $FILTERED_PA -t $TRAIT_FILE -o $OUT_DIR --p_value_cutoff 0.05 --correction BH -s 4

echo ">>> P-GWAS complete. Check results in $OUT_DIR"
