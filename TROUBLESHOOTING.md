# Technical Troubleshooting & Metrics Log: blongpi

This document records the specific technical challenges, bug fixes, and data metrics encountered during the analysis of 28 *B. longum* MAGs.

---

## 1. Pangenome Metrics & Thresholding (Step 1)

### 1.1 The "0 Core Gene" Paradox
- **Initial Result**: `summary_statistics.txt` reported **0 Core genes** and **0 Soft-core genes**.
- **The Metrics**:
    - **Total Genes**: 4,064 clusters.
    - **Shell Genes (15-95%)**: 2,104 clusters.
    - **Cloud Genes (0-15%)**: 1,960 clusters.
- **The Issue**: Panaroo defaults to a 95% threshold for tree building. Because MAGs were fragmented (one as small as 928KB), no genes were present in 27/28 genomes.
- **The Fix**: Redefined `--core_threshold` to `0.5`.
- **Tree Metrics**: 
    - **Genes Meeting 50% Threshold**: 1,550 clusters.
    - **Final Alignment**: `core_gene_alignment.aln` was built using these **1,550 genes**.
    - **Note**: The summary statistics still show 0 core genes because they use fixed reporting bins, but the tree-building logic correctly used the 1,550 genes.

### 1.2 Biopython Deprecation Warnings
- **Bug**: Thousands of warnings: `BiopythonDeprecationWarning: Previously, the FASTA parser silently ignored comments...`
- **Root Cause**: Modern Biopython versions (3.9+) are stricter about leading comments in FASTA files generated during the pangenome graph construction.
- **Fix**: Suppressed via environment variable in the shell script:
  ```bash
  export PYTHONWARNINGS="ignore::BiopythonDeprecationWarning"
  ```

---

## 2. Functional HMM & Data Collation (Step 2)

### 2.1 Pfam ID Verification
- **Strategy Mismatch**: The original document listed `PF01391` for TadA.
- **Verification**: `PF01391` is actually **Collagen**.
- **Correction Table**:
    - **TadA**: `PF00437` (T2SS_ATPase)
    - **TadE/F**: `PF07811`
    - **Flp**: `PF04964` (Flp_Fap)
    - **Sortase**: `PF04203`
    - **GH95**: `PF22124` (Corrected from KaiC `PF06745`)

### 2.2 Collation Logic Bug
- **Issue**: The original `scripts/02.search_functions.sh` used `awk` to grab hits but didn't include the **Genome ID**, making it impossible for R to link functions to specific infants.
- **Fix**: Modified the loop to extract the filename and prepend it to every hit in `combined_markers.tsv`.

---

## 3. R Visualization & Statistics (Step 3)

### 3.1 Metadata ID Alignment
- **Issue**: The metadata used `user_genome` with `.fa.gz` extensions, while the pangenome matrix used filenames with `.fa.fna`.
- **Fix**: Implemented a universal regex cleaner in R:
  ```r
  gsub(".fa.gz|.fa.fna|.fa|.fasta", "", basename(id))
  ```

### 3.2 ggtree Metadata Join
- **Bug**: The phylogenetic tree in the PDF had no color for Exposure/Control groups.
- **Root Cause**: The `ggtree` `%<+%` operator requires the ID column to be the **first column** in the metadata dataframe.
- **Fix**: Reordered columns using `select(mag_id_clean, everything())`.

### 3.3 ComplexHeatmap "Break Values" Error
- **Bug**: `Error: You should have at least two distinct break values.`
- **Root Cause**: The "Efficiency" metadata column had too many `NA` values or lacked variance among the 28 MAGs, causing the color scale generator to fail.
- **Fix**: Removed the continuous "Efficiency" annotation from the heatmap to ensure stable rendering of the primary "Group" categories.

---

## 4. Environment & Paths
- **Path Inconsistency**: The pipeline was moved from `/home/zhujie/projects/ohanalysis/HOLA` to `/mnt/store/users/zhujie/HOLA/assay/blongpi`.
- **Fix**: Updated all scripts with absolute paths or used the `BASE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"` logic for portability.
- **Conda**: Tools are executed within `env-blongpi`.
