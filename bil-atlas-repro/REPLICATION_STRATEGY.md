# Replication Strategy: Bifidobacterium Global Atlas

This folder contains the framework to replicate the analysis logic from **Shao et al. (2026) Cell**.

## Analysis Logic (The "Atlas" Workflow)

1.  **Species Delineation**:
    *   The paper argues *B. infantis* and *B. longum* are distinct species.
    *   Run `scripts/02.run_ani_matrix.sh` to build the ANI matrix and confirm the >95% clusters.

2.  **Biogeographic Stratification**:
    *   The authors identified 36 clades tied to geography.
    *   Run the Pangenome (`scripts/03.run_pangenome.sh`) and correlate the accessory gene patterns with the geographic labels in Supplementary Table S1.

3.  **Specialist vs. Generalist**:
    *   *B. infantis* is a "Hyper-specialist" (HMO focus).
    *   *B. longum* is a "Generalist" (Plant glycan focus).
    *   Run `scripts/04.metabolic_mapping.sh` to compare the H1/H2/H3 loci (Infantis) vs. diverse Carbohydrate-Active Enzymes (Longum).

4.  **Probiotic Design Logic**:
    *   Safety first: Run `scripts/05.probiotic_screening.sh` for AMR/Virulence.
    *   Efficacy: Prioritize clades that utilize both Urea and HMOs (LMIC-optimized) or focus on industrial-diet markers (HIC-optimized).

## How to use this framework
1.  **Reference Genomes**: The 4,098 genomes are already symlinked in `data/raw_genomes/`.
2.  **Metadata**: Download Table S1 from the Cell paper and place it in `data/metadata/`.
3.  **Scripts**: The scripts in `scripts/` provide the technical "skeletons" for each major step.
