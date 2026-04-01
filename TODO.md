# bifidpi - Progress & Next Steps

## Completed (Stage 1: Foundation)
- [x] `scripts/00.prepare_hmms.sh` - HMM database prepared (10 targeted Pfam markers)
- [x] `scripts/01.run_genomics.sh` - Bakta annotation + Panaroo pangenome + IQ-TREE phylogeny (local MAGs only)
- [x] `scripts/02.search_functions.sh` - HMMER functional search completed
- [x] `scripts/03.analysis_viz.R` - Statistical analysis + visualization completed
- [x] Downloaded 4,098 global *B. longum* genomes from Shao et al. 2026 (Cell)

---

## Phase 2: Contextualized Genomics (Stage 2) - COMPLETE (2026-04-01)
- [x] `scripts/07.run_skani.sh` - ANI screening completed (114,744 comparisons)
- [x] `scripts/06.select_references.R` - Selected 40 top ANI matches (>99%)
- [x] Re-run `scripts/01.run_genomics.sh` with integrated local + reference genomes (68 genomes)
  - Total genes: 4,968 clusters
  - Shell genes (15-95%): 2,098
  - Soft core (95-99%): 6
- [x] IQ-TREE phylogeny completed with 68 genomes
- [x] `scripts/04.run_pgwas.sh` - P-GWAS completed (no BH-significant genes, n=28 underpowered)
- [x] `scripts/05.check_operons.py` - MtrAB-LpqB intact (1.0), FructosePTS intact (1.0), TadPili partial
- [x] `scripts/08.run_dbcan.sh` - CAZyme profiling completed (57 genomes, see REPORT.md)

---

## Phase 3: Downstream Analysis - COMPLETE (2026-04-01)
- [x] `scripts/09.phase3_analysis.R` - Motor & Shield prevalence + CAZyme integration
- [x] `scripts/10.better_viz.R` - Publication-ready heatmap + barplot
- [x] Motor & Shield enrichment: 100% Exposure vs 19% Control (p<0.0001)
- [x] Updated REPORT.md with Phase 3 findings

## Phase 4: Advanced Analyses (Future)
- [ ] Mobile Genetic Element (MGE) detection (VIBRANT/IslandPath)
- [ ] Synteny/operon conservation in 40 global reference genomes
- [ ] Population structure analysis (PCA on gene presence/absence)
- [ ] Genome-scale metabolic modeling (GEMs)
- [ ] Expand cohort size for statistical power
