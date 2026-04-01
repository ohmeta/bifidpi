#!/usr/bin/env Rscript
# bifidpi - Phase 3: Integrated Analysis
# - Phylogenetic tree with Exposure/Control + Reference annotations
# - Motor & Shield gene prevalence across groups
# - CAZyme heatmap for key families
# - Multi-site lineage analysis (inspired by Shao et al. 2026)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtree)
  library(ggtreeExtra)
  library(ComplexHeatmap)
  library(circlize)
})

BASE_DIR <- "/mnt/store/users/zhujie/HOLA/assay/bifidpi"
OUT_DIR  <- file.path(BASE_DIR, "results/phase3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat(">>> Phase 3 Analysis Starting...\n")

# ──────────────────────────────────────────────
# 1. Load data
# ──────────────────────────────────────────────

# 1a. Phylogenetic tree (68 genomes)
tree_file <- file.path(BASE_DIR, "results/phylogeny/blongum_core.treefile")
tree <- read.tree(tree_file)
cat(sprintf("  Tree: %d tips\n", length(tree$tip.label)))

# 1b. Metadata (group assignments for local MAGs)
meta_file <- file.path(BASE_DIR, "data/metadata_traits.csv")
meta <- read.csv(meta_file, stringsAsFactors = FALSE) %>%
  rename(mag_id = Name, is_exposure = is_exposure) %>%
  mutate(Group = ifelse(is_exposure == 1, "Exposure", "Control"))

# 1c. Identify genome types
all_tips <- data.frame(mag_id = tree$tip.label, stringsAsFactors = FALSE) %>%
  left_join(meta, by = "mag_id") %>%
  mutate(
    Group = ifelse(is.na(Group), "Reference", Group),
    Source = case_when(
      grepl("^HOLA\\.", mag_id) ~ "Local MAG",
      grepl("^MGYG", mag_id) ~ "Global Ref (Shao)",
      grepl("^GCA_", mag_id) ~ "Global Ref (NCBI)",
      grepl("^[0-9]+_", mag_id) ~ "Other MAG",
      TRUE ~ "Unknown"
    )
  )

cat(sprintf("  Exposure: %d, Control: %d, Reference: %d\n",
    sum(all_tips$Group == "Exposure"),
    sum(all_tips$Group == "Control"),
    sum(all_tips$Group == "Reference")))

# 1d. Gene presence/absence matrix
pa_file <- file.path(BASE_DIR, "results/pangenome/gene_presence_absence.csv")
pa <- read.csv(pa_file, stringsAsFactors = FALSE, check.names = FALSE)

# 1e. CAZyme summaries (key families per genome)
cazyme_dir <- file.path(BASE_DIR, "results/dbcan")
cazyme_summary <- list()
target_families <- c("GH29", "GH95", "GH33", "GH13", "GH3", "GH43", "GT2", "CE9", "GH25", "GH51")

for (d in list.dirs(cazyme_dir, recursive = FALSE)) {
  id <- basename(d)
  if (id == "test") next
  ov_file <- file.path(d, "overview.txt")
  if (!file.exists(ov_file)) next
  ov <- tryCatch(read.delim(ov_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(ov)) next
  
  # Count each target family
  counts <- sapply(target_families, function(fam) {
    sum(grepl(fam, ov$HMMER, fixed = TRUE))
  })
  cazyme_summary[[id]] <- counts
}

cazyme_mat <- do.call(rbind, cazyme_summary)
cazyme_df <- as.data.frame(cazyme_mat)
cazyme_df$mag_id <- rownames(cazyme_df)
rownames(cazyme_df) <- NULL

cat(sprintf("  CAZyme data: %d genomes × %d families\n", nrow(cazyme_df), ncol(cazyme_mat)))

# ──────────────────────────────────────────────
# 2. Motor & Shield gene prevalence
# ──────────────────────────────────────────────

cat(">>> Checking Motor & Shield gene prevalence...\n")

# Genes of interest
motor_genes <- c("fruG", "fruF", "fruK", "fruE")
shield_genes <- c("mtrA", "mtrB", "lpqB")
adhesion_genes <- c("tadA")

# Extract presence/absence for these genes
extract_gene_pa <- function(pa, gene_names) {
  result <- data.frame(mag_id = character(), gene = character(), present = integer(),
                       stringsAsFactors = FALSE)
  for (g in gene_names) {
    # Find rows matching this gene (exact match on Gene column or annotation)
    idx <- which(pa$Gene == g | grepl(paste0("^", g, "$"), pa$Gene))
    if (length(idx) == 0) {
      # Try partial match
      idx <- which(grepl(g, pa$Gene, ignore.case = TRUE))
    }
    if (length(idx) == 0) next
    
    for (i in idx) {
      gene_row <- as.character(pa[i, -(1:3)])
      present <- !is.na(gene_row) & gene_row != "" & gene_row != "0"
      mags <- colnames(pa)[-(1:3)]
      result <- rbind(result, data.frame(
        mag_id = mags,
        gene = g,
        present = as.integer(present),
        stringsAsFactors = FALSE
      ))
    }
  }
  # Keep only first match per gene (deduplicate)
  result <- result %>%
    group_by(mag_id, gene) %>%
    summarise(present = max(present), .groups = "drop")
  return(result)
}

motor_pa <- extract_gene_pa(pa, motor_genes)
shield_pa <- extract_gene_pa(pa, shield_genes)

# Combine into prevalence matrix
all_gene_pa <- rbind(motor_pa, shield_pa)
gene_wide <- all_gene_pa %>%
  pivot_wider(names_from = gene, values_from = present, values_fill = 0)

# Join with metadata
gene_wide <- gene_wide %>%
  left_join(all_tips %>% select(mag_id, Group, Source), by = "mag_id")

# Motor score: how many fru genes present (0-4)
# Shield score: how many mtr/lpq genes present (0-3)
gene_wide <- gene_wide %>%
  mutate(
    Motor_score = rowSums(across(all_of(motor_genes))),
    Shield_score = rowSums(across(all_of(shield_genes)))
  )

# Summary by group
cat("\n  Motor & Shield prevalence by group:\n")
group_summary <- gene_wide %>%
  filter(Group %in% c("Exposure", "Control")) %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    Motor_complete = sum(Motor_score == 4),
    Shield_complete = sum(Shield_score == 3),
    Both_complete = sum(Motor_score == 4 & Shield_score == 3),
    Motor_mean = mean(Motor_score),
    Shield_mean = mean(Shield_score),
    .groups = "drop"
  )
print(as.data.frame(group_summary))

# Fisher's test for enrichment
if (nrow(group_summary) == 2) {
  # Motor complete
  motor_tab <- matrix(c(
    group_summary$Motor_complete[1], group_summary$n[1] - group_summary$Motor_complete[1],
    group_summary$Motor_complete[2], group_summary$n[2] - group_summary$Motor_complete[2]
  ), nrow = 2)
  motor_ft <- fisher.test(motor_tab)
  cat(sprintf("\n  Motor complete enrichment: OR=%.2f, p=%.4f\n", motor_ft$estimate, motor_ft$p.value))
  
  # Shield complete
  shield_tab <- matrix(c(
    group_summary$Shield_complete[1], group_summary$n[1] - group_summary$Shield_complete[1],
    group_summary$Shield_complete[2], group_summary$n[2] - group_summary$Shield_complete[2]
  ), nrow = 2)
  shield_ft <- fisher.test(shield_tab)
  cat(sprintf("  Shield complete enrichment: OR=%.2f, p=%.4f\n", shield_ft$estimate, shield_ft$p.value))
  
  # Both complete
  both_tab <- matrix(c(
    group_summary$Both_complete[1], group_summary$n[1] - group_summary$Both_complete[1],
    group_summary$Both_complete[2], group_summary$n[2] - group_summary$Both_complete[2]
  ), nrow = 2)
  both_ft <- fisher.test(both_tab)
  cat(sprintf("  Both complete enrichment: OR=%.2f, p=%.4f\n", both_ft$estimate, both_ft$p.value))
}

# Save prevalence table
write.csv(gene_wide, file.path(OUT_DIR, "motor_shield_prevalence.csv"), row.names = FALSE)

# ──────────────────────────────────────────────
# 3. Visualization: Tree + Annotations + CAZymes
# ──────────────────────────────────────────────

cat(">>> Generating integrated visualization...\n")

# Color palette
group_colors <- c("Exposure" = "#E64B35", "Control" = "#4DBBD5", "Reference" = "#999999")
source_colors <- c("Local MAG" = "#E64B35", "Global Ref (Shao)" = "#00A087",
                    "Global Ref (NCBI)" = "#3C5488", "Other MAG" = "#F39B7F")

# 3a. Build tree with group annotation
p1 <- ggtree(tree, layout = "rectangular", size = 0.3) %<+% all_tips +
  geom_tippoint(aes(color = Group), size = 1.5, alpha = 0.8) +
  scale_color_manual(values = group_colors) +
  theme(legend.position = "bottom") +
  labs(title = "68-genome B. longum phylogeny",
       subtitle = "Core gene ML tree (IQ-TREE, 1550 genes)")

# 3b. Add Motor/Shield heatmap strip
gene_tip_data <- gene_wide %>%
  filter(mag_id %in% tree$tip.label) %>%
  select(mag_id, all_of(c(motor_genes, shield_genes))) %>%
  column_to_rownames("mag_id")

# Ensure order matches tree tips
gene_tip_data <- gene_tip_data[tree$tip.label[tree$tip.label %in% rownames(gene_tip_data)], , drop = FALSE]

# Create heatmap annotation
if (nrow(gene_tip_data) > 0) {
  # Convert to matrix for heatmap
  hm_mat <- as.matrix(gene_tip_data)
  colnames(hm_mat) <- gsub("fru", "fru ", colnames(hm_mat))
  colnames(hm_mat) <- gsub("mtr", "mtr ", colnames(hm_mat))
  colnames(hm_mat) <- gsub("lpq", "lpq ", colnames(hm_mat))
  
  ht <- Heatmap(hm_mat,
    name = "Present",
    col = c("0" = "white", "1" = "#E64B35"),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    column_title = "Motor & Shield Genes",
    column_title_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 7),
    border = TRUE,
    width = unit(4, "cm")
  )
} else {
  ht <- NULL
}

# 3c. CAZyme heatmap for local MAGs only
local_mags <- all_tips$mag_id[all_tips$Group %in% c("Exposure", "Control")]
cazyme_local <- cazyme_df %>%
  filter(mag_id %in% local_mags) %>%
  column_to_rownames("mag_id")

# Order by group
cazyme_local <- cazyme_local[order(
  all_tips$Group[match(rownames(cazyme_local), all_tips$mag_id)]
), , drop = FALSE]

if (nrow(cazyme_local) > 0) {
  group_annotation <- rowAnnotation(
    Group = all_tips$Group[match(rownames(cazyme_local), all_tips$mag_id)],
    col = list(Group = group_colors),
    show_legend = TRUE
  )
  
  ht_cazyme <- Heatmap(as.matrix(cazyme_local),
    name = "Count",
    col = colorRamp2(c(0, 5, 10, 20), c("white", "#FEE08B", "#FDAE61", "#E64B35")),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_gp = gpar(fontsize = 5),
    column_title = "CAZyme Families (Local MAGs)",
    column_title_gp = gpar(fontsize = 9),
    column_names_gp = gpar(fontsize = 8),
    left_annotation = group_annotation,
    border = TRUE
  )
} else {
  ht_cazyme <- NULL
}

# ──────────────────────────────────────────────
# 4. Save outputs
# ──────────────────────────────────────────────

cat(">>> Saving outputs...\n")

# 4a. Tree plot
pdf(file.path(OUT_DIR, "phylogeny_annotated.pdf"), width = 12, height = 16)
print(p1)
dev.off()

# 4b. Motor/Shield heatmap
if (!is.null(ht)) {
  pdf(file.path(OUT_DIR, "motor_shield_heatmap.pdf"), width = 8, height = 16)
  draw(ht)
  dev.off()
}

# 4c. CAZyme heatmap
if (!is.null(ht_cazyme)) {
  pdf(file.path(OUT_DIR, "cazyme_heatmap_local.pdf"), width = 8, height = 12)
  draw(ht_cazyme)
  dev.off()
}

# 4d. Summary statistics
write.csv(group_summary, file.path(OUT_DIR, "group_summary.csv"), row.names = FALSE)

# 4e. CAZyme summary table
write.csv(cazyme_df, file.path(OUT_DIR, "cazyme_all_genomes.csv"), row.names = FALSE)

cat(sprintf("\n>>> Phase 3 complete! Outputs in %s\n", OUT_DIR))
cat("    - phylogeny_annotated.pdf: Tree with group colors\n")
cat("    - motor_shield_heatmap.pdf: Gene presence heatmap on tree\n")
cat("    - cazyme_heatmap_local.pdf: CAZyme profiles by group\n")
cat("    - motor_shield_prevalence.csv: Per-genome gene scores\n")
cat("    - group_summary.csv: Enrichment statistics\n")
cat("    - cazyme_all_genomes.csv: Full CAZyme matrix\n")
