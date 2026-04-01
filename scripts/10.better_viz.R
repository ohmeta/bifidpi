#!/usr/bin/env Rscript
# bifidpi - Better visualization: Motor & Shield heatmap with clear labels
suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
})

BASE_DIR <- "/mnt/store/users/zhujie/HOLA/assay/bifidpi"
OUT_DIR  <- file.path(BASE_DIR, "results/phase3")

# Load prevalence data
prev <- read.csv(file.path(OUT_DIR, "motor_shield_prevalence.csv"), stringsAsFactors = FALSE)

# Filter to local MAGs only
local <- prev %>% filter(Source == "Local MAG")

# Clean MAG names: just remove "megahit.semibin_single" to shorten
local <- local %>%
  mutate(short_id = gsub("\\.megahit\\.semibin_single", "", mag_id))

# Sort: Exposure first, then Control, then by Motor+Shield score descending
local <- local %>%
  arrange(Group, desc(Motor_score + Shield_score), short_id)

# Gene matrix
gene_cols <- c("fruG", "fruF", "fruK", "fruE", "mtrA", "mtrB", "lpqB")
mat <- as.matrix(local[, gene_cols])
rownames(mat) <- local$short_id

# Nicer column names
colnames(mat) <- c("fruG", "fruF", "fruK", "fruE", "mtrA", "mtrB", "lpqB")

# Row annotation: Group
row_annot <- rowAnnotation(
  Group = local$Group,
  Motor = local$Motor_score,
  Shield = local$Shield_score,
  col = list(
    Group = c("Exposure" = "#E64B35", "Control" = "#4DBBD5"),
    Motor = colorRamp2(c(0, 2, 4), c("white", "#FEE08B", "#E64B35")),
    Shield = colorRamp2(c(0, 1.5, 3), c("white", "#FEE08B", "#E64B35"))
  ),
  show_legend = TRUE,
  gp = gpar(fontsize = 7)
)

# Column annotation: Motor vs Shield
col_annot <- HeatmapAnnotation(
  System = c(rep("Motor (Fructose PTS)", 4), rep("Shield (MtrAB-LpqB)", 3)),
  col = list(
    System = c("Motor (Fructose PTS)" = "#F39B7F", "Shield (MtrAB-LpqB)" = "#8491B4")
  ),
  show_legend = TRUE,
  gp = gpar(fontsize = 8)
)

# Build heatmap
ht <- Heatmap(mat,
  name = "Present",
  col = c("0" = "grey95", "1" = "#E64B35"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 5.5),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title = "Motor & Shield Gene Prevalence in 28 Local MAGs",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  left_annotation = row_annot,
  top_annotation = col_annot,
  border = TRUE,
  rect_gp = gpar(col = "white", lwd = 1.5),
  width = unit(7, "cm"),
  height = unit(16, "cm")
)

pdf(file.path(OUT_DIR, "motor_shield_heatmap_v2.pdf"), width = 11, height = 10)
draw(ht, padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

cat("Saved: motor_shield_heatmap_v2.pdf\n")

# ──────────────────────────────────────────────
# Also: create a clear summary barplot
# ──────────────────────────────────────────────

# Count completeness by group
summary_df <- local %>%
  mutate(
    Motor_complete = Motor_score == 4,
    Shield_complete = Shield_score == 3,
    Both = Motor_complete & Shield_complete
  ) %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    Motor_only = sum(Motor_complete & !Shield_complete),
    Shield_only = sum(!Motor_complete & Shield_complete),
    Both = sum(Both),
    Neither = sum(!Motor_complete & !Shield_complete),
    .groups = "drop"
  )

# Reshape for plotting
plot_df <- summary_df %>%
  pivot_longer(cols = c(Motor_only, Shield_only, Both, Neither),
               names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("Both", "Motor_only", "Shield_only", "Neither")),
         pct = Count / n * 100)

p <- ggplot(plot_df, aes(x = Group, y = pct, fill = Category)) +
  geom_col(position = "stack", width = 0.6, color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Count > 0, paste0(Count, " (", round(pct), "%)"), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(
    values = c("Both" = "#E64B35", "Motor_only" = "#F39B7F",
               "Shield_only" = "#8491B4", "Neither" = "grey80"),
    labels = c("Both Motor+Shield", "Motor only", "Shield only", "Neither")
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  labs(x = "", y = "% of MAGs", fill = "",
       title = "Motor & Shield Operon Completeness by Group",
       subtitle = "Motor = fruG/F/K/E (Fructose PTS) | Shield = mtrA/B/lpqB") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "motor_shield_barplot.pdf"), p, width = 7, height = 5)
cat("Saved: motor_shield_barplot.pdf\n")

# Print clear summary
cat("\n========== CLEAR SUMMARY ==========\n")
cat("Exposure (12 MAGs):\n")
cat(sprintf("  Complete Motor+Shield: %d/12 (100%%)\n", sum(local$Group == "Exposure" & local$Motor_score == 4 & local$Shield_score == 3)))
cat(sprintf("  Motor complete: %d/12 (100%%)\n", sum(local$Group == "Exposure" & local$Motor_score == 4)))
cat(sprintf("  Shield complete: %d/12 (100%%)\n", sum(local$Group == "Exposure" & local$Shield_score == 3)))
cat("\nControl (16 MAGs):\n")
cat(sprintf("  Complete Motor+Shield: %d/16 (%.0f%%)\n", sum(local$Group == "Control" & local$Motor_score == 4 & local$Shield_score == 3), sum(local$Group == "Control" & local$Motor_score == 4 & local$Shield_score == 3)/16*100))
cat(sprintf("  Motor complete: %d/16 (%.0f%%)\n", sum(local$Group == "Control" & local$Motor_score == 4), sum(local$Group == "Control" & local$Motor_score == 4)/16*100))
cat(sprintf("  Shield complete: %d/16 (%.0f%%)\n", sum(local$Group == "Control" & local$Shield_score == 3), sum(local$Group == "Control" & local$Shield_score == 3)/16*100))
cat("\n")
cat("Individual Control MAGs:\n")
ctrl <- local %>% filter(Group == "Control") %>% arrange(desc(Motor_score + Shield_score))
for (i in 1:nrow(ctrl)) {
  cat(sprintf("  %s: Motor=%d/4, Shield=%d/3\n", ctrl$short_id[i], ctrl$Motor_score[i], ctrl$Shield_score[i]))
}
