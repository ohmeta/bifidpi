#!/usr/bin/env Rscript
# blongpi - Step 6: Select best global references from skani output

library(tidyverse)
library(here)

base_dir <- here("blongpi")
ani_path <- file.path(base_dir, "results/ani/mag_vs_bil_ani.txt")
out_path <- file.path(base_dir, "data/selected_references.txt")

if (!file.exists(ani_path)) {
  stop("ANI results not found. Run Step 7 (skani) first.")
}

# Load skani results
# Columns: Query, Ref, ANI, Align_fraction_query, Align_fraction_ref, Ref_name, Query_name
ani_data <- read_tsv(ani_path)

echo ">>> Selecting top global references for each MAG..."

selected <- ani_data %>%
  # Filter for high quality matches
  filter(ANI > 99) %>%
  group_by(Query_Name) %>%
  # Pick top 2 most similar global genomes for each MAG
  slice_max(order_by = ANI, n = 2) %>%
  ungroup() %>%
  distinct(Ref_Name) %>%
  pull(Ref_Name)

# Save the full paths of the selected genomes
write_lines(selected, out_path)

cat(paste(">>> Selected", length(selected), "reference genomes from the Cell atlas.\n"))
cat(paste(">>> List saved to:", out_path, "\n"))
