library(tidyverse)

# Read the TSV file
res_data <- read_tsv("./enrichment.all.tsv")

# ==========================================
# Extract Top Enriched Pathways
# ==========================================

# Define thresholds
es_threshold <- 2.0
fdr_threshold <- 0.05

# Filter the dataset and extract the requested columns
# Note: Backticks (`) are required around these specific column names
top_enrichment <- res_data %>%
  filter(
    `enrichment score` > es_threshold,
    `false discovery rate` < fdr_threshold
  ) %>%
  select(
    `term description`,
    `matching proteins in your input (IDs)`
  )

# View the results
cat("Number of pathways meeting criteria:", nrow(top_enrichment), "\n\n")
print(top_enrichment)
