# ----------------------------
# Libraries needed for analysis
# ----------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  library(ggrepel)
  library(ggvenn)
  library(readxl)
})

# ==============================================================================
#  SECTION 1: DATA LOADING & TIDYING
# ==============================================================================
# 1. Read the raw RNA differential expression data from the "Results" sheet
rna_raw <- read_excel("data/aneuploid_differential_expression.xlsx", sheet = "Results")

# 2. Tidy the results table using the standard pipeline
rna_tidy <- rna_raw |>
  pivot_longer(
    cols = c(starts_with("logFC"), starts_with("P.Value")),
    names_to = c(".value", "cell_line"),
    names_pattern = "(.*)_(.*)"
  ) |>
  mutate(
    cell_line = factor(cell_line, levels = c("Htr5", "Hte5")),
    chromosome_name = as.character(chromosome_name),
    # Flag genes encoded on Chromosome 5
    chr_group = if_else(chromosome_name == "5", "Chromosome 5", "Other Chromosomes")
  )

# Calculate FDR per cell line
rna_sig <- rna_tidy |>
  group_by(cell_line) |>
  mutate(
    FDR = p.adjust(P.Value, method = "BH"),
    Significance = case_when(
      FDR < 0.05 & logFC > 0 ~ "Upregulated",
      FDR < 0.05 & logFC < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  ) |>
  ungroup()

# ==============================================================================
#  SECTION 2: CHROMOSOME 5 DOSAGE ANALYSIS
# ==============================================================================

# Expected logFC values based on theoretical copy number additions
expected_fc <- data.frame(
  cell_line = c("Htr5", "Hte5"),
  # Trisomy 5: 3 copies vs 2 copies
  # Tetrasomy 5: 4 copies vs 2 copies
  yint = c(log2(3/2), log2(4/2))
)

# Plot 1: Expression of Chr 5 vs Other Chromosomes
plot_chr5_expr <- ggplot(rna_sig, aes(x = chr_group, y = logFC, fill = chr_group)) +
  geom_violin(alpha = 0.6, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(data = expected_fc, aes(yintercept = yint), color = "red", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ cell_line) +
  scale_fill_manual(values = c("Chromosome 5" = "#D95F02", "Other Chromosomes" = "gray70")) +
  labs(
    title = "Gene Expression: Chromosome 5 vs Rest of Genome",
    x = "Genomic Location",
    y = "Log2 Fold Change"
  ) +
  theme_classic() +
  theme(legend.position = "none", text = element_text(face = "bold"))

ggsave("./results/chr5_expression_violin.pdf", plot = plot_chr5_expr, width = 8, height = 6)

# Plot 2: DNA Copy Number vs RNA Expression (Chromosome 5 Only)

# Re-import and tidy the DNA copy number data
genomic_raw <- read_excel("data/aneuploidy_copy_number.xlsx")

clean_genomic_data <- genomic_raw |>
  pivot_longer(
    cols = c("H", "Htr5", "Hte5"),
    names_to = "cell_line",
    values_to = "value",
    values_drop_na = TRUE
  ) |>
  mutate(
    cell_line = factor(cell_line, levels = c("H", "Htr5", "Hte5")),
    chr = as.character(chr)
  )

# Extract DNA values for Chromosome 5
dna_chr5 <- clean_genomic_data |>
  filter(chr == "5", cell_line %in% c("Htr5", "Hte5")) |>
  select(cell_line, metric_value = value) |>
  mutate(Assay = "DNA (Copy Number)")

# Extract RNA values for Chromosome 5
rna_chr5 <- rna_sig |>
  filter(chromosome_name == "5") |>
  select(cell_line, metric_value = logFC) |>
  mutate(Assay = "RNA (Expression)")

combined_chr5 <- bind_rows(dna_chr5, rna_chr5)

# Generate Plot
plot_dna_rna <- ggplot(combined_chr5, aes(x = Assay, y = metric_value, fill = Assay)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8) +
  facet_wrap(~ cell_line) +
  stat_compare_means(method = "t.test", comparisons = list(c("DNA (Copy Number)", "RNA (Expression)"))) +
  scale_fill_manual(values = c("DNA (Copy Number)" = "#1B9E77", "RNA (Expression)" = "#7570B3")) +
  labs(
    title = "Dosage Compensation: DNA Ratio vs RNA Fold Change (Chr 5)",
    x = "Assay Type",
    y = "Log2 Ratio / Log2 Fold Change"
  ) +
  theme_classic() +
  theme(legend.position = "none", text = element_text(face = "bold"))

ggsave("./results/chr5_dna_vs_rna.pdf", plot = plot_dna_rna, width = 8, height = 6)
# ==============================================================================
#  SECTION 3: DIFFERENTIAL EXPRESSION VISUALIZATION
# ==============================================================================

# 1. Venn Diagrams
up_genes <- list(
  Htr5 = rna_sig |> filter(cell_line == "Htr5", Significance == "Upregulated") |> pull(hgnc_symbol) |> na.omit(),
  Hte5 = rna_sig |> filter(cell_line == "Hte5", Significance == "Upregulated") |> pull(hgnc_symbol) |> na.omit()
)

down_genes <- list(
  Htr5 = rna_sig |> filter(cell_line == "Htr5", Significance == "Downregulated") |> pull(hgnc_symbol) |> na.omit(),
  Hte5 = rna_sig |> filter(cell_line == "Hte5", Significance == "Downregulated") |> pull(hgnc_symbol) |> na.omit()
)

venn_up <- ggvenn(up_genes, fill_color = c("#1B9E77", "#D95F02"), set_name_size = 5) + ggtitle("Upregulated Genes (FDR < 0.05)")
venn_down <- ggvenn(down_genes, fill_color = c("#1B9E77", "#D95F02"), set_name_size = 5) + ggtitle("Downregulated Genes (FDR < 0.05)")

ggsave("./results/venn_upregulated.pdf", plot = venn_up, width = 6, height = 6)
ggsave("./results/venn_downregulated.pdf", plot = venn_down, width = 6, height = 6)

# Extract top common genes for Volcano labeling
common_up <- intersect(up_genes$Htr5, up_genes$Hte5)
common_down <- intersect(down_genes$Htr5, down_genes$Hte5)

# Identify the top 10 most significant common genes across both lines for labeling
top_common_genes <- rna_sig |>
  filter(hgnc_symbol %in% c(common_up, common_down)) |>
  group_by(hgnc_symbol) |>
  summarize(mean_p = mean(P.Value)) |>
  slice_min(mean_p, n = 10) |>
  pull(hgnc_symbol)

# 2. Volcano Plots
plot_volcano <- function(df, c_line, top_genes) {

  df_plot <- df |>
    filter(cell_line == c_line) |>
    mutate(
      log10_fdr = -log10(FDR),
      Highlight = case_when(
        Significance != "Not Significant" & chr_group == "Chromosome 5" ~ "Significant (Chr 5)",
        Significance != "Not Significant" & chr_group == "Other Chromosomes" ~ "Significant (Other)",
        TRUE ~ "Not Significant"
      ),
      Label = if_else(hgnc_symbol %in% top_genes, hgnc_symbol, NA_character_)
    )

  p <- ggplot(df_plot, aes(x = logFC, y = log10_fdr, color = Highlight, label = Label)) +
    geom_vline(xintercept = 0, color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "dimgrey") +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_label_repel(
      fontface = "bold",
      max.overlaps = 20,
      size = 3.5,
      show.legend = FALSE,
      color = "black",
      fill = "white"
    ) +
    scale_color_manual(values = c(
      "Significant (Chr 5)" = "#D95F02",
      "Significant (Other)" = "#1B9E77",
      "Not Significant" = "gray80"
    )) +
    labs(
      title = paste("Volcano Plot:", c_line),
      x = "Log2 Fold Change",
      y = "-Log10(FDR)"
    ) +
    theme_classic() +
    theme(text = element_text(face = "bold"), legend.position = "bottom")

  ggsave(paste0("./results/volcano_", c_line, ".pdf"), plot = p, width = 8, height = 7)
  return(p)
}

plot_volcano(rna_sig, "Htr5", top_common_genes)
plot_volcano(rna_sig, "Hte5", top_common_genes)

message("Differential Expression Pipeline Complete.")
