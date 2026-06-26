library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rCGH)
library(DNAcopy)

#=========================================
# 1. Reading and tidying up the data
#=========================================
raw_data <- read_excel("data/aneuploidy_copy_number.xlsx")

tidy_data <- raw_data %>%
  pivot_longer(
    cols = c("H", "Htr5", "Hte5"),
    names_to = "cell_line",
    values_to = "value",
    values_drop_na = TRUE
  )

clean_genomic_data <- tidy_data %>%
  mutate(
    cell_line = factor(cell_line, levels = c("H", "Htr5", "Hte5")),
    chr = factor(chr, levels = c(as.character(1:22), "X", "Y"))
  ) %>%
  filter(!is.na(chr))

#=========================================
# Global Plot Variables
#=========================================
val_monosomic  <- log2(1/2)  # -1.00
val_disomic    <- log2(2/2)  #  0.00
val_trisomic   <- log2(3/2)  #  0.58
val_tetrasomic <- log2(4/2)  #  1.00

data(hg19)

centromere_all <- hg19 %>%
  rename(chr = chrom) %>%
  mutate(
    chr = factor(chr, levels = c(as.character(1:22), "X", "Y")),
    centromere_mid = (centromerStart + centromerEnd) / 2
  ) %>%
  filter(!is.na(chr))

#=========================================
# 2 & 3. Location plot with centromeres (Chr 3 & 5)
#=========================================
chr3_5_data <- clean_genomic_data %>%
  filter(chr %in% c("3", "5"))

centromere_chr3_5 <- centromere_all %>%
  filter(chr %in% c("3", "5"))

location_plot <- ggplot(chr3_5_data, aes(x = start, y = value)) +
  geom_point(alpha = 0.8, size = 0.5, color = "gray30") +

  geom_hline(yintercept = val_monosomic, color = "magenta", linetype = "dashed") +
  geom_hline(yintercept = val_disomic, color = "black", linetype = "solid", linewidth = 0.8) +
  geom_hline(yintercept = val_trisomic, color = "red", linetype = "dashed") +
  geom_hline(yintercept = val_tetrasomic, color = "green", linetype = "dashed") +

  geom_vline(
    data = centromere_chr3_5,
    aes(xintercept = centromere_mid),
    color = "gray50",
    linetype = "dotted",
    linewidth = 0.8
  ) +

  facet_grid(cell_line ~ chr, scales = "free_x") +
  labs(
    title = "Copy Number Ratios across Chromosomes 3 and 5",
    x = "Genomic Position (bp)",
    y = "Log2 Copy Number Ratio"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    text = element_text(face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave("./results/chr3_5_location_plot_centromeres.pdf", plot = location_plot, width = 12, height = 6, device = "pdf")

#=========================================
# 4. Visualize copy number segments (Chr 3 & 5)
#=========================================
cell_lines = colnames(raw_data)[1:3]

CNA.object = CNA(genomdat = raw_data[,cell_lines],
                 sampleid = cell_lines,
                 chrom = raw_data$chr,
                 maploc = raw_data$start,
                 data.type = "logratio")

CNA.object.smooth = smooth.CNA(CNA.object)
segment.CNA = segment(CNA.object.smooth)

segment_data <- segment.CNA$output %>%
  rename(cell_line = ID, chr = chrom) %>%
  mutate(
    cell_line = factor(cell_line, levels = c("H", "Htr5", "Hte5")),
    chr = factor(chr, levels = c(as.character(1:22), "X", "Y"))
  ) %>%
  filter(chr %in% c("3", "5"))

segmented_plot <- location_plot +
  geom_segment(
    data = segment_data,
    aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean),
    color = "darkorange",
    linewidth = 1.2,
    inherit.aes = FALSE
  )

ggsave("./results/chr3_5_segmented_plot.pdf", plot = segmented_plot, width = 12, height = 6, device = "pdf")

#=========================================
# 5. Visualize Segmented Karyotypes
#=========================================
segment_all <- segment.CNA$output %>%
  rename(cell_line = ID, chr = chrom) %>%
  mutate(
    cell_line = factor(cell_line, levels = c("H", "Htr5", "Hte5")),
    chr = factor(chr, levels = c(as.character(1:22), "X", "Y"))
  ) %>%
  filter(!is.na(chr))

karyotype_plot <- ggplot(clean_genomic_data, aes(x = start, y = value)) +
  geom_point(alpha = 0.3, size = 0.1, color = "gray30") +

  geom_hline(yintercept = val_monosomic, color = "magenta", linetype = "dashed") +
  geom_hline(yintercept = val_disomic, color = "black", linetype = "solid", linewidth = 0.8) +
  geom_hline(yintercept = val_trisomic, color = "red", linetype = "dashed") +
  geom_hline(yintercept = val_tetrasomic, color = "green", linetype = "dashed") +

  geom_vline(
    data = centromere_all,
    aes(xintercept = centromere_mid),
    color = "gray50",
    linetype = "dotted",
    linewidth = 0.5
  ) +

  geom_segment(
    data = segment_all,
    aes(x = loc.start, xend = loc.end, y = seg.mean, yend = seg.mean),
    color = "darkorange",
    linewidth = 1,
    inherit.aes = FALSE
  ) +

  facet_grid(cell_line ~ chr, scales = "free_x", space = "free_x") +
  labs(
    title = "Whole Genome Copy Number Profile",
    x = "Chromosome",
    y = "Log2 Copy Number Ratio"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    strip.background = element_rect(fill = "gray90", color = "black"),

    text = element_text(face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.1, "lines")
  )

ggsave("./results/whole_genome_karyotype.pdf", plot = karyotype_plot, width = 20, height = 8, device = "pdf")
