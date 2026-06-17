library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(rCGH)

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
  )

head(clean_genomic_data)

#=========================================
# 2. Generating a location plot
#=========================================

# a. Filter for chromosomes 3 and 5
chr3_5_data <- clean_genomic_data %>%
  filter(chr %in% c("3", "5"))

# Define expected copy number ratio baselines (assuming log2 ratios)
# If your data represents raw ratios instead, change these to 0.5, 1, 1.5, and 2.0.
val_monosomic  <- log2(1/2)  # -1.00
val_disomic    <- log2(2/2)  #  0.00
val_trisomic   <- log2(3/2)  #  0.58
val_tetrasomic <- log2(4/2)  #  1.00

# b. Generate the faceted point plot
location_plot <- ggplot(chr3_5_data, aes(x = start, y = value)) +

  # Point plot for genomic positions
  geom_point(alpha = 0.4, size = 0.5, color = "gray30") +

  # Add expected copy number lines with clear visual distinctions
  geom_hline(yintercept = val_monosomic, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = val_disomic, color = "black", linetype = "solid", linewidth = 0.8) +
  geom_hline(yintercept = val_trisomic, color = "red", linetype = "dashed") +
  geom_hline(yintercept = val_tetrasomic, color = "purple", linetype = "dashed") +

  # Facet rows by cell line, columns by chromosome
  # scales = "free_x" is required so each chromosome scales to its own physical length
  facet_grid(cell_line ~ chr, scales = "free_x") +

  labs(
    title = "Copy Number Ratios across Chromosomes 3 and 5",
    x = "Genomic Position (bp)",
    y = "Log2 Copy Number Ratio",
    caption = "Lines: Blue (Monosomic), Black (Disomic), Red (Trisomic), Purple (Tetrasomic)"
  ) +
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold")
  )

# View the plot
print(location_plot)

# Export the plot
ggsave("./results/chr3_5_location_plot.pdf", plot = location_plot, width = 12, height = 6, device = "pdf")