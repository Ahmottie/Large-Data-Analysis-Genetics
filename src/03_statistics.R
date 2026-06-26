# Load packages
library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)

df <- read_excel("data/aneuploid_micronuclei.xlsx")

df <- df %>%
  mutate(
    condition = gsub(" ", "", condition),
    perc_micronuclei = (micronuclei / nuclei) * 100,
    condition = factor(condition, levels = c("Htr5", "Hte5", "HCT116", "HCT116+APH"))
  )

df_summary <- df %>%
  group_by(condition) %>%
  summarise(
    N = n(),
    mean_perc = mean(perc_micronuclei),
    se_perc = sd(perc_micronuclei) / sqrt(N),
    .groups = "drop"
  )

print(df_summary)

comparisons <- list(
  c("HCT116+APH", "HCT116"),
  c("Htr5", "HCT116"),
  c("Hte5", "HCT116"),
  c("HCT116+APH", "Hte5")
)

custom_colors <- c("Htr5" = "#1B9E77", "Hte5" = "#7570B3", "HCT116" = "#D95F02", "HCT116+APH" = "#D95F02")
custom_lines <- c("Htr5" = "solid", "Hte5" = "solid", "HCT116" = "solid", "HCT116+APH" = "dashed")

#=========================================
# Bar Plot
#=========================================
bar_plot <- ggplot(df, aes(x = condition, y = perc_micronuclei)) +
  geom_col(
    data = df_summary,
    aes(x = condition, y = mean_perc, fill = condition, linetype = condition),
    inherit.aes = FALSE,
    color = "black",
    width = 0.5
  ) +
  geom_errorbar(
    data = df_summary,
    aes(x = condition, ymin = mean_perc - se_perc, ymax = mean_perc + se_perc),
    inherit.aes = FALSE,
    width = 0.2
  ) +
  stat_compare_means(
    comparisons = comparisons,
    method = "t.test",
    step.increase = 0.08
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_linetype_manual(values = custom_lines) +
  labs(
    title = "Average Percentage of Micronuclei per Condition",
    x = "Condition",
    y = "Average Micronuclei %"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17),
    legend.position = "none"
  )

print(bar_plot)
ggsave("./results/bar_plot.pdf", plot = bar_plot, width = 10, height = 5, device = "pdf")

#=========================================
# Box Plot
#=========================================
box_plot <- ggplot(df, aes(x = condition, y = perc_micronuclei)) +
  stat_boxplot(aes(linetype = condition), geom = "errorbar", width = 0.2, coef = Inf) +
  geom_boxplot(aes(fill = condition, linetype = condition), color = "black", coef = Inf, width = 0.4) +
  stat_compare_means(comparisons = comparisons, method = "t.test", step.increase = 0.08) +
  scale_fill_manual(values = custom_colors) +
  scale_linetype_manual(values = custom_lines) +
  labs(
    title = "Distribution of Micronuclei Percentages",
    x = "Cell Line Condition",
    y = "Micronuclei per Observation (%)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    legend.position = "none"
  )

print(box_plot)
ggsave("./results/box_plot.pdf", plot = box_plot, width = 7, height = 6, device = "pdf")

# Tetraploidy has potentially a larger effect on the formation of micronucleai
#=========================================
# Explicit T-Tests
#=========================================
test_a <- t.test(perc_micronuclei ~ condition, data = subset(df, condition %in% c("HCT116+APH", "HCT116")), var.equal = FALSE)
test_b <- t.test(perc_micronuclei ~ condition, data = subset(df, condition %in% c("Htr5", "HCT116")), var.equal = FALSE)
test_c <- t.test(perc_micronuclei ~ condition, data = subset(df, condition %in% c("Hte5", "HCT116")), var.equal = FALSE)
test_d <- t.test(perc_micronuclei ~ condition, data = subset(df, condition %in% c("HCT116+APH", "Hte5")), var.equal = FALSE)

print("HCT116+APH vs HCT116 (Positive vs Negative Control)")
print(test_a)
print("Htr5 vs HCT116 (Aneuploid Model 1 vs Negative Control)")
print(test_b)
print("Hte5 vs HCT116 (Aneuploid Model 2 vs Negative Control)")
print(test_c)
print("HCT116+APH vs Hte5 (Positive Control vs Aneuploid Model 2)")
print(test_d)

# Hte5 has significant increase in micronuclei. We can interperate that tetraploidi
# has a more drastic effect on  the formation of micronucli in comparison to triploidi
# Hte5 has a significant difference campred to negative
# control by has the same approximate percentage of micronuclei
# The Htr5 condition does not significantly increase the percentage of micronuclei

