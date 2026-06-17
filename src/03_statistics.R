library(readxl)
library(dplyr)
library(tidyverse)

df <- read_excel("data/aneuploid_micronuclei.xlsx")

df <- df %>%
  mutate(perc_micronuclei = (micronuclei / nuclei) * 100)

# Group by condition and calculate summary statistics
df_summary <- df %>%
  group_by(condition) %>%
  summarise(
    N = n(),
    mean_perc = mean(perc_micronuclei),
    se_perc = sd(perc_micronuclei) / sqrt(N)

  )

print(df_summary)
library(ggplot2)

bar_plot <- ggplot(df_summary, aes(x = condition, y = mean_perc)) +
  geom_col(fill = "gray70", color = "black", width = 0.5) +
  geom_errorbar(aes(ymin = mean_perc - se_perc, ymax = mean_perc + se_perc),
                width = 0.2) +
  labs(title = "Average Percentage of Micronuclei per Condition",
       x = "Condition",
       y = "Average Micronuclei %") +
  theme(
    plot.title = element_text(size = 20),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 17),
    axis.title.y = element_text(size = 17)
  )

print(bar_plot)
ggsave("./results/bar_plot.pdf", plot = bar_plot, width = 10, height = 5, device = "pdf")

library(ggpubr)

comparisons <- list(
  c("HCT116+APH", "HCT116"),
  c("Htr5", "HCT116"),
  c("Hte5", "HCT116"),
  c("HCT116+APH", "Hte5")
)
# box plot
box_plot <- ggplot(df, aes(x = condition, y = perc_micronuclei)) +
  stat_boxplot(geom = "errorbar", width = 0.2, coef = Inf) +
  geom_boxplot(fill = "lightblue", color = "black", coef = Inf, width = 0.4) +
  # t.test
  stat_compare_means(comparisons = comparisons, method = "t.test", step.increase = 0.08) +

  labs(title = "Distribution of Micronuclei Percentages",
       x = "Cell Line Condition",
       y = "Micronuclei per Observation (%)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14)
  )

print(box_plot)

ggsave("./results/box_plot.pdf", plot = box_plot, width = 7, height = 6, device = "pdf")

# log scale
box_plot_log <- box_plot +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 2.5, 5, 7.5, 10)) +
  labs(y = "Micronuclei per Observation (%) [Log Scale]")

print(box_plot_log)
ggsave("./results/box_plot_log.pdf", plot = box_plot_log, width = 7, height = 6, device = "pdf")

# a. HCT116+APH vs HCT116 (Positive vs Negative Control)
test_a <- t.test(perc_micronuclei ~ condition,
                 data = subset(df, condition %in% c("HCT116+APH", "HCT116")),
                 var.equal = FALSE)

# b. Htr5 vs HCT116 (Aneuploid Model 1 vs Negative Control)
test_b <- t.test(perc_micronuclei ~ condition,
                 data = subset(df, condition %in% c("Htr5", "HCT116")),
                 var.equal = FALSE)

# c. Hte5 vs HCT116 (Aneuploid Model 2 vs Negative Control)
test_c <- t.test(perc_micronuclei ~ condition,
                 data = subset(df, condition %in% c("Hte5", "HCT116")),
                 var.equal = FALSE)

# d. HCT116+APH vs Hte5 (Positive Control vs Aneuploid Model 2)
test_d <- t.test(perc_micronuclei ~ condition,
                 data = subset(df, condition %in% c("HCT116+APH", "Hte5")),
                 var.equal = FALSE)

# Print the results
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
