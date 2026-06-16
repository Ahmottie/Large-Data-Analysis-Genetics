# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyverse)

# Read the data into a data frame
df <- read_excel("data/aneuploid_micronuclei.xlsx")

# Add the percentage of micronuclei per observation
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

# View the results
print(df_summary)


# Load ggplot2 if not already loaded
library(ggplot2)

# Create the bar plot with error bars
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


# Create the box plot
box_plot <- ggplot(df, aes(x = condition, y = perc_micronuclei)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Distribution of Micronuclei Percentages",
       x = "Cell Line Condition",
       y = "Micronuclei per Observation (%)") +
  theme_classic()

print(box_plot)
ggsave("./results/box_plot.pdf", plot = box_plot, width = 10, height = 5, device = "pdf")


# Tetraploidy has potentially a larger effect on the formation of micronucleai

# Updated box plot with min/max whiskers and error bar caps
box_plot <- ggplot(df, aes(x = condition, y = perc_micronuclei)) +
  # 1. Add the error bar caps first so they render behind the box
  # coef = Inf forces them to the absolute minimum and maximum values
  stat_boxplot(geom = "errorbar", width = 0.2, coef = Inf) +

  # 2. Add the box plot itself
  # coef = Inf prevents individual outlier points from being drawn
  geom_boxplot(fill = "lightblue", color = "black", coef = Inf) +

  labs(title = "Distribution of Micronuclei Percentages",
       x = "Cell Line Condition",
       y = "Micronuclei per Observation (%)") +
  theme_classic() +
  # Applying the larger font sizes for consistency
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 14)
  )

print(box_plot)
ggsave("./results/box2_plot.pdf", plot = box_plot, width = 10, height = 5, device = "pdf")


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

# Print the results to view the p-values and test statistics
print(test_a)
print(test_b)
print(test_c)
print(test_d)