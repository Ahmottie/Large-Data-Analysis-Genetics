# Load packages
library(ggplot2)

normalized_stats <- readRDS("./bin/normalized.rds")

#=========================================
# 3.2 Plotting the data
#=========================================
# x, y, and grouping by cell_line
growth_plot <- ggplot(normalized_stats, aes(x = time, y = mean, color = cell_line)) +
  # a. Line plot styled by cell line (defaults to solid when no linetype is mapped)
  geom_line() +
  # c. Add mean RLUs as points
  geom_point() +
  # e. Add error bars for standard error
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  # b. Separate measurements from different experiments into separate panels
  facet_grid(. ~ experiment, labeller = as_labeller(function(x) paste("Experiment", x))) +
  # d. Apply proper axis and legend labels
  labs(
    title = "Normalized Cell Growth Over Time",
    x = "Time (days)",
    y = "Average RLU (millions)",
    color = "Cell Line"
  ) +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")

#=========================================
# f. Export figure as PDF
#=========================================
ggsave("./results/plot.pdf", plot = growth_plot, width = 10, height = 5, device = "pdf")
