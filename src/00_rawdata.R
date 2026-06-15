# Load packages
library(dplyr)
library(tidyr)

#=========================================
# 1. Raw Data Processing
#=========================================
# 1. preprocessing
# a. read data
data1 <- read.csv("./data/aneuploid_growth_data1.csv")
data2 <- read.csv("./data/aneuploid_growth_data2.csv")
data3 <- read.csv("./data/aneuploid_growth_data3.csv")

head(data1)
head(data2)
head(data3)


# ----------------------------------
# Tidy data and transform data
# ----------------------------------
# b. join and tidy
tidy_data <- bind_rows(data1, data2, data3) %>%
  drop_na()

# c. Transform the time points from hours to days and the RLU into RLU per million
# d. Specify the cell line names as a factor with levels ordered "H", "Htr5" and "Hte5"
clean_data <- tidy_data %>%
  mutate(
    time = time / 24,
    RLU = RLU / 1000000,
    cell_line = factor(cell_line, levels = c("H", "Htr5", "Hte5"))
  )

# ----------------------------------
# Save clean data
# ----------------------------------
saveRDS(clean_data, "./bin/cleaned_data.rds")
