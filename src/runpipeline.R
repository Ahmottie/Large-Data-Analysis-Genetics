# Load and tidy raw data
source("src/00_rawdata.R")

# Process data and calculate statistics
source("src/01_processing.R")

# Generate and export visualizations
source("src/02_visualization.R")

print("Pipeline execution complete.")