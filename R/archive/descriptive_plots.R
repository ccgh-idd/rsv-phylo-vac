# Load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")

library(ggplot2)
library(dplyr)
library(readr)

# --- File Paths ---
base_dir <- "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/phlo/rsv_vac/rsvphylovac"
output_dir <- file.path(base_dir, "figures", "descriptive")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

data_file <- file.path(base_dir, "summaries", "all_clean_meta.csv")

if (!file.exists(data_file)) {
  stop("The file all_clean_meta.csv does not exist. Please run prepare_final_data.R first.")
}

data <- read_csv(data_file, show_col_types = FALSE)

# --- Plot 1: Samples per year by subtype ---
p1 <- ggplot(data, aes(x = year, fill = subtype)) +
  geom_bar(position = "dodge") +
  theme_minimal() +
  labs(title = "Temporal Distribution of RSV Samples (>= 2010)", x = "Year", y = "Count")

ggsave(file.path(output_dir, "temporal_distribution.png"), plot = p1, width = 10, height = 6, dpi = 300)
cat("Saved temporal_distribution.png\n")

# --- Plot 2: Top 20 countries by sample count ---
top_countries <- data %>% 
  count(geoLocCountry, sort = TRUE) %>% 
  head(20)

p2 <- data %>% 
  filter(geoLocCountry %in% top_countries$geoLocCountry) %>%
  ggplot(aes(x = reorder(geoLocCountry, geoLocCountry, function(x) -length(x)), fill = subtype)) +
  geom_bar() +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 20 Countries by Sample Count", x = "Country", y = "Count")

ggsave(file.path(output_dir, "geographical_distribution.png"), plot = p2, width = 10, height = 8, dpi = 300)
cat("Saved geographical_distribution.png\n")

# --- Plot 3: Samples by year for Top 12 countries ---
top_12_countries <- head(top_countries$geoLocCountry, 12)

p3 <- data %>%
  filter(geoLocCountry %in% top_12_countries) %>%
  ggplot(aes(x = year, fill = subtype)) +
  geom_bar(position = "stack") +
  facet_wrap(~geoLocCountry, scales = "free_y") +
  theme_minimal() +
  labs(title = "RSV Samples by Year (Top 12 Countries)", x = "Year", y = "Count")

ggsave(file.path(output_dir, "rsv_samples_by_country_year.png"), plot = p3, width = 12, height = 8, dpi = 300)
cat("Saved rsv_samples_by_country_year.png\n")

cat("\nAll descriptive plots successfully generated in figures/descriptive/\n")