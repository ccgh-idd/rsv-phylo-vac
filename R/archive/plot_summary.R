# Load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
library(ggplot2)
library(dplyr)
library(readr)

# --- File Paths ---
base_dir <- "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/phlo/rsv_vac/rsvphylovac"
summaries_dir <- file.path(base_dir, "summaries")
output_dir <- file.path(base_dir, "figures")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# --- Read and Combine Data ---
rsv_a_summary <- read_csv(file.path(summaries_dir, "rsv-a_country_year_summary.csv")) %>%
  mutate(type = "RSV-A")
rsv_b_summary <- read_csv(file.path(summaries_dir, "rsv-b_country_year_summary.csv")) %>%
  mutate(type = "RSV-B")

combined_summary <- bind_rows(rsv_a_summary, rsv_b_summary) %>%
  rename(country = geoLocCountry)

# Filter for top N countries to keep the plot readable
top_countries <- combined_summary %>%
  group_by(country) %>%
  summarise(total = sum(count)) %>%
  top_n(20, total) %>%
  pull(country)

plot_data <- combined_summary %>%
  filter(country %in% top_countries)

# --- Create Plot ---
sample_plot <- ggplot(plot_data, aes(x = reorder(country, -count), y = count, fill = as.factor(year))) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~type, scales = "free_x") +
  labs(
    title = "Number of RSV Samples by Country and Year",
    x = "Country",
    y = "Number of Samples",
    fill = "Year"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Save Plot ---
output_plot_file <- file.path(output_dir, "rsv_samples_by_country_year.png")
ggsave(output_plot_file, plot = sample_plot, width = 14, height = 8, dpi = 300)

print(paste("Plot saved to", output_plot_file))
