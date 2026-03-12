# Load necessary packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", repos = "https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")

library(ggplot2)
library(dplyr)
library(readr)
library(here)

# --- File Paths ---
meta_file <- here("summaries", "all_clean_meta.csv")
vacc_file <- here("data", "country_vaccination_status.csv")
output_dir <- here("figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Read files
if (!file.exists(meta_file) || !file.exists(vacc_file)) {
  stop("Metadata or vaccination status file not found.")
}

meta <- read_csv(meta_file, show_col_types = FALSE)
vacc_status <- read_csv(vacc_file, show_col_types = FALSE)

# 2. Join and determine vaccination setting
plot_data <- meta %>%
  left_join(vacc_status, by = c("geoLocCountry" = "country")) %>%
  mutate(
    # A sample is in a "Vaccinated Setting" if the country has vaccination AND the sample year is >= the introduction year
    is_vaccinated_setting = case_when(
      has_rsv_vaccination == 1 & !is.na(introduction_year) & year >= introduction_year ~ "Vaccinated Setting",
      TRUE ~ "Unvaccinated Setting"
    )
  )

# 3. Create the plot
p <- ggplot(plot_data, aes(x = year, fill = is_vaccinated_setting)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c("Vaccinated Setting" = "#00BFC4", "Unvaccinated Setting" = "#F8766D")) +
  theme_minimal() +
  labs(
    title = "RSV Samples per Year by Vaccination Setting",
    x = "Year",
    y = "Number of Samples",
    fill = "Status"
  ) +
  theme(legend.position = "bottom")

# 4. Save the plot
output_plot_file <- here("figures", "vaccination_samples_per_year.png")
ggsave(output_plot_file, plot = p, width = 10, height = 6, dpi = 300)

cat("Successfully generated plot at:", output_plot_file, "\n")
