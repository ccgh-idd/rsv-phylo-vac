# Load necessary packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")

library(dplyr)
library(readr)
library(here)

# --- File Paths ---
meta_file <- here("summaries", "all_clean_meta.csv")
output_file <- here("data", "country_vaccination_status.csv")

# 1. Read clean metadata to get the exact list of countries
meta <- read_csv(meta_file, show_col_types = FALSE)
unique_countries <- sort(unique(meta$geoLocCountry))

# 2. Create a basic template
vacc_data <- data.frame(
  country = unique_countries,
  has_rsv_vaccination = 0,    # 1 for Yes, 0 for No
  introduction_year = NA      # Year of introduction
)

# 3. Pre-populate known adopters based on national regulator approvals / programme launches
# Sources: EMA (Arexvy/Abrysvo approved June/Aug 2023), FDA, Health Canada, PMDA, ANMAT,
#          Swissmedic, NIAC, Folkhälsomyndigheten, Medsafe, Wikipedia RSV vaccine article.
# High confidence (formal national programme confirmed):
adopters_2023 <- c(
  "USA",          # FDA approved May/Jul 2023; ACIP recommended
  "Canada",       # Health Canada approved Aug 2023
  "France",       # HAS recommended 2023 season
  "Spain",        # CISNS programme 2023
  "Germany",      # STIKO recommended 2023
  "Argentina",    # ANMAT approved Abrysvo 2023
  "Japan",        # PMDA approved Arexvy Sep 2023
  "Israel",       # MoH recommendation 2023
  "Italy",        # EMA + regional programme adoption 2023
  "Ireland",      # NIAC recommended 2023
  "Switzerland",  # Swissmedic approved 2023
  "Austria",      # EMA-approved; national programme 2023
  "Finland",      # THL recommendation 2023
  "Hungary",      # EMA-approved; national adoption 2023
  "Slovakia"      # EMA-approved; national adoption 2023
)
adopters_2024 <- c(
  "United Kingdom",  # JCVI recommended; NHSE programme autumn 2024
  "Australia",       # TGA approved; NIP listed 2024
  "Netherlands",     # RIVM programme 2024
  "Belgium",         # Belgian vaccine committee 2024
  "Sweden",          # Folkhälsomyndigheten recommended 75+ in 2024
  "New Zealand"      # Medsafe approved Arexvy 2024
)

vacc_data <- vacc_data %>%
  mutate(
    has_rsv_vaccination = case_when(
      country %in% c(adopters_2023, adopters_2024) ~ 1,
      TRUE ~ 0
    ),
    introduction_year = case_when(
      country %in% adopters_2023 ~ 2023,
      country %in% adopters_2024 ~ 2024,
      TRUE ~ NA_real_
    )
  )

# 4. Write to CSV (which can be natively opened/edited in Excel)
write_csv(vacc_data, output_file)

cat("Vaccination status spreadsheet generated at:", output_file, "\n")
