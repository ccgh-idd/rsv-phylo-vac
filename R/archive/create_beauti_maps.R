suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(lubridate)
  library(here)
})

# 1. Load data
meta <- read_csv(here("data", "beast_input", "subsampled_metadata.csv"), show_col_types = FALSE)
vacc <- read_csv(here("data", "country_vaccination_status.csv"), show_col_types = FALSE)

# 2. Add the Vaccination Trait
plot_data <- meta %>%
  left_join(vacc, by = c("geoLocCountry" = "country")) %>%
  mutate(
    # Format date as YYYY-MM-DD for BEAST interpretation (defaulting to middle of year/month if missing)
    formatted_date = format(parsed_date, "%Y-%m-%d"),
    trait_status = case_when(
      has_rsv_vaccination == 1 & !is.na(introduction_year) & year >= introduction_year ~ "Vaccinated",
      TRUE ~ "Unvaccinated"
    )
  )

# 3. Create BEAUti mapping files for RSV-A
data_a <- plot_data %>% filter(subtype == "RSV-A")
# BEAUti expects a simple 2-column TSV with no header for imports
write.table(data_a %>% select(accessionVersion, formatted_date), 
            file = here("data", "beast_input", "beauti_dates_A.tsv"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data_a %>% select(accessionVersion, trait_status), 
            file = here("data", "beast_input", "beauti_traits_A.tsv"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


# 4. Create BEAUti mapping files for RSV-B
data_b <- plot_data %>% filter(subtype == "RSV-B")
write.table(data_b %>% select(accessionVersion, formatted_date), 
            file = here("data", "beast_input", "beauti_dates_B.tsv"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(data_b %>% select(accessionVersion, trait_status), 
            file = here("data", "beast_input", "beauti_traits_B.tsv"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Successfully created BEAUti mapping files in data/beast_input/\n")
