# Load necessary packages
if (!requireNamespace("seqinr", quietly = TRUE)) {
  options(timeout = 300, repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages("seqinr")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("lubridate", quietly = TRUE)) {
  install.packages("lubridate")
}
library(seqinr)
library(dplyr)
library(lubridate)

# --- Function to process and summarize data ---
process_rsv_data <- function(fasta_file, metadata_file, output_summary_file, output_country_year_file) {

  # --- Part 1: Cleaning Summary ---
  sequences <- read.fasta(file = fasta_file, seqtype = "AA")
  initial_count <- length(sequences)

  # Filter out sequences with more than 50% 'X'
  cleaned_sequences_list <- sequences[sapply(sequences, function(seq) {
    seq_string <- getSequence(seq, as.string = TRUE)
    x_count <- nchar(gsub("[^X]", "", seq_string))
    x_count / nchar(seq_string) < 0.5
  })]
  final_count <- length(cleaned_sequences_list)
  removed_count <- initial_count - final_count

  # Write cleaning summary
  summary_text <- paste(
    "Cleaning Summary for:", basename(fasta_file), "\n",
    "--------------------------------------------------\n",
    "Initial number of sequences:", initial_count, "\n",
    "Number of sequences removed (more than 50% 'X'):", removed_count, "\n",
    "Final number of sequences:", final_count, "\n"
  )
  writeLines(summary_text, output_summary_file)
  print(paste("Cleaning summary saved to", output_summary_file))

  # --- Part 2: Metadata Analysis & Filtering ---
  metadata <- read.delim(metadata_file, sep = "\t")

  # Get the names of the cleaned sequences
  cleaned_ids <- names(cleaned_sequences_list)

  # Filter metadata to include only cleaned sequences and parse dates
  cleaned_metadata <- metadata %>%
    filter(accessionVersion %in% cleaned_ids) %>%
    mutate(
      parsed_date = parse_date_time(sampleCollectionDate, orders = c("d.m.y", "dmy", "Y-m", "Y")),
      year = year(parsed_date)
    ) %>%
    filter(!is.na(year) & year >= 2010) # Keep dates from 2010 onwards

  # --- Part 3: Create new filtered FASTA files ---
  # Get the IDs of the sequences to keep
  final_ids_to_keep <- cleaned_metadata$accessionVersion
  final_sequences_to_keep <- cleaned_sequences_list[names(cleaned_sequences_list) %in% final_ids_to_keep]

  # Write the final filtered sequences to a new FASTA file
  write.fasta(sequences = final_sequences_to_keep, names = names(final_sequences_to_keep), file.out = output_filtered_fasta_file)
  print(paste("Final filtered FASTA file saved to", output_filtered_fasta_file))


  # --- Part 4: Summarize by country and year ---
  country_year_summary <- cleaned_metadata %>%
    group_by(geoLocCountry, year) %>%
    summarise(count = n(), .groups = 'drop')

  # Write country-year summary to a CSV file
  write.csv(country_year_summary, output_country_year_file, row.names = FALSE)
  print(paste("Country-year summary saved to", output_country_year_file))
}

# --- File Paths ---
base_dir <- "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/phlo/rsv_vac/rsvphylovac"
raw_data_dir <- file.path(base_dir, "raw_data")
output_dir <- file.path(base_dir, "summaries")

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# --- Process RSV-A ---
process_rsv_data(
  fasta_file = file.path(raw_data_dir, "rsv-a_aligned-aa-F_2026-03-05T1242.fasta"),
  metadata_file = file.path(raw_data_dir, "rsv-a_metadata_2026-03-05T1245.tsv"),
  output_summary_file = file.path(output_dir, "rsv-a_cleaning_summary.txt"),
  output_country_year_file = file.path(output_dir, "rsv-a_country_year_summary.csv")
)

# --- Process RSV-B ---
process_rsv_data(
  fasta_file = file.path(raw_data_dir, "rsv-b_aligned-aa-F_2026-03-05T1244.fasta"),
  metadata_file = file.path(raw_data_dir, "rsv-b_metadata_2026-03-05T1245.tsv"),
  output_summary_file = file.path(output_dir, "rsv-b_cleaning_summary.txt"),
  output_country_year_file = file.path(output_dir, "rsv-b_country_year_summary.csv")
)
