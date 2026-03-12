# Load required packages
if (!requireNamespace("seqinr", quietly = TRUE)) install.packages("seqinr", repos = "https://cloud.r-project.org")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate", repos = "https://cloud.r-project.org")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")

library(seqinr)
library(dplyr)
library(lubridate)
library(readr)
library(here)

prepare_data <- function(fasta_path, meta_path, out_fasta, out_meta_csv, subtype_name) {
  cat(sprintf("\n--- Processing %s ---\n", subtype_name))
  
  # 1. Read FASTA
  cat("Reading fasta...\n")
  seqs <- read.fasta(file = fasta_path, seqtype = "AA")
  
  # Filter < 10% X
  valid_seqs <- seqs[sapply(seqs, function(s) {
    seq_str <- getSequence(s, as.string = TRUE)
    x_count <- nchar(gsub("[^X]", "", seq_str))
    (x_count / nchar(seq_str)) <= 0.10
  })]
  
  cat(sprintf("FASTA valid sequences (<=10%% X): %d / %d\n", length(valid_seqs), length(seqs)))
  
  # 2. Read Metadata
  cat("Reading metadata...\n")
  meta <- read.delim(meta_path, sep = "\t", stringsAsFactors = FALSE)
  
  # 3. Filter Metadata
  clean_meta <- meta %>%
    filter(accessionVersion %in% names(valid_seqs)) %>%
    mutate(
      parsed_date = parse_date_time(sampleCollectionDate, orders = c("d.m.y", "dmy", "Y-m", "Y", "y-m-d", "ymd")),
      year = year(parsed_date)
    ) %>%
    filter(!is.na(year) & year >= 2010) %>%
    filter(!is.na(geoLocCountry) & geoLocCountry != "")
    
  cat(sprintf("Metadata valid rows (>=2010, Has Date, Has Country): %d\n", nrow(clean_meta)))
  
  # 4. Intersect
  final_ids <- clean_meta$accessionVersion
  final_seqs <- valid_seqs[final_ids]
  cat(sprintf("Final intersection count: %d\n", length(final_seqs)))
  
  # 5. Write outputs
  cat("Writing clean outputs...\n")
  write.fasta(sequences = final_seqs, names = names(final_seqs), file.out = out_fasta)
  write_csv(clean_meta, out_meta_csv)
  
  return(clean_meta %>% mutate(subtype = subtype_name))
}

dir.create(here("data", "sequences"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("summaries"), recursive = TRUE, showWarnings = FALSE)

meta_a <- prepare_data(
  here("raw_data", "rsv-a_aligned-aa-F_2026-03-05T1242.fasta"),
  here("raw_data", "rsv-a_metadata_2026-03-05T1245.tsv"),
  here("data", "sequences", "rsv-a_analysis_ready.fasta"),
  here("summaries", "rsv-a_clean_meta.csv"),
  "RSV-A"
)

meta_b <- prepare_data(
  here("raw_data", "rsv-b_aligned-aa-F_2026-03-05T1244.fasta"),
  here("raw_data", "rsv-b_metadata_2026-03-05T1245.tsv"),
  here("data", "sequences", "rsv-b_analysis_ready.fasta"),
  here("summaries", "rsv-b_clean_meta.csv"),
  "RSV-B"
)

all_meta <- bind_rows(meta_a, meta_b)
write_csv(all_meta, here("summaries", "all_clean_meta.csv"))
cat("\nPipeline complete! Final summary data saved to summaries/all_clean_meta.csv\n")