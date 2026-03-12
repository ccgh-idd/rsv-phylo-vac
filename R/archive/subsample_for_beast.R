# Load necessary packages
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = "https://cloud.r-project.org")
if (!requireNamespace("seqinr", quietly = TRUE)) install.packages("seqinr", repos = "https://cloud.r-project.org")
if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = "https://cloud.r-project.org")
if (!requireNamespace("here", quietly = TRUE)) install.packages("here", repos = "https://cloud.r-project.org")

library(dplyr)
library(seqinr)
library(readr)
library(here)

set.seed(42) # For reproducible subsampling

target_n <- 800 # Target number of sequences per subtype suitable for BEAST 2

# 1. Read metadata
meta_file <- here("summaries", "all_clean_meta.csv")
meta <- read_csv(meta_file, show_col_types = FALSE)

# 2. Subsample intelligently
# We want to preserve geographic and temporal diversity without being totally 
# swamped by massively overrepresented countries/years (e.g., USA 2023).
# We calculate a weight based on the square root of the group count.
cat("Calculating optimal subsample representation...\n")

subsample_meta <- meta %>%
  group_by(subtype, geoLocCountry, year) %>%
  mutate(group_count = n()) %>%
  ungroup() %>%
  mutate(sampling_weight = 1 / sqrt(group_count)) %>% # Underrepresented groups get higher relative weight
  group_by(subtype) %>%
  slice_sample(n = target_n, weight_by = sampling_weight) %>%
  ungroup()

dir.create(here("data", "beast_input"), recursive = TRUE, showWarnings = FALSE)
write_csv(subsample_meta, here("data", "beast_input", "subsampled_metadata.csv"))

# 3. Read and subsample FASTA for RSV-A
cat("Subsampling RSV-A FASTA...\n")
fasta_a_path <- here("data", "sequences", "rsv-a_analysis_ready.fasta")
seqs_a <- read.fasta(file = fasta_a_path, seqtype = "AA")

a_ids <- subsample_meta %>% filter(subtype == "RSV-A") %>% pull(accessionVersion)
seqs_a_sub <- seqs_a[names(seqs_a) %in% a_ids]
write.fasta(sequences = seqs_a_sub, names = names(seqs_a_sub), file.out = here("data", "beast_input", "rsv-a_subsampled.fasta"))

# 4. Read and subsample FASTA for RSV-B
cat("Subsampling RSV-B FASTA...\n")
fasta_b_path <- here("data", "sequences", "rsv-b_analysis_ready.fasta")
seqs_b <- read.fasta(file = fasta_b_path, seqtype = "AA")

b_ids <- subsample_meta %>% filter(subtype == "RSV-B") %>% pull(accessionVersion)
seqs_b_sub <- seqs_b[names(seqs_b) %in% b_ids]
write.fasta(sequences = seqs_b_sub, names = names(seqs_b_sub), file.out = here("data", "beast_input", "rsv-b_subsampled.fasta"))

cat(sprintf("\nSubsampling complete! %d sequences per subtype saved to data/beast_input/\n", target_n))