# Load the seqinr package
if (!requireNamespace("seqinr", quietly = TRUE)) {
  options(timeout = 300)
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages("seqinr")
}
library(seqinr)

# Function to clean a FASTA file
clean_fasta <- function(input_file, output_file) {
  # Read the FASTA file
  sequences <- read.fasta(file = input_file, seqtype = "AA")

  # Filter out sequences with more than 50% 'X'
  cleaned_sequences <- sequences[sapply(sequences, function(seq) {
    seq_string <- getSequence(seq, as.string = TRUE)
    # Count the number of 'X' characters
    x_count <- nchar(gsub("[^X]", "", seq_string))
    # Keep the sequence if less than 50% are 'X'
    x_count / nchar(seq_string) < 0.5
  })]

  # Check if the sequences are aligned
  if (length(cleaned_sequences) > 1) {
    seq_lengths <- sapply(cleaned_sequences, function(seq) nchar(getSequence(seq, as.string = TRUE)))
    if (length(unique(seq_lengths)) != 1) {
      warning(paste("Sequences in", input_file, "are not all of the same length after cleaning. They may not be aligned."))
    } else {
      print(paste("Cleaned sequences in", input_file, "are aligned."))
    }
  }

  # Write the cleaned sequences to a new FASTA file
  write.fasta(sequences = cleaned_sequences, names = names(cleaned_sequences), file.out = output_file)
  print(paste("Cleaned file saved to", output_file))
}

# Create output directory if it doesn't exist
output_dir <- "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/phlo/rsv_vac/rsvphylovac/data/sequences"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Clean the RSV-A file
clean_fasta(
  input_file = "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/phlo/rsv_vac/rsvphylovac/raw_data/rsv-a_aligned-aa-F_2026-03-05T1242.fasta",
  output_file = file.path(output_dir, "rsv-a_cleaned.fasta")
)

# Clean the RSV-B file
clean_fasta(
  input_file = "/Users/davidhodgson/Dropbox/Mac (3)/Documents/research/phlo/rsv_vac/rsvphylovac/raw_data/rsv-b_aligned-aa-F_2026-03-05T1244.fasta",
  output_file = file.path(output_dir, "rsv-b_cleaned.fasta")
)
