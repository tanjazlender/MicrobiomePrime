# A SCRIPT TO VERIFY WHETHER THE INPUTS ARE CORRECTLY FORMATTED

# THE INPUT FILES MUST MEET THE FOLLOWING REQUIREMENTS:
# All tables (relabund_tab, metadata, taxonomy) must be in .tsv or .csv format
# Column names in relabund_tab must be the same as sequence IDs in the FASTA file
# Column names in relabund_tab must be the same as sequence IDs in the taxonomy file
# Row names in relabund_tab must be the same as sample names (column "Sample") in metadata
# The row sums of relabund_tab must be 1

################################################################################
############################# Read input files #################################
error_messages <- list()

# Read relabund_tab (in .tsv or .csv format)
if (file.exists("data/input_files/relabund_tab.tsv")) {
  relabund_tab <- read.delim("data/input_files/relabund_tab.tsv",
                             header = TRUE,
                             row.names = 1)
} else if (file.exists("data/input_files/relabund_tab.csv")) {
  relabund_tab <- read.csv("data/input_files/relabund_tab.csv",
                           header = TRUE,
                           row.names = 1)
} else {
  error_messages <- c(error_messages, "Error: No valid relabund_tab file found in the input_files folder.\n")
  # Convert the list to a character vector
  error_messages <- unlist(error_messages)
  # Write the errors to the file
  writeLines(error_messages, "scripts/input_formatting_verification/input_formatting_verification.log")
}

# Read metadata (in .tsv or .csv format)
if (file.exists("data/input_files/metadata.tsv")) {
  metadata <- read.delim("data/input_files/metadata.tsv",
                         header = TRUE)
} else if (file.exists("data/input_files/metadata.csv")) {
  metadata <- read.csv("data/input_files/metadata.csv",
                       header = TRUE)
} else {
  error_messages <- c(error_messages, "Error: No valid metadata file found in the input_files folder.\n")
  # Convert the list to a character vector
  error_messages <- unlist(error_messages)
  # Write the errors to the file
  writeLines(error_messages, "scripts/input_formatting_verification/input_formatting_verification.log")
}

# Read the taxonomy file (in .tsv or .csv format)
if (file.exists("data/input_files/taxonomy.tsv")) {
  taxonomy <- read.delim("data/input_files/taxonomy.tsv",
                         header = TRUE)
} else if (file.exists("data/input_files/taxonomy.csv")) {
  taxonomy <- read.csv("data/input_files/taxonomy.csv",
                       header = TRUE)
} else {
  error_messages <- c(error_messages, "Error: No valid taxonomy file found in the input_files folder.\n")
  # Convert the list to a character vector
  error_messages <- unlist(error_messages)
  # Write the errors to the file
  writeLines(error_messages, "scripts/input_formatting_verification/input_formatting_verification.log")
}

# Read the FASTA file
library(seqinr)
if (file.exists("data/input_files/sequences.fa")) {
  fasta_file <- read.fasta("data/input_files/sequences.fa",
                           forceDNAtolower = FALSE)
} else {
  error_messages <- c(error_messages, "Error: FASTA file not found or could not be read in the input_files folder.")
  # Convert the list to a character vector
  error_messages <- unlist(error_messages)
  # Write the errors to the file
  writeLines(error_messages, "scripts/input_formatting_verification/input_formatting_verification.log")
}

################################################################################
############################ Verify input formatting ###########################

# Column names in relabund_tab must be the same as sequence IDs in the FASTA file
fasta_seqIDs <- sapply(fasta_file, function(x) attr(x, "name"))
fasta_seqIDs <- sort(as.character(fasta_seqIDs))
relabund_tab_seqIDs <- sort(colnames(relabund_tab))

if (length(relabund_tab_seqIDs) != length(fasta_seqIDs)) {
  # Print an error if they do not match
  error_messages <- c(error_messages, "Error: The number of column names in relabund_tab does not match the number of sequence IDs in the FASTA file.")
} else {
  # Check if sequence IDs of the FASTA file match those in the relabund table
  if (!identical(relabund_tab_seqIDs, fasta_seqIDs)) {
    # Print an error if they do not match
    error_messages <- c(error_messages, "Error: Column names of the relabund_tab do not match sequence IDs of the FASTA file.")
  }
}

# Column names in relabund_tab must be the same as sequence IDs in the taxonomy file
taxonomy_seqIDs <- sort(taxonomy$SeqID)

if (length(relabund_tab_seqIDs) != length(taxonomy_seqIDs)) {
  # Print an error if they do not match
  error_messages <- c(error_messages, "Error: The number of columns of the relabund_tab does not match the number of sequence IDs in the taxonomy table.")
} else {
  # Check if sequence IDs of the relabund table match those in the taxonomy file
  if (!identical(relabund_tab_seqIDs, taxonomy_seqIDs)) {
    # Print an error if they do not match
    error_messages <- c(error_messages, "Error: Column names of the relabund_tab do not match the sequence IDs in the taxonomy table.")
  }
}

# Row names in relabund_tab must be the same as the samples in metadata
relabund_tab_samples <- sort(rownames(relabund_tab))
metadata_samples <- sort(metadata$Sample)

if (length(relabund_tab_samples) != length(metadata_samples)) {
  # Print an error if they do not match
  error_messages <- c(error_messages, "Error: The number of rows in the relabund_tab does not match the number of samples in the metadata table.")
} else {
  # Check if row names in the relabund_tab match sample names in the metadata table
  if (!identical(relabund_tab_samples, metadata_samples)) {
    # Print an error if they do not match
    error_messages <- c(error_messages, "Error: Row names in the relabund_tab do not match sample names in the metadata table.")
  }
}

# The row sums of relabund_tab must be 1
tolerance <- 0.05
if (any(abs(rowSums(relabund_tab) - 1) > tolerance)) {
  # Print an error if any row sum is not approximately 1
  error_messages <- c(error_messages, "Error: The sum of each row in the relabund_tab is not 1.")
}

# Write error messages to input_formatting_verification.log
if (length(error_messages) > 0) {
  # Convert the list to a character vector
  error_messages <- unlist(error_messages)
  # Write the errors to the file
  writeLines(error_messages, "scripts/input_formatting_verification/input_formatting_verification.log")
} else {
  # Write success message to the file
  writeLines("The input files are formatted correctly.", "scripts/input_formatting_verification/input_formatting_verification.log")
}

