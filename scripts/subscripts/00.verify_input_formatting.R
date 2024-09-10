# A SCRIPT TO VERIFY WHETHER THE INPUTS AND VARIABLES ARE CORRECTLY FORMATTED

# THE INPUT FILES MUST MEET THE FOLLOWING REQUIREMENTS:
# All tables (relabund_tab, metadata, taxonomy) must be in .tsv or .csv format
# Column names in relabund_tab must be the same as sequence IDs in the FASTA file
# Column names in relabund_tab must be the same as sequence IDs in the taxonomy file
# Row names in relabund_tab must be the same as sample names (column "Sample") in metadata
# The row sums of relabund_tab must be 1

# Set the working directory to the project directory
setwd("../")

library(dplyr)

################################################################################
############################# Read input files #################################
error_messages <- list()
warning_messages <- list()

log_file <- "scripts/current_logs/00.verify_input_formatting.log"

# Check if the log directory exists; if not, create it
log_dir <- dirname(log_file)
if (!dir.exists(log_dir)) {
  dir.create(log_dir, recursive = TRUE)
}

# Delete the log file if it exists
if (file.exists(log_file)) {
  file.remove(log_file)
}

# Open a connection to the log file
con <- file(log_file, "w")
on.exit(close(con), add = TRUE)  # Ensure the connection is closed on exit


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
  writeLines(error_messages, log_file)
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
  writeLines(error_messages, log_file)
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
  writeLines(error_messages, log_file)
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
  writeLines(error_messages, log_file)
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

################################################################################
########################## Verify variables formatting #########################

# Check whether are key variables are not empty

library(config)
library(ini)

# Read parameters from variables.ini
variables <- read.ini("scripts/variables.ini")

# Check if target is NULL and print an error message if necessary
target_raw <- variables$settings$target

if (is.null(target_raw)) {
  error_messages <- c(error_messages, "Error: target value is missing in the variables.ini file.")
} else {
  target_list <- strsplit(target_raw, ",")
  target <- trimws(unlist(target_list))
  
  # Check if any targets are not found in metadata and print an error if necessary
  targets_not_in_metadata <- character()
  for (t in target) {
    print(t)
    if (!(t %in% metadata$Source)) {
      targets_not_in_metadata <- c(targets_not_in_metadata, t)
    }
    
    # Print an error if any of the defined target sources was not found in metadata
    if (length(targets_not_in_metadata) > 0) {
      error_messages <- c(error_messages, paste("Error: One or more defined targets (", paste(targets_not_in_metadata, collapse = ", "), ") in the variables.ini file is not present in the `Source` column of metadata.", sep = ""))
    }
  }
}

# If multiple targets are defined and no target_group_name is set, print a warning
target_group_name <- variables$settings$target_group_name

if (length(target) > 1 & is.null(target_group_name)) {
  warning_messages <- c(warning_messages, paste("Warning: setting target_group_name in the variables.ini file is recommended when more than one target source is defined."))
}

# Check whether kmer_size is set
kmer_size <- as.numeric(variables$settings$kmer_size)

if (length(kmer_size) == 0) {
  error_messages <- c(error_messages, "Error: kmer_size value is missing in the variables.ini file.")
}

# Check whether sensitivity and specificity criteria are set
kmer_sensitivity_cutoff <- as.numeric(variables$settings$kmer_sensitivity_cutoff)
kmer_specificity_cutoff <- as.numeric(variables$settings$kmer_specificity_cutoff)
marker_sensitivity_cutoff <- as.numeric(variables$settings$marker_sensitivity_cutoff)
marker_specificity_cutoff <- as.numeric(variables$settings$marker_specificity_cutoff)

if (length(kmer_sensitivity_cutoff) == 0) {
  error_messages <- c(error_messages, "Error: kmer_sensitivity_cutoff value is missing in the variables.ini file.")
}

if (length(kmer_specificity_cutoff) == 0) {
  error_messages <- c(error_messages, "Error: kmer_specificity_cutoff value is missing in the variables.ini file.")
}

if (length(marker_sensitivity_cutoff) == 0) {
  error_messages <- c(error_messages, "Error: marker_sensitivity_cutoff value is missing in the variables.ini file.")
}

if (length(marker_specificity_cutoff) == 0) {
  error_messages <- c(error_messages, "Error: marker_specificity_cutoff value is missing in the variables.ini file.")
}


# Check whether mismatch or primer-binding stability variables are set
max_mismatch <- as.numeric(variables$settings$max_mismatch)
max_primer_delta <- as.numeric(variables$settings$max_primer_delta)

if (length(max_mismatch) == 0 & length(max_primer_delta) == 0) {
  error_messages <- c(error_messages, "Error: setting max_mismatch in the variables.ini file is recommended. Alternatively you can set max_primer_delta.")
}


# Write error messages to input_formatting_verification.log
success_message <- "The input files are formatted correctly.\nDONE: The script has completed successfully."

if (length(error_messages) > 0) {
  # Convert the list to a character vector
  error_messages <- unlist(error_messages)
  warning_messages <- unlist(warning_messages)
  # Write the errors and warnings to the log file
  writeLines(c(error_messages, warning_messages), log_file)
} else {
  # Convert the list to a character vector
  warning_messages <- unlist(warning_messages)
  # Write the warnings and a success message to the log file
  if (length(warning_messages) > 0) {
    writeLines(c(warning_messages, success_message), log_file)
  } else {
    writeLines(success_message, log_file)
  }
  
}





