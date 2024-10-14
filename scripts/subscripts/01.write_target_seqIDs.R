# This script creates a list of Sequence IDs found in target samples.
# This list will be used in another script to extract FASTA sequences.

# Set the working directory to the project directory
setwd("../")

library(dplyr)
library(tidyr)

################################################################################
############################# Read input files #################################
cat("Reading input files.\n")

# Read relabund_tab
if (file.exists("data/input_files/relabund_tab.tsv")) {
  relabund_tab <- read.delim("data/input_files/relabund_tab.tsv", 
                                header = TRUE, 
                                row.names = 1)
} else if (file.exists("data/input_files/relabund_tab.csv")) {
  relabund_tab <- read.csv("data/input_files/relabund_tab.csv", 
                              header = TRUE, 
                              row.names = 1)
} else {
  stop("No relabund_tab found in the input_files folder.\n")
}

# Read metadata
if (file.exists("data/input_files/metadata.tsv")) {
  metadata <- read.delim("data/input_files/metadata.tsv", 
                         header = TRUE)
} else if (file.exists("data/input_files/metadata.csv")) {
  metadata <- read.csv("data/input_files/metadata.csv", 
                       header = TRUE)
} else {
  stop("No metadata found in the input_files folder.\n")
}

################################################################################
############################### Read variables #################################

# Read parameters from variables.ini
cat("Reading variables.\n")
variables <- ini::read.ini("scripts/variables.ini")
target_list <- strsplit(variables$settings$target, ",")
target <- trimws(unlist(target_list))
target_combined <- paste(target, collapse = " ")
target_group_name <- variables$settings$target_group_name

# If target group name is not defined, set automatic target group name
if (is.null(target_group_name)) {
  target_group_ID <- gsub(" ", "-", fixed=TRUE, target_combined)
} else {
  target_group_ID <- gsub(" ", "-", fixed=TRUE, target_group_name)
}

# Check if all targets are found in metadata and print an error if necessary
cat("Verifying if all specified targets are present in the metadata file.\n")
missing_targets <- target[!target %in% metadata$Source]

if (length(missing_targets) > 0) {
  error_message <- paste("ERROR: these targets were not found in the metadata file:", paste(missing_targets, collapse = ", "))
  stop(error_message, "\n")
} else {
  cat("All defined targets are found in metadata.\n")
}

################################################################################
######################## Write sequence ID lists ###############################
cat("Writing sequence IDs.\n")

# Create a subdirectory for K-mers
kmer_dir <- file.path("data/generated_files/kmers")

if (!file.exists(kmer_dir)) {
  dir.create(kmer_dir, showWarnings = FALSE)
}

# Create subdirectories for sequence ID list files and FASTA files to be generated
sequences_dir <- file.path("data/generated_files/sequences")
fasta_dir <- file.path("data/generated_files/sequences/fasta_files")
seqID_dir <- file.path("data/generated_files/sequences/seqID_lists")

if (!file.exists(sequences_dir)) {
  dir.create(sequences_dir, showWarnings = FALSE)
}

if (!file.exists(fasta_dir)) {
  dir.create(fasta_dir, showWarnings = FALSE)
}

if (!file.exists(seqID_dir)) {
  dir.create(seqID_dir, showWarnings = FALSE)
}


# List to keep track of all generated file paths
file_paths <- vector("list", length(target))
names(file_paths) <- sapply(target, function(t) gsub(" ", "-", fixed=TRUE, t))

# Generate a list of sequence IDs for each target source
for (t in target) {
  t_ID <- gsub(" ", "-", fixed=TRUE, t)
  
  # Filter metadata (keep only target t samples)
  metadata_target <- metadata %>%
    filter(Source == t) 
  
  # Filter relabund_tab (keep only target samples)
  relabund_tab_target <- filter(relabund_tab, 
                                rownames(relabund_tab) %in% unique(metadata_target$Sample))
  
  # Find out, which sequence IDs are detected in target t samples
  target_seqIDs <- relabund_tab_target %>%
    mutate(Sample = rownames(.)) %>%
    pivot_longer(names_to = "seqID", 
                 values_to = "Relabund", 
                 -Sample) %>%
    filter(Relabund>0)
  
  # List unique sequence IDs found in target t samples
  target_seqIDs_unique <- unique(target_seqIDs$seqID)
  
  # Write a list of sequence IDs found in target samples
  output_file_path <- paste("data/generated_files/sequences/seqID_lists/", t_ID, ".txt", sep="")
  
  write.table(target_seqIDs_unique, 
              output_file_path, 
              quote=F, 
              col.names=F,
              row.names=F)
  
  # Add file path to the list
  file_paths[[t_ID]] <- output_file_path
}

# Check if all output files were created
missing_files <- file_paths[!sapply(file_paths, file.exists)]

if (length(missing_files) > 0) {
  cat("\nERROR: The following files were not created:\n")
  cat(paste(names(missing_files), collapse = ", "), "\n")
} else {
  cat("\nDONE: The script has completed successfully.\n")
}

