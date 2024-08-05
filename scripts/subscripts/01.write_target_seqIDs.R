# This script creates a list of Sequence IDs found in target samples.
# This list will be used in another script to extract FASTA sequences.

# Set the working directory to the parent folder (useful if running the script
# from its originating directory).
setwd("../../")

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
  cat("\nERROR: No relabund_tab found in the input_files folder.\n")
}

# Read metadata
if (file.exists("data/input_files/metadata.tsv")) {
  metadata <- read.delim("data/input_files/metadata.tsv", 
                         header = TRUE)
} else if (file.exists("data/input_files/metadata.csv")) {
  metadata <- read.csv("data/input_files/metadata.csv", 
                       header = TRUE)
} else {
  cat("\n: No metadata found in the input_files folder.\n")
}

################################################################################
############################### Read variables #################################
cat("Reading variables.\n")

library(config)
library(ini)

# Read parameters from variables.ini
variables <- read.ini("scripts/variables.ini")

# Access parameters
# Add a main source (or a group of main sources)
target1 <- variables$settings$target1
if (is.null(target1)) {stop("ERROR: target1 needs to be defined!\n")}

target2 <- variables$settings$target2
if (is.null(target2)) {target2 <- "Not specified"}

target3 <- variables$settings$target3
if (is.null(target3)) {target3 <- "Not specified"}

target4 <- variables$settings$target4
if (is.null(target4)) {target4 <- "Not specified"}

target5 <- variables$settings$target5
if (is.null(target5)) {target5 <- "Not specified"}

# Add group name; this is necessary if you are looking for primer pairs associated with a group of sources
target_group_name <- variables$settings$target_group_name
target_group_ID <- gsub(" ", "-", fixed=TRUE, target_group_name)

# If you are looking for primer pairs associated with a single source, the group name will be the name of that source
# If you have multiple target sources, you MUST define a group name!
if (target2 == "Not specified" && target3 == "Not specified" && target4 == "Not specified" && target5 == "Not specified") {
  target_group_name <- target1
  target_group_ID <- gsub(" ", "-", fixed=TRUE, target_group_name)
}

detach("package:config", unload = TRUE)
detach("package:ini", unload = TRUE)

################################################################################
######################## Write sequence ID lists ###############################
cat("Writing sequence IDs.\n")

# Create subdirectories for sequence ID list files and FASTA files to be generated
fasta_dir <- file.path("data/generated_files/sequences/fasta_files")
seqID_dir <- file.path("data/generated_files/sequences/seqID_lists")

if (!file.exists(fasta_dir)) {
  dir.create(fasta_dir, showWarnings = FALSE)
}

if (!file.exists(seqID_dir)) {
  dir.create(seqID_dir, showWarnings = FALSE)
}

# Filter metadata (keep only target samples)
metadata_target <- metadata %>%
  filter(Source == target1 | 
           Source == target2 | 
           Source == target3 | 
           Source == target4 | 
           Source == target5) 

# Filter relabund_tab (keep only target samples)
relabund_tab_target <- filter(relabund_tab, 
                              rownames(relabund_tab) %in% unique(metadata_target$Sample))

# Find out, which sequence IDs are detected in target samples
target_seqIDs <- relabund_tab_target %>%
  mutate(Sample = rownames(.)) %>%
  pivot_longer(names_to = "seqID", 
               values_to = "Relabund", 
               -Sample) %>%
  filter(Relabund>0)

# List unique sequence IDs found in target samples
target_seqIDs_unique <- unique(target_seqIDs$seqID)

# Write a list of sequence IDs found in target samples
output_file_path <- paste("data/generated_files/sequences/seqID_lists/", target_group_ID, ".txt", sep="")

write.table(target_seqIDs_unique, 
            output_file_path, 
            quote=F, 
            col.names=F,
            row.names=F)

# Check if the final output was written
if (file.exists(output_file_path)) {
  cat("\nDONE: The script has completed successfully.\n")
} else {
  cat("\n: The file", output_file_path, "was not created.")
}
