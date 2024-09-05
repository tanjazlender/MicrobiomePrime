# This script generates primers.
# Primers are actually K-mers that meet our sensitivity or both sensitivity and specificity criteria.

# Set the working directory to the parent folder (useful if running the script
# from its originating directory).
setwd("../")

################################################################################
############################### Read variables #################################
library(config)
library(ini)

# Read parameters from variables.ini
cat("Reading variables.\n")
variables <- read.ini("scripts/variables.ini")

kmer_sensitivity_cutoff <- as.numeric(variables$settings$kmer_sensitivity_cutoff)
kmer_specificity_cutoff <- as.numeric(variables$settings$kmer_specificity_cutoff)
target_list <- strsplit(variables$settings$target, ",")
target <- trimws(unlist(target_list))
target_combined <- paste(target, collapse = " ")
target_group_name <- variables$settings$target_group_name

# Set target_group_ID
if (is.null(target_group_name)  || !nzchar(target_group_name)) {
  target_group_ID <- gsub(" ", "-", fixed=TRUE, target_combined)
} else {
  target_group_ID <- gsub(" ", "-", fixed=TRUE, target_group_name)
}

# Extract specificity_exception_raw, defaulting to an empty string if it does not exist or is NULL
specificity_exception_raw <- if (!is.null(variables$settings$specificity_exception)) {
  variables$settings$specificity_exception
} else {
  # Default to an empty string if not available
  ""
}

# Check if the string specificity_exceptionis not empty
if (nchar(specificity_exception_raw) > 0) {
  # Split by comma if there are any commas
  specificity_exception <- unlist(strsplit(specificity_exception_raw, ",\\s*"))
} else {
  # If the string is empty, set an empty character vector
  specificity_exception <- character(0)
}

detach("package:config", unload = TRUE)
detach("package:ini", unload = TRUE)

################################################################################
############################# Read input files #################################

library(data.table)
library(tidyr)
library(dplyr)

cat("Reading input files.\n")

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

# Check if all specificity exceptions are found in metadata and print an error if necessary
cat("Verifying if all specified specificity exceptions are present in the metadata file.\n")
exceptions_not_found <- specificity_exception[!specificity_exception %in% metadata$Source]

if (length(exceptions_not_found) > 0) {
  stop(paste0("These targets were not found in the metadata file: ", paste(exceptions_not_found, collapse = ", "), "\n"))
} else {
  cat("All defined targets are found in metadata.\n")
}

# Create an empty data.table to store the combined data
kmer_seqIDs <- data.table()

# Read target K-mers and their associated sequence IDs from all sources
for (target_source_x in target) {
  # Replace spaces with "-"
  target_source_x_ID <- gsub(" ", "-", target_source_x)
  
  # Define input file path
  input_file_path <- paste0("data/generated_files/kmers/", target_source_x_ID, "_kmers.csv")
  
  if (file.exists(input_file_path)) {
    # Read the input .csv file using data.table::fread for efficiency
    kmer_seqIDs_partial <- fread(input_file_path, sep = ";", header = TRUE)
    
    # Append the data table to the list
    kmer_seqIDs <- rbind(kmer_seqIDs, kmer_seqIDs_partial)
    
  } else {
    stop(paste0("A list of kmers for target ", target_source_x, 
                " (", input_file_path, ") does not exist. ",
                "The list has to be created and saved in file path: ",
                input_file_path))
  }
}

# Remove duplicated rows
kmer_seqIDs <- unique(kmer_seqIDs)


################################################################################
################# Calculate K-mer sensitivity and specificity ##################

# Separate seqIDs
kmer_seqID <- separate_rows(kmer_seqIDs, seqIDs, sep = ",") %>%
  rename("seqID" = "seqIDs") %>%
  as.data.table

# Group by seqIDs and collapse K-mers into a single string
seqID_kmers <- kmer_seqID[, .(kmers = paste(unique(kmer), collapse = ", ")), by = seqID]

# Create a long table with relative abundances
relabund_tab_long <- data.frame(relabund_tab, Sample = rownames(relabund_tab)) %>%
  pivot_longer(names_to = "seqID", values_to = "relabund", -Sample) %>%
  left_join(select(metadata, c("Sample", "Source")))



############################## K-mer sensitivity ###############################
cat("\nCalculating K-mer sensitivity.\n")

metadata_target <- metadata %>%
  filter(Source %in% target)

no_target_samples <- length(unique(metadata_target$Sample))

relabund_tab_long_target <- relabund_tab_long %>%
  filter(., Source %in% metadata_target$Source)

seqID_kmers_target <- filter(seqID_kmers, seqID %in% unique(relabund_tab_long_target$seqID))

sensitivity <- relabund_tab_long_target %>%
  filter(relabund>0) %>%
  left_join(seqID_kmers_target) %>%
  filter(kmers != "NA") %>%
  tidyr::separate_rows(kmers, sep = ", ") %>%
  rename("kmer" = "kmers") %>%
  group_by(kmer) %>%
  summarise(TP = length(unique(Sample)),
            FN = no_target_samples - TP,
            sensitivity = TP/(TP+FN)*100) %>%
  filter(sensitivity >= kmer_sensitivity_cutoff) %>%
  arrange(desc(sensitivity))

# Check if the sensitivity data frame contains data and print an error if necessary
if (!is.null(sensitivity) && nrow(sensitivity) > 0) {
  sensitivity <- sensitivity %>%
    mutate(ID = row_number())
  cat("K-mer sensitivity calculation completed.\n")
} else {
  stop("K-mer sensitivity calculation failed or resulted in an empty data frame.\n")
}

############################# K-mer specificity ################################
cat("\nCalculating K-mer specificity.\n")

metadata_nontarget <- metadata %>%
  filter(!Source %in% target &
           !Source %in% specificity_exception)

no_nontarget_samples <- length(unique(metadata_nontarget$Sample))

relabund_tab_long_nontarget <- relabund_tab_long %>%
  filter(., Source %in% metadata_nontarget$Source)

highly_sensitive_nontarget_kmers <- kmer_seqID %>%
  filter(seqID %in% unique(relabund_tab_long_nontarget$seqID),
         kmer %in% sensitivity$kmer) %>%
  group_by(seqID) %>%
  summarise(kmers = paste0(unique(kmer), collapse = ", "))

highly_sensitive_relabund <- highly_sensitive_nontarget_kmers %>%
  left_join(relabund_tab_long_nontarget) %>%
  tidyr::separate_rows(kmers, sep = ", ") 

library(data.table)

# Convert to data.table
setDT(highly_sensitive_relabund)

# Perform the operations
nontarget_PA <- highly_sensitive_relabund[, .(relabund = sum(relabund)), by = .(Sample, kmers, Source)
                                          ][, PA := as.integer(relabund > 0)][
                                            , .(kmer = kmers, relabund, PA)]

specificity <- nontarget_PA %>%
  group_by(kmer) %>%
  summarise(FP = sum(PA),
            TN = no_nontarget_samples - FP,
            specificity = TN/(TN+FP)*100) %>%
  filter(., specificity >= kmer_specificity_cutoff) %>%
  left_join(sensitivity)

# Check if the specificity data frame contains data and print an error if necessary
if (!is.null(specificity) && nrow(specificity) > 0) {
  cat("K-mer specificity calculation completed.\n")
} else {
  stop("K-mer specificity calculation failed or resulted in an empty data frame.\n")
}

################################################################################
#################### Generate forward and reverse primers ######################
# Generate forward and reverse primers. Forward primers are created based on 
# K-mer sensitivity or/and specificity thresholds, while reverse primers are 
# the reverse complements of the forward primers.

library(seqinr)

############################# Forward primers ##################################
cat("Generating forward primers.\n")
# Create a list of forward primers that satisfy the K-mer sensitivity threshold
sensitive_fw <- sensitivity$kmer

sensitive_fw_IDs <- paste(target_group_ID, 
                          sensitivity$ID, "Fw", 
                          sep = "")

# Create a list of forward primers that satisfy the K-mer sensitivity and specificity thresholds
sensitive_specific_fw <- specificity$kmer

sensitive_specific_fw_IDs <- paste(target_group_ID, 
                                   specificity$ID, "Fw", 
                                   sep = "")

############################# Reverse primers ##################################
cat("Generating reverse primers.\n")
# Create reverse primers (reverse complements of the forward primers)
reverse_complement <- function(seq) {
  complement <- chartr("ATGC", "TACG", seq)
  reverse_complement <- rev(strsplit(complement, NULL)[[1]])
  return(paste(reverse_complement, collapse = ""))
}

# Create a list of reverse primers that satisfy the K-mer sensitivity threshold
sensitive_rv <- sapply(sensitive_fw, reverse_complement, USE.NAMES = FALSE)

sensitive_rv_IDs <- paste(target_group_ID, 
                          sensitivity$ID, "Rv", 
                          sep = "")

# Create a list of reverse primers that satisfy the K-mer sensitivity and specificity thresholds
sensitive_specific_rv <- sapply(sensitive_specific_fw, reverse_complement, USE.NAMES = FALSE)

sensitive_specific_rv_IDs <- paste(target_group_ID, 
                                   specificity$ID, "Rv", 
                                   sep = "")

###################### Write primers into fasta files ##########################
cat("Writing the generated forward and reverse primers into FASTA files.\n")
fasta_path <- paste0("out/", target_group_ID, 
                     "/sens", kmer_sensitivity_cutoff, "_spec", 
                     kmer_specificity_cutoff, 
                     "/primers/original/")

write.fasta(sequences = as.list(sensitive_fw), 
            names = sensitive_fw_IDs, 
            open = "w", 
            file.out = paste0(fasta_path, 
                              "primers_sens", kmer_sensitivity_cutoff, 
                              "_fw.fa"))
                              
write.fasta(sequences = as.list(sensitive_specific_fw), 
            names = sensitive_specific_fw_IDs, 
            open = "w", 
            file.out = paste0(fasta_path, 
                              "primers_sens", kmer_sensitivity_cutoff,
                              "_spec", kmer_specificity_cutoff, 
                              "_fw.fa"))
                              
write.fasta(sequences = as.list(sensitive_rv), 
            names = sensitive_rv_IDs, 
            open = "w", 
            file.out = paste0(fasta_path, 
                              "primers_sens", kmer_sensitivity_cutoff,
                              "_rv.fa"))
                              
write.fasta(sequences = as.list(sensitive_specific_rv), 
            names = sensitive_specific_rv_IDs, 
            open = "w", 
            file.out = paste0(fasta_path, 
                              "primers_sens", kmer_sensitivity_cutoff,
                              "_spec", kmer_specificity_cutoff, 
                              "_rv.fa"))

# Check if the FASTA files were written and write a log file
file_fw <- paste0(fasta_path, 
                  "primers_sens", kmer_sensitivity_cutoff, 
                  "_fw.fa")

file_spec_fw <- paste0(fasta_path, 
                       "primers_sens", kmer_sensitivity_cutoff,
                       "_spec", kmer_specificity_cutoff, 
                       "_fw.fa")

file_rv <- paste0(fasta_path, 
                  "primers_sens", kmer_sensitivity_cutoff,
                  "_rv.fa")

file_spec_rv <- paste0(fasta_path, 
                       "primers_sens", kmer_sensitivity_cutoff,
                       "_spec", kmer_specificity_cutoff, 
                       "_rv.fa")

# Check if the FASTA files exist and print an error if necessary
if (all(file.exists(c(file_fw, file_spec_fw, file_rv, file_spec_rv)))) {
  cat("Primer FASTA files were written.\n")
  cat("\nDONE: The script has completed successfully.")
} else {
  stop("Primer FASTA files were NOT written.") 
}
