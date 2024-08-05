# This script generates primers.
# Primers are actually K-mers that meet our sensitivity or both sensitivity and specificity criteria.

# Set the working directory to the parent folder (useful if running the script
# from its originating directory).
setwd("../../")

#### Script je potrebno še enkrat optimizirati, ko uredimo Lucijin del. Ideja je, da se k-meri generirajo le na sekvencah tarčne živali, prisotnost teh k-mer pa se preveri na celotni zbirki OTU-jev

# Ta del generira output, ki bi ga sicer želela od Lucije !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
library(dplyr)

all_kmers <- read.csv2("../../projects/2024-05-02.methodological.article/data/kmers/kmers_all/01-All-merged-DS10-noth-ALL-ss-40rule.csv", 
                       header = TRUE) %>%
  select(c("kmer", "zotus")) %>%
  rename(c("SeqIDs" = "zotus"))

library(tidyr)
all_kmers_split <- all_kmers %>%
  separate_rows(., SeqIDs, sep = ", ")

taxonomy <- read.csv("data/input_files/taxonomy.csv")

all_kmers_split_filt <- filter(all_kmers_split, SeqIDs %in% taxonomy$SeqID)

kmer_seqIDs <- all_kmers_split_filt %>%
  group_by(kmer) %>%
  summarise(seqID = paste0(SeqIDs, collapse = ", "))

################################################################################
######################## Load variables and parameters #########################
library(config)
library(ini)

# Read parameters from variables.ini
cat("Reading variables.\n")

variables <- read.ini("scripts/variables.ini")

# Add a target source (or a group of target sources)
target1 <- variables$settings$target1
if (is.null(target1)) {cat("ERROR: target1 needs to be defined!\n")}

target2 <- variables$settings$target2
if (is.null(target2)) {target2 <- "Not specified"}

target3 <- variables$settings$target3
if (is.null(target3)) {target3 <- "Not specified"}

target4 <- variables$settings$target4
if (is.null(target4)) {target4 <- "Not specified"}

target5 <- variables$settings$target5
if (is.null(target5)) {target5 <- "Not specified"}

# Add an additional source/sources you don't want to incorporate in your specificity calculations
specificity_exception1 <- variables$settings$specificity_exception1
if (is.null(specificity_exception1)) {specificity_exception1 <- "Not specified"}

specificity_exception2 <- variables$settings$specificity_exception2
if (is.null(specificity_exception2)) {specificity_exception2 <- "Not specified"}

specificity_exception3 <- variables$settings$specificity_exception3
if (is.null(specificity_exception3)) {specificity_exception3 <- "Not specified"}

specificity_exception4 <- variables$settings$specificity_exception4
if (is.null(specificity_exception4)) {specificity_exception4 <- "Not specified"}

specificity_exception5 <- variables$settings$specificity_exception5
if (is.null(specificity_exception5)) {specificity_exception5 <- "Not specified"}

# Add group name; this is necessary if you are looking for primer pairs associated with a group of sources
target_group_name <- variables$settings$target_group_name
target_group_ID <- gsub(" ", "-", fixed=TRUE, target_group_name)

# If you are looking for primer pairs associated with a single source, the group name will be the name of that source
# If you have multiple target sources, you MUST define a group name!
if (target2 == "Not specified" && 
    target3 == "Not specified" && 
    target4 == "Not specified" && 
    target5 == "Not specified") {
  target_group_name <- target1
  target_group_ID <- gsub(" ", "-", fixed=TRUE, target_group_name)
}

# Add other parameters
kmer_sensitivity_cutoff <- as.numeric(variables$settings$kmer_sensitivity_cutoff)
kmer_specificity_cutoff <- as.numeric(variables$settings$kmer_specificity_cutoff)

detach("package:config", unload = TRUE)
detach("package:ini", unload = TRUE)

################################################################################
################# Calculate K-mer sensitivity and specificity ##################
library(data.table)
library(tidyr)

cat("Reading input files.\n")

# Read metadata
if (file.exists("data/input_files/metadata.tsv")) {
  metadata <- read.delim("data/input_files/metadata.tsv", 
                         header = TRUE)
} else if (file.exists("data/input_files/metadata.csv")) {
  metadata <- read.csv("data/input_files/metadata.csv", 
                       header = TRUE)
} else {
  cat("ERROR: No metadata found in the input_files folder.\n")
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
  cat("ERROR: No relabund_tab found in the input_files folder.\n")
}

# Split and unnest the seqIDs column
kmer_seqIDs_dt <- as.data.table(kmer_seqIDs)

kmer_seqID <- kmer_seqIDs_dt[, .(kmers = rep(kmer, sapply(strsplit(seqID, ", "), length)), 
                                 seqID = unlist(strsplit(seqID, ", "))), by = 1:nrow(kmer_seqIDs_dt)]  %>%
  rename("kmer" = "kmers")

# Create a long table with relative abundances
relabund_tab_long <- data.frame(relabund_tab, Sample = rownames(relabund_tab)) %>%
  pivot_longer(names_to = "seqID", values_to = "relabund", -Sample) %>%
  left_join(select(metadata, c("Sample", "Source")))

# Group by seqIDs and collapse K-mers into a single string
seqID_kmers <- kmer_seqID[, .(kmers = paste(unique(kmer), collapse = ", ")), by = seqID]

############################## K-mer sensitivity ###############################
cat("\nCalculating K-mer sensitivity.\n")

metadata_target <- metadata %>%
  filter(Source == target1 | 
           Source == target2 | 
           Source == target3 | 
           Source == target4 | 
           Source == target5)

no_target_samples <- length(unique(metadata_target$Sample))

relabund_tab_long_target <- relabund_tab_long %>%
  filter(., Source %in% metadata_target$Source)

seqID_kmers_target <- filter(seqID_kmers, seqID %in% unique(relabund_tab_long_target$seqID))

sensitivity <- relabund_tab_long_target %>%
  filter(relabund>0) %>%
  left_join(seqID_kmers_target) %>%
  tidyr::separate_rows(kmers, sep = ", ") %>%
  rename("kmer" = "kmers") %>%
  group_by(kmer) %>%
  summarise(TP = length(unique(Sample)),
            FN = no_target_samples - TP,
            sensitivity = TP/(TP+FN)*100) %>%
  filter(sensitivity >= kmer_sensitivity_cutoff) %>%
  dplyr::arrange(desc(sensitivity)) %>%
  mutate(ID = 1:nrow(.))

# Check if the sensitivity data frame contains data and print an error if necessary
if (!is.null(sensitivity) && nrow(sensitivity) > 0) {
  cat("K-mer sensitivity calculation completed.\n")
} else {
  cat("ERROR: K-mer sensitivity calculation failed or resulted in an empty data frame.\n")
}

############################# K-mer specificity ################################
cat("\nCalculating K-mer specificity.\n")

metadata_nontarget <- metadata %>%
  filter(Source != target1 & 
           Source != target2 & 
           Source != target3 & 
           Source != target4 & 
           Source != target5 &
           Source != specificity_exception1 & 
           Source != specificity_exception2 & 
           Source != specificity_exception3 &
           Source != specificity_exception4 & 
           Source != specificity_exception5)

no_nontarget_samples <- length(unique(metadata_nontarget$Sample))

relabund_tab_long_nontarget <- relabund_tab_long %>%
  filter(., Source %in% metadata_nontarget$Source)

highly_sensitive_nontarget_kmers <- kmer_seqID %>%
  filter(seqID %in% unique(relabund_tab_long_nontarget$seqID),
         kmer %in% sensitivity$kmer) %>%
  group_by(seqID) %>%
  summarise(kmers = paste0(unique(kmer), collapse = ", "))

nontarget_PA <- highly_sensitive_nontarget_kmers %>%
  left_join(relabund_tab_long_nontarget) %>%
  tidyr::separate_rows(kmers, sep = ", ") %>%
  group_by(Sample, kmers, Source) %>%
  summarise(relabund = sum(relabund)) %>%
  rename("kmer" = "kmers") %>%
  mutate(PA = if_else(relabund>0, 1, 0))

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
  cat("ERROR: K-mer specificity calculation failed or resulted in an empty data frame.\n")
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
  cat("\nERROR: Primer FASTA files were NOT written.") 
}
