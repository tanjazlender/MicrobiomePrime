# This script computes sensitivity, specificity, and additional metrics
# for each primer pair.
# Tables including all relevant metrics are generated for primers that meet our 
# sensitivity and specificity criteria.

# Set the working directory to the parent folder (useful if running the script
# from its originating directory).
setwd("../../")

################################################################################
######################## Load variables and parameters #########################
cat("Reading variables.\n")

library(config)
library(ini)

# Read parameters from variables.ini
variables <- read.ini("scripts/variables.ini")
settings <- read.ini("scripts/settings.ini")

# Access parameters
# Add a target source (or a group of target sources)
target1 <- variables$settings$target1
if (is.null(target1)) {
  stop("\nERROR: target1 needs to be defined!\n")
  }

target2 <- variables$settings$target2
if (is.null(target2)) {
  target2 <- "Not specified"
  }

target3 <- variables$settings$target3
if (is.null(target3)) {
  target3 <- "Not specified"
  }

target4 <- variables$settings$target4
if (is.null(target4)) {
  target4 <- "Not specified"
  }

target5 <- variables$settings$target5
if (is.null(target5)) {
  target5 <- "Not specified"
  }

# Add an additional source/sources you don't want to incorporate in your specificity calculations
specificity_exception1 <- variables$settings$specificity_exception1
if (is.null(specificity_exception1)) {
  specificity_exception1 <- "Not specified"
  }

specificity_exception2 <- variables$settings$specificity_exception2
if (is.null(specificity_exception2)) {
  specificity_exception2 <- "Not specified"
  }

specificity_exception3 <- variables$settings$specificity_exception3
if (is.null(specificity_exception3)) {
  specificity_exception3 <- "Not specified"
  }

specificity_exception4 <- variables$settings$specificity_exception4
if (is.null(specificity_exception4)) {
  specificity_exception4 <- "Not specified"
  }

specificity_exception5 <- variables$settings$specificity_exception5
if (is.null(specificity_exception5)) {
  specificity_exception5 <- "Not specified"
  }

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
marker_sensitivity_cutoff <- as.numeric(variables$settings$marker_sensitivity_cutoff)
marker_specificity_cutoff <- as.numeric(variables$settings$marker_specificity_cutoff)
min_amplicon_length <- as.numeric(settings$settings$min_amplicon_length)

detach("package:config", unload = TRUE)
detach("package:ini", unload = TRUE)

################################################################################
############################## Read input files ################################
cat("Reading input files.\n")

library(dplyr)
library(tidyr)
library(stringr)
library(tibble)


# Read relabund_tab and pivot into relabund_tab_long
if (file.exists("data/input_files/relabund_tab.tsv")) {
  relabund_tab <- read.delim("data/input_files/relabund_tab.tsv", 
                                header = TRUE, 
                                row.names = 1)
} else if (file.exists("data/input_files/relabund_tab.csv")) {
  relabund_tab <- read.csv("data/input_files/relabund_tab.csv", 
                              header = TRUE, 
                              row.names = 1)
} else {
  stop("\nERROR: No relabund_tab found in the input_files folder.\n")
}

relabund_tab_long <- pivot_longer(data.frame(Sample = rownames(relabund_tab), relabund_tab), 
                            names_to = "SeqID", 
                            values_to = "Relabund", 
                            -Sample) %>%
  mutate(Percent_abundance = Relabund*100) %>%
  select(-Relabund)

# Read metadata
if (file.exists("data/input_files/metadata.tsv")) {
  metadata <- read.delim("data/input_files/metadata.tsv", 
                         header = TRUE)
} else if (file.exists("data/input_files/metadata.csv")) {
  metadata <- read.csv("data/input_files/metadata.csv", 
                       header = TRUE)
} else {
  stop("\nERROR: No metadata found in the input_files folder.\n")
}

cat("Calculating the number of samples.\n")

# Calculate the number of samples per source
nsamples_per_source <- metadata %>%
  group_by(Source) %>%
  summarise(Nsamples = length(unique(Sample)))

# Calculate the number of target samples
nsamples_target <- nsamples_per_source %>%
  filter(Source == target1 |
           Source == target2 |
           Source == target3 |
           Source == target4 |
           Source == target5) %>%
  summarise(sum(Nsamples)) %>%
  pull()

# Calculate the number of non-target samples
nontarget_samples <- filter(metadata, 
                            Source != target1 &
                              Source != target2 &
                              Source != target3 &
                              Source != target4 &
                              Source != target5 &
                              Source != specificity_exception1 &
                              Source != specificity_exception2 &
                              Source != specificity_exception3 &
                              Source != specificity_exception4 &
                              Source != specificity_exception5)

nsamples_nontarget <- length(unique(nontarget_samples$Sample))

cat("Preparing data for sensitivity and specificity calculations.\n")

# Join relabund table and metadata
relabund_tab_long_metadata <-  relabund_tab_long %>%
  left_join(select(metadata, c("Sample", "Source")))

# Find target sequence IDs and the samples in which they were found
seqIDs_samples_target <- relabund_tab_long_metadata %>%
  filter(Source == target1 | 
           Source == target2 | 
           Source == target3 | 
           Source == target4 | 
           Source == target5) %>%
  group_by(SeqID, Source) %>%
  summarise(samples = paste0(Sample, collapse = ", "))

# Find nontarget sequences and the samples in which they were found
seqIDs_samples_nontarget <- relabund_tab_long_metadata %>%
  filter(Source != target1 & 
           Source != target2 & 
           Source != target3 & 
           Source != target4 & 
           Source != target5 &
           Source != specificity_exception1 & 
           Source != specificity_exception2 & 
           Source != specificity_exception3 &
           Source != specificity_exception4 & 
           Source != specificity_exception5) %>%
  group_by(SeqID, Source) %>%
  summarise(Samples = paste0(Sample, collapse = ", "))

# Find sequence IDs of samples excluded from specificity calculations
seqIDs_samples_exceptions <- relabund_tab_long_metadata %>%
  filter(Source == specificity_exception1 |
           Source == specificity_exception2 | 
           Source == specificity_exception3 |
           Source == specificity_exception4 | 
           Source == specificity_exception5) %>%
  group_by(SeqID, Source) %>%
  summarise(Samples = paste0(Sample, collapse = ", "))

################################################################################
##################### Load and rearrange taxonomy info #########################
# Read taxonomy data
if (file.exists("data/input_files/taxonomy.tsv")) {
  taxonomy <- read.delim("data/input_files/taxonomy.tsv", 
                                header = TRUE)
} else if (file.exists("data/input_files/taxonomy.csv")) {
  taxonomy <- read.csv("data/input_files/taxonomy.csv", 
                              header = TRUE)
} else {
  stop("\nERROR: Taxonomy file was not found in the input_files folder.\n")
}

# Function to find the last non-empty value (taxa) in each row
get_last_taxon <- function(row) {
  non_empty_values <- row[row != ""]
  if (length(non_empty_values) > 0) {
    last_value <- tail(non_empty_values, 1)
    last_value_col <- names(row)[tail(which(row == last_value), 1)]
    
    if (last_value_col != "Genus") {
      last_value <- paste0(last_value, " (unknown)")
    }
    
    return(last_value)
  } else {
    return(NA)
  }
}

# Apply the function to each row and remove unnecesary columns
taxonomy$Last_taxon <- apply(taxonomy, 1, get_last_taxon)
taxonomy_last_taxon <- taxonomy[, c("SeqID", "Last_taxon")]

################################################################################
################# Loop through partial tntblast_results files ##################
cat("Computing sensitivity, specificity, and other key metrics for each primer pair.\n")

# Set the path to the directory of the input files
input_directory_path <- paste0("out/", target_group_ID, 
                               "/sens", kmer_sensitivity_cutoff, 
                               "_spec", kmer_specificity_cutoff, 
                               "/tntblast_tables/", 
                               sep = "")

# Set input file name pattern
input_name_pattern <- paste0(target_group_ID, 
                        "_sens", kmer_sensitivity_cutoff, 
                        "_spec", kmer_specificity_cutoff, 
                        "_results_table", 
                        "\\d+.txt", 
                        sep = "")

# Start_index 
# Change only, if you want to start the loop from a certain file in the file_list, not from the beginning
start_index <- 1

# Get the list of files in the directory
file_list <- list.files(input_directory_path, pattern = input_name_pattern)

file_list_length <- length(file_list)

if (start_index < 1 || start_index > file_list_length) {
  stop("Invalid start index.")
}

file_list <- file_list[start_index:file_list_length]

counter <- 0

# Loop through the files and print the corresponding numbers
for (filename in file_list){
  
  counter <- counter+1
  
  cat(paste0("\nFile ", counter, "/", file_list_length, "\n"))
  
  # Extract numbers between "table" and ".txt" to get the file number
  file_number <- as.numeric(sub(".*table(\\d+)\\.txt", "\\1", filename))
  cat(paste0("Analysing file ", filename, "\n"))
  
  ############ Load and reorganize ThermonucleotideBLAST results ##############
  tntblast_results <- read.table(paste0(input_directory_path, target_group_ID, 
                                        "_sens", kmer_sensitivity_cutoff, 
                                        "_spec", kmer_specificity_cutoff, 
                                        "_results_table", 
                                        file_number,".txt"), 
                                 sep ="\t", header = TRUE) %>%
    filter(Amplicon_size >= min_amplicon_length)
  
  # Skip this iteration if the tntblast_results data frame is empty
  if (nrow(tntblast_results) == 0) {
    cat(paste0(filename, "does not contain primer pairs meeting the specificied criteria.\n"))
    next  # Skip to the next iteration of the loop
  }
  
  ####################### Rearrange primer information #########################
  
  # Presence of sequences in target and nontarget samples
  seqID_presence <- relabund_tab_long_metadata %>%
    mutate(Target_nontarget = if_else(Source == target1 |
                                        Source == target2 |
                                        Source == target3 |
                                        Source == target4 |
                                        Source == target5, "T", "NT")) %>%
    group_by(SeqID) %>%
    summarise(Target_nontarget = paste0(unique(Target_nontarget), 
                                        collapse = ",")) %>%
    mutate(Target_nontarget = str_replace(Target_nontarget, "NT,T", "T,NT"))
  
  # Check whether given primer pairs amplify target and/or nontarget samples
  pp_presence <- tntblast_results %>%
    left_join(seqID_presence) %>%
    separate_rows(Target_nontarget, sep = ",") %>%
    group_by(PP_ID, PrimerF, PrimerR) %>%
    summarise(Target_nontarget = paste0(unique(Target_nontarget), 
                                        collapse = ",")) %>%
    mutate(Target_nontarget = str_replace(Target_nontarget, "NT,T", "T,NT"))
  
  # Check amplicon sizes for each primer pair
  amplicon_sizes <- select(tntblast_results, c("PP_ID", "Amplicon_size", "SeqID")) %>%
    left_join(pp_presence) %>%
    separate_rows(Target_nontarget, sep = ",") %>%
    group_by(PP_ID, Target_nontarget) %>%
    summarise(Amplicon_sizes = paste0(sort(unique(Amplicon_size)), 
                                      collapse = ", ")) %>%
    pivot_wider(names_from = Target_nontarget, values_from = Amplicon_sizes)
  
  # Rename columns
  if ("T" %in% colnames(amplicon_sizes)) {
    colnames(amplicon_sizes)[colnames(amplicon_sizes) == "T"] <- "Amplicon_sizes_target"
  } else {
    amplicon_sizes$Amplicon_sizes_target <- NA
  }
  
  if ("NT" %in% colnames(amplicon_sizes)) {
    colnames(amplicon_sizes)[colnames(amplicon_sizes) == "NT"] <- "Amplicon_sizes_nontarget"
  } else {
    amplicon_sizes$Amplicon_sizes_nontarget <- NA
  }

  # Write the range of Tm for each primer (each primer can anneal to its template with a non-perfect match, but has lower Tm)
  primersF_target <- select(tntblast_results, c("PrimerF", 
                                                "TmF", 
                                                "MismatchF", 
                                                "HeuristicsF", 
                                                "Heterodimer_Tm", 
                                                "PP_ID", 
                                                "SeqID")) %>%
    mutate(TmF = round(TmF, digits = 2)) %>%
    left_join(pp_presence) %>%
    separate_rows(Target_nontarget, sep = ",") %>%
    group_by(PP_ID, PrimerF, Target_nontarget) %>%
    summarise(TmF = paste0(sort(unique(TmF)), collapse = ", "),
              MismatchF = paste0("[", min(MismatchF), 
                                 "-", max(MismatchF), "]", 
                                 collapse = "")) %>%
    pivot_wider(names_from = Target_nontarget, 
                values_from = c("TmF", "MismatchF")) 
  
  # Rename columns
  if ("TmF_T" %in% colnames(primersF_target)) {
    colnames(primersF_target)[colnames(primersF_target) == "TmF_T"] <- "TmF_target"
    colnames(primersF_target)[colnames(primersF_target) == "MismatchF_T"] <- "MismatchF_target"
  } else {
    primersF_target$TmF_target <- NA
    primersF_target$MismatchF_target <- NA
  }
  
  
  if ("TmF_NT" %in% colnames(primersF_target)) {
    colnames(primersF_target)[colnames(primersF_target) == "TmF_NT"] <- "TmF_nontarget"
    colnames(primersF_target)[colnames(primersF_target) == "MismatchF_NT"] <- "MismatchF_nontarget"
  } else {
    primersF_target$TmF_nontarget <- NA
    primersF_target$MismatchF_nontarget <- NA
  }
  
  # Find primer melting temperatures and heuristics
  primersF <- select(tntblast_results, c("PrimerF", 
                                         "TmF", 
                                         "MismatchF", 
                                         "HeuristicsF", 
                                         "Heterodimer_Tm")) %>%
    group_by(PrimerF) %>%
    summarise(TmF_max = round(max(TmF), digits=2),
              HeuristicsF = paste0(unique(HeuristicsF), 
                                   collapse = "/")) %>%
    left_join(primersF_target)
  
  primersR_target <- select(tntblast_results, c("PrimerR", 
                                                "TmR", 
                                                "MismatchR", 
                                                "HeuristicsR", 
                                                "Heterodimer_Tm", 
                                                "PP_ID", 
                                                "SeqID")) %>%
    mutate(TmR = round(TmR, digits = 2)) %>%
    left_join(pp_presence) %>%
    separate_rows(Target_nontarget, sep = ",") %>%
    group_by(PP_ID, PrimerR, Target_nontarget) %>%
    summarise(TmR = paste0(sort(unique(TmR)), 
                           collapse = ", "),
              MismatchR = paste0("[", min(MismatchR), 
                                 "-", max(MismatchR), 
                                 collapse = "", "]")) %>%
    pivot_wider(names_from = Target_nontarget, values_from = c("TmR", "MismatchR")) 
    
    
  # Rename columns
    if ("TmR_T" %in% colnames(primersR_target)) {
    colnames(primersR_target)[colnames(primersR_target) == "TmR_T"] <- "TmR_target"
    colnames(primersR_target)[colnames(primersR_target) == "MismatchR_T"] <- "MismatchR_target"
  } else {
    primersR_target$TmR_target <- NA
    primersR_target$MismatchR_target <- NA
  }
  
  if ("TmR_NT" %in% colnames(primersR_target)) {
    colnames(primersR_target)[colnames(primersR_target) == "TmR_NT"] <- "TmR_nontarget"
    colnames(primersR_target)[colnames(primersR_target) == "MismatchR_NT"] <- "MismatchR_nontarget"
  } else {
    primersR_target$TmR_nontarget <- NA
    primersR_target$MismatchR_nontarget <- NA
  }
  
  primersR <- select(tntblast_results, c("PrimerR", 
                                         "TmR", 
                                         "MismatchR", 
                                         "HeuristicsR", 
                                         "Heterodimer_Tm")) %>%
    group_by(PrimerR) %>%
    summarise(TmR_max = round(max(TmR), digits=2),
              HeuristicsR = paste0(unique(HeuristicsR), 
                                   collapse = "/")) %>%
    left_join(primersR_target)
  
  # Create the final data frame with primer information
  primers_info <- select(tntblast_results, c("PP_ID", "PrimerF", "PrimerR")) %>%
    distinct() %>%
    left_join(primersF) %>%
    left_join(primersR)  %>%
    # Remove rows where either forward or reverse primer have mismatches to all target samples
    filter(grepl("^\\[0", MismatchF_target)) %>%
    filter(grepl("^\\[0", MismatchR_target))
  
  # Skip this iteration if the primers_info data frame is empty
  if (nrow(primers_info) == 0) {
    cat(paste0(filename, "does not contain primer pairs meeting the specificied criteria.\n"))
    next  # Skip to the next iteration of the loop
  }
  
  # Filter tntblast_results table - remove primer pairs that have mismatches on target DNA
  tntblast_results_filt <- filter(tntblast_results, PP_ID %in% primers_info$PP_ID)
  
  ##################### Calculate primer pair sensitivity ######################
  # Calculate the relative abundance of each marker within each target sample
  abund_long_target <- select(tntblast_results_filt, c("SeqID", PP_ID)) %>%
    group_by(SeqID) %>%
    summarise(PP_ID = paste0(PP_ID, collapse = ", ")) %>%
    filter(., SeqID %in% seqIDs_samples_target$SeqID) %>%
    left_join(seqIDs_samples_target) %>%
    separate_rows(PP_ID, sep = ", ") %>%
    separate_rows(samples, sep = ", ") %>%
    rename(Sample = samples) %>%
    left_join(relabund_tab_long_metadata) %>% 
    left_join(taxonomy_last_taxon) %>%
    left_join(nsamples_per_source) %>%
    group_by(PP_ID, Source, Nsamples, Sample) %>%
    summarise(SeqIDs = paste0(SeqID, collapse = ", "), 
              Percent_abundance = round(sum(Percent_abundance), digits=4),
              Taxonomy_target = paste0(unique(Last_taxon), 
                                       collapse = ", ")) %>%
    ungroup()
  
  # Skip this iteration if the abund_long_target data frame is empty
  if (nrow(abund_long_target) == 0) {
    cat(paste0(filename, "does not contain primer pairs meeting the specificied criteria.\n"))
    next  # Skip to the next iteration of the loop
  }
  
  # Calculate the sensitivity of each primer pair
  pp_sensitivity <- abund_long_target %>%
    filter(Percent_abundance != 0) %>%
    group_by(PP_ID) %>%
    summarise(N_positive_target_samples = length(unique(Sample))) %>%
    mutate(Sensitivity = round(N_positive_target_samples/nsamples_target*100, digits = 2),
           Sensitivity2 = paste0(N_positive_target_samples, "/", nsamples_target, sep = "")) %>%
    select(-c(N_positive_target_samples)) %>%
    filter(Sensitivity>=marker_sensitivity_cutoff)
  
  # Skip this iteration if the pp_sensitivity data frame is empty
  if (nrow(pp_sensitivity) == 0) {
    cat(paste0(filename, "does not contain primer pairs meeting the specificied criteria.\n"))
    next  # Skip to the next iteration of the loop
  }
  
  ##################### Calculate primer pair specificity ######################
  # Calculate the relative abundance of each chosen marker within each nontarget sample
  abund_long_nontarget <- select(tntblast_results_filt, c("SeqID", PP_ID)) %>% 
    filter(PP_ID %in% pp_sensitivity$PP_ID) %>% # Keep only primer pairs that are sensitive enough in target samples
    group_by(SeqID) %>%
    summarise(PP_ID = paste0(PP_ID, collapse = ", ")) %>%
    filter(., SeqID %in% seqIDs_samples_nontarget$SeqID) %>% # Keep only sequence IDs found in nontarget samples
    left_join(seqIDs_samples_nontarget) %>%
    separate_rows(PP_ID, sep = ", ") %>%
    separate_rows(Samples, sep = ", ") %>%
    rename(Sample = Samples) %>%
    left_join(relabund_tab_long_metadata) %>% # Add number of reads per sample per sequence
    left_join(taxonomy_last_taxon) %>%
    group_by(PP_ID, Source, Sample) %>%
    summarise(SeqIDs = paste0(SeqID, collapse = ", "), 
              Percent_abundance = round(sum(Percent_abundance), digits=4),
              Taxonomy_nontarget = paste0(unique(Last_taxon), 
                                          collapse = ", ")) %>%
    ungroup()
  
  # Prepare data for specificity calculations
  pp_samples_nontarget <- abund_long_nontarget %>%
    group_by(PP_ID) %>%
    summarise(Sample = paste0(unique(Sample), collapse = ", "))
  
  pp_nontarget_df <- pp_sensitivity %>%
    left_join(pp_samples_nontarget) %>%
    separate_rows(Sample, sep = ", ") %>%
    left_join(abund_long_nontarget)
  
  # Calculate primer pair specificity
  pp_specificity_sensitivity <- pp_nontarget_df %>%
    mutate(PA = if_else(!is.na(Percent_abundance) & Percent_abundance != 0, 1, 0)) %>%
    group_by(PP_ID, Sensitivity, Sensitivity2) %>%
    summarise(N_positive_nontarget_samples = sum(PA),
              TN = nsamples_nontarget-N_positive_nontarget_samples,
              Specificity = round((TN)/nsamples_nontarget*100, digits = 2),
              Specificity2 = paste0(TN, "/", nsamples_nontarget)) %>%
    select(-c("TN")) %>%
    filter(Specificity >= marker_specificity_cutoff) %>%
    ungroup()
  
  # Skip this iteration if the pp_sensitivity data frame is empty
  if (nrow(pp_specificity_sensitivity) == 0) {
    cat(paste0(filename, " does not contain primer pairs meeting the specificied criteria.\n"))
    next  # Skip to the next iteration of the loop
  }
  
  ##################### Calculate and write out other parameters ######################
  
  pp_sensitivity_detailed <- abund_long_target %>%
    filter(Percent_abundance != 0) %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    group_by(Source, PP_ID, Nsamples) %>%
    summarise(N_positive_target_samples = length(unique(Sample)),
              Positive_target_samples = paste0(unique(Sample), 
                                               collapse = ", ")) %>%
    mutate(Sensitivity2_detailed = paste0(Source, " (", 
                                          N_positive_target_samples, "/", 
                                          Nsamples, ")")) %>%
    group_by(PP_ID) %>%
    summarise(Sensitivity2_detailed = paste0(Sensitivity2_detailed, 
                                             collapse = ", "),
              Positive_target_samples = paste0(Source, " (", 
                                               Positive_target_samples, ")", 
                                               collapse = ", "))
  
  # Calculate relative abundances of markers
  pp_abundance_target <- abund_long_target %>%
    filter(Percent_abundance > 0) %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    group_by(PP_ID) %>%
    summarise(Percent_abundance2 = round(mean(Percent_abundance), digits = 4), 
              Percent_abundance2_SD = round(sd(Percent_abundance), digits = 4),
              Percent_abundance_target = paste0(Percent_abundance2, "+/-", 
                                                Percent_abundance2_SD)) %>%
    select(-c("Percent_abundance2", "Percent_abundance2_SD"))
  
  pp_abundance_target_detailed <- abund_long_target %>%
    filter(Percent_abundance > 0) %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    group_by(PP_ID, Source, Nsamples) %>%
    summarise(Percent_abundance2 = round(mean(Percent_abundance), digits = 4), 
              Percent_abundance2_SD = round(sd(Percent_abundance), digits = 4)) %>%
    group_by(PP_ID) %>%
    summarise(Percent_abundance_target_detailed = paste0(Source, " (", Percent_abundance2, "+/-", 
                                                         Percent_abundance2_SD, ")", 
                                                         collapse = ", "))
  
  # Assign taxonomy to each marker
  pp_taxonomy_target <- abund_long_target %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    separate_rows(Taxonomy_target, sep = ", ") %>%
    group_by(PP_ID) %>%
    summarise(Taxonomy_target = paste0(unique(sort(Taxonomy_target)), 
                                       collapse = ", "))
  
  # Find sequences not detected by the given primer pairs
  pp_negative_target_samples <- abund_long_target %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    filter(Percent_abundance == 0) %>%
    group_by(PP_ID, Source) %>%
    summarise(Negative_target_samples = paste0(Sample, collapse = ", ")) %>%
    group_by(PP_ID) %>%
    summarise(Negative_target_samples = paste0(Source, " (", Negative_target_samples, ")", collapse = ", "))
  
  # Call the garbage collector to free up unused memory
  gc()
  
  # Calculate the number of samples detected within using the given primer pairs
  pp_nontarget_df_filt <- pp_nontarget_df %>%
    filter(Percent_abundance>0) %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID)
  
  pp_presence_nontarget <- pp_nontarget_df_filt %>%
    left_join(nsamples_per_source) %>%
    group_by(PP_ID, Source, Nsamples) %>%
    summarise(N_positive_nontarget_samples = length(unique(Sample)),
              Positive_nontarget_samples = paste0(unique(Sample), collapse = ", ")) %>%
    group_by(PP_ID) %>%
    summarise(Presence_nontarget_samples = paste0(Source, " (", 
                                                  N_positive_nontarget_samples, "/", 
                                                  Nsamples, ")", 
                                                  collapse = ", "),
              Positive_nontarget_samples = paste0(Source, " (", 
                                                  Positive_nontarget_samples, ")", 
                                                  collapse = ", "))   
  
  # Calculate the relative abundance [%] of markers detected using given primer pairs
  pp_abundance_nontarget <- pp_nontarget_df_filt %>%
    group_by(PP_ID) %>%
    summarise(Percent_abund = round(mean(Percent_abundance), digits = 4), 
              Percent_abund_SD = round(sd(Percent_abundance), digits = 4),
              Percent_abundance_nontarget = paste0(Percent_abund, "+/-", 
                                                   Percent_abund_SD)) %>%
    select(-c("Percent_abund", "Percent_abund_SD"))
  
  
  pp_abundance_nontarget_detailed <- pp_nontarget_df_filt %>%
    group_by(PP_ID, Source) %>%
    summarise(Percent_abund_nontarget = round(mean(Percent_abundance), digits = 4), 
              Percent_abund_nontarget_SD = round(sd(Percent_abundance), digits = 4)) %>%
    group_by(PP_ID) %>%
    summarise(Percent_abundance_nontarget_detailed = paste0(Source, " (", 
                                                            Percent_abund_nontarget, "+/-", 
                                                            Percent_abund_nontarget_SD, ")", 
                                                            collapse = ", "))

  # Taxonomy of detected nontarget sequences
  pp_taxonomy_nontarget <- select(pp_nontarget_df_filt, c("PP_ID", "Taxonomy_nontarget")) %>%
    separate_rows(Taxonomy_nontarget, sep = ", ") %>%
    group_by(PP_ID) %>%
    summarise(Taxonomy_nontarget = paste0(unique(sort(Taxonomy_nontarget)), 
                                          collapse = ", "))
  
  # Join data frames
  pp_other_parameters_join <- select(pp_specificity_sensitivity, c("PP_ID")) %>%
    left_join(pp_sensitivity_detailed) %>%
    left_join(pp_abundance_target) %>%
    left_join(pp_abundance_target_detailed) %>%
    left_join(pp_negative_target_samples) %>%
    left_join(pp_taxonomy_target) %>%
    left_join(pp_presence_nontarget) %>%
    left_join(pp_abundance_nontarget) %>%
    left_join(pp_abundance_nontarget_detailed) %>%
    left_join(pp_taxonomy_nontarget) 
  
  #################### Specificity exceptions calculations #####################
  # Specificity exeptions are sources (e.g. species or groups of samples) that were excluded from specificity calculations
  # This is useful for example if you are trying to find swan-associated markers, but some of your samples are of unidentified birds
  # In this case you would treat unidentified birds as specificity exceptions (because they might or might not be swans)

  # Write all specificity exceptions if there are any
  specificity_exception_values <- c(specificity_exception1, 
                                    specificity_exception2, 
                                    specificity_exception3, 
                                    specificity_exception4, 
                                    specificity_exception5)
  
  if (any(specificity_exception_values != "Not specified")) {
    specificity_exceptions <- specificity_exception_values[specificity_exception_values != "Not specified"]
  } else {
    specificity_exceptions <- NA
  }
  
  # Calculate marker sensitivity for exceptions (excluded from specificity calculations) and add the info to the joined data frame
  if (any(specificity_exception_values != "Not specified")){
    
    # Calculate the relative abundance of each marker within each sample excluded from specificity calculations
    abund_long_exceptions <- select(tntblast_results_filt, c("SeqID", PP_ID)) %>%
      filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
      group_by(SeqID) %>%
      summarise(PP_ID = paste0(PP_ID, collapse = ", ")) %>%
      filter(., SeqID %in% seqIDs_samples_exceptions$SeqID) %>%
      left_join(seqIDs_samples_exceptions) %>%
      separate_rows(PP_ID, sep = ", ") %>%
      separate_rows(Samples, sep = ", ") %>%
      rename(Sample = Samples) %>%
      left_join(relabund_tab_long_metadata) %>% 
      left_join(taxonomy_last_taxon) %>%
      left_join(nsamples_per_source) %>%
      group_by(PP_ID, Source, Nsamples, Sample) %>%
      summarise(SeqIDs = paste0(SeqID, collapse = ", "), 
                Percent_abundance = round(sum(Percent_abundance), digits=4),
                Taxonomy_exceptions = paste0(unique(Last_taxon), 
                                             collapse = ", ")) %>%
      ungroup()
    
    # See which samples excluded from specificity calculations were detected using the given primer pairs
    pp_presence_exceptions <- abund_long_exceptions %>%
      filter(Percent_abundance > 0) %>%
      filter(Source %in% specificity_exceptions) %>%
      left_join(nsamples_per_source) %>%
      group_by(Source, PP_ID, Nsamples) %>%
      summarise(N_positive_samples = length(unique(Sample)),
                Positive_exceptions_samples = paste0(Sample, collapse = ", ")) %>%
      mutate(Sensitivity2_detailed = paste0(Source, " (", N_positive_samples, "/", 
                                            Nsamples, ")")) %>%
      group_by(PP_ID) %>%
      summarise(Presence_exceptions_samples = paste0(Sensitivity2_detailed, collapse = ", "),
                Positive_exceptions_samples = paste0(Source, " (", 
                                                     Positive_exceptions_samples, ")", 
                                                     collapse = ", ")) 
    
    # Calculate relative abundances of markers within samples excluded from specificity calculations
    pp_abundance_exceptions <- abund_long_exceptions %>%
      filter(Percent_abundance > 0) %>%
      filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
      group_by(PP_ID, Source, Nsamples) %>%
      summarise(Percent_abundance_exceptions = round(mean(Percent_abundance), digits = 4), 
                Percent_abundance_exceptions_SD = round(sd(Percent_abundance), digits = 4)) %>%
      group_by(PP_ID) %>%
      summarise(Percent_abundance_exceptions_detailed = paste0(Source, " (", 
                                                               Percent_abundance_exceptions, "+", 
                                                               Percent_abundance_exceptions_SD, ")", 
                                                               collapse = ", "))
    
    # Assign taxonomy to each marker found in specificity exception samples
    pp_taxonomy_exceptions <- abund_long_exceptions %>%
      filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
      separate_rows(Taxonomy_exceptions, sep = ", ") %>%
      group_by(PP_ID) %>%
      summarise(Taxonomy_exceptions = paste0(unique(sort(Taxonomy_exceptions)), collapse = ", "))
    
    # Create a joined table
    pp_exceptions_join <- select(pp_specificity_sensitivity, c("PP_ID")) %>%
      mutate(Exceptions = paste0(specificity_exceptions)) %>%
      left_join(pp_presence_exceptions) %>%
      left_join(pp_abundance_exceptions) %>%
      left_join(pp_taxonomy_exceptions)
    
  } else {
    
    pp_exceptions_join <- select(pp_specificity_sensitivity, c("PP_ID")) %>%
      mutate(Exceptions = NA,
             Presence_exceptions_samples = NA,
             Positive_exceptions_samples = NA,
             Percent_abundance_exceptions_detailed = NA,
             Taxonomy_exceptions = NA)
  }
  
  # Call the garbage collector to free up unused memory
  gc()
  
  ############# Write sequence IDs detected by best primer pairs ###############
  # Write out all detected sequences
  pp_all_detected_seqIDs <- tntblast_results_filt %>%
    filter(., PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    group_by(PP_ID) %>%
    summarise(SeqIDs_all = paste0(SeqID, collapse = ", "))
  
  # Write sequence IDs that are detected by best primer pairs in target samples
  pp_seqIDs_target_separated <- select(abund_long_target, c("PP_ID", "SeqIDs")) %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    separate_rows(SeqIDs, sep = ", ") %>%
    distinct(PP_ID, SeqIDs)
  
  pp_seqIDs_target <- pp_seqIDs_target_separated %>% 
    group_by(PP_ID) %>%
    summarise(SeqIDs_target = paste(SeqIDs[order(as.numeric(str_extract(SeqIDs, "\\d+")))], collapse = ", "),
              N_seqIDs_target = length(unique(SeqIDs)))
  
  pp_seqIDs_nontarget_separated <- select(abund_long_nontarget, c("PP_ID", "SeqIDs")) %>%
    filter(PP_ID %in% pp_specificity_sensitivity$PP_ID) %>%
    separate_rows(SeqIDs, sep = ", ")
  
  pp_seqIDs_nontarget <- pp_seqIDs_nontarget_separated %>% 
    distinct(PP_ID, SeqIDs) %>%
    group_by(PP_ID) %>%
    summarise(SeqIDs_nontarget = paste(SeqIDs[order(as.numeric(str_extract(SeqIDs, "\\d+")))], collapse = ", "),
              N_seqIDs_nontarget = length(unique(SeqIDs)))
  
  # Perform anti join to find sequence IDs detected in target but not in nontarget samples
  pp_seqIDs_target_only <- pp_seqIDs_target_separated %>%
    anti_join(pp_seqIDs_nontarget_separated, by = c("PP_ID", "SeqIDs")) %>%
    group_by(PP_ID) %>%
    summarise(seqIDs_target_only = paste(SeqIDs, collapse = ", "))
  
  # Perform anti join to find sequence IDs detected in nontarget but not in target samples
  pp_seqIDs_nontarget_only <- pp_seqIDs_nontarget_separated %>%
    anti_join(pp_seqIDs_target_separated, by = c("PP_ID", "SeqIDs")) %>%
    group_by(PP_ID) %>%
    summarise(seqIDs_nontarget_only = paste(SeqIDs, collapse = ", "))
  
  # Join all info
  pp_seqIDs_all <- select(pp_specificity_sensitivity, c("PP_ID")) %>%
    left_join(pp_all_detected_seqIDs) %>%
    left_join(pp_seqIDs_target) %>%
    left_join(pp_seqIDs_nontarget) %>%
    left_join(pp_seqIDs_target_only) %>%
    left_join(pp_seqIDs_nontarget_only) %>%
    select("PP_ID", "SeqIDs_all", "seqIDs_target_only", "seqIDs_nontarget_only", "N_seqIDs_target", "N_seqIDs_nontarget") %>%
    mutate(File_number = file_number)
  
  ######################### Join all info and save it ##########################
  
  # Final results table
  full_table <- pp_specificity_sensitivity %>%
    left_join(pp_other_parameters_join) %>%
    left_join(pp_exceptions_join) %>%
    left_join(primers_info) %>%
    left_join(amplicon_sizes) %>%
    left_join(select(pp_seqIDs_all, c("PP_ID", "N_seqIDs_target"))) %>%
    select("PP_ID", "Specificity", "Specificity2", 
           "Sensitivity", "Sensitivity2", "Sensitivity2_detailed", 
           "Presence_nontarget_samples",  
           "Percent_abundance_target", "Percent_abundance_nontarget", 
           "Percent_abundance_target_detailed", "Percent_abundance_nontarget_detailed", 
           "Taxonomy_target", "Taxonomy_nontarget", "N_seqIDs_target",
           "Positive_target_samples", "Positive_nontarget_samples", "Negative_target_samples",
           "Exceptions", "Presence_exceptions_samples", "Percent_abundance_exceptions_detailed", 
           "Taxonomy_exceptions", "Positive_exceptions_samples", 
           "PrimerF", "PrimerR", 
           "TmF_max", "TmR_max", 
           "TmF_target", "TmF_nontarget", 
           "TmR_target", "TmR_nontarget", 
           "MismatchF_target", "MismatchF_nontarget", 
           "MismatchR_target", "MismatchR_nontarget", 
           "Amplicon_sizes_target", "Amplicon_sizes_nontarget", 
           "HeuristicsF", "HeuristicsR")
  
  # Write the full_table
  output_directory_path <- paste0("out/", target_group_ID, 
                                  "/sens", kmer_sensitivity_cutoff, 
                                  "_spec", kmer_specificity_cutoff, 
                                  "/")
  
  marker_table_output_name_prefix <- paste0(target_group_ID,
                               "_msens", marker_sensitivity_cutoff, 
                               "_mspec", marker_specificity_cutoff)
  
  write.table(full_table, file = paste0(output_directory_path, "marker_tables/", 
                                        marker_table_output_name_prefix, "_markers", file_number, ".tsv", sep=""), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE,
              fileEncoding = "UTF-8")
  
  # Write a table with all sequence IDs detected by selected primer pairs
  seqID_table_output_name_prefix <- paste0(target_group_ID, "_", 
                                            "msens", marker_sensitivity_cutoff, 
                                            "_mspec", marker_specificity_cutoff)
  
  write.table(pp_seqIDs_all, file = paste0(output_directory_path, "detected_sequences/", 
                                         seqID_table_output_name_prefix, "_seqIDs", file_number, ".tsv", sep=""), 
              quote = FALSE, 
              sep = "\t", 
              row.names = FALSE)
  
  # If the full table is not empty, print a success message, otherwise print a warning
  if (nrow(full_table) > 0) {
    # If full_table is not empty, print success message
    cat(paste0("Tables ", marker_table_output_name_prefix, " and ", seqID_table_output_name_prefix, "were written.\n"))
  } else {
    # If full_table is empty, print warning message
    cat(paste0(filename, " does not contain primer pairs meeting the specified criteria.\n"))
  }
  
  # Call the garbage collector to free up unused memory
  gc()
}

# List the contents of the folder
folder_contents <- list.files(paste0(output_directory_path, "marker_tables/"))

# Check if the folder is empty and print a message
if(length(folder_contents) == 0) {
  # Print a message if no markers were found
  cat("\nThe script has completed. No markers were found.\n")
} else {
  # Print the contents of the folder
  cat("\nDONE: The script has completed successfully.")
}

