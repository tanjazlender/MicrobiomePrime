# This file joins all results of a particular type into a single file

# Set the working directory to the parent folder (useful if running the script
# from its originating directory).
setwd("../../")

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

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
  cat("\nERROR: target1 needs to be defined!\n")
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

detach("package:config", unload = TRUE)
detach("package:ini", unload = TRUE)

################################################################################
#################### Loop through partial markers files ########################

host_directory_path <-  paste0("out/", target_group_ID, 
                               "/sens", kmer_sensitivity_cutoff, 
                               "_spec", kmer_specificity_cutoff, 
                               "/")

input_directory_path <- paste0(host_directory_path, 
                         "marker_tables", 
                         sep = "")

# Construct the regular expression pattern
pattern <- paste0(target_group_ID, 
                  "_msens", marker_sensitivity_cutoff, 
                  "_mspec", marker_specificity_cutoff, 
                  "_markers", "\\d+.tsv")

# Write out a list of files
file_list <- list.files(input_directory_path, 
                        pattern = pattern, 
                        full.names = TRUE)

# start_index (change only, if you want to start the loop from a certain file in the file_list, not from the beginning)
start_index <- 1

if (start_index < 1 || start_index > length(file_list)) {
  stop("Invalid start index.")
}

file_list <- file_list[start_index:length(file_list)]

# Set an empty log file
log_file <- paste0(host_directory_path, "log_files/", 
                   target_group_ID, 
                   "_msens", marker_sensitivity_cutoff, 
                   "_mspec", marker_specificity_cutoff, 
                   "_final_table.log")


# Initialize an empty data frame to store filtered results
final_results <- data.frame()

counter <- 0

# Loop through the marker files
for (filename in file_list){
  counter <- (counter+1)
  
  # Extract file numbers
  file_number <- as.numeric(sub(".*markers(\\d+)\\.tsv", "\\1", filename))
  
  # Read results
  results_table <- read.table(filename, header = TRUE, sep = "\t")
  
  detected_sequences <- read.table(paste0("out/", target_group_ID, 
                                          "/sens", kmer_sensitivity_cutoff, 
                                          "_spec", kmer_specificity_cutoff, 
                                          "/detected_sequences/",
                                          target_group_ID, 
                                          "_msens", marker_sensitivity_cutoff, 
                                          "_mspec", marker_specificity_cutoff, 
                                          "_seqIDs", file_number, ".tsv"), 
                                   header = TRUE, 
                                   sep = "\t")
  
  results_table$TmF_target <- as.character(results_table$TmF_target)
  results_table$TmR_target <- as.character(results_table$TmR_target)
  results_table$TmF_nontarget <- as.character(results_table$TmF_nontarget)
  results_table$TmR_nontarget <- as.character(results_table$TmR_nontarget)
  results_table$Amplicon_sizes_target <- as.character(results_table$Amplicon_sizes_target)
  results_table$Amplicon_sizes_nontarget <- as.character(results_table$Amplicon_sizes_nontarget)
  results_table$File_number <- file_number
  
  
  # Append filtered results to the final_results data frame
  final_results <- dplyr::bind_rows(final_results, results_table)
  
  cat(paste0("Finished processing file number ", file_number, 
               "(", round(counter/length(file_list)*100, digits = 2), 
               "% finished)\n"))
}  



# Define the output file path
output_file_path <- paste0(host_directory_path, 
                    "final_results/", target_group_ID, 
                    "_msens", marker_sensitivity_cutoff, 
                    "_mspec", marker_specificity_cutoff, 
                    "_final_results.tsv")

# Write the table
write.table(final_results, 
            output_file_path, 
            row.names = FALSE, 
            quote = FALSE,
            sep = "\t")


# Check if the final results table was written
if (file.exists(output_file_path)) {
  cat(paste0("The final results have been saved to: ", output_file_path, "\n"))
  cat("\nDONE: The script has completed successfully.")
} else {
  print(cat("\nERROR: The final results could not be written."))
}

