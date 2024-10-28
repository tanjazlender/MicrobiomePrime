# This file joins all results of a particular type into a single file

# Set the working directory to the parent folder (useful if running the script
# from its originating directory).
setwd("../")

library(dplyr)
library(tibble)
library(tidyr)
library(stringr)

################################################################################
######################## Load variables and parameters #########################
cat("Reading variables.\n")

# Read parameters from variables.ini
cat("Reading variables.\n")
variables <- ini::read.ini("scripts/variables.ini")

kmer_sensitivity_cutoff <- as.numeric(variables$settings$kmer_sensitivity_cutoff)
kmer_specificity_cutoff <- as.numeric(variables$settings$kmer_specificity_cutoff)
marker_sensitivity_cutoff <- as.numeric(variables$settings$marker_sensitivity_cutoff)
marker_specificity_cutoff <- as.numeric(variables$settings$marker_specificity_cutoff)
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

################################################################################
################################ Join results ##################################

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

# Check if there are any input files and print an error if necessary
input_files <- list.files(path = paste0(input_directory_path), full.names = TRUE)

if (length(input_files) == 0) {
  stop("No primer pairs meeting the defined criteria found.")
}

# List existing files in the directory and delete them
existing_files <- list.files(path = paste0(host_directory_path, "/final_results"), full.names = TRUE)

cat("Checking for and removing any existing files in the output directory.\n")
if (length(existing_files) > 0) {
  # Delete files
  file.remove(existing_files)
}


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
file_list_length <- length(file_list)

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
  
  # Convert to character to ensure compatibility when joining tables
  results_table$TmF_target <- as.character(results_table$TmF_target)
  results_table$TmR_target <- as.character(results_table$TmR_target)
  results_table$TmF_nontarget <- as.character(results_table$TmF_nontarget)
  results_table$TmR_nontarget <- as.character(results_table$TmR_nontarget)
  results_table$Amplicon_sizes_target <- as.character(results_table$Amplicon_sizes_target)
  results_table$Amplicon_sizes_nontarget <- as.character(results_table$Amplicon_sizes_nontarget)
  results_table$Percent_abundance_target <- as.character(results_table$Percent_abundance_target)
  results_table$Percent_abundance_nontarget <- as.character(results_table$Percent_abundance_nontarget)
  results_table$Percent_abundance_target_detailed <- as.character(results_table$Percent_abundance_target_detailed)
  results_table$Percent_abundance_nontarget_detailed <- as.character(results_table$Percent_abundance_nontarget_detailed)
  
  # Append filtered results to the final_results data frame
  final_results <- dplyr::bind_rows(final_results, results_table)
  
  cat(paste0("Finished processing file number ", file_number, 
               " (progress:", counter, "/", file_list_length, ")\n"))
}  

# Define a function to replace empty strings with NA only for character columns
replace_empty_with_na <- function(df) {
  df %>%
    mutate(across(where(is.character), ~ifelse(.x == "", NA, .x)))
}

# Apply the function to your data frame
full_table <- replace_empty_with_na(final_results)


# Define the output file path
output_file_path <- paste0(host_directory_path, 
                    "final_results/", target_group_ID, 
                    "_msens", marker_sensitivity_cutoff, 
                    "_mspec", marker_specificity_cutoff, 
                    "_final_results.tsv")

# Write the table
write.table(full_table, 
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

