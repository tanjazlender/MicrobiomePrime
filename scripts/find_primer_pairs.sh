#!/bin/bash

# Exit script if any command fails and handle errors in pipelines
set -e
set -o pipefail

# Define log file paths
current_log_path="current_logs"
default_log_path="${current_log_path}/find_primer_pairs.log"

# Create log directory if it does not exist
mkdir -p "$current_log_path"

# Redirect all output (stdout and stderr) to the default log file
exec 2> >(tee -a "$default_log_path")

# Define a function to run commands and check their exit status
run_and_check() {
    local cmd="$1"
    local log_file="$2"

    echo "Running: $cmd"
    stdbuf -oL -eL $cmd 2>&1 | tee "$log_file"
    if [ ${PIPESTATUS[0]} -ne 0 ]; then
        echo "Error occurred while executing: $cmd"
        exit 1
    fi
}

# Verify input formatting
run_and_check "Rscript subscripts/00.verify_input_formatting.R" "${current_log_path}/00.verify_input_formatting.log"

# Create a list of sequence IDs found in target samples
run_and_check "Rscript subscripts/01.write_target_seqIDs.R" "${current_log_path}/01.write_target_seqIDs.log"

# Extract FASTA sequences based on the list of sequence IDs
run_and_check "python -u subscripts/02.extract_fasta_files.py" "${current_log_path}/02.extract_fasta_files.log"

# Organize directories
run_and_check "python -u subscripts/03.organize_directories.py" "${current_log_path}/03.organize_directories.log"

# Extract valid kmers
run_and_check "python -u subscripts/04.extract_valid_kmers.py" "${current_log_path}/04.extract_valid_kmers.log"

# Generate primers
run_and_check "Rscript subscripts/05.generate_primers.R" "${current_log_path}/05.generate_primers.log"

# Generate primer pairs
run_and_check "python -u subscripts/06.generate_PPs.py" "${current_log_path}/06.generate_PPs.log"

# Split the list of primer pairs into multiple files (otherwise the tntblast output files are too large)
run_and_check "python -u subscripts/07.split_PP_lists.py" "${current_log_path}/07.split_PP_lists.log"

# Run ThermonucleotideBLAST
run_and_check "python -u subscripts/08.run_tntblast.py" "${current_log_path}/08.run_tntblast.log"

# Generate a table out of ThermonucleotideBLAST output
run_and_check "python -u subscripts/09.rearrange_tntblast_output.py" "${current_log_path}/09.rearrange_tntblast_output.log"

# Calculate the sensitivity and specificity of primer pairs (PPs)
run_and_check "python -u subscripts/10.calculate_sensitivity_specificity.py" "${current_log_path}/10.calculate_sensitivity_specificity.log"

# Generate the results table
run_and_check "Rscript subscripts/11.join_all_results.R" "${current_log_path}/11.join_all_results.log"

# Read the variables
config_file="variables.ini"

sensitivity=$(awk -F= '/kmer_sensitivity_cutoff/ {print $2}' "$config_file")
specificity=$(awk -F= '/kmer_specificity_cutoff/ {print $2}' "$config_file")
target_raw=$(awk -F= '/target/ {print $2}' "$config_file")
target_group_name=$(awk -F= '/target_group_name/ {gsub(/ /, "-", $2); print $2}' "$config_file")

# Trim leading and trailing whitespace from target_raw
target_raw=$(echo "$target_raw" | xargs)

# Process target
IFS=',' read -r -a target_list <<< "$target_raw"
target_combined=$(printf '%s ' "${target_list[@]}" | xargs)  # Join array elements with space and trim

# Debug: Print processed target
echo "Processed target list: ${target_list[@]}"
echo "Target combined: '$target_combined'"

# Set target_group_ID
if [[ -n "$target_group_name" ]]; then
    target_group_ID=$(echo "$target_group_name" | tr -s ' ' '-' | xargs)  # Replace whitespace with dashes and trim
else
    target_group_ID=$(echo "$target_combined" | tr -s ' ' '-' | xargs)  # Replace whitespace with dashes and trim
fi

echo "Target group ID: '$target_group_ID'"

# Define log file path and create the main output folder for a given target
final_log_path="../out/${target_group_ID}/sens${sensitivity}_spec${specificity}/log_files"

echo "$final_log_path"

# Ensure the final log directory exists
mkdir -p "$final_log_path"

# Move all files from current_log_path to final_log_path
for file in "$current_log_path"/*; do
    mv "$file" "$final_log_path"
    echo "Moved: $file to $final_log_path"
done
