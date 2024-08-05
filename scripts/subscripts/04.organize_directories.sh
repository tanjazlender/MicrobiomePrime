# This script creates directories for storing analysis outputs and results.

#!/bin/bash

echo "Reading variables."

# Source variables
target_group_name=$(awk -F= '/target_group_name/ {gsub(/ /, "-", $2); print $2}' ../variables.ini)
sensitivity=$(awk -F= '/kmer_sensitivity_cutoff/ {print $2}' ../variables.ini)
specificity=$(awk -F= '/kmer_specificity_cutoff/ {print $2}' ../variables.ini)

target1=$(awk -F= '/target1/ {print $2}' ../variables.ini)
target2=$(awk -F= '/target2/ {print $2}' ../variables.ini)
target3=$(awk -F= '/target3/ {print $2}' ../variables.ini)
target4=$(awk -F= '/target4/ {print $2}' ../variables.ini)
target5=$(awk -F= '/target5/ {print $2}' ../variables.ini)

# Check if target_group_name is empty
if [ -z "$target_group_name" ]; then
  if [ -z "$target2" ] && [ -z "$target3" ] && [ -z "$target4" ] && [ -z "$target5" ]; then
  # Set target_group_name to target1
  target_group_name=$target1
  else
  echo "Error: Target group name not found. For analyzing multiple target sources, you need to set target_group_name in the scripts/variables.ini file."
  exit 1
  fi
fi

insilico_directory="../../out/${target_group_name}/sens${sensitivity}_spec${specificity}"

echo "Creating directories."

if [ -n "$target_group_name" ]; then
    # Create directories
    mkdir -p "$insilico_directory/marker_tables"
    mkdir -p "$insilico_directory/tntblast"
    mkdir -p "$insilico_directory/final_results"
    mkdir -p "$insilico_directory/tntblast_tables"
    mkdir -p "$insilico_directory/detected_sequences"
    mkdir -p "$insilico_directory/primers/split_primer_pairs"
    mkdir -p "$insilico_directory/primers/original"
    mkdir -p "$insilico_directory/log_files"
    mkdir -p "$insilico_directory/variables_settings"

    # Copy variables.ini to the primers/original directory
    cp ../variables.ini "$insilico_directory/variables_settings/variables_copy.ini"
    cp ../settings.ini "$insilico_directory/variables_settings/settings_copy.ini"
    
    # Check if all directories and the file exist
    if [ -d "$insilico_directory/marker_tables" ] &&
       [ -d "$insilico_directory/tntblast" ] &&
       [ -d "$insilico_directory/final_results" ] &&
       [ -d "$insilico_directory/tntblast_tables" ] &&
       [ -d "$insilico_directory/detected_sequences" ] &&
       [ -d "$insilico_directory/primers/split_primer_pairs" ] &&
       [ -d "$insilico_directory/primers/original" ] &&
       [ -d "$insilico_directory/log_files" ]; then
        echo "Directories for the in-silico PCR targeting $target_group_name created."
        printf "\nDONE: The script has completed successfully.\n"
    else
        printf "\nERROR: Failed to create all necessary directories.\n"
    fi
fi
