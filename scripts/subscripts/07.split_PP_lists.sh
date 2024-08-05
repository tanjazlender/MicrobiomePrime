# This script divides the primer pairs file into multiple smaller files
# to simplify analysis and reduce the size of the output files.

#!/bin/bash

# Source variables
echo "Reading variables."

sensitivity=$(awk -F= '/kmer_sensitivity_cutoff/ {print $2}' ../variables.ini)
specificity=$(awk -F= '/kmer_specificity_cutoff/ {print $2}' ../variables.ini)

target_group_name=$(awk -F= '/target_group_name/ {gsub(/ /, "-", $2); print $2}' ../variables.ini)
target1=$(awk -F= '/target1/ {print $2}' ../variables.ini)
target2=$(awk -F= '/target2/ {print $2}' ../variables.ini)
target3=$(awk -F= '/target3/ {print $2}' ../variables.ini)
target4=$(awk -F= '/target4/ {print $2}' ../variables.ini)
target5=$(awk -F= '/target5/ {print $2}' ../variables.ini)

# Check if target_group_name is empty; set the target_group_name to target1 if necessary
if [ -z "$target_group_name" ] && [ -z "$target2" ] && [ -z "$target3" ] && [ -z "$target4" ] && [ -z "$target5" ]; then
  # Set target_group_name to target1
  target_group_name=$target1  
fi

# Define the input file and output prefix
echo "Reading the input file"

primers_folder="../../out/${target_group_name}/sens${sensitivity}_spec${specificity}/primers"
input_file="$primers_folder/original/primer_pairs_sens${sensitivity}_spec${specificity}.tsv"
output_prefix="$primers_folder/split_primer_pairs/${target_group_name}_sens${sensitivity}_spec${specificity}_primers"

# Check if the input file exists
if [ ! -f "$input_file" ]; then
  echo "\nERROR: Input file not found: $input_file"
  exit 1
fi

# Delete existing files in the output directory
for file in "${output_prefix}"*.txt; do
  rm -f "$file"
done

# Set the maximum number of lines per split file
max_lines_per_file=3000

# Split the file into smaller parts using awk
echo "Splitting primer pairs (${max_lines_per_file} lines/file)."

split -l $max_lines_per_file --additional-suffix=.txt "$input_file" "$output_prefix"

# Rename the split files with the desired numbering convention
counter=1
for split_file in "${output_prefix}"*.txt
do
  if [ "$split_file" != "$input_file" ]; then
    new_name="${output_prefix}${counter}.txt"
    mv "$split_file" "$new_name"
    ((counter++))
  fi
done

echo "Primer pairs were split into multiple files."

# Check if any split files were created
if ls "${output_prefix}"*.txt 1> /dev/null 2>&1; then
  echo -e "\nDONE: The script has completed successfully."
else
  echo -e "\nERROR: The primer pair file was not split into multiple smaller files."
fi
