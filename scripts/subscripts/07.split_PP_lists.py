# This script divides the primer pairs file into multiple smaller files
# to simplify analysis and reduce the size of the output files.

import os
import configparser
import re
from pathlib import Path

# Change to the parent directory
os.chdir('..')

################################################################################
############################### Read variables #################################

# Read variables from the config file
config = configparser.ConfigParser()
config.read('scripts/variables.ini')

# Extract variables
sensitivity = config.get('settings', 'kmer_sensitivity_cutoff')
specificity = config.get('settings', 'kmer_specificity_cutoff')
target_group_name = config.get('settings', 'target_group_name')
target = config.get('settings', 'target')

# Process target
target_list = [t.strip() for t in target.split(',')]
target_combined = ' '.join(target_list).strip()

# Set target_group_ID
if target_group_name:
    target_group_ID = re.sub(r'\s+', '-', target_group_name.strip())
else:
    target_group_ID = re.sub(r'\s+', '-', target_combined)

################################################################################
########################## Divide the primer pair file #########################
print("Setting file paths.")

# Define the input file and output prefix
print("Reading the input file")

primers_folder = f"out/{target_group_ID}/sens{sensitivity}_spec{specificity}/primers"
input_file = f"{primers_folder}/original/primer_pairs_sens{sensitivity}_spec{specificity}.tsv"
output_path = f"{primers_folder}/split_primer_pairs"
output_file = f"{output_path}/{target_group_ID}_sens{sensitivity}_spec{specificity}_primers"

# Clear existing files in the output directory
print("Checking for and removing any existing files in the output directory.")
if os.path.isdir(output_path):
    for file_name in os.listdir(output_path):
        file_path = os.path.join(output_path, file_name)
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Warning: Could not delete file {file_path}: {e}")
else:
    print(f"Directory {output_path} does not exist.")

# Check if the input file exists
if not Path(input_file).is_file():
    print(f"\nERROR: Input file not found: {input_file}")
    exit(1)

# Delete existing files in the output directory
for file in Path(primers_folder).glob(f"{output_file}*.txt"):
    file.unlink()

# Set the maximum number of lines per split file
max_lines_per_file = 3000

# Split the file into smaller parts
print(f"Splitting primer pairs ({max_lines_per_file} lines/file).")

output_dir = Path(output_file).parent
output_dir.mkdir(parents=True, exist_ok=True)

def split_file(input_file, output_file, max_lines):
    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        total_lines = len(lines)
        num_files = (total_lines // max_lines) + (1 if total_lines % max_lines else 0)
        print(f"Total lines in input file: {total_lines}")
        for i in range(num_files):
            output_file_x = f"{output_file}{i+1}.txt"
            with open(output_file_x, 'w') as outfile:
                start = i * max_lines
                end = start + max_lines
                outfile.writelines(lines[start:end])

split_file(input_file, output_file, max_lines_per_file)

# Check if any split files were created
split_files = list(Path(output_dir).glob(f"{target_group_ID}_sens{sensitivity}_spec{specificity}_primers*.txt"))

if split_files:
    print(f"\nPrimer pairs were split into {len(split_files)} files.")
    print("\nDONE: The script has completed successfully.")
else:
    print("\nERROR: The primer pair file was not split into multiple smaller files.")
    exit(1)

