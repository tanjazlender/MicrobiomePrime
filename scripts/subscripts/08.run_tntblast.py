# This script performs in silico PCR on sequences from the 'sequence_singleline.fa' file
# using the primer pairs generated by our primer design process.

import os
import re
import configparser
import subprocess
from pathlib import Path

# Change to the parent directory
os.chdir('..') 

################################################################################
############################### Read variables #################################
print("Reading variables.")

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

max_amplicon_length = config.get('settings', 'max_amplicon_length')
min_primer_tm = config.get('settings', 'min_primer_tm')
max_primer_tm = config.get('settings', 'max_primer_tm')
max_primer_delta = config.get('settings', 'max_primer_delta')
min_primer_delta = config.get('settings', 'min_primer_delta')
primer_clamp = config.get('settings', 'primer_clamp')
max_mismatch = config.get('settings', 'max_mismatch')

config.read('scripts/settings.ini')
tntblast_path = config.get('settings', 'tntblast_path')
cpus = config.get('settings', 'cpus')

# Ensure the TNTBLAST executable exists
if not Path(tntblast_path).is_file():
    print(f"ERROR: ThermonucleotideBLAST executable not found at {tntblast_path}")
    exit(1)

print(f"Target group ID = {target_group_ID}")


################################################################################
############################ Perform in-silico PCR #############################

# Set file paths and read the FASTA file
print("Setting file paths.")

file_path = f"out/{target_group_ID}/sens{sensitivity}_spec{specificity}"
singleline_fasta = "data/generated_files/sequences/fasta_files/sequences_singleline.fa"

# Set a prefix for the output file
prefix = f"{target_group_ID}_sens{sensitivity}_spec{specificity}"

# Create the output directory if it does not exist
output_dir = Path(file_path) / 'tntblast'
try:
    output_dir.mkdir(parents=True, exist_ok=True)
except Exception as e:
    print(f"ERROR: Failed to create output directory {output_dir}. Details: {e}")
    exit(1)
    
# Remove any existing files in the output directory
print("Clearing existing files in the output directory.")
if output_dir.is_dir():
    for file in output_dir.iterdir():
        if file.is_file():
            try:
                file.unlink()
            except Exception as e:
                print(f"Error deleting {file}: {e}")

# Construct optional parameters
print("Reading optional PCR parameters.")

optional_parameters = ""
if max_amplicon_length:
    optional_parameters += f" -l {max_amplicon_length}"
if min_primer_tm:
    optional_parameters += f" -e {min_primer_tm}"
if max_primer_tm:
    optional_parameters += f" -x {max_primer_tm}"
if max_primer_delta:
    optional_parameters += f" -g {max_primer_delta}"
if min_primer_delta:
    optional_parameters += f" -z {min_primer_delta}"
if primer_clamp:
    optional_parameters += f" --primer-clamp {primer_clamp}"
if max_mismatch:
    optional_parameters += f" --max-mismatch {max_mismatch}"

# Iterate over input files
input_pattern = Path(file_path) / 'primers' / 'split_primer_pairs' / f"{prefix}_primers*.txt"
split_files = list((Path(file_path) / 'primers' / 'split_primer_pairs').glob('*'))

# Check if any primer files were found
if len(split_files) == 0:
    print(f"ERROR: No primer files found in {input_pattern.parent}")
    exit(1)

print(f"Found {len(split_files)} primer files.")

# Initialize the success flag
successfully_written = False

# Run ThermonucleotideBLAST on the split primer pair files
for input_file in split_files:
    number = re.search(r'\d+', input_file.stem.split('_primers')[-1]).group()
    output_file = output_dir / f"{prefix}_results{number}.txt"

    # Construct the command
    cmd = [
        'mpirun', '-np', cpus,
        str(tntblast_path),
        '-i', str(input_file),
        '-o', str(output_file),
        '-D', singleline_fasta
    ]
    if optional_parameters:
        cmd.extend(optional_parameters.split())
        
    
    # Print command for debugging
    print(f"Running command: {' '.join(cmd)}")
    
    # Run the command
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout.decode())
        print(f"Finished processing {input_file}. Results saved to {output_file}.")
        successfully_written = True  # Set the flag if at least one file was processed successfully

    except subprocess.CalledProcessError as e:
        print(f"ERROR: Running in silico PCR failed on input file {input_file}")
        print(f"Details: {e.stderr.decode()}")
        exit(1)


if successfully_written:
  print("\nDONE: The script has completed successfully.")