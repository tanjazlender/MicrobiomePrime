# This script generates primer pairs.
# One primer in each pair must meet both sensitivity and specificity criteria.
# The other primer must meet only the sensitivity criteria.

from Bio import SeqIO
import configparser
import os
import re

os.chdir('..')

################################################################################
################################ Read variables ################################
print("Reading variables.")

# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the variables
config.read('scripts/variables.ini')
sensitivity = config.get('settings', 'kmer_sensitivity_cutoff')
specificity = config.get('settings', 'kmer_specificity_cutoff')
target_group_name = config.get('settings', 'target_group_name')
target = config.get('settings', 'target')

target_list = [t.strip() for t in target.split(',')]
target_combined = ' '.join(target_list)

# Set target_group_ID
if target_group_name !='':
    target_group_ID = re.sub(r'\s+', '-', target_group_name)
else:
    target_group_ID = re.sub(r'\s+', '-', target_combined)

# Define paths
primers_folder = f"out/{target_group_ID}/sens{sensitivity}_spec{specificity}/primers/original"
combined_output_file_path = f"{primers_folder}/primer_pairs_sens{sensitivity}_spec{specificity}.tsv"

# Remove any existing files in the output folder
for file_name in os.listdir(primers_folder):
    if file_name.startswith("primer_pairs"):
        file_path = os.path.join(primers_folder, file_name)
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
                print(f"Deleted {file_path}")
            except Exception as e:
                print(f"Error deleting {file_path}: {e}")

################################################################################
############################## Create primer pairs #############################
print("Creating primer pairs.")

# Combine primers into primer pairs
def read_fasta(file_path):
    primers = {}
    with open(file_path, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            primers[record.id] = str(record.seq)
    return primers

def generate_combined_table(primers_fw, primers_rv, output_file):
    with open(output_file, "a") as output:
        for id_fw, seq_fw in primers_fw.items():
            for id_rv, seq_rv in primers_rv.items():
                pair = f"{id_fw}:{id_rv}"
                output.write(f"{pair}\t{seq_fw}\t{seq_rv}\n")

# Remove duplicates from the combined output file
def remove_duplicates(file_path):
    seen = set()
    with open(file_path, "r") as file:
        lines = file.readlines()
    
    with open(file_path, "w") as file:
        for line in lines:
            if line not in seen:
                file.write(line)
                seen.add(line)

# Write primer pairs into files
fw1_file_path = f"{primers_folder}/primers_sens{sensitivity}_spec{specificity}_fw.fa"
fw2_file_path = f"{primers_folder}/primers_sens{sensitivity}_fw.fa"
rv1_file_path = f"{primers_folder}/primers_sens{sensitivity}_rv.fa"
rv2_file_path = f"{primers_folder}/primers_sens{sensitivity}_spec{specificity}_rv.fa"
    
primers_fw1 = read_fasta(fw1_file_path)
primers_fw2 = read_fasta(fw2_file_path)
primers_rv1 = read_fasta(rv1_file_path)
primers_rv2 = read_fasta(rv2_file_path)

generate_combined_table(primers_fw1, primers_rv1, combined_output_file_path)
generate_combined_table(primers_fw2, primers_rv2, combined_output_file_path)

# Remove duplicates from the combined output file
print("Removing duplicate rows from the combined output file.")

remove_duplicates(combined_output_file_path)
    
if os.path.exists(combined_output_file_path) and os.path.getsize(combined_output_file_path) > 0:
        print(f"Primers for target source(s) {target_group_ID} were combined into primer pairs (fw1-rv1 and fw2-rv2) and written into file {combined_output_file_path}.")
        print("\nDONE: The script has completed successfully.")
else:
        print("\nERROR: The output file was not created or is empty.")
