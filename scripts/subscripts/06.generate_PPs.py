# This script generates primer pairs.
# One primer in each pair must meet both sensitivity and specificity criteria.
# The other primer must meet only the sensitivity criteria.

from Bio import SeqIO
import configparser
import os

print("Reading variables.")

# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the variables
config.read('../variables.ini')
sensitivity = config.get('settings', 'kmer_sensitivity_cutoff')
specificity = config.get('settings', 'kmer_specificity_cutoff')

# Check if target_group_name is empty; if yes, set target name to target1
target_group_name = config.get('settings', 'target_group_name')

if target_group_name != '':
    target = target_group_name.replace(' ', '-')
else:
    target = config.get('settings', 'target1').replace(' ', '-')

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

# Write primer pairs into files
primers_folder = f"../../out/{target}/sens{sensitivity}_spec{specificity}/primers/original"

fw1_file_path = f"{primers_folder}/primers_sens{sensitivity}_spec{specificity}_fw.fa"
fw2_file_path = f"{primers_folder}/primers_sens{sensitivity}_fw.fa"
rv1_file_path = f"{primers_folder}/primers_sens{sensitivity}_rv.fa"
rv2_file_path = f"{primers_folder}/primers_sens{sensitivity}_spec{specificity}_rv.fa"
    
primers_fw1 = read_fasta(fw1_file_path)
primers_fw2 = read_fasta(fw2_file_path)
primers_rv1 = read_fasta(rv1_file_path)
primers_rv2 = read_fasta(rv2_file_path)

combined_output_file_path = f"{primers_folder}/primer_pairs_sens{sensitivity}_spec{specificity}.tsv"

generate_combined_table(primers_fw1, primers_rv1, combined_output_file_path)
generate_combined_table(primers_fw2, primers_rv2, combined_output_file_path)
    
if os.path.exists(combined_output_file_path) and os.path.getsize(combined_output_file_path) > 0:
        print(f"Primers for target source {target} were combined into primer pairs (fw1-rv1 and fw2-rv2) and written into file {combined_output_file_path}.")
        print("\nDONE: The script has completed successfully.")
else:
        print("\nERROR: The output file was not created or is empty.")
