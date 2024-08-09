# This script extracts target sequences from a FASTA file.
# The target sequences are extracted based on the list of sequence IDs.

from Bio import SeqIO
import configparser
import os

os.chdir('..')

#######################################################################
############ Create a FASTA file in a single line format ##############

# Define a function to convert a multiline FASTA to a singleline
def convert_fasta_to_singleline(input_fasta, output_fasta):
    # Check if the output file does not exist or if file sizes differ
    if not os.path.exists(output_fasta) or (os.path.exists(input_fasta) and os.path.getsize(input_fasta) != os.path.getsize(output_fasta)):
        # Open the output file in write mode
        with open(output_fasta, "w") as output_handle:
            # Iterate over each record in the input FASTA file
            for record in SeqIO.parse(input_fasta, "fasta"):
                # Write the sequence ID
                output_handle.write(f">{record.id}\n")
                # Write the sequence in a single line
                output_handle.write(str(record.seq).replace('\n', '') + "\n")
        # Check if the output file was written successfully and write the completion or error message
        if os.path.exists(output_fasta):
            print(f"Conversion of multiline to singleline FASTA completed.")
        else:
            print("\nERROR: Failed to write the singleline FASTA output file.")
            exit(1)
    else:
        print("Singleline FASTA file already exists.")

# Define input and output files
input_fasta = "data/input_files/sequences.fa"
output_fasta = "data/generated_files/sequences/fasta_files/sequences_singleline.fa"

# Convert to a singleline FASTA
convert_fasta_to_singleline(input_fasta, output_fasta)


#######################################################################
################### Extract source FASTA file #########################

# Read the variables
config = configparser.ConfigParser()
config.read('scripts/variables.ini')

target = config.get('settings', 'target')
target_list = [t.strip() for t in target.split(',')]

# Define fasta file path
singleline_fasta = "data/generated_files/sequences/fasta_files/sequences_singleline.fa"

# Generate a FASTA file for each target
for t in target_list:
    seqID_list = f"data/generated_files/sequences/seqID_lists/{t}.txt"

    # Check if the list of sequence IDs exists before proceeding
    if os.path.exists(seqID_list):
        # Define the output file path
        output_target_fasta = f"data/generated_files/sequences/fasta_files/{t}.fa"
        
        # Read the sequence IDs found in target samples
        with open(seqID_list, "r") as list_file:
            sequence_IDs = [line.strip() for line in list_file]
        
        # Extract sequences based on sequence IDs
        with open(output_target_fasta, "w") as output_file:
            for record in SeqIO.parse(singleline_fasta, "fasta"):
                if record.id in sequence_IDs:
                    # Write the sequence ID and sequence in a single line to the output file
                    output_file.write(f">{record.id}\n")
                    output_file.write(str(record.seq).replace('\n', '') + "\n")
        
        # Check if the output file was written successfully
        if os.path.exists(output_target_fasta) and os.path.getsize(output_target_fasta) > 0:
            print(f"FASTA file for source {t} written.")
        else:
            print(f"\nERROR: The output file for source {t} was not created or is empty.")
            exit(1)

    else:
        print(f"\nERROR: Sequence ID list for source {t} is not found in the data/generated_files/sequences/seqID_lists folder.")
        exit(1)
        
print("\nDONE: The script has completed successfully.")
