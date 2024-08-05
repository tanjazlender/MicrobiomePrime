# This script formats the output from ThermonucleotideBLAST, which performs in silico PCR,
# into a structured table format for easier analysis and interpretation.

import re
import pandas as pd
import configparser
import glob
import sys

print('Reading variables.')

# Create a ConfigParser object
config = configparser.ConfigParser()

# Read the variables
config.read('../variables.ini')
sensitivity = config.get('settings', 'kmer_sensitivity_cutoff')
specificity = config.get('settings', 'kmer_specificity_cutoff')

# Check if target_group_name is empty; if yes, set target group name to target1
target_group_name = config.get('settings', 'target_group_name')
sys.stdout.flush()  # Force flush

if target_group_name != '':
    target = target_group_name.replace(' ', '-')
else:
    target = config.get('settings', 'target1').replace(' ', '-')

print('Rearranging ThermonucleotideBLAST output files into formatted tables.')
print(f'Target = {target}, kmer_sensitivity_cutoff = {sensitivity}, kmer_specificity_cutoff = {specificity}\n')
sys.stdout.flush()  # Force flush

# Function to extract information from the text
def extract_info(text):
    # Use regular expressions to extract relevant information
    SeqID = re.search(r'>([A-Za-z]+[0-9]+)', text)
    Amplicon_size = re.search(r'amplicon length = (\d+)', text)
    Amplicon_range = re.search(r'amplicon range = (-?\d+) \.\. (\d+)', text)
    PP_ID = re.search(r'name = (.+)', text)
    PrimerF = re.search(r'forward primer = 5\' (.+) 3\'', text).group(1)
    PrimerR = re.search(r'reverse primer = 5\' (.+) 3\'', text).group(1)
    TmF = re.search(r'forward primer tm = (.+)', text).group(1)
    TmR = re.search(r'reverse primer tm = (.+)', text).group(1)
    GC_F = re.search(r'forward primer %GC = (.+)', text).group(1)
    GC_R = re.search(r'reverse primer %GC = (.+)', text).group(1)
    MismatchF = int(re.search(r'forward primer mismatches = (\d+)', text).group(1))
    MismatchR = int(re.search(r'reverse primer mismatches = (\d+)', text).group(1))
    Clamp3F = int(re.search(r'min 3\' clamp = (\d+)', text).group(1))
    Clamp3R = int(re.search(r'max 3\' clamp = (\d+)', text).group(1))
    HeuristicsF = re.search(r'forward primer heuristics = (.+)', text).group(1)
    HeuristicsR = re.search(r'reverse primer heuristics = (.+)', text).group(1)
    Heterodimer_Tm = re.search(r'heterodimer tm = (\d+(\.\d+)?)', text).group(1)

    if all(match is not None for match in [SeqID, Amplicon_size, Amplicon_range, PP_ID, PrimerF, PrimerR, TmF, TmR, GC_F, GC_R, MismatchF, MismatchR, Clamp3F, Clamp3R, HeuristicsF, HeuristicsR, Heterodimer_Tm]):
        return {
            'SeqID': SeqID.group(1),
            'Amplicon_size': int(Amplicon_size.group(1)),
            'Amplicon_range': f"{Amplicon_range.group(1)}-{Amplicon_range.group(2)}",
            'PP_ID': PP_ID.group(1),
            'PrimerF': PrimerF,
            'PrimerR': PrimerR,
            'TmF': TmF,
            'TmR': TmR,
            'GC_F': GC_F,
            'GC_R': GC_R,
            'MismatchF': MismatchF,
            'MismatchR': MismatchR,
            'Clamp3F': Clamp3F,
            'Clamp3R': Clamp3R,
            'HeuristicsF': HeuristicsF,
            'HeuristicsR': HeuristicsR,
            'Heterodimer_Tm': Heterodimer_Tm
        }
    else:
        print("\nERROR: Failed to extract information from section:")
        sys.stdout.flush()  # Force flush
        print(text)
        return None

base_path = f'../../out/{target}/sens{sensitivity}_spec{specificity}'

file_prefix = f'{target}_sens{sensitivity}_spec{specificity}_results'

input_files = glob.glob(f'{base_path}/tntblast/{file_prefix}[0-9]*.txt')

print('Rearranging ThermonucleotideBLAST output into a table format.')

# Iterate through the list of input files
for input_file in input_files:
    
    # Extract the number from the input file name
    number_match = re.search(r'results(\d+)\.txt', input_file)
    if number_match:
        i = int(number_match.group(1))
    else:
        print(f"\nERROR: Unable to extract number from file name {file_prefix}{i}.txt\n")
        continue

    # Read the text file content
    with open(input_file, 'r') as file:
        # Initialize a list to store the extracted information for each primer pair
        primer_info_list = []

        # Initialize a flag to check if the current line is part of the section
        in_section = False

        # Iterate through each line in the file
        for line in file:
            if line.startswith("name ="):
                in_section = True
                section_text = line
            elif in_section:
                if line.strip() == "":
                    in_section = False
                    primer_info = extract_info(section_text)
                    if primer_info is not None:
                        primer_info_list.append(primer_info)
                else:
                    section_text += line

        # Create a DataFrame from the extracted information
        df = pd.DataFrame(primer_info_list)

        if not df.empty:
            # Generate output file name
            output_file = f'{base_path}/tntblast_tables/{file_prefix}_table{i}.txt'

            # Write the DataFrame to a text file
            df.to_csv(output_file, index=False, sep='\t')

            print(f'Formatted table for file no.{i} written to {file_prefix}_table{i}.txt')
            sys.stdout.flush()  # Force flush
        else:
            print(f'{file_prefix}{i}.txt is empty - table not written.')
            sys.stdout.flush()  # Force flush

print("\nDONE: The script has completed successfully.")
