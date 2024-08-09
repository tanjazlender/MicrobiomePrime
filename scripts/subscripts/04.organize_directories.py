# This script creates directories for storing analysis outputs and results.

import os
import shutil
import configparser
import re

os.chdir('..')

################################################################################
############################### Read variables #################################

config = configparser.ConfigParser()
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
  
################################################################################
############################# Create directories ###############################

# Define a function of creating directories
def create_directories(base_path, dirs):
    # Create directories if they do not exist.
    # Ensure the base directory is created
    os.makedirs(base_path, exist_ok=True)
    print(f"Created base directory: {base_path}")
    
    for directory in dirs:
        path = os.path.join(base_path, directory)
        os.makedirs(path, exist_ok=True)
        print(f"Created directory: {path}")

# Define directories to create
insilico_directory = f"out/{target_group_ID}/sens{sensitivity}_spec{specificity}"
directories = [
  'marker_tables',
  'tntblast',
  'final_results',
  'tntblast_tables',
  'detected_sequences',
  'primers/split_primer_pairs',
  'primers/original',
  'log_files',
  'variables']    

print("Creating directories.")

# Create the directories
create_directories(insilico_directory, directories)

# Create a copy of ini files
shutil.copy('scripts/variables.ini', os.path.join(insilico_directory, 'variables/variables_copy.ini'))

# Verify that all directories were created
all_directories_exist = all(
    os.path.isdir(os.path.join(insilico_directory, directory))
    for directory in directories
)

if all_directories_exist:
    print("\nDONE: The script has completed successfully.")
else:
    print("\nERROR: Failed to create all necessary directories.")
    exit(1)

