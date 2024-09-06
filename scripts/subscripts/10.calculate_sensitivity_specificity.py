import os
import gc
import sys
import pandas as pd
import psutil
from concurrent.futures import ProcessPoolExecutor

# Set up the environment
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, project_root)

# Ensure the current working directory is the project root
os.chdir(project_root)

# Import project-specific modules
from src.config import Config
import src.spec_sens_analysis as analysis
import src.data_processing as dp 
from src.utils import (
    clean_directory,
    check_if_empty,
    estimate_nworkers)

################################################################################
######################## Load variables and parameters #########################

print("Reading variables.")

# Create an instance of the Config class
config = Config()

kmer_sensitivity_cutoff = config.kmer_sensitivity_cutoff
kmer_specificity_cutoff = config.kmer_specificity_cutoff 
marker_sensitivity_cutoff = config.marker_sensitivity_cutoff 
marker_specificity_cutoff = config.marker_specificity_cutoff 
min_amplicon_length = config.min_amplicon_length 
target = config.target
target_group_ID = config.target_group_ID
specificity_exception = config.specificity_exception

################################################################################
############################## Read input files ################################

# Read and process relabund_tab
print("Reading relabund_tab.")
relabund_tab = dp.read_relabund_tab()
relabund_tab_long = dp.transform_relabund_tab(relabund_tab)

# Read metadata
print("Reading metadata.")
metadata = dp.read_metadata()

# Read taxonomy
print("Reading taxonomy.")
taxonomy = dp.read_taxonomy()

# Calculate the number of samples per source, number of target samples,
# nontarget samples and exception samples
print("Calculating the number of samples.")
nsamples_per_source, nsamples_target, nsamples_specificity_exception, nsamples_nontarget, target_samples, nontarget_samples, specificity_exception_samples = dp.calculate_nsamples(metadata, target, specificity_exception)

# Join relabund table and metadata
print("Joining relabund_table and metadata.")
relabund_tab_long_metadata = dp.relabund_metadata(relabund_tab_long, metadata)

# Check sequence IDs found in target samples, nontarget samples and specificity exceptions
seqIDs_samples_target, seqIDs_samples_nontarget, seqIDs_samples_exception = dp.seqIDs_samples(relabund_tab_long_metadata, target_samples, nontarget_samples, specificity_exception_samples)

################################################################################
####################### Set paths and read input files #########################
# Set paths
main_directory = f"out/{target_group_ID}/sens{kmer_sensitivity_cutoff}_spec{kmer_specificity_cutoff}"
input_directory = f"{main_directory}/tntblast_tables/"
output_directory_joined = f"{main_directory}/marker_tables/"
output_directory_seqs = f"{main_directory}/detected_sequences/"

# Remove existing files in the output directories
clean_directory(output_directory_joined)
clean_directory(output_directory_seqs)

#################################################################################
####################### Define the processing function #########################

def process_file(filename):
  
    try:  
    
        # Extract file number
        file_number = dp.extract_file_number(filename)
    
        # Read tntblast_results and remove amplicons with amplicon size < min_amplicon_length
        tntblast_results = pd.read_csv(os.path.join(input_directory, filename), sep='\t', on_bad_lines='warn')
        tntblast_results = tntblast_results[tntblast_results['Amplicon_size'] >= min_amplicon_length]
        
        # Skip to next iteration if tntblast_results is empty
        if check_if_empty(tntblast_results, filename):
            return  # Exit function if no data to process
    
        # Check whether given primer pairs amplify target and/or nontarget samples
        pp_presence = analysis.get_pp_presence(relabund_tab_long_metadata, tntblast_results, target)
    
        # Check amplicon sizes for each primer pair
        amplicon_sizes = analysis.get_amplicon_sizes(tntblast_results, pp_presence)
    
        # Rearrange primer pair information
        primers_info = analysis.get_primers_info(tntblast_results, pp_presence)
        
        # Skip to next iteration if primers_info is empty
        if check_if_empty(primers_info, filename):
            return  # Exit function if no data to process
        
        ########### Calculate sensitivity, specificity, and other metrics ###########
        # Filter tntblast_results
        tntblast_results_filt = analysis.filter_tntblast_results(tntblast_results, primers_info)
    
        # Create target, nontarget and exception abundance dataframes
        target_abundance = analysis.get_abundance_df(
            tntblast_results_filt, seqIDs_samples_target, relabund_tab_long_metadata, taxonomy, nsamples_per_source, sample_type='target')
        
        # Skip to next iteration if target_abundance is empty
        if check_if_empty(target_abundance, filename):
            return  # Exit function if no data to process
        
        nontarget_abundance = analysis.get_abundance_df(
            tntblast_results_filt, seqIDs_samples_nontarget, relabund_tab_long_metadata, taxonomy, nsamples_per_source, sample_type='nontarget')
        exception_abundance = analysis.get_abundance_df(
            tntblast_results_filt, seqIDs_samples_exception, relabund_tab_long_metadata, taxonomy, nsamples_per_source, sample_type='exception')
    
        # Calculate sensitivity and specificity of each primer pair        
        sensitivity_specificity = analysis.calculate_sensitivity_specificity(
            target_abundance, nontarget_abundance, nsamples_target, nsamples_nontarget, marker_sensitivity_cutoff, marker_specificity_cutoff)
        
        # Skip to next iteration if sensitivity_specificity is empty
        if check_if_empty(sensitivity_specificity, filename):
            return  # Exit function if no data to process
    
        sensitivity_detailed = analysis.get_sensitivity_detailed(target_abundance, sensitivity_specificity)
    
        # Calculate the percent abundance of each marker in target source
        abundance_target = analysis.calculate_abundance_target(target_abundance, sensitivity_specificity)
    
        # Taxonomy of target amplicons
        taxonomy_target = analysis.get_taxonomy_target(target_abundance, sensitivity_specificity)
    
        # Negative target samples
        negative_target_samples = analysis.process_negative_target_samples(target_abundance, sensitivity_specificity)
    
        # Detected nontarget samples
        detected_nontarget = analysis.process_detected_nontarget_samples(nontarget_abundance, sensitivity_specificity)
    
        # Get marker presence, abundance
        exceptions_info = analysis.process_exception_samples(exception_abundance, sensitivity_specificity, taxonomy, specificity_exception)
    
        pp_seqIDs = analysis.process_sequence_ids(tntblast_results_filt, sensitivity_specificity, target_abundance, nontarget_abundance, file_number)
    
        joined = analysis.join_info(
            sensitivity_specificity, sensitivity_detailed, amplicon_sizes, primers_info, abundance_target, taxonomy_target,
            negative_target_samples, detected_nontarget, exceptions_info, pp_seqIDs, file_number)
    
        ############################# Write output tables ########################## 
        # Define output file names
        output_name_prefix_joined = f"{target_group_ID}_msens{marker_sensitivity_cutoff}_mspec{marker_specificity_cutoff}"
        output_name_prefix_seqs = f"{target_group_ID}_msens{marker_sensitivity_cutoff}_mspec{marker_specificity_cutoff}"
    
        # Write the joined table to a TSV file
        joined_table_path = os.path.join(output_directory_joined, f"{output_name_prefix_joined}_markers{file_number}.tsv")
        joined.to_csv(joined_table_path, sep='\t', index=False, quoting=False, encoding='utf-8')
    
        # Write the pp_seqIDs to a TSV file
        seqID_table_path = os.path.join(output_directory_seqs, f"{output_name_prefix_seqs}_markers{file_number}.tsv")
        pp_seqIDs.to_csv(seqID_table_path, sep='\t', index=False, quoting=False)
    
        # Memory usage after processing
        mem_info_after = psutil.Process(os.getpid()).memory_info()
        
        print(f"Finished analysing file {filename}."
              f" Memory after: {mem_info_after.rss / (1024 * 1024):.2f} MB"
              f" Virtual memory: {mem_info_after.vms / (1024 * 1024):.2f} MB")
                     
    except Exception as e:
        print(f"An unexpected error accured when processing file {filename}: {str(e)}")
        # Re-raise the exception to stop further processing
        sys.exit(1)

#################################################################################
############################## Parallel Execution ##############################

# Estimate the number of workers (processes) based on available memory
nworkers = estimate_nworkers(relabund_tab, complexity_factor=1)

if __name__ == '__main__':
    file_list = sorted(os.listdir(input_directory), key=extract_file_number)
    
    # Set the number of workers, or let it default to the number of CPUs
    with ProcessPoolExecutor(nworkers) as executor:
        # Submit all tasks for parallel execution
        executor.map(process_file, file_list)
