import os
import pandas as pd
    
# Function to read input files
def read_relabund_tab():
    file_paths = [
        "data/input_files/relabund_tab.tsv",
        "data/input_files/relabund_tab.csv"]
    for file_path in file_paths:
        if os.path.exists(file_path):
            return pd.read_csv(file_path, sep='\t' if file_path.endswith('.tsv') else ',', index_col=0)
    raise FileNotFoundError(f"No valid file found in the paths: {file_paths}")
    
def read_metadata():
    # Define paths inside the function
    file_paths = [
        "data/input_files/metadata.tsv",
        "data/input_files/metadata.csv"]
    for file_path in file_paths:
        if os.path.exists(file_path):
            # Read the file with appropriate separator based on its extension
            return pd.read_csv(file_path, sep='\t' if file_path.endswith('.tsv') else ',', index_col=None)
    
    # Raise an error if no valid file is found
    raise FileNotFoundError(f"No valid file found in the paths: {file_paths}")
    
def read_taxonomy():
    # Define file paths
    file_paths = [
        "data/input_files/taxonomy.tsv",
        "data/input_files/taxonomy.csv"
    ]
    
    # Read the taxonomy file
    taxonomy = None
    for file_path in file_paths:
        if os.path.exists(file_path):
            taxonomy = pd.read_csv(file_path, sep='\t' if file_path.endswith('.tsv') else ',', index_col=None)
            break
    
    if taxonomy is None:
        raise FileNotFoundError(f"No valid file found in the paths: {file_paths}")
    
    # Function to find the last non-empty value (taxa) in each row
    def get_last_taxon(row):
        non_empty_values = row.dropna()
        if not non_empty_values.empty:
            # Get the last non-empty value
            last_value = non_empty_values.iloc[-1]
            # Get the column of the last non-empty value
            last_value_col = non_empty_values.index[-1]
            
            if last_value_col != "Genus":
                last_value = f"{last_value} (unknown)"
            
            return last_value
        else:
            return None
    
    # Apply the function to each row
    taxonomy['Last_taxon'] = taxonomy.apply(get_last_taxon, axis=1)
    
    # Select only the necessary columns
    taxonomy_last_taxon = taxonomy[['SeqID', 'Last_taxon']]
    
    return taxonomy_last_taxon

# Function to transform relabund_tab
def transform_relabund_tab(relabund_tab):
    relabund_tab_long = relabund_tab.reset_index().melt(id_vars=['index'], var_name='SeqID', value_name='Relabund')
    relabund_tab_long.rename(columns={'index': 'Sample'}, inplace=True)
    relabund_tab_long['Percent_abundance'] = relabund_tab_long['Relabund'] * 100
    return relabund_tab_long.drop(columns=['Relabund'])
    
# Function to calculate the number of samples based on metadata
def calculate_nsamples(metadata, target, specificity_exception):
        
    # Find target and nontarget samples
    target_samples = metadata.loc[metadata['Source'].isin(target), 'Sample'].tolist()
    nontarget_samples = metadata.loc[~metadata['Source'].isin(target + specificity_exception), 'Sample'].tolist()
    specificity_exception_samples = metadata.loc[metadata['Source'].isin(specificity_exception), 'Sample'].tolist()
    
    # Calculate number of samples per source
    nsamples_per_source = metadata.groupby('Source')['Sample'].nunique().reset_index(name='Nsamples')
    
    # Calculate number of target and non-target samples
    nsamples_target = len(target_samples)
    nsamples_nontarget = len(nontarget_samples)
    nsamples_specificity_exception = len(specificity_exception_samples)
    
    return nsamples_per_source, nsamples_target, nsamples_specificity_exception, nsamples_nontarget, target_samples, nontarget_samples, specificity_exception_samples
    
# Function to join and group data
def relabund_metadata(relabund_tab_long, metadata):
    
    # Select only 'Sample' and 'Source' columns from metadata
    metadata_subset = metadata[['Sample', 'Source']]
    
    # Join relabund table and metadata
    relabund_tab_long_metadata = pd.merge(relabund_tab_long, metadata_subset, on='Sample', how='left')
    
    return relabund_tab_long_metadata

# Function to create seqIDs_samples dataframes
def seqIDs_samples(relabund_tab_long_metadata, target_samples, nontarget_samples, specificity_exception_samples):
    def group_and_filter(data, samples):
        # Filter data for the provided samples
        filtered_data = data[data['Sample'].isin(samples)]
        
        # Check if the filtered data is empty
        if filtered_data.empty:
            # Return an empty DataFrame with the desired columns
            return pd.DataFrame(columns=['SeqID', 'Source', 'Samples'])
        
        # If not empty, group by SeqID and Source and concatenate Sample values
        grouped = filtered_data.groupby(['SeqID', 'Source']).agg({
            'Sample': lambda x: ", ".join(x)
        }).reset_index()
        
        # Rename the aggregated column to 'Samples'
        grouped.rename(columns={'Sample': 'Samples'}, inplace=True)
        
        return grouped

    # Process each of the sample lists
    seqIDs_samples_target = group_and_filter(relabund_tab_long_metadata, target_samples)
    seqIDs_samples_nontarget = group_and_filter(relabund_tab_long_metadata, nontarget_samples)
    seqIDs_samples_exceptions = group_and_filter(relabund_tab_long_metadata, specificity_exception_samples)
    
    return seqIDs_samples_target, seqIDs_samples_nontarget, seqIDs_samples_exceptions

def extract_file_number(filename):
    # Split the filename based on non-digit characters
    parts = filename.split('_')[-1].split('table')
    number_part = parts[-1].replace('.txt', '')  # Remove file extension
    
    return int(number_part) if number_part.isdigit() else 0
