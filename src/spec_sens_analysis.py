import pandas as pd
import gc
import re

def get_pp_presence(relabund_tab_long_metadata, tntblast_results, target):
    # Create seqID_presence DataFrame
    relabund_tab_long_metadata['Percent_abundance'] = pd.to_numeric(relabund_tab_long_metadata['Percent_abundance'], errors='coerce')
    # Filter rows where 'Percent_abundance' is greater than 0
    filtered_metadata = relabund_tab_long_metadata[relabund_tab_long_metadata['Percent_abundance'] > 0]
    # Create seqID_presence DataFrame
    seqID_presence = (
        filtered_metadata
        .assign(Target_nontarget=lambda df: df['Source'].apply(lambda x: 'T' if x in target else 'NT'))
        .groupby('SeqID')
        .agg({'Target_nontarget': lambda x: ','.join(x.unique())})
        .reset_index()
        .assign(Target_nontarget=lambda df: df['Target_nontarget'].str.replace('NT,T', 'T,NT'))
    )
    # Create pp_presence DataFrame
    pp_presence = (
        tntblast_results
        .merge(seqID_presence, on='SeqID', how='left')
        .assign(Target_nontarget=lambda df: df['Target_nontarget'].str.split(','))
        .explode('Target_nontarget')
        .groupby(['PP_ID', 'PrimerF', 'PrimerR'])
        .agg({'Target_nontarget': lambda x: ','.join(x.unique())})
        .reset_index()
        .assign(Target_nontarget=lambda df: df['Target_nontarget'].str.replace('NT,T', 'T,NT'))
    )
    return pp_presence


def get_amplicon_sizes(tntblast_results, pp_presence):
    # Merge with pp_presence DataFrame
    merged_df = (tntblast_results[['PP_ID', 'Amplicon_size', 'SeqID']]
                 .merge(pp_presence, on='PP_ID', how='left')
                )
    # Process the Target_nontarget column and aggregate amplicon sizes
    amplicon_sizes_df = (merged_df
                         .assign(Target_nontarget=lambda df: df['Target_nontarget'].str.split(','))
                         .explode('Target_nontarget')
                         .groupby(['PP_ID', 'Target_nontarget'])
                         .agg({'Amplicon_size': lambda x: ', '.join(map(str, sorted(x.unique())))})
                         .reset_index()
                         .pivot(index='PP_ID', columns='Target_nontarget', values='Amplicon_size')
                         .rename(columns={'T': 'Amplicon_sizes_target', 'NT': 'Amplicon_sizes_nontarget'})
                         .fillna({'Amplicon_sizes_target': pd.NA, 'Amplicon_sizes_nontarget': pd.NA})
                         .reset_index()
                        )
    return amplicon_sizes_df
    
def get_primers_info(tntblast_results, pp_presence):
    """
    Process primer information to get detailed target-specific and general primer info for forward and reverse primers.
    
    Parameters:
    - tntblast_results (pd.DataFrame): DataFrame containing primer results with columns like PrimerF, TmF, etc.
    - pp_presence (pd.DataFrame): DataFrame containing primer presence information with columns like PP_ID, PrimerF, etc.
    
    Returns:
    - pd.DataFrame: Merged DataFrame with target-specific and general primer information for both forward and reverse primers.
    """
    
    primer_types = ['F', 'R']
    all_primers_info = []
    
    for primer_type in primer_types:
        suffix = primer_type
        
        # Define column renaming mapping based on primer type
        column_names = {
            f'Tm{suffix}_T': f'Tm{suffix}_target',
            f'Mismatch{suffix}_T': f'Mismatch{suffix}_target',
            f'Tm{suffix}_NT': f'Tm{suffix}_nontarget',
            f'Mismatch{suffix}_NT': f'Mismatch{suffix}_nontarget',
        }
        
        # Process primer target-specific information
        primers_target_info = (
            tntblast_results[[f'Primer{suffix}', f'Tm{suffix}', f'Mismatch{suffix}', 'PP_ID', 'SeqID']]
            .assign(**{f'Tm{suffix}': lambda df: df[f'Tm{suffix}'].round(2)})
            .merge(pp_presence, on=['PP_ID', f'Primer{suffix}'], how='left')
            .assign(Target_nontarget=lambda df: df['Target_nontarget'].str.split(','))
            .explode('Target_nontarget')
            .groupby([f'Primer{suffix}', 'Target_nontarget', 'PP_ID'])
            .agg({
                f'Tm{suffix}': lambda x: ', '.join(map(str, sorted(x.unique()))),
                f'Mismatch{suffix}': lambda x: f"[{x.min()}-{x.max()}]"
            })
            .reset_index()
            .pivot(index=f'PP_ID', columns='Target_nontarget', values=[f'Tm{suffix}', f'Mismatch{suffix}'])
            .reset_index()
        )
        
        
        # Flatten MultiIndex columns and rename
        primers_target_info.columns = [
            column_names.get('_'.join(col).strip('_'), '_'.join(col).strip('_'))
            for col in primers_target_info.columns.values
        ]
        
        # Check for missing columns and add them if necessary
        missing_columns = set(column_names.values()) - set(primers_target_info.columns)
        for col in missing_columns:
            primers_target_info[col] = pd.NA
        
        # Process general primer information
        primers_general_info = (
            tntblast_results[[f'Primer{suffix}', f'Tm{suffix}', f'Mismatch{suffix}', f'Heuristics{suffix}', 
                f'Hairpin_Tm_{suffix}', f'Homodimer_Tm_{suffix}', 'Heterodimer_Tm']]
            .groupby(f'Primer{suffix}')
            .agg(
                Tm_max_temp=(f'Tm{suffix}', 'max'),
                Heuristics_temp=(f'Heuristics{suffix}', lambda x: '/'.join(x.unique())),
                Hairpin_Tm=(f'Hairpin_Tm_{suffix}', 'first'),
                Homodimer_Tm=(f'Homodimer_Tm_{suffix}', 'first'),
                Heterodimer_Tm=('Heterodimer_Tm', 'first')
            )
            .assign(Tm_max_temp=lambda df: df['Tm_max_temp'].round(2))
            .rename(columns={'Heuristics_temp': f'Heuristics{suffix}', 'Tm_max_temp': f'Tm{suffix}_max',
                'Hairpin_Tm': f'Hairpin_Tm_{suffix}', 'Homodimer_Tm': f'Homodimer_Tm_{suffix}'})
        )
    
        # Merge the general and target-specific primer information
        merged_primers_info = (
            tntblast_results[['PP_ID', f'Primer{suffix}']]
            .drop_duplicates()
            .merge(primers_target_info, on = f'PP_ID', how='left')
            .merge(primers_general_info, on = f'Primer{suffix}', how='left')
        )
        
        # Store results based on primer type
        if suffix == 'F':
            merged_primers_infoF = merged_primers_info
        else:
            merged_primers_infoR = merged_primers_info

    # Select and merge primer information
    primers_info = (
        merged_primers_infoF
        .merge(merged_primers_infoR.drop(columns=['Heterodimer_Tm']), 
            left_on='PP_ID', right_on='PP_ID', how='left')
    )

    # Filter rows where mismatches are not allowed
    primers_info_final = (
        primers_info
        .loc[primers_info['MismatchF_target'].str.contains(r'^\[0', na=False)]
        .loc[primers_info['MismatchR_target'].str.contains(r'^\[0', na=False)]
    )
 
    
    # If primers_info_final does not exist, create an empty DataFrame
    try:
        primers_info_final
    except NameError:
        primers_info_final = pd.DataFrame()
        print("Primers_info file is empty")
    
    return primers_info_final

    
def filter_tntblast_results(tntblast_results, primers_info):
    if primers_info.empty:
        return None
    else:
        # Filter tntblast_results to include only rows with PP_IDs present in primers_info
        tntblast_results_filt = tntblast_results[tntblast_results['PP_ID'].isin(primers_info['PP_ID'])]
        return tntblast_results_filt


def get_abundance_df(tntblast_results_filt, seqIDs_samples, relabund_tab_long_metadata, taxonomy, nsamples_per_source, sample_type):
    """
    Calculate non-target abundance DataFrame based on sensitivity information.
    
    Parameters:
    - tntblast_results_filt (pd.DataFrame): Filtered TNTBLAST results containing SeqID and PP_ID.
    - pp_sensitivity (pd.DataFrame): DataFrame with primer pair sensitivity information.
    - seqIDs_samples_nontarget (pd.DataFrame): DataFrame mapping SeqID to non-target sample information.
    - relabund_tab_long_metadata (pd.DataFrame): Metadata containing percent abundance and other details.
    - taxonomy (pd.DataFrame): Taxonomy DataFrame with taxon information.
    - sample_type: Type of samples ('target', 'nontarget' or 'exception').
    
    Returns:
    - pd.DataFrame: Non-target abundance DataFrame.
    """
    
    suffix = sample_type
    
    abundance_df = (
         tntblast_results_filt[['SeqID', 'PP_ID']]
        .groupby('SeqID')
        .agg({'PP_ID': lambda x: ', '.join(x)})
        .reset_index()
        .merge(seqIDs_samples, on='SeqID')
        .assign(PP_ID=lambda df: df['PP_ID'].str.split(', '))
        .explode('PP_ID')
        .assign(Sample=lambda df: df['Samples'].str.split(', '))
        .explode('Sample')
        .merge(relabund_tab_long_metadata, on=['Sample', 'SeqID', 'Source'])
        .merge(taxonomy, on='SeqID')
        .merge(nsamples_per_source, on='Source')
        .groupby(['PP_ID', 'Source', 'Nsamples', 'Sample'])
        .agg({
            'SeqID': lambda x: ', '.join(x),
            'Percent_abundance': lambda x: round(sum(x), 4),
            'Last_taxon': lambda x: ', '.join(x.unique())
        })
        .rename(columns={'Last_taxon': f'Taxonomy_{suffix}'})
        .reset_index()
        )
    return abundance_df
    
def calculate_sensitivity(target_abundance_df, nsamples_target, marker_sensitivity_cutoff):
    """
    Calculate sensitivity from the target abundance DataFrame.

    Parameters:
    - target_abundance_df (pd.DataFrame): DataFrame containing target abundances.
    - nsamples_target (int): Number of target samples.
    - marker_sensitivity_cutoff (float): Sensitivity cutoff value.

    Returns:
    - pd.DataFrame: DataFrame containing sensitivity metrics.
    """
    return (
        target_abundance_df[target_abundance_df['Percent_abundance'] != 0]
        .groupby('PP_ID')
        .agg(N_positive_target_samples=('Sample', lambda x: len(x.unique())))
        .assign(
            Sensitivity=lambda df: round(df['N_positive_target_samples'] / nsamples_target * 100, 2),
            Sensitivity2=lambda df: df['N_positive_target_samples'].astype(str) + '/' + str(nsamples_target)
        )
        .drop(columns=['N_positive_target_samples'])
        .query('Sensitivity >= @marker_sensitivity_cutoff')
        .reset_index()
    )

def calculate_specificity(pp_sensitivity, nontarget_abundance, nsamples_nontarget, marker_specificity_cutoff):
    """
    Calculate specificity from the non-target abundance DataFrame.

    Parameters:
    - pp_sensitivity (pd.DataFrame): DataFrame containing sensitivity metrics.
    - nontarget_abundance (pd.DataFrame): DataFrame containing non-target abundances.
    - nsamples_nontarget (int): Number of non-target samples.
    - marker_specificity_cutoff (float): Specificity cutoff value.

    Returns:
    - pd.DataFrame: DataFrame containing specificity metrics.
    """
    nontarget_abundance = nontarget_abundance.drop(columns=['Nsamples'])
    
    pp_samples_nontarget = (
        nontarget_abundance
        .groupby('PP_ID')
        .agg(Sample=('Sample', lambda x: ', '.join(x.unique())))
        .reset_index()
    )
    
    pp_nontarget = (
        pp_sensitivity
        .merge(pp_samples_nontarget, on='PP_ID', how='left')
        .assign(Sample=lambda df: df['Sample'].str.split(', '))
        .explode('Sample')
        .merge(nontarget_abundance, on=['PP_ID', 'Sample'], how='left')
    )
    
    
    pp_nontarget['PA'] = pp_nontarget.apply(
        lambda row: 1 if not pd.isna(row['Percent_abundance']) and row['Percent_abundance'] != 0 else 0, axis=1
    )    
    

    specificity_sensitivity_df = (
        pp_nontarget
        .groupby(['PP_ID', 'Sensitivity', 'Sensitivity2'])
        .agg(
            N_positive_nontarget_samples=('PA', 'sum'),
            TN=('PA', lambda x: nsamples_nontarget - x.sum())
        )
        .assign(
            Specificity=lambda df: round(df['TN'] / nsamples_nontarget * 100, 2),
            Specificity2=lambda df: df['TN'].astype(str) + '/' + str(nsamples_nontarget)
        )
        .drop(columns=['TN'])
        .query('Specificity >= @marker_specificity_cutoff')
        .reset_index()
    )
    
    return specificity_sensitivity_df

def calculate_sensitivity_specificity(target_abundance_df, nontarget_abundance, nsamples_target, nsamples_nontarget, marker_sensitivity_cutoff, marker_specificity_cutoff):
    """
    Calculate both sensitivity and specificity for primer pairs based on target and non-target samples.

    Parameters:
    - target_abundance_df (pd.DataFrame): DataFrame containing target abundances.
    - nontarget_abundance (pd.DataFrame): DataFrame containing non-target abundances.
    - nsamples_target (int): Number of target samples.
    - nsamples_nontarget (int): Number of non-target samples.
    - marker_sensitivity_cutoff (float): Sensitivity cutoff value.
    - marker_specificity_cutoff (float): Specificity cutoff value.

    Returns: pp_sensitivity_specificity
    """

    # Calculate sensitivity
    pp_sensitivity = calculate_sensitivity(target_abundance_df, nsamples_target, marker_sensitivity_cutoff)
    
    # Calculate specificity
    sensitivity_specificity = calculate_specificity(pp_sensitivity, nontarget_abundance, nsamples_nontarget, marker_specificity_cutoff)

    return sensitivity_specificity
    
    
def get_sensitivity_detailed(target_abundance, sensitivity_specificity):
    """
    Calculate detailed sensitivity metrics for primer pairs based on target abundances.

    Parameters:
    - target_abundance (pd.DataFrame): DataFrame containing target abundances and other details.
    - sensitivity_specificity (pd.DataFrame): DataFrame containing specificity and sensitivity information.

    Returns:
    - pd.DataFrame: DataFrame with detailed sensitivity information.
    """
    # Filter based on 'Percent_abundance' not equal to 0
    filtered_target_abundance = target_abundance[target_abundance['Percent_abundance'] != 0]

    # Filter based on 'PP_ID' present in 'sensitivity_specificity'
    filtered_target_abundance = filtered_target_abundance[
        filtered_target_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID'])]

    # Calculate detailed sensitivity
    target_abundance_grouped = (
        filtered_target_abundance
        .groupby(['Source', 'PP_ID', 'Nsamples'])
        .agg(
            N_positive_target_samples=('Sample', 'nunique'),
            Positive_target_samples=('Sample', lambda x: ', '.join(x.unique())))
        .reset_index()
    )
    
    target_abundance_grouped['Sensitivity2_detailed'] = (
        target_abundance_grouped.apply(lambda row: f"{row['Source']} ({row['N_positive_target_samples']}/{row['Nsamples']})", axis=1)
        )

    sensitivity_detailed = (
        target_abundance_grouped
        .groupby('PP_ID')
        .apply(lambda df: pd.Series({
            'Sensitivity2_detailed': ', '.join(
                f"{row['Source']} ({row['N_positive_target_samples']}/{row['Nsamples']})"
                for _, row in df.iterrows()
            ),
            'Positive_target_samples': ', '.join(
                f"{row['Source']} ({row['Positive_target_samples']})"
                for _, row in df.iterrows()
            )
        }))
        .reset_index()
    )

    return sensitivity_detailed
    
def calculate_abundance_target(target_abundance, sensitivity_specificity):
    """
    Calculate detailed abundance metrics for primer pairs based on target samples.

    Parameters:
    - target_abundance (pd.DataFrame): DataFrame containing target abundances and other details.
    - sensitivity_specificity (pd.DataFrame): DataFrame containing specificity and sensitivity information.

    Returns:
    - pd.DataFrame: DataFrame with detailed abundance information for primer pairs.
    """
    # Filter based on 'Percent_abundance' > 0
    filtered_target_abundance = target_abundance[target_abundance['Percent_abundance'] > 0]

    # Filter based on 'PP_ID' present in 'sensitivity_specificity'
    filtered_target_abundance = filtered_target_abundance[
        filtered_target_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID'])
    ]

    # Calculate relative abundances of markers in target samples
    pp_abundance_target = (
        filtered_target_abundance
        .groupby('PP_ID')
        .agg(
            Percent_abundance2=('Percent_abundance', 'mean'),
            Percent_abundance2_SD=('Percent_abundance', 'std')
        )
        .assign(
            Percent_abundance_target=lambda df: df.apply(
                lambda row: f"{row['Percent_abundance2']:.4f}, SD={row['Percent_abundance2_SD']:.4f}", 
                axis=1
            )
        )
        .reset_index()
        .drop(columns=['Percent_abundance2', 'Percent_abundance2_SD'])
    )

    # Detailed relative abundances in target samples
    pp_abundance_target_detailed = (
        filtered_target_abundance
        .groupby(['PP_ID', 'Source', 'Nsamples'])
        .agg(
            Percent_abundance2=('Percent_abundance', 'mean'),
            Percent_abundance2_SD=('Percent_abundance', 'std')
        )
        .reset_index()
        .assign(
            Percent_abundance_target_detailed=lambda df: df.apply(
                lambda row: f"{row['Source']} ({row['Percent_abundance2']:.4f}" +
                    (f", SD={row['Percent_abundance2_SD']:.4f}" if pd.notna(row['Percent_abundance2_SD']) else "") + ")", 
                    axis=1
            )
        )
        .groupby('PP_ID')
        .agg(
            Percent_abundance_target_detailed=('Percent_abundance_target_detailed', ', '.join)
        )
        .reset_index()
        .merge(pp_abundance_target, on='PP_ID', how='left')
    )

    return pp_abundance_target_detailed
    

def get_taxonomy_target(target_abundance, sensitivity_specificity):
    """
    Calculate taxonomy information for both target and non-target primer pairs.

    Parameters:
    - target_abundance (pd.DataFrame): DataFrame containing target abundances and taxonomy information.
    - sensitivity_specificity (pd.DataFrame): DataFrame containing specificity and sensitivity information.

    Returns: pp_taxonomy_target (pd.DataFrame): DataFrame with aggregated taxonomy for target primer pairs.
    """
    
    # Calculate taxonomy for target primer pairs
    pp_taxonomy_target = (
        target_abundance[target_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID'])]
        .assign(Taxonomy_target=lambda df: df['Taxonomy_target'].str.split(', '))
        .explode('Taxonomy_target')
        .groupby('PP_ID')
        .agg(Taxonomy_target=('Taxonomy_target', lambda x: ', '.join(sorted(set(x)))))
        .reset_index()
    )
    return pp_taxonomy_target
    
def process_negative_target_samples(target_abundance, sensitivity_specificity):
    """
    Calculate negative target samples based on the given DataFrame.

    Parameters:
    - target_abundance (pd.DataFrame): DataFrame containing abundance data with 'PP_ID', 'Source', 'Sample', and 'Percent_abundance'.
    - sensitivity_specificity (pd.DataFrame): DataFrame containing PP_IDs that need to be filtered.

    Returns:
    - pd.DataFrame: DataFrame with aggregated negative target samples information, or an empty DataFrame with specified columns if no data is found.
    """
    # Filter the DataFrame
    filtered_df = target_abundance[
        (target_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID'])) & 
        (target_abundance['Percent_abundance'] == 0)
    ]

    # If no negative target samples are found, return an empty DataFrame with specified columns
    if filtered_df.empty:
        return pd.DataFrame(columns=['PP_ID', 'Negative_target_samples'])

    # Group by PP_ID and Source, then aggregate Sample
    negative_samples_grouped = (
        filtered_df
        .groupby(['PP_ID', 'Source'])
        .agg(Negative_target_samples=('Sample', lambda x: ', '.join(x)))
        .reset_index()
    )

    # Group by PP_ID and aggregate the Source and Negative_target_samples
    pp_negative_target_samples = (
        negative_samples_grouped
        .groupby('PP_ID')
        .apply(lambda x: ', '.join(
            f"{row['Source']} ({row['Negative_target_samples']})"
             for _, row in x.iterrows()
        ))
        .reset_index(name='Negative_target_samples')
    )

    return pp_negative_target_samples

def process_detected_nontarget_samples(nontarget_abundance, sensitivity_specificity):
    """
    Calculate the number of detected non-target samples and relative abundances of the given markers.

    Parameters:
    - nontarget_abundance (pd.DataFrame): DataFrame containing non-target abundances and other details.
    - sensitivity_specificity (pd.DataFrame): DataFrame containing specificity and sensitivity information.

    Returns:
    - pd.DataFrame: DataFrame with summarized marker detection information for non-target samples.
    """

    # Filter pp_nontarget_df for Percent_abundance > 0 and PP_ID present in sensitivity_specificity
    filtered_df = nontarget_abundance[
        (nontarget_abundance['Percent_abundance'] > 0) &
        (nontarget_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID']))
    ]
    
    if filtered_df.empty:
        return pd.DataFrame(columns=['PP_ID', 'Sensitivity', 'Sensitivity2', 'Sample', 'Source', 'SeqIDs', 'Percent_abundance', 'Taxonomy_nontarget'])
    else:

        positive_samples = (
            filtered_df
            .groupby(['PP_ID', 'Nsamples', 'Source'])
            .agg(
                Percent_abundance_nontarget=('Percent_abundance', 'sum'),
                N_positive_nontarget_samples=('Sample', lambda x: len(set(x))),
                Positive_nontarget_samples=('Sample', lambda x: ', '.join(sorted(set(x))))
            )
            .reset_index()
        )
        
        # Group by PP_ID to aggregate the information
        detected_nontarget = (
            positive_samples
            .groupby('PP_ID')
            .apply(lambda df: pd.Series({
                'Presence_nontarget_samples': ', '.join(
                    f"{row['Source']} ({row['N_positive_nontarget_samples']}/{row['Nsamples']})"
                    for _, row in df.iterrows()
                ),
                'Positive_nontarget_samples': ', '.join(
                    f"{row['Source']} ({row['Positive_nontarget_samples']})"
                    for _, row in df.iterrows()
                )
            }))
            .reset_index()
        )
        
        # Calculate the relative abundance [%] of markers detected
        pp_abundance_nontarget_sum = (
            filtered_df
            .groupby('PP_ID')
            .agg(
                Percent_abund=('Percent_abundance', 'mean'),
                Percent_abund_SD=('Percent_abundance', 'std'),
            )
            .reset_index()
            .assign(
                Percent_abundance_nontarget=lambda df: df.apply(
                    lambda row: (
                        f"{round(row['Percent_abund'], 4)}" +
                        (f", SD={round(row['Percent_abund_SD'], 4)}" if pd.notna(row['Percent_abund_SD']) else "")),
                    axis=1
                )
            )
        )
        
        pp_abundance_nontarget_sum = pp_abundance_nontarget_sum[['PP_ID', 'Percent_abundance_nontarget']]
    
        # Calculate detailed abundance
        pp_abundance_nontarget_detailed = (
            filtered_df
            .groupby(['PP_ID', 'Source'])
            .agg(
                Percent_abund_nontarget=('Percent_abundance', 'mean'),
                Percent_abund_nontarget_SD=('Percent_abundance', 'std')
            )
            .reset_index()
            .groupby('PP_ID')
            .apply(lambda df: pd.Series({
                'Percent_abundance_nontarget_detailed': ', '.join(
                    f"{row['Source']} ({round(row['Percent_abund_nontarget'], 4)}" +
                    (f", SD={round(row['Percent_abund_nontarget_SD'], 4)}" if pd.notna(row['Percent_abund_nontarget_SD']) else "") +
                    ")"
                    #f"{row['Source']} ({round(row['Percent_abund_nontarget'], 4)} SD={round(row['Percent_abund_nontarget_SD'], 4)})"
                    for _, row in df.iterrows()
                )
            }))
            .reset_index()
        )
    
        # Taxonomy of detected nontarget sequences
        pp_taxonomy_nontarget = (
            filtered_df
            .assign(Taxonomy_nontarget=filtered_df['Taxonomy_nontarget'].str.split(', '))
            .explode('Taxonomy_nontarget')
            .groupby('PP_ID')
            .agg(
                Taxonomy_nontarget=('Taxonomy_nontarget', lambda x: ', '.join(sorted(x.unique())))
            )
            .reset_index()
        )

        # Merge all results into a single DataFrame
        detected_nontarget_summary = (
            detected_nontarget
            .merge(pp_abundance_nontarget_sum, on='PP_ID')
            .merge(pp_abundance_nontarget_detailed, on='PP_ID')
            .merge(pp_taxonomy_nontarget, on='PP_ID')
        )     
        
        return detected_nontarget_summary


def process_exception_samples(exception_abundance, sensitivity_specificity, taxonomy, specificity_exception):
    
    filtered_exception_abundance = exception_abundance[
        (exception_abundance['Percent_abundance'] > 0) &
        (exception_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID']))
    ]
    
    if len(filtered_exception_abundance) > 0:
        # Calculate detailed sensitivity
        exception_abundance_grouped = (
            filtered_exception_abundance 
            .groupby(['Source', 'PP_ID', 'Nsamples'])
            .agg(
                N_positive_exceptions_samples=('Sample', 'nunique'),
                Positive_exceptions_samples=('Sample', lambda x: ', '.join(x.unique()))
            )
            .reset_index()
        )
    
        exception_abundance_grouped['Presence_exceptions_samples'] = (
            exception_abundance_grouped
            .apply(lambda row: f"{row['Source']} ({row['N_positive_exceptions_samples']}/{row['Nsamples']})", axis=1)
        )

        positive_exceptions = (
            exception_abundance_grouped
            .groupby('PP_ID')
            .apply(lambda df: pd.Series({
                'Presence_exceptions_samples': ', '.join(
                    f"{row['Source']} ({row['N_positive_exceptions_samples']}/{row['Nsamples']})"
                    for _, row in df.iterrows()
                ),
                'Positive_exceptions_samples': ', '.join(
                    f"{row['Source']} ({row['Positive_exceptions_samples']})"
                    for _, row in df.iterrows()
                )
            }))
            .reset_index()
        )
        
        
        # Calculate relative abundances of markers within samples excluded from specificity calculations
        exception_abundance_mean = (
            filtered_exception_abundance
            .groupby(['PP_ID', 'Source', 'Nsamples'], as_index=False)
            .agg({
                'Percent_abundance': ['mean', 'std']
            })
        )
        
        # Flatten the MultiIndex columns
        exception_abundance_mean.columns = ['_'.join(col).strip().rstrip('_') for col in exception_abundance_mean.columns]

        pp_abundance_exceptions = (
            exception_abundance_mean
            .rename(columns={
                'Percent_abundance_mean': 'Percent_abundance_exceptions',
                'Percent_abundance_std': 'Percent_abundance_exceptions_SD'
            })
            .groupby('PP_ID', as_index=False)
            .agg({
                'Source': 'first',
                'Percent_abundance_exceptions': lambda x: ', '.join(map(str, x.fillna(''))),
                'Percent_abundance_exceptions_SD': lambda x: ', '.join(map(str, x.fillna('')))
            })
            .assign(
                Percent_abundance_exceptions_detailed=lambda df: df.apply(
                    lambda row: (
                        f"{row['Source']} ({round(float(row['Percent_abundance_exceptions']), 4)}"
                        + (f" SD={round(float(row['Percent_abundance_exceptions_SD']), 4)}" 
                           if pd.notna(row['Percent_abundance_exceptions_SD']) and row['Percent_abundance_exceptions_SD'].strip() != ''
                           else "") +
                        ")"
                    ),
                    axis=1
                )
            )
            .drop(columns=['Percent_abundance_exceptions', 'Percent_abundance_exceptions_SD'])
            .rename(columns={'Source': 'Exceptions'})
        )
        
        # Assign taxonomy to each marker found in specificity exception samples
        pp_taxonomy_exceptions = (
            filtered_exception_abundance 
            .assign(Taxonomy_exceptions=filtered_exception_abundance['Taxonomy_exception'].str.split(', ')) 
            .explode('Taxonomy_exceptions')
            .groupby('PP_ID', as_index=False)
            .agg({'Taxonomy_exceptions': lambda x: ', '.join(sorted(set(x)))})
        )
        
        # Create a joined table        
        pp_exceptions_join = (
            sensitivity_specificity[['PP_ID']]
            .merge(positive_exceptions, on='PP_ID', how='left')
            .merge(pp_abundance_exceptions, on='PP_ID', how='left')
            .merge(pp_taxonomy_exceptions, on='PP_ID', how='left')   
        )
    
    else:
        
        pp_exceptions_join = (
            sensitivity_specificity[['PP_ID']]
            .assign(
                Presence_exceptions_samples=None,
                Positive_exceptions_samples=None,
                Percent_abundance_exceptions_detailed=None,
                Taxonomy_exceptions=None
            )
        )
    
    return pp_exceptions_join

def process_sequence_ids(tntblast_results_filt, sensitivity_specificity, target_abundance, nontarget_abundance, file_number):
    
    # Write out all detected sequences
    tntblast_results_subset = tntblast_results_filt[tntblast_results_filt['PP_ID'].isin(sensitivity_specificity['PP_ID'])]
    
    pp_all_detected_seqIDs = (
        tntblast_results_subset
        .groupby('PP_ID', as_index=False)
        .agg({'SeqID': lambda x: ', '.join(sorted(set(x), key=lambda y: int(re.search(r'\d+', y).group())))})
        .rename(columns={'SeqID': 'SeqIDs_all'})
    )
    
    # Write sequence IDs that are detected by best primer pairs in target samples
    target_abundance_subset = (
        target_abundance[target_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID'])]
        .query('Percent_abundance > 0')
        [['PP_ID', 'SeqID']]
    )
    
    pp_seqIDs_target_separated = (
        target_abundance_subset
        .assign(SeqID=lambda df: df['SeqID'].str.split(', '))
        .explode('SeqID')
        .drop_duplicates(subset=['PP_ID', 'SeqID'])
    )
    
    
    pp_seqIDs_target = (
        pp_seqIDs_target_separated
        .groupby('PP_ID', as_index=False)
        .agg(
            SeqIDs_target=pd.NamedAgg(
                column='SeqID',
                aggfunc=lambda x: ', '.join(sorted(set(x), key=lambda y: int(re.search(r'\d+', y).group())))
            ),
            N_seqIDs_target=pd.NamedAgg(
                column='SeqID',
                aggfunc=lambda x: len(set(x))
            )
        )
    )
    
    nontarget_abundance_subset = (
        nontarget_abundance[nontarget_abundance['PP_ID'].isin(sensitivity_specificity['PP_ID'])]
        .query('Percent_abundance > 0')
        [['PP_ID', 'SeqID']]
    )
    
    if nontarget_abundance_subset.empty:
        pp_seqIDs_nontarget_separated = pd.DataFrame(columns=['PP_ID', 'SeqID'])
        pp_seqIDs_nontarget = pd.DataFrame(columns=['PP_ID', 'SeqIDs_nontarget', 'N_seqIDs_nontarget'])   
    else:    
        pp_seqIDs_nontarget_separated = (
            nontarget_abundance_subset
            .assign(SeqID=lambda df: df['SeqID'].str.split(', '))
            .explode('SeqID')
            .drop_duplicates(subset=['PP_ID', 'SeqID'])
        )
        pp_seqIDs_nontarget = (
            pp_seqIDs_nontarget_separated
            .groupby('PP_ID', as_index=False)
            .agg(
                SeqIDs_nontarget=pd.NamedAgg(
                    column='SeqID',
                    aggfunc=lambda x: ', '.join(sorted(set(x), key=lambda y: int(re.search(r'\d+', y).group())))
                ),
                N_seqIDs_nontarget=pd.NamedAgg(
                    column='SeqID',
                    aggfunc=lambda x: len(set(x))
                )
            )
        )
    
    # Perform anti join to find sequence IDs detected in target but not in nontarget samples
    joins = (
        pp_seqIDs_target_separated
        .merge(pp_seqIDs_nontarget_separated, on=['PP_ID', 'SeqID'], how = 'left', indicator=True)
    )
        
    pp_seqIDs_target_only = (
        joins
        .query('_merge == "left_only"')
        .groupby('PP_ID', as_index=False)
        .agg({'SeqID': lambda x: ', '.join(sorted(set(x), key=lambda y: int(re.search(r'\d+', y).group())))})
        .rename(columns={'SeqID': 'SeqIDs_target_only'})
    )
    
    # Perform anti join to find sequence IDs detected in nontarget but not in target samples
    nontarget_joins = (
        joins
        .query('_merge == "right_only"')
    )
    
    if nontarget_joins.empty:
        pp_seqIDs_nontarget_only = pd.DataFrame(columns=['PP_ID', 'SeqIDs_nontarget_only'])
    else:
        pp_seqIDs_nontarget_only = (
            nontarget_joins
            .groupby('PP_ID', as_index=False)
            .agg({'SeqID': lambda x: ', '.join(x)})
            .rename(columns={'SeqID': 'SeqIDs_nontarget_only'})
        )
    
    # Join all info
    pp_seqIDs_all = (
        sensitivity_specificity[['PP_ID']]
        .merge(pp_all_detected_seqIDs, on='PP_ID', how='left')
        .merge(pp_seqIDs_target, on='PP_ID', how='left')
        .merge(pp_seqIDs_nontarget, on='PP_ID', how='left')
        .merge(pp_seqIDs_target_only, on='PP_ID', how='left')
        .merge(pp_seqIDs_nontarget_only, on='PP_ID', how='left')
        .assign(File_number=file_number)
        .reindex(columns=['PP_ID', 'SeqIDs_all', 'seqIDs_target_only', 'seqIDs_nontarget_only', 'N_seqIDs_target', 'N_seqIDs_nontarget', 'File_number'])
    )
    
    return pp_seqIDs_all
    
def join_info(sensitivity_specificity, sensitivity_detailed, amplicon_sizes, primers_info, abundance_target, taxonomy_target, negative_target_samples, detected_nontarget, exceptions_info, pp_seqIDs, file_number, specificity_exception):
    
    if specificity_exception:
        exception_string = ', '.join(specificity_exception)
    else:
        exception_string = 'NA'
    
    # Perform left joins
    merged_df = (
        sensitivity_specificity
        .merge(sensitivity_detailed, how='left')
        .merge(amplicon_sizes, how='left')
        .merge(primers_info, how='left')
        .merge(abundance_target, how='left')
        .merge(taxonomy_target, how='left')
        .merge(negative_target_samples, how='left')
        .merge(detected_nontarget, how='left')
        .merge(exceptions_info, how='left')
        .merge(pp_seqIDs[['PP_ID', 'N_seqIDs_target']], how='left', on='PP_ID')
        .assign(
            File_number=file_number,
            Exceptions=exception_string
        )
    )
    
    # Change columns order
    desired_columns = [
        "PP_ID", "Specificity", "Specificity2", 
        "Sensitivity", "Sensitivity2", "Sensitivity2_detailed", 
        "Presence_nontarget_samples",  
        "Percent_abundance_target", "Percent_abundance_nontarget", 
        "Percent_abundance_target_detailed", "Percent_abundance_nontarget_detailed", 
        "Taxonomy_target", "Taxonomy_nontarget", "N_seqIDs_target",
        "Positive_target_samples", "Positive_nontarget_samples", "Negative_target_samples",
        "Exceptions", "Presence_exceptions_samples", "Percent_abundance_exceptions_detailed", 
        "Taxonomy_exceptions", "Positive_exceptions_samples", 
        "PrimerF", "PrimerR", 
        "TmF_max", "TmR_max", 
        "TmF_target", "TmF_nontarget", 
        "TmR_target", "TmR_nontarget", 
        "MismatchF_target", "MismatchF_nontarget", 
        "MismatchR_target", "MismatchR_nontarget", 
        "Amplicon_sizes_target", "Amplicon_sizes_nontarget", 
        "HeuristicsF", "HeuristicsR", 
        "Hairpin_Tm_F", "Hairpin_Tm_R", 
        "Homodimer_Tm_F", "Homodimer_Tm_R",
        "Heterodimer_Tm",
        "File_number"
    ]
    
    # Reindex the DataFrame to include only the desired columns in the specified order
    full_table = merged_df.reindex(columns=desired_columns)
    
    return full_table
