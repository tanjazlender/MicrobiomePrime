import os
import psutil
import pandas as pd

def clean_directory(path):
    if os.path.exists(path):
        existing_files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
        for f in existing_files:
            os.remove(os.path.join(path, f))
            

def check_if_empty(df, filename):
    """
    Checks if the DataFrame df has any rows.
    If it does not, prints a message and returns True to indicate processing should be skipped.
    If it has rows, returns False to indicate processing should continue.

    Parameters:
    df (DataFrame): The DataFrame to check.
    filename (str): File name associated with the DataFrame.
    
    Returns:
    bool: True if processing should be skipped (i.e., DataFrame is empty), False if it should continue (i.e., DataFrame is not empty).
    """
    if df.shape[0] == 0:
        print(f"{filename} does not contain primer pairs meeting the specified criteria.")
        return True  # Return True to indicate skipping
    return False  # Return False to indicate continuing


def get_memory_info():
    """Read the memory setting from the given INI file, or get available system memory if empty."""
    import configparser
    import psutil
    
    config = configparser.ConfigParser()
    config.read('scripts/settings.ini')

    try:
        # Get the 'memory' setting from the 'settings' section
        memory_str = config.get('settings', 'memory').strip()
        
        if memory_str:
            # Convert to integer if not empty
            memory = int(memory_str)
            print(f"Memory from settings: {memory} GB")
            return memory
        else:
            # If memory is empty, use available system memory
            available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)  # Convert bytes to GB
            print(f"Available system memory: {available_memory_gb:.2f} GB")
            return available_memory_gb
    
    except configparser.NoSectionError:
        # Handle the case where the section is missing
        print("Error: 'settings' section is missing in the settings file.")
        available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)  # Convert bytes to GB
        print(f"Available system memory: {available_memory_gb:.2f} GB")
        return available_memory_gb
    
    except configparser.NoOptionError:
        # Handle the case where the option is missing
        print("Error: 'memory' option is missing in the settings file.")
        available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)  # Convert bytes to GB
        print(f"Available system memory: {available_memory_gb:.2f} GB")
        return available_memory_gb
    
    except ValueError as e:
        # Handle the case where the value is not an integer
        print(f"Error converting memory value: {e}")
        available_memory_gb = psutil.virtual_memory().available / (1024 ** 3)  # Convert bytes to GB
        print(f"Available system memory: {available_memory_gb:.2f} GB")
        return available_memory_gb


def estimate_nworkers(df, complexity_factor=1):
    """
    Estimate the number of workers based on the number of rows and columns.

    Parameters:
    - df (DataFrame): The DataFrame to check.
    - complexity_factor (int or float): Adjusts the estimate based on processing complexity.

    Returns:
    - float: Estimated memory usage per worker in GB.
    """
    
    nrows, ncols = df.shape
    
    # Calculate factor based on nrows, ncols
    factor = (nrows + 0.1 * ncols) /8277
    
    # Estimate total memory based on rows, columns, and average memory per cell
    estimated_memory_per_worker_gb = 6**(factor)+1
    
    # Check allocated memory
    allocated_memory = get_memory_info()
    
    # Check if estimated memory per worker is larger than allocated memory
    if estimated_memory_per_worker_gb > allocated_memory:
        estimated_memory_per_worker_gb = allocated_memory
    
    # Estimate the number of workers
    nworkers = int(allocated_memory/estimated_memory_per_worker_gb * complexity_factor)
    
    print(f"Number of processes based on input data and available memory: {nworkers}")
    
    return nworkers




