#!/bin/bash

# Exit script if any command fails and handle errors in pipelines
set -e
set -o pipefail

# Load settings from settings.ini
CPUS=$(awk -F '=' '/cpus/ {print $2}' settings.ini)
MEMORY=$(awk -F '=' '/memory/ {print $2}' settings.ini)
TIME=$(awk -F '=' '/time/ {print $2}' settings.ini)

# Define log file paths
current_log_path="current_logs"
default_log_path="${current_log_path}/find_primer_pairs.log"

# Create log directory if it does not exist
mkdir -p "$current_log_path"

# Submit the job using sbatch with dynamically allocated resources
sbatch --job-name=MicrobiomePrime \
       --output=/dev/null \
       --error=/dev/null \
       --mem=${MEMORY}G \
       --ntasks=${CPUS} \
       --cpus-per-task=1 \
       --time=${TIME} \
       --output=/dev/null \
       --error=/dev/null \
       find_primer_pairs.sh
