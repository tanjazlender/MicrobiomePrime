# MicrobiomePrime
MicrobiomePrime is a tool for identifying primer pairs with high specificity and sensitivity for a particular source of microbiota by analysing amplicon sequences. The pipeline was originally developed for use in Microbial Source Tracking (MST).
You can find more about MicrobiomePrime at: (link available soon)

## Contents
- Installation and requirements
- Data preprocessing
- Inputs
- Variables and settings
- Code overview
- Running the code
- Progress and error monitoring
- Outputs
- Definitions

## Installation and requirements
MicrobiomePrime is intended to be run in a x86-64 Linux OS (tested on Ubuntu). 

**1. Create a conda environment**

Use the provided environment.yml file in the `requirements` folder to create a Conda environment with all necessary dependencies:

```
conda env create -f environment.yml
```
*If Conda is not installed on your system, you can find installation instructions [here](https://conda.io/projects/conda/en/latest/index.html).*

**2. Activate the environment**

Once the environment named MicrobiomePrime is created, activate it using:
```
conda activate MicrobiomePrime
```

**3. Install ThermonucleotideBLAST**

ThermonucleotideBLAST is a tool for *in silico* PCR that must be installed separately. Follow the installation instructions on the [official ThermonucleotideBLAST page](https://public.lanl.gov/jgans/tntblast/tntblast_doc.html).

**4. Install OpenMPI (Optional but highly recommended)**

OpenMPI is used for parallelization to speed up the analysis. While ThermonucleotideBLAST can run without OpenMPI, it will take significantly longer. To install OpenMPI, follow the instructions on the [OpenMPI website](https://www.open-mpi.org/).

>Note: Due to the complexity and high computational demands of sequencing data, using a High Performance Computing (HPC) system may be necessary to ensure efficient processing and analysis of large datasets.

## Data preprocessing
Before using MicrobiomePrime, preprocess your raw amplicon sequences by removing primer sequences, quality filtering, detecting and removing chimeras, clustering or denoising sequences, and normalizing the data. Ensure each sequence is assigned a unique sequence ID (SeqID) with letters followed by numbers (e.g. OTU001, OTU002, OTU003).

Suitable software options data preprocessing include [Usearch](https://www.drive5.com/usearch/), [Qiime2](https://qiime2.org/), [DADA2](https://benjjneb.github.io/dada2/), and [Mothur](https://mothur.org/).

In the final step of data preprocessing, calculate the relative abundances of each sequence within each sample. Ensure that the preprocessed file outputs adhere to the required format specified in the "Inputs" section.

## Inputs

For the analysis, you will need the following four files:

<details>
<summary> Metadata file</summary>

File name: `metadata.csv` or `metadata.tsv`

This file contains the information about the samples being analysed. 

The metadata file should include columns labeled `Sample`, which represents the names of the samples, and `Source`, which indicates the source of the microbiota for each sample. 

Example metadata table:
Sample | Source 
--- | --- 
Sample1 | Human feces 
Sample2 | Cattle feces
Sample3 | Pig feces
Sample4 | Human feces

</details>


<details>
<summary> Relative abundance table</summary>

File name: `relabund_tab.csv` or `relabund_tab.tsv`

This table provides the relative abundances of various microbial sequences (identified as OTUs, ASVs, or ZOTUs or any other sequence ID) across different samples. 
Each row corresponds to a specific sample, and each column represents a unique sequence ID. 

The sum of each row should be 1.

Example table with relative abundances:

|        | Otu1   | Otu2   | Otu3   | Otu4   | Otu5   | Otu6   | Otu7   |
|--------|--------|--------|--------|--------|--------|--------|--------|
| Sample1| 0.1000 | 0.0200 | 0.8800 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample2| 0.0200 | 0.3500 | 0.0000 | 0.0000 | 0.2550 | 0.3750 | 0.0000 |
| Sample3| 0.0000 | 0.0920 | 0.4050 | 0.3253 | 0.0000 | 0.0000 | 0.1777 |
| Sample4| 0.4250 | 0.1005 | 0.0000 | 0.0000 | 0.2145 | 0.0600 | 0.2000 |

</details>



<details>
<summary> FASTA file</summary>

File name: `sequences.fa`

This file contains the DNA sequences detected in the samples in a FASTA format. 

Each sequence is associated with a unique identifier or sequence ID (e.g. OTU001, OTU002, OTU003,...).

Example metadata file:

```plaintext
>Otu1
TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGG
CCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTC
ATTGACG
>Otu2
TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAGCGATGAAGG
CCTTCGGGTCGTAAAGCTCTGTCCTCAAGGAAGATAATGACGGTACTTGAGGAGGAAGCC
>Otu3
TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGATGACGG
CCTTCGGGTTGTAAAGCTCTGTCTTCAGGGACGATAATGACGGTACCTGAGGAGGAAGCC
ACGG
>Otu4
TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGATGAAGG
TTTTCGGATCGTAAAGCTCTGTCTTTGGGGAAGATAATGACGGTACCCAAGGA
>Otu5
TAGGGAATCTTCGGCAATGGGGGCAACCCTGACCGAGCAACGCCGCGTGAGTGAAGAAGG
TTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGTGTGTGAGAGTGGAAAGTTCACA
>Otu6
TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCAACGCCGCGTGAGCGATGAAGG
CCTTCGGGTCGTAAAGCTCTGTCCTCAAGGAAGATAATGACGGTACTTGAGGAGGAA
>Otu7
TAGGGAATCTTCGGCAATGGGGGCAACCCTGACCGAGCAACGCCGCGTGAGTGAAGAAGG
TTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGTGTGTG
```

The sequence IDs in the FASTA file should match the sequence IDs in the relative abundance and taxonomy tables.

</details>


<details>
<summary> Taxonomy file</summary>

File name: `taxonomy.csv` or `taxonomy.tsv`

This file contains the taxonomic classification of the sequences identified in the samples. 
Each sequence ID is assigned to different taxonomic ranks, such as domain, phylum, class, order, family, and genus.

The table should contain the following column names: `SeqID`, `Domain`, `Phylum`, `Class`, `Order`, `Family` and `Genus`. If a sequence cannot be identified at a particular taxonomic level, leave the corresponding cell empty. 

SeqID | Domain | Phylum | Class | Order | Family | Genus  
--- | --- | --- | --- | --- | --- | ---
otu1 | Bacteria	| Firmicutes	| Clostridia	| Clostridiales	| Peptostreptococcaceae	| Romboutsia
otu2 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Carnobacteriaceae	| Catellicoccus
otu3 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Streptococcaceae	| Streptococcus
otu4 | Bacteria	| Bacteroidetes	| Bacteroidia	| Bacteroidales	| Prevotellaceae	| Prevotella
otu5 | Bacteria	| Fusobacteria	| Fusobacteriia	| Fusobacteriales	| Fusobacteriaceae | 
otu6 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Enterococcaceae	| Enterococcus
otu7 | Bacteria	| Bacteroidetes	| Bacteroidia	| Bacteroidales	|	| 

The sequence IDs in this file should match those in the relative abundances table and FASTA file.

</details>

Ensure that the input files are formatted according to our specifications. Carefully verify the file names, row names, column names, and table formatting.

## Variables and settings

To adjust the MicrobiomePrime analysis to your specific needs, you can adjust variables and settings in the `variables.ini` and `settings.ini` files located in the `scripts` folder.

#### Variables
Key variables that need to be defined are: 

- `target`
- `kmer_size`
- `kmer_sensitivity_cutoff`
- `kmer_specificity_cutoff`
- `marker_sensitivity_cutoff`
- `marker_specificity_cutoff`
- `max_mismatch` and/or `max_primer_delta`

The rest of the variables can stay undefined.

<details>
  <summary>See all MicrobiomePrime variables and their descriptions</summary>

| Variable                 | Explanation                                                                                               | Example                              | Default Value |
|--------------------------|-----------------------------------------------------------------------------------------------------------|--------------------------------------|---------------|
| `target_group_name`      | Name of the target source group to be analyzed. Only needs to be set if you define more than one target sources. | `source_group_name=Anatids`          | /             |
| `target`                 | Target source(s). Multiple sources must be separated with commas.                                         | `target=Stork, Duck, Goose`          | /             |
| `specificity_exception`  | Source(s) excluded from specificity calculations. Multiple sources must be separated with commas.         | `specificity_exception=Bird_unknown` | /             |
| `kmer_size`              | Size of K-mers (and primers) to be generated.                                                             | `kmer_size=22`                       | `22`          |
| `kmer_sensitivity_cutoff`| Minimum sensitivity value for a K-mer to be used as a primer in in-silico PCR.                            | `kmer_sensitivity_cutoff=50`         | /             |
| `kmer_specificity_cutoff`| Minimum specificity value for a K-mer to be used as a primer. At least one primer (forward or reverse) must meet this criteria. | `kmer_specificity_cutoff=70`         | /             |
| `marker_sensitivity_cutoff` | Minimum sensitivity of markers amplified with a given primer pair.  | `marker_sensitivity_cutoff=60`       | /             |
| `marker_specificity_cutoff` | Minimum specificity of markers amplified with a given primer pair.  | `marker_specificity_cutoff=95`       | /             |
| `max_amplicon_length`    | Maximum length of the amplicons.                                                                         | `max_amplicon_length=150`            | `2000`        |
| `min_amplicon_length`    | Minimum length of the amplicons.                                                        | `min_amplicon_length=70`             | `0`           |
| `max_primer_tm`          | Maximum allowed melting temperature (in °C) for a primer to bind to the target sequence.                   | `max_primer_tm=70`                   | `9999`        |
| `min_primer_tm`          | Minimum allowed melting temperature (in °C) for a primer to bind to the target sequence.                   | `min_primer_tm=50`                   | `50`          |
| `max_primer_delta`       | Maximum allowed delta G (in kcal/mol) for a primer to bind to the target sequence.                        | `max_primer_delta=-1`                | `9999`        |
| `min_primer_delta`       | Minimum allowed delta G (in kcal/mol) for a primer to bind to the target sequence.                        | `min_primer_delta=-10`               | `-9999`       |
| `max_mismatch`           | Maximum number of mismatches allowed in an oligonucleotide match.                                        | `max_mismatch=2`                     | `999`         |
| `primer_clamp`           | Number of bases at the 3' end of each primer that must perfectly match the target sequence.               | `primer_clamp=2`                     | `0`           |

*The examples in this table do not correspond with the example dataset, where there is only one target source - Pig feces.
</details>

>Note: The `kmer_sensitivity_cutoff` and `marker_sensitivity_cutoff` are related but serve different purposes. The `kmer_sensitivity_cutoff` is used prior to creating primer pairs, determining the minimum proportion of target samples that must contain the K-mer for it to be considered a valid primer. The `marker_sensitivity_cutoff` applies to the final marker and reflects the sensitivity of the entire marker amplification, which is influenced by the performance of both (forward and reverse) primers.

>Similarly, the `kmer_specificity_cutoff` and `marker_sensitivity_cutoff` are related but serve different purposes. The `kmer_specificity_cutoff` is used prior to creating primer pairs and determines the minimum proporion of non-target samples that do not contain the K-mer for it to be considered a valid primer. The `marker_specificity_cutoff` applies to the final marker and reflects the specificity of the entire marker amplification, which is influenced by the performance of both (forward and reverse) primers.

#### Settings

In the settings file:

- Define `tntblast_path` unless it is accessible via system's `PATH`.
- Define `cpus` and `memory` only when running the code via Slurm.

<details>
  <summary>See all MicrobiomePrime settings and their descriptions.</summary>

MicrobiomePrime settings:
Setting | Explanation | Example 
--- | --- | --- 
`tntblast_path` | The path of ThermonucleotideBLAST on your computer. Leave empty if the executable is accessible via system's `PATH`.  | `tntblast_path=/usr/bin/thermonucleotideBLAST/tntblast` 
`cpus` | Number of CPUs allocated for running the program. | `cpus=250`
`memory` | Amount of RAM (in GB) available for this analysis. Ensure that you allocate sufficient memory based on your input data requirements. | `memory=500` 

</details>


## Code overview

The code consists of four main sections:

<details>
<summary> Section 1: Preparing the data</summary>

In this section, the data is prepared for K-mer generation.

It consists of the following subscripts found in the `scripts/subscripts` folder:
- 00.verify_input_formatting.R: verifies whether the input is formatted correctly.
- 01.write_target_seqIDs.R: creates a list of sequence IDs found in target samples.
- 02.extract_fasta_files.py: extract FASTA sequences based on the list of sequence IDs found in target samples.
- 03.organize_directories.py: creates and organizes directories for storing analysis outputs.

</details>


<details>
<summary> Section 2: Generating K-mers</summary>

In the second section, amplicon sequences are split into K-mers that are length of a primer. The sequences are split in a one bp window slide approach as shown on the picture below.
<p align="center">
  <img src="https://github.com/tanjazlender/MicrobiomePrime/assets/100705053/0300193e-dc1b-44b1-bc9f-6231b781fafb" alt="splitting kmers_small">
</p>

This section encompasses the following scripts found in the `scripts/subscripts` folder.

</details>

#### Section 2: Creating primer pairs
Primers are essentially K-mers produced in Section 1. 
When designing primer pairs, we use two key cutoffs:
- kmer_sensitivity_cutoff
- kmer_specificity_cutoff
  
For a primer pair to be considered valid, both primers must meet or exceed the kmer_sensitivity_cutoff. Additionally, at least one of the primers must meet or exceed the kmer_specificity_cutoff.

#### Section 3: Assessing the sensitivity and specificity of primer pairs in an *in silico* PCR analysis
In the final section, we conduct an in silico PCR analysis using the primer pairs generated in Section 2. The main two parameters we calculate here are source sensitivity and specificity.
Source sensitivity measures how effectively the primer pair detects samples from the target source. Specificity, on the other hand, evaluates whether the primers also recognize sequences from nontarget microbiotas, ensuring they are not falsely detected in unrelated samples.
Although 100% sensitivity and 100% specificity would be ideal, it is often challenging to achieve in practice.

## Running the code
**1. Navigate to the `scripts` directory:**
```
cd scripts/
```
**2. Run the script:**

- To run the script **using Slurm**, execute the following command::  
     ```bash
     bash submit_slurm_job.sh
     ```

- **(For small datasets)** To run the script **directly in a Linux environment** without Slurm, use:  
     ```bash
     bash find_primer_pairs.sh
     ```

## Progress and error monitoring



## Outputs
MicrobiomePrime generates a primary results table detailing the performance of primer pairs, along with several additional tables containing information about the detected sequence IDs.

<details>
  <summary>Primary results table</summary>

The primary results table provides detailed metrics on primer pair performance, including specificity, sensitivity, and abundance measures. 

Description of column names in the main results table detailing the performance of primer pairs:
| **Column name**                      | **Definition**                                                                                                    | **Example**                                                                                     |
|--------------------------------------|-------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| **PP_ID**                            | The identifier (name) of a primer pair.                                                                            | Anatids101Fw:Anatids117Rv                                                                       |
| **Specificity**                      | Primer pair specificity [%].                                                                                       | 97.92                                                                                           |
| **Specificity2**                     | Primer pair specificity [number of non-detected non-target samples/number of nontarget samples].                   | 470/480                                                                                         |
| **Sensitivity**                      | Primer pair sensitivity [%].                                                                                       | 94.44                                                                                           |
| **Sensitivity2**                     | Primer pair sensitivity [number of detected target samples/number of target samples].                              | 58/60                                                                                           |
| **Sensitivity2_detailed**            | Primer pair sensitivity of each target host. Relevant when analyzing multiple target hosts [source (number of detected target samples from source/number of all samples from source)]. | Duck (30/30), Swan (29/30), Goose (26/30) |
| **N_detected_nontarget_samples**     | Number of detected nontarget samples [source (number of detected samples from source/number of all samples from source)]. | Chicken (3/30), Pigeon (7/30) |
| **Percent_abundance_target**         | The average relative abundance of target sequences amplified by the primer pair [% abundance, SD].                 | 2.2333 SD=1.0550                                                                 |
| **Percent_abundance_nontarget**      | The average relative abundance of nontarget sequences amplified by the primer pair [% abundance, SD]               | 0.5051 SD=0.0530                                                                 |
| **Percent_abundance_target_detailed**| The average relative abundance of target sequences amplified by the primer pair for each source [source (% abundance, SD)]. | Duck (2.8500, SD=1.024), Swan (2.0325, SD=1.1203), Goose (1.7029, SD=0.9588) |
| **Percent_abundance_nontarget_detailed** | The average relative abundance of nontarget sequences amplified by the primer pair for each source [source (% abundance, SD)]. | Chicken (0.1923, SD=0.5207), Pigeon (0.7942, SD=0.6032) |
| **Taxonomy_target**                  | The taxonomy of target sequences amplified by the primer pair, including multiple taxa if applicable. If the taxonomic classification does not reach the genus level, "(unknown)" is appended to the taxa name. | Catellicoccus                                                                                  |
| **Taxonomy_nontarget**               | The taxonomy of nontarget sequences amplified by the primer pair, including multiple taxa if applicable. If the taxonomic classification does not reach the genus level, "(unknown)" is appended to the taxa name. | Catellicoccus, Bacteria (unknown)                                                              |
| **N_seqIDs_target**                  | The number of unique target sequences amplified by the primer pair.                                                | 5                                                                                               |
| **Positive_target_samples**          | The list of target sequence IDs that were amplified by the primer pair [source (target sequence IDs)].             | Duck (AF001, AF002, AF003,…), Swan (AF004, AF006, AF007,...), Goose (AF008, AF010, AF014,...)   |
| **Positive_nontarget_samples**       | The list of non-target sequence IDs that were amplified by the primer pair [source (target sequence IDs)].         | Chicken (AF095, AF096, AF100), Pigeon (AF121, AF122, AF123, AF127, AF152, AF158, AF165)         |
| **Negative_target_samples**          | The list of target sequence IDs that were not amplified by the primer pair [source (target sequence IDs)].         | Swan (AF005), Goose (AF009, AF011, AF012, AF013)                                                |
| **Exceptions**                       | Sources that are neither classified as target nor non-target. Therefore, they are excluded from sensitivity, specificity and mean abundance calculations. | Bird unknown                                                                                   |
| **N_detected_exception_samples**     | Number of detected samples that are neither classified as target nor non-target [source (number of detected samples from source/number of all samples from source)]. | Bird unknown (1/7)                                                                             |
| **Percent_abundance_exceptions_detailed** | The average relative abundance of exception sequences (classified neither as target nor non-target) amplified by the primer pair for each source [source (% abundance, SD)]. | Bird unknown (2.5783) |

The primary results table can be found in path: `MicrobiomePrime/out/{Source}/sens{kmer_sensitivity_cutoff}_spec{kmer_specificity_cutoff}/final_table/`.

</details>

<details>
  <summary>Detected sequences tables</summary>

The output tables provide detailed information about detected sequences for each primer pair. 

Description of column names in the output tables containing the information about detected sequences. 
| **Column name**          | **Definition**                                                                                      | **Example**                                                                                      |
|--------------------------|-----------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------|
| **PP_ID**                | The identifier (name) of a primer pair.                                                             | Anatids101Fw:Anatids117Rv                                                                        |
| **SeqIDs_all**           | List of all sequence IDs that were amplified by the primer pair in both target and non-target samples. | Otu0001, Otu0002, Otu0003, Otu0004, Otu0005, Otu0006, Otu0007, Otu0008, Otu0009                  |
| **SeqIDs_target_only**   | List of sequence IDs that were amplified by the primer pair in target samples only.                 | Otu001, Otu0002, Otu0003, Otu0004, Otu005, Otu006                                                |
| **SeqIDs_nontarget_only**| List of sequence IDs that were amplified by the primer pair in non-target samples only.             | Otu0008                                                                                          |
| **N_seqIDs_target**      | The number of sequence IDs from non-target sequences that were amplified by the primer pair.        | 8                                                                                                |
| **N_seqIDs_nontarget**   | The number of sequence IDs from non-target sequences that were amplified by the primer pair.        | 3                                                                                                |
| **File_number**          | The specific number assigned to each file in the analysis. This number helps locate intermediate results generated during the analysis. | 6 |

The detected sequences tables can be found in path: `MicrobiomePrime/out/{Source}/sens{kmer_sensitivity_cutoff}_spec{kmer_specificity_cutoff}/detected_sequences/`.

</details>

>**Note:** The primer pairs should always be validated in the laboratory with both target and non-target samples. Test the sensitivity and specificity of the primer pairs using methods such as conventional PCR, real-time PCR, or digital PCR.

## Definitions

- **K-mer** - A short nucleotide sequence, typically between 18-24 bases long (the length of an optimal primer pair).
- **Marker** - any genetic sequence or a group of sequences that are amplified by a specific primer pair and are highly associated with a particular microbiota source.
- **Microbiota source** - the environment from which a community of microorganisms originates, such as pig gut, human skin or soil.
- **Non-target samples** - samples from non-target sources.
- **Nontarget source(s)** - the environment(s) or host(s) that are not of primary interest and are used to test the specificity of primer pairs, ensuring that they do not amplify markers from these sources.
- **Primer pair** - two short nucleotide sequences designed to bind to specific regions of a DNA sequence to facilitate the amplification of target genetic material during PCR.
- **Specificity** - the proportion non-target samples that was not detected using the primer pair, indicating the primer's ability to avoid false positives.
- **Specificity exception(s)** - source(s) not included in specificity calculations.
- **Sensitivity** - the proportion of target samples that was detected using the primer pair, indicating the primer's ability to identify true positives.
- **Target samples** - samples from target sources.
- **Target source(s)** - the environment(s) or host(s) of interest for which we aim to develop highly specific and sensitive primer pairs.


## License and third-party software
MicrobiomePrime is distributed under a ??? licence. Additionally, it redistributes the following third party software:
- [ThermonucleotideBLAST](https://public.lanl.gov/jgans/tntblast/tntblast_doc.html)
- [OpenMPI](https://www.open-mpi.org/)

The licenses for all dependencies used in this pipeline are detailed in the NOTICE file.

MicrobiomePrime is developed by Tanja Zlender, Lucija Brezocnik, Vili Podgorelec and Maja Rupnik.
For support, please contact tanja.zlender@nlzoh.si.


## Other - add to the readme
- log files
