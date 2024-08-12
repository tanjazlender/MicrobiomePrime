# MicrobiomePrime
MicrobiomePrime is a tool for identifying primer pairs with high specificity and sensitivity for a particular source of microbiota by analysing amplicon sequences. The pipeline was originally developed for use in Microbial Source Tracking (MST).
You can find more about MicrobiomePrime at: (link available soon)

## Contents
- Installation
- Data preprocessing
- Inputs
- Variables and settings
- Code overview
- Progress and error monitoring
- Definitions

## Installation
MicrobiomePrime is intended to be run in a x86-64 Linux OS (tested on Ubuntu). The best way to start is to create a conda environment with all the necessary dependencies using the provided environment.yml file in the requirements folder:
```
conda env create -f environment.yml
```
*If Conda is not installed on your system, you can find installation instructions [here](https://conda.io/projects/conda/en/latest/index.html).*

Once the environment named MicrobiomePrime is created, activate it using:
```
conda activate MicrobiomePrime
```
You also need to install ThermonucleotideBLAST, a program for *in silico* PCR. It can be installed following the [instructions](https://public.lanl.gov/jgans/tntblast/tntblast_doc.html) on their official page. For this program to work across multiple CPUs, we installed OpenMPI.

## Data preprocessing
Before searching for primers using MicrobiomePrime, you need need to preprocess raw amplicon sequences. This typically involves steps such as removal of primer sequences, quality filtering, chimera detection and removal, sequence clustering or denoising and some form of normalization. A sequence ID (SeqID) must be assigned to each sequence and must consist of letters first and then numbers.

Suitable software options data preprocessing include [Usearch](https://www.drive5.com/usearch/), [Qiime2](https://qiime2.org/), [DADA2](https://benjjneb.github.io/dada2/), and [Mothur](https://mothur.org/).

In the final step of data preprocessing, relative abundances of each sequence within each sample must be calculated. For details on the required format of the preprocessed file outputs (which serve as inputs for MicrobiomePrime analysis), please refer to the section titled "Inputs."

## Inputs
Make sure that the MicrobiomePrime input files are formatted to match our example. Do not forget to check file names, row names, column names and table formatting.

For the analysis itself, you will need the following four files:
1. Metadata file
2. Relative abundances table
3. FASTA file
4. Taxonomy file

The specific structure and format of each file are detailed below.

**1. Metadata**

*File name: metadata.csv or metadata.tsv*
This is a file that contains the information about the samples being analysed. To this analysis the key information is the source of microbiota.
If trying to find primers for tracking fecal pollution deriving from pigs, these sources can be for example human, cattle and pig feces.

Sample | Source 
--- | --- 
Sample1 | Human feces 
Sample2 | Cattle feces
Sample3 | Pig feces
Sample4 | Human feces

**2. Relative abundances table**

A table with relative abundances. Samples must be in rows and sequence IDs (e.g. names of OTUs, ZOTUs or ASVs) must be in columns. 
The sum of each row should be 1.

*File name: relabund_tab.csv or relabund_tab.tsv*

|        | Otu1   | Otu2   | Otu3   | Otu4   | Otu5   | Otu6   | Otu7   |
|--------|--------|--------|--------|--------|--------|--------|--------|
| Sample1| 0.1000 | 0.0200 | 0.8800 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample2| 0.0200 | 0.3500 | 0.0000 | 0.0000 | 0.2550 | 0.3750 | 0.0000 |
| Sample3| 0.0000 | 0.0920 | 0.4050 | 0.3253 | 0.0000 | 0.0000 | 0.1777 |
| Sample4| 0.4250 | 0.1005 | 0.0000 | 0.0000 | 0.2145 | 0.0600 | 0.2000 |

**3. FASTA file**

*File name: sequences.fa*
```plaintext
>otu1
TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGG
CCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTC
ATTGACG
>otu2
TGGGGAATATTGCACAATGGGCGAAAGCCTGATGCAGCAACGCCGCGTGAGCGATGAAGG
CCTTCGGGTCGTAAAGCTCTGTCCTCAAGGAAGATAATGACGGTACTTGAGGAGGAAGCC
>otu3
TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGATGACGG
CCTTCGGGTTGTAAAGCTCTGTCTTCAGGGACGATAATGACGGTACCTGAGGAGGAAGCC
ACGG
>otu4
TGGGGAATATTGCACAATGGGGGAAACCCTGATGCAGCAACGCCGCGTGAGTGATGAAGG
TTTTCGGATCGTAAAGCTCTGTCTTTGGGGAAGATAATGACGGTACCCAAGGA
>otu5
TAGGGAATCTTCGGCAATGGGGGCAACCCTGACCGAGCAACGCCGCGTGAGTGAAGAAGG
TTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGTGTGTGAGAGTGGAAAGTTCACA
>otu6
TGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCAACGCCGCGTGAGCGATGAAGG
CCTTCGGGTCGTAAAGCTCTGTCCTCAAGGAAGATAATGACGGTACTTGAGGAGGAA
>otu7
TAGGGAATCTTCGGCAATGGGGGCAACCCTGACCGAGCAACGCCGCGTGAGTGAAGAAGG
TTTTCGGATCGTAAAGCTCTGTTGTAAGAGAAGAACGTGTGTG
```

**4. Taxonomy**

*File name: taxonomy.csv or taxonomy.tsv*

SeqID | Domain | Phylum | Class | Order | Family | Genus  
--- | --- | --- | --- | --- | --- | ---
otu1 | Bacteria	| Firmicutes	| Clostridia	| Clostridiales	| Peptostreptococcaceae	| Romboutsia
otu2 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Carnobacteriaceae	| Catellicoccus
otu3 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Streptococcaceae	| Streptococcus
otu4 | Bacteria	| Bacteroidetes	| Bacteroidia	| Bacteroidales	| Prevotellaceae	| Prevotella
otu5 | Bacteria	| Fusobacteria	| Fusobacteriia	| Fusobacteriales	| Fusobacteriaceae | 
otu6 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Enterococcaceae	| Enterococcus
otu7 | Bacteria	| Bacteroidetes	| Bacteroidia	| Bacteroidales	| Bacteroidaceae	| 

Ensure that the column names in your file match our format. The first column should be named "SeqID".

## Variables and settings
You can change variables and settings in the variables.ini and settings.ini files found in scripts folder.

### Variables
Variable | Explanation | Example | Default value
--- | --- | --- | ---
target_group_name | The name of the target source group you are analysing (trying to design primers for). The source group name only needs to be set if you define more than one target sources. | source_group_name=Anatids|
target | Target source(s). If you define multiple target sources, separate them with commas. | target=Stork, Duck, Goose |
specificity_exception=Bird_unknown | Source(s) excluded from the specificity calculations. This is useful if you have samples of unknown origin that could include target samples. If you define multiple specificity exceptions, separate them with commas. | specificity_exception=Bird_unknown |
kmer_size | The size of K-mers (and primers) to be generated. | kmer_size=22 | kmer_size=22
kmer_sensitivity_cutoff | A minimum value for sensitivity of a K-mer to be used as a primer in in-silico PCR | kmer_sensitivity_cutoff=50 |
kmer_specificity_cutoff | A minimum value for specificity of a K-mer to be used as a primer. Only one of the primers (forward OR reverse) has to match specificity criteria. | kmer_specificity_cutoff=70 |
marker_sensitivity_cutoff | Minimum sensitivity of markers amplified with a given primer pair. This value is different from the kmer_sensitivity_cutoff and can be either the same or higher | marker_sensitivity_cutoff=60 |
marker_specificity_cutoff | Minimum specificity of markers amplified with a given primer pair. This value is different from the kmer_specificity_cutoff and can be either the same or higher | marker_specificity_cutoff=95 |
max_amplicon_length | Maximum length of the amplicons. | max_amplicon_length=150 | max_amplicon_length=2000
min_amplicon_length | Minimum length of the amplicons in the *in silico* PCR. | min_amplicon_length=70 | min_amplicon_length=0
max_primer_tm | The maximum allowed temperature (in °C) for a primer oligo to bind to a target sequence. | max_primer_tm=70 |  max_primer_tm=9999
min_primer_tm | The minimum allowed temperature (in °C) for a primer oligo to bind to a target sequence | min_primer_tm=50 | min_primer_tm=50
max_primer_delta | The maximum allowed delta G (in Kcal/Mole) for a primer oligo to bind a target sequence |max_primer_delta=-1 | max_primer_delta=9999
min_primer_delta | The minimum allowed delta G (in Kcal/Mole) for a primer oligo to bind a target sequence | min_primer_delta=-10 | min_primer_delta=-9999
max_mismatch | The maximum number of mismatches allowed in an oligonucleotide match. | max_mismatch=2 | max_mismatch=999
primer_clamp | Specifies the number of bases at the 3' end of each primer that must perfectly match the target sequence. | primer_clamp=2 | primer_clamp=0

*The examples in this table do not correspond with the example dataset, where there is only one target source - Pig feces.

### Settings
Setting | Explanation | Example | Default value
--- | --- | --- | ---
tntblast_path | The path of ThermonucleotideBLAST on your computer. | tntblast_path=/usr/bin/thermonucleotideBLAST/tntblast | 
cpus | Number of CPUs allocated for running the program. | cpus=200 | 
memory | Amount of RAM (in GB) available for this analysis. Ensure that you allocate sufficient memory based on your input data requirements. | memory=100 | 

## Code overview
T
The code consists of three main parts:
1. Generating K-mers
2. Creating primer pairs
3. Assessing the sensitivity and specificity of primer pairs in an *in silico* PCR analysis

### Section 1: Generating K-mers
In the first part, amplicon sequences are split into K-mers that are length of a primer. The sequences are split in a one bp window slide approach.
<p align="center">
  <img src="https://github.com/tanjazlender/MicrobiomePrime/assets/100705053/0300193e-dc1b-44b1-bc9f-6231b781fafb" alt="splitting kmers_small">
</p>

### Section 2: Creating primer pairs
Primers are essentially K-mers produced in Section 1. 
When designing primer pairs, we use two key thresholds:
- kmer_sensitivity_threshold
- kmer_specificity_threshold
For a primer pair to be considered valid, both primers must meet or exceed the kmer_sensitivity_threshold. Additionally, at least one of the primers must meet or exceed the kmer_specificity_threshold.

Here are some important factors to consider when setting thresholds for primer pairs:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




### Section 3: Assessing the sensitivity and specificity of primer pairs in an *in silico* PCR analysis
In the final section, we conduct an in silico PCR analysis using the primer pairs generated in Section 2. The main two parameters we calculate here are source sensitivity and specificity.
Source sensitivity measures how effectively the primer pair detects samples from the target source. Specificity, on the other hand, evaluates whether the primers also recognize sequences from nontarget microbiotas, ensuring they are not falsely detected in unrelated samples.
Although 100% sensitivity and 100% specificity would be ideal, it is often challenging to achieve in practice.

Write about marker_sensitivity_threshold and marker_specificity_threshold

**Important:** The primer pairs should always be validated in a laboratory on multiple target (where you want the primers to amplify) and non-target samples. If you are testing sensitivity and specificity of primer pairs, you can use e.g. end-point PCR, real-time PCR or digital PCR.

## Progress and error monitoring






## Definitions
Marker
Primer pair
Target samples
Non-target samples
Specificity
Sensitivity
K-mer

## License and third-party software
MicrobiomePrime is distributed under a ??? licence. Additionally, it redistributes the following third party software:
- [ThermonucleotideBLAST](https://public.lanl.gov/jgans/tntblast/tntblast_doc.html)

The licenses for all dependencies used in this pipeline are detailed in the NOTICE file.

MicrobiomePrime is developed by Tanja Zlender, Lucija Brezocnik and Vili Podgorelec.
For support, please contact tanja.zlender@nlzoh.si.


## Other - add to the readme
- log files
