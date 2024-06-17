# MicrobiomePrime
A tool for identifying primer pairs with high sensitivity and specificity for a particular source of microbiota. It is designed for analysing amplicon sequencing data.
You can find the MicrobiomePrime paper at: (link)

## Contents
- About
- Installation
- Variables
- Code overview
- Primer pair validation
- Definitions

## About
MicrobiomePrime is a tool for identifying identifying microbiome source-associated markers that can be detected using Polymerase Chain Reaction.
The pipeline was originally developed for use in Microbial Source Tracking (MST).

## Installation
MicrobiomePrime is intended to be run in a x86-64 Linux OS (tested on Ubuntu). The best way to start is to create a conda environment with all the necessary dependencies using the provided environment.yml file:
```
conda env create -f environment.yml
```
*If Conda is not installed on your system, you can find installation instructions [here](https://conda.io/projects/conda/en/latest/index.html). Alternatively, you can manually install all the dependencies in the environment.yml file.*

Once the environment named MicrobiomePrime is created, activate it using:
```
conda activate MicrobiomePrime
```
You also need to install ThermonucleotideBLAST, a program for *in silico* PCR. It can be installed following the [instructions](https://public.lanl.gov/jgans/tntblast/tntblast_doc.html) on their official page. For this program to work across multiple CPUs, we installed OpenMPI.

## Inputs
To begin, you need to generate operational taxonomic units or amplicon sequence variants. Suitable software options include [Usearch](https://www.drive5.com/usearch/), [Qiime2](https://qiime2.org/), [DADA2](https://benjjneb.github.io/dada2/), and [Mothur](https://mothur.org/).

For the analysis itself, you will need the following four files:
- Metadata file
- OTU or ASV table
- Taxonomy file
- FASTA file

The specific structure and format of each file are detailed below.

OTU or ASV table with relative abundance values. Samples must be in rows and OTUs/ASVs in columns. The sum of each row should be 1.
metadata.tsv
Sample | Source 
--- | --- 
Sample1 | Human feces 
Sample2 | Cattle feces
Sample3 | Pig feces
Sample4 | Human feces

otutab_relabund.tsv
|        | Otu1   | Otu2   | Otu3   | Otu4   | Otu5   | Otu6   | Otu7   |
|--------|--------|--------|--------|--------|--------|--------|--------|
| Sample1| 0.1000 | 0.0200 | 0.8800 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Sample2| 0.0200 | 0.3500 | 0.0000 | 0.0000 | 0.2550 | 0.3750 | 0.0000 |
| Sample3| 0.0000 | 0.0920 | 0.4050 | 0.3253 | 0.0000 | 0.0000 | 0.1777 |
| Sample4| 0.4250 | 0.1005 | 0.0000 | 0.0000 | 0.2145 | 0.0600 | 0.2000 |

taxonomy.tsv
Otu | Domain | Phylum | Class | Order | Family | Genus  
--- | --- | --- | --- | --- | --- | ---
otu1 | Bacteria	| Firmicutes	| Clostridia	| Clostridiales	| Peptostreptococcaceae	| Romboutsia
otu2 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Carnobacteriaceae	| Catellicoccus
otu3 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Streptococcaceae	| Streptococcus
otu4 | Bacteria	| Bacteroidetes	| Bacteroidia	| Bacteroidales	| Prevotellaceae	| Prevotella
otu5 | Bacteria	| Fusobacteria	| Fusobacteriia	| Fusobacteriales	| Fusobacteriaceae | 
otu6 | Bacteria	| Firmicutes	| Bacilli	| Lactobacillales	| Enterococcaceae	| Enterococcus
otu7 | Bacteria	| Bacteroidetes	| Bacteroidia	| Bacteroidales	| Bacteroidaceae	| 

Ensure that the column names in your file match our format. The first column should be named "Otu" regardless of whether you are analysing ASVs or OTUs.

otus.fa
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

## Variables
You can change the variables in the variables.ini file found in scripts folder.
Variable | Explanation | Example | Default value
--- | --- | --- | --- 
source1, source2, source3, source4, source5 | Target sources. | source1=Swan, source2=Wild duck, source3=Domestic duck, source4=Goose | /
source_group_name | The name of the target source group you are analysing (trying to design primers for). The source group name only needs to be set if you define more than one target sources. | source_group_name=Anatids | /
specificity_exception1, specificity_exception2, specificity_exception3, specificity_exception4, specificity_exception5 | Sources excluded from the specificity calculations | specificity_exception1=Bird unknown | /
kmer_sensitivity_cutoff | A minimum value for sensitivity of a K-mer to be used as a primer in in-silico PCR | kmer_sensitivity_cutoff=50 | /
kmer_specificity_cutoff | A minimum value for specificity of a K-mer to be used as a primer. Only one of the primers (forward OR reverse) has to match specificity criteria. | kmer_specificity_cutoff=70 | /
marker_sensitivity_cutoff | Minimum sensitivity of markers amplified with a given primer pair. This value is different from the kmer_sensitivity_cutoff and can be either the same or higher | marker_sensitivity_cutoff=60 | /
marker_specificity_cutoff | Minimum specificity of markers amplified with a given primer pair. This value is different from the kmer_specificity_cutoff and can be either the same or higher | marker_specificity_cutoff=95 | /
minimum_amplicon_length | Minimum length of an amplicon. Deafult value is set to 70 | min_amplicon_length=70 | /
cpus | The number of CPUs to run the program on. | /

## Code overview
T
The code is split into three main sections. In the first section, we analyse raw sequencing data and generate taxonomic units or sequence variants. In the second part, the filtered amplicon sequences are split into K-mers that are length of a primer. The sequences are split in a one bp window slide approach. Finally, selected K-mers linked to a specific microbiota source function as primers in an in silico PCR.  Primer pairs with highest sensitivity and specificity (or differential abundance) can potentially be used to amplify microbiota source-associated markers.

### Input
???

### Section 1: Generation of taxonomic units or sequence variants
In this section, we eliminate primer sequences and filter out low-quality reads before generating Zero Radius Operational Taxonomic Units (ZOTUs) using Usearch software. Alternatively you can use other bioinformatic pipelines for the analysis of amplicon sequence data such as Mothur, Qiime2 and DADA2 to generate Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs). If you are not using Usearch, make sure that the output files are formatted to match our example (check file names, row names and table formatting).

### Section 2: Generation of K-mers

<p align="center">
  <img src="https://github.com/tanjazlender/MicrobiomePrime/assets/100705053/0300193e-dc1b-44b1-bc9f-6231b781fafb" alt="splitting kmers_small">
</p>

### Section 3: In silico PCR and selection of best primer pairs
In the last section, we perform an in silico PCR using selected K-mers as primers. When looking for primer pairs that detect markers found in a target microbiota source and not in other sources, we calculate:
- Source sensitivity: Evaluating how many samples from a specific source can be detected using the given primer pair.
- Specificity to a Particular Source: Assessing whether the primers also detect microbiota from other sources.
This is useful in 
- Differential Abundance: Determining the strength of the association with the target source (is a marker more abundant in target than in non-target samples?)

## Validation
The primer pairs should always be validated in a laboratory on multiple target (where you want the primers to amplify) and non-target samples (the number of samples in the validation process depends on the study itself). If you are testing sensitivity and specificity of primer pairs, you can use e.g. end-point PCR, real-time PCR or digital PCR. If the differentiation of target and non-target samples is based on differential abundance of a particular marker, the amplification should be quantified (whether using real-time or digital PCR).



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
