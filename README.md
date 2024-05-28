# MicrobiomePrime
A tool for identifying primer pairs with high sensitivity and specificity for a particular source of microbiota. It is designed for analysing amplicon sequencing data, however it could (with some modifications) also be used on shotgun sequencing data.
You can find the MicrobiomePrime paper at: (link)

## Contents
- About
- Installation
- 

## About
MicrobiomePrime is a tool for identifying identifying microbiome-associated markers that can be detected using Polymerase Chain Reaction.
The pipeline was originally developed for use in Microbial Source Tracking (MST), but can also be used for designing primers in medicine or environmental ecology.

## Installation
MicrobiomePrime is intended to be run in a x86-64 Linux OS (tested on Ubuntu). The best way to start is to create a conda environment with all the necessary dependencies using the provided environment.yml file:
```
conda env create -f environment.yml
```
*If Conda is not installed on your system, you can find installation instructions [here](https://conda.io/projects/conda/en/latest/index.html). Alternatively, you can manually install all the dependencies in the environment.yml file.*

Once the environment named MicrobiomePrime is created, activate it using:
Now you created a new conda environment named MicrobiomePrime, which must always be activated before using the pipline:
```
conda activate MicrobiomePrime
```
You also need to install two programs:
- **ThermonucleotideBLAST**, a program for *in silico* PCR. It can be installed following the [instructions](https://public.lanl.gov/jgans/tntblast/tntblast_doc.html) on their official page. For this program to work across multiple CPUs, we installed OpenMPI.
- a **software for generation of taxonomic units** (OTUs or ZOTUs) **or amplicon sequence variants** (ASVs). We used Usearch for which you need to purchase a licence (there is a free version for very small data sets). If you would like to install Usearch, follow the instructions on their [official page](https://www.drive5.com/usearch/). Alternatively, you can use other software, such as Qiime2, DADA2 or Mothur.

## Code overview
The code is split into three main sections. In the first section, we analyse raw sequencing data and generate taxonomic units or sequence variants. In the second part, the filtered amplicon sequences are split into K-mers that are length of a primer. The sequences are split in a one bp window slide approach. Finally, selected K-mers linked to a specific microbiota source function as primers in an in silico PCR.  Primer pairs with highest sensitivity and specificity (or differential abundance) can potentially be used to amplify microbiota source-associated markers.

### Section 1: Generation of taxonomic units or sequence variants
In this section, we eliminate primer sequences and filter out low-quality reads before generating Zero Radius Operational Taxonomic Units (ZOTUs) using Usearch software. Alternatively you can use other bioinformatic pipelines for the analysis of amplicon sequence data such as Mothur, Qiime2 and DADA2 to generate Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs). If you are not using Usearch, make sure that the output files are formatted to match our example (check file names, row names and table formatting).

### Section 2: Generation of K-mers
![splitting kmers_small](https://github.com/tanjazlender/MicrobiomePrime/assets/100705053/0300193e-dc1b-44b1-bc9f-6231b781fafb)

### Section 3: In silico PCR and selection of best primer pairs
In the last section, we perform an in silico PCR using selected K-mers as primers. When looking for primer pairs that detect markers found in a target microbiota source and not in other sources, we calculate:
- Source sensitivity: Evaluating how many samples from a specific source can be detected using the given primer pair.
- Specificity to a Particular Source: Assessing whether the primers also detect microbiota from other sources.
This is useful in 
- Differential Abundance: Determining the strength of the association with the target source.
- differential abundance (based on how strong the association is)





Copyright (c) 2024 Zlender T. tanja.zlender@nlzoh.si (see LICENSE)
