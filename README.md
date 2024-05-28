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

## About the code
The code is split into three main parts.

### Part 1: Generation of taxonomic units or sequence variants
In this part we generate Zero Radius Operational Taxonomic Units (ZOTUs) using Usearch-UNOISE3 software. Alternatively you can use other bioinformatic pipelines for the analysis of amplicon sequence data such as Mothur, Qiime2 and DADA2 to generate Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs). Just be careful with the formatting of the output files ( check the row names, is the table transposed differently? etc.) as it has to be the same as in our example.

Copyright (c) 2024 Zlender T. tanja.zlender@nlzoh.si (see LICENSE)
