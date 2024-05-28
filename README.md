# MicrobiomePrime
A tool for identifying primer pairs with high sensitivity and specificity for a particular source of microbiota
You can find the MicrobiomePrime paper at: (link)

## Contents
- About
- Installation

## About MicrobiomePrime
MicrobiomePrime is a tool for identifying identifying microbiome-associated markers that can be detected using Polymerase Chain Reaction.
The pipeline was originally developed for use in Microbial Source Tracking (MST), but can also be used for designing primers in medicine or environmental ecology.

## Installation
MicrobiomePrime is intended to be run in a x86-64 Linux OS (tested on Ubuntu). The easiest way to install all the necessary libraries is by using conda. 

This will create a new conda environment named MicrobiomePrime, which must always be activated before using the pipline:
conda activate MicrobiomePrime

## About the code
The code is split into three main parts.

### Part 1: Generation of taxonomic units
In this part we generate Zero Radius Operational Taxonomic Units (ZOTUs) using Usearch-UNOISE3 software. Alternatively you can use other bioinformatic pipelines for the analysis of amplicon sequence data such as Mothur, Qiime2 and DADA2 to generate Operational Taxonomic Units (OTUs) or Amplicon Sequence Variants (ASVs). Just be careful with the formatting of the output files ( check the row names, is the table transposed differently? etc.) as it has to be the same as in our example.

