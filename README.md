# MicrobiomePrime
A tool for identifying primer pairs with high sensitivity and specificity for a particular source of microbiota. It is designed for analysing amplicon sequencing data, however it could (with some modifications) also be used on shotgun sequencing data.
You can find the MicrobiomePrime paper at: (link)

## Contents
- About
- Installation
- Variables
- Code overview
- Primer pair validation
- Definitions

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
- a **software for generation of taxonomic units** (OTUs or ZOTUs) **or amplicon sequence variants** (ASVs). For this you can use [Usearch](https://www.drive5.com/usearch/) or other similar software such as [Qiime2](https://qiime2.org/), [DADA2](https://benjjneb.github.io/dada2/) and [Mothur](https://mothur.org/). We recommend to use Usearch because the outputs of this program are correctly formatted for our analysis.

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
rarefy_cutoff | Rarefaction depth. | rarefy_cutoff=10000 | rarefy_cutoff=10000
zotu_prevalence_cutoff | ZOTUs prevalent in less than this threshold [%] within a single source, it is treated as a low prevalent ZOTU | zotu_prevalence_cutoff=30 | zotu_prevalence_cutoff=30
rare_zotu_cutoff | If a ZOTU is low prevalent and has a relative abundance <= rare_zotu_cutoff, it is removed from the given sample | rare_zotu_cutoff=0.0001 | rare_zotu_cutoff=0.0001




## Code overview
The code is split into three main sections. In the first section, we analyse raw sequencing data and generate taxonomic units or sequence variants. In the second part, the filtered amplicon sequences are split into K-mers that are length of a primer. The sequences are split in a one bp window slide approach. Finally, selected K-mers linked to a specific microbiota source function as primers in an in silico PCR.  Primer pairs with highest sensitivity and specificity (or differential abundance) can potentially be used to amplify microbiota source-associated markers.

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
- [Usearch](https://www.drive5.com/usearch/)
The licenses for all dependencies used in this pipeline are detailed in the NOTICE file.

MicrobiomePrime is developed by Tanja Zlender, Lucija Brezocnik and Vili Podgorelec.
For support, please contact tanja.zlender@nlzoh.si.
