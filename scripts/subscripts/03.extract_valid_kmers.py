# This script extracts valid kmers that satisfy the requirements
# and provide kmers corresponding ZOTUs

import pandas as pd
import configparser
import os
from Bio import SeqIO


#######################################################################
########### Retrieve variables and define files location ##############


config = configparser.ConfigParser()
config.read('../variables.ini')

target_group_name = config.get('settings','target_group_name')

if target_group_name !='':
    target_group_ID = target_group_name
else:
    target_group_ID = config.get('settings', 'target1')

# convert the comma-separated string to a list
target_group_ID = [name.strip() for name in target_group_ID.split(',')]

# define the length of the kmer. If nothing is set in the variables,
# set the length of kmer to 20
kmer_size = config.get('settings', 'kmer_size')
if kmer_size !='':
    kmer_size = int(config.get('settings', 'kmer_size'))
else:
    kmer_size = 20

# provide location of fasta files
fasta_files_location = '../../data/generated_files/sequences/fasta_files/'

# provide location of the output file(s)
output_files_location = '../../data/generated_files/kmers/'

# check if the folder exists, and if not, create it
if not os.path.exists(output_files_location):
    os.makedirs(output_files_location)




#######################################################################
#### Defined functions for kmers and corresponding ZOTUs extraction ###


# load content from FASTA file
def loadFasta(file):
    sample=[]
    for zotu in SeqIO.parse(file, "fasta"):
        sample.append(zotu)
    return sample

# check if kmer has repeats
def hasRepeats(kmer, size=4):
    for i in range(0,len(kmer)-size-size+1):
        if (kmer[i:i+size] in kmer[i+size:]):
            return True
    return False

# provide inverse kmer
def getInverseKmer(kmer):
    s=''
    for i in range(len(kmer)):
        if (kmer[i]=='C'):
            s='G'+s
        if (kmer[i]=='G'):
            s='C'+s
        if (kmer[i]=='A'):
            s='T'+s
        if (kmer[i]=='T'):
            s='A'+s
    return s

# check if kmer has repeats
def hasInverseRepeats(kmer, size=4):
    for i in range(0,len(kmer)-size-size+1):
        invkmer=getInverseKmer(kmer[i:i+size])
        if (invkmer in kmer[i+size:]):
            return True
    return False

# check if kmer is valid based on defined requirements
def isValidKmer(kmer):
    if ('CCCC' in kmer):
        return False
    if ('GGGG' in kmer):
        return False
    if ('AAAA') in kmer:
        return False
    if ('TTTT') in kmer:
        return False
    countGC=kmer.count('G')+kmer.count('C')
    if (countGC<9):
        return False
    if (countGC>13):
        return False
    if hasRepeats(kmer):
        return False
    if hasInverseRepeats(kmer):
        return False
    return True

# extract valid kmers from FASTA file
def extractValidKmers(sample, ksize=kmer_size):
    allKmers=0
    repeatedKmers=0
    nonValidKmers=0
    kmers=[]
    for s in range(len(sample)):
        seq=str(sample[s].seq)
        nKmers=len(seq)-ksize+1
        allKmers=allKmers+nKmers
        for i in range(nKmers):
            kmer=seq[i:i+ksize]
            if (isValidKmer(kmer)):
                if (not(kmer in kmers)):
                    kmers.append(kmer)
                else:
                    repeatedKmers=repeatedKmers+1
            else:
                nonValidKmers=nonValidKmers+1
    print('From all', allKmers, 'k-mers, there were', nonValidKmers, 'non-valid k-mers and',
          repeatedKmers, 'valid k-mers were repeated.')
    return kmers, allKmers

# main function of processing fasta files
def processSample(folder, species, verbose=False):
    file=folder+species
    print(f"Processing file'{file}'.fa", end=' \n')
    sample=loadFasta(file+'.fa')
    kmers,allKmers=extractValidKmers(sample, kmer_size)
    if verbose:
        print('Sample:')
        print('  zotus:',len(sample))
        print('  possible kmers:',allKmers)
        print('  extracted kmers:',len(kmers))
    sampleItem={'species':species, 'sample':sample, 'kmers':kmers, 'numKmers':len(kmers), 'allKmers':allKmers}
    return kmers, sampleItem

# prepare valid kmers
def prepareValidKmers(folder, species):
    kmers, item = processSample(folder, species, verbose=True)
    return kmers

# add ZOTU(s) to corresponding kmer
def prepareKmersWithZotus(kmers):
    df = pd.DataFrame(dict.fromkeys(['kmer', 'Zotus'], []))
    for kmer in kmers:
        zotus_list = []
        for zotu in singleline_fasta:
            if kmer in zotu.seq:
                zotus_list.append(zotu.name)
        df = pd.concat([df, pd.DataFrame({'kmer': [kmer], 'Zotus': [', '.join(zotus_list)]})], ignore_index=True)
    return df



#######################################################################
######### Main program that extracts kmers and prepares output ########


singleline_fasta = loadFasta(fasta_files_location+'sequences_singleline.fa')

# go through all fasta files (animals) and
# 1.) prepare valid kmers
# 2.) add ZOTUs to corresponding kmer
# 3.) save outputs in folder "kmers" for each animal
for s in target_group_ID:
    kmers = prepareValidKmers(fasta_files_location, s)
    df = prepareKmersWithZotus(kmers)
    df.to_csv(output_files_location+s+'_kmers.csv', sep=';', decimal=',', index=False)



