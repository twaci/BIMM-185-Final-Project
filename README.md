# BIMM-185-Final-Project

Introduction
------------
This pipeline generates possible sgRNAs for a user-inputted E. coli biological pathway and sgRNA length (from 17-22 bp). The generated sgRNAs are unique in the E. coli genome, contain a PAM site (NGG), and has an ideal GC content of 40-80%.

Requirements
------------
This tool requires the following:
* Python
* Biopython (for use of Entrez)

Files:
------
* ecoli_**mers.txt: Text files containing all the unique k-mers (of length 17-22) in the reference E. coli genome.
* ecoli_pathways.txt: Text file containing E. coli biological pathways and their corresponding genes.
* only_sequence.fasta: FASTA file with the sequence of the forward strand of the reference E. coli genome.
* only_sequence_comp.fasta: FASTA file with the sequence of the reverse strand of the reference E. coli genome.
* kmerCount.py: Program used to generate the ecoli_**mers.txt files.
* pathwayCRISPR.py: Main pipeline that generates sgRNAs based on user-inputted biological pathway and desired sgRNA length.

How to Use
----------
Run pathwayCRISPR.py with the command-line arguments to get the sgRNAs associated with a biological pathway:
python pathwayCRISPR.py -l LENGTH_OF_SGRNA -p BIOLOGICAL_PATHWAY_NAME

-l : length of sgRNA, from 17 to 22 inclusive
-p : name of biological pathway

To generate the list of biological pathways to choose from:
python pathwayCRISPR.py -lp

Additional help:
python pathwayCRISPR.py [-h | --help]
