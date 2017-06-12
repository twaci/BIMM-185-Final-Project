# Main pipeline that generates a list of sgRNAs based on user-inputted E. coli
# biological pathway and sgRNA length (between 17 and 22 bp inclusive)

import sys
import argparse
from Bio import Entrez
Entrez.email = 'y8qin@ucsd.edu'

# Function to get the reverse complement of a sequence
def reverseComp(sequence):
	seq = sequence[::-1]
	forward_seq = ""
	for i in seq:
		if i == "A": forward_seq += "T"
		elif i == "T": forward_seq += "A"
		elif i == "C": forward_seq += "G"
		elif i == "G": forward_seq += "C"
	return forward_seq

# Parse arguments kmer length and pathway
parser = argparse.ArgumentParser()
parser.add_argument('-l', action = "store", type = int, dest = "kmer_length", help = "Length of sgRNA")
parser.add_argument('-p', action = "store", dest = "pathway", help = "Biological Pathway Code")
parser.add_argument('-lp', action = "store_true", dest= "list_pathways", help = "List Biological Pathways")
args = parser.parse_args()
kmer_length = args.kmer_length
pathway = args.pathway

# Correct usage of program
if(len(sys.argv) != 5 and len(sys.argv) != 2):
	print "Usage: python pathwayCRISPR.py -l <length of sgRNA> -p <Biological Pathway Code>"
	print "Usage: python pathwayCRISPR.py -h"
	sys.exit()

# Status to help user if -lp entered
helpUser = False
if args.list_pathways:
	helpUser = True

num_tab = 0
with open("ecoli_pathways.txt", "r") as input:
	txt = input.read().split('\n')

# Extract genes from pathways
geneinPath = dict()
pathName = dict()
for i in txt:
	cur = i.split('\t')
	pathCode = cur[0]
	name = cur[1]
	pathName[pathCode] = name
	cur = cur[2::]
	geneinPath[pathCode] = cur
	if cur[len(cur)-1] is '':
		cur = cur[:-1]
	geneinPath[pathCode] = cur

# List pathways if user enters -lp
if helpUser is True:
	for p in pathName:
		print p + "\t" + pathName[p]
	sys.exit()

# Check kmer_length
allowed = [17, 18, 19, 20, 21, 22]
if kmer_length not in allowed:
	print "ERROR: Invalid sgRNA length."
	sys.exit()

# Check whether user entered pathway code is valid or not.
# If not valid, print error message and exit program.
# Else, get the genes inside the pathway.
if pathway in geneinPath:
	genes = geneinPath[pathway]
else:
	print "ERROR: Invalid Pathway."
	sys.exit()

count = 0 # count of sequences to return
all_pos = [] # positions of gene sequences
genome_length = 4641652 # length of E. coli genome
comp_or_forward = {} # store if gene position is from forward or reverse strand
comp = [] # keep track of whether positions are complement
geneName = [] # list of gene names

# Print pathway information and genes inside
print pathway + ": " + pathName[pathway]
print "Genes in the pathway: " + ', '.join(str(x) for x in genes)
print "======== Gene Annotation ========"

# Parse through genes and get their corresponding positions & sequences
for i in genes:
	item = i
	animal = 'Escherichia coli'
	search_string = item+"[All Fields] AND (alive[prop]) AND Escherichia coli[Organism]"

	handle = Entrez.esearch(db="gene", term=search_string)
	record = Entrez.read(handle)
	ids = record['IdList']
	if len(ids) == 0:
		continue
	seq_id = ids[0]
	handle = Entrez.efetch(db="gene", id=seq_id, rettype="fasta", retmode="txt")
	records = handle.read()
	r = records.split('\n')
	# Record positions of genes in genome and if from forward or reverse strand
	for entry in r:
		if "Annotation" in entry:
			annotation = entry.split()
			anno = annotation[1:]
			if "complement)" in annotation:
				comp_or_forward[i] = True
			else:
				comp_or_forward[i] = False
			print i + "\t" + "\t".join(str(x) for x in anno)
			positions = annotation[2].strip('(').strip(',').strip(')').split('..')
			if comp_or_forward[i]:
				positions[0] = genome_length - int(positions[0])
				positions[1] = genome_length - int(positions[1])
				comp.append(True)
				comp.append(True)
				geneName.append(i)
				all_pos.append(positions[1])
				all_pos.append(positions[0])
			else:
				positions[0] = int(positions[0]) - 1
				positions[1] = int(positions[1]) - 1
				comp.append(False)
				comp.append(False)
				geneName.append(i)
				all_pos.append(positions[0])
				all_pos.append(positions[1])

with open("only_sequence_comp.fasta", "r") as ref:
	reference_comp = ref.read().strip()

with open("only_sequence.fasta", "r") as ref:
	reference = ref.read().strip()

# Find corresponding gene sequences from reference genome
gene_sequences = []
i = 0
while i < len(all_pos)-1:
	if comp[i] and comp[i+1]:
		forward_seq = reverseComp(reference_comp[all_pos[i]:all_pos[i+1]+1])
		gene_sequences.append(forward_seq)
	else:
		forward_seq = reference[all_pos[i]:all_pos[i+1]+1]
		gene_sequences.append(forward_seq)
	i += 2

# Filter sequences based on uniqueness and GC content
k = kmer_length+3 # length of kmer
kmerlist = [] # list of kmers to return
gene_kmers = [] # unique kmers
for seq in gene_sequences:
	for i in range(0, len(seq) - k + 1):
		GC = 0
		kmer = seq[i:(i+k)]
		for index in range(0, len(kmer)):
			if kmer[index] == "G" or kmer[index] == "C":
				GC += 1
		if kmer not in gene_kmers and kmer[len(kmer)-2:len(kmer)] == "GG":
			if 1.0*GC / kmer_length <= 0.80 and 1.0*GC / kmer_length >= 0.40:
				gene_kmers.append(kmer)
	kmerlist.append(gene_kmers)
	gene_kmers = []

# Parse through unique kmers
if kmer_length == 17:
	with open("ecoli_17mers.txt", "r") as ref:
		reference_kmers = ref.read().split()
elif kmer_length == 18:
	with open("ecoli_18mers.txt", "r") as ref:
		reference_kmers = ref.read().split()
elif kmer_length == 19:
	with open("ecoli_19mers.txt", "r") as ref:
		reference_kmers = ref.read().split()
elif kmer_length == 20:
	with open("ecoli_20mers.txt", "r") as ref:
		reference_kmers = ref.read().split()
elif kmer_length == 21:
	with open("ecoli_21mers.txt", "r") as ref:
		reference_kmers = ref.read().split()
else:
	with open("ecoli_22mers.txt", "r") as ref:
		reference_kmers = ref.read().split()

print "======== sgRNA Design ========"
check = False
count = 0
for k in range(len(kmerlist)):
	print "*** " + geneName[k] + " unique sgRNAs: ***"
	for gene in kmerlist[k]:
		if gene in reference_kmers:
			print gene[:-3] + " " + gene[-3:]
			check = True
			count += 1

if check is False:
	print "No unique sgRNA found."

else:
	print "Number of sgRNAs found: " + str(count)
