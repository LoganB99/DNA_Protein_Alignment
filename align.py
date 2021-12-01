import sys
import re
from convert_to_amino import *
from get_score import *
from affine import * 

#Open dna sequence file
dna_seq_file = sys.argv[1]
f = open(dna_seq_file, "r")
#Open protein sequence file
protein_seq_file = sys.argv[2]
prot_file = open(protein_seq_file, "r")

#Open codon to amino acid letter translation and Blosum62 scoring matrix
codons = open("codon.txt", "r")
blosum = open("blosum62.txt", "r")

#create scoring table
table = list()
for line in blosum:
    line = line.strip()
    table.append(line.split(","))

#dictionary of codons for dna conversion
codon_dict = dict()
for line in codons:
    cols = line.split("\t")
    if (cols[2] == "O"):
        cols[2] = "."
    codon_dict[cols[0]] = cols[2]
#Read DNA and Protein Fafsas
count = 0
dna_name = ""
dna_seq = ""

for line in f:
    if(not count):
        dna_name = line
    else:
        seq = line.strip()
        dna_seq += seq
    count += 1
count = 0
prot_name = ""
prot_seq = ""
for line in prot_file:
    if(not count):
        prot_name = line
    else:
        seq = line.strip()
        prot_seq += seq
    count += 1
reading_frames = convert_to_amino(dna_seq, codon_dict)
smith_waterman_gotoh(reading_frames[0],reading_frames[1], reading_frames[2], prot_seq, table, 1, 1, 1)

