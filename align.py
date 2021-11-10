import sys
import re
from convert_to_amino import *
from get_score import *

dna_seq_file = sys.argv[1]
f = open(dna_seq_file, "r")
codons = open("codon.txt", "r")
blosum = open("blosum62.txt", "r")
table = list()
for line in blosum:
    line = line.strip()
    table.append(line.split(","))
codon_dict = dict()
for line in codons:
    cols = line.split("\t")
    if (cols[2] == "O"):
        cols[2] == "."
    codon_dict[cols[0]] = cols[2]
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
print(convert_to_amino(dna_seq, codon_dict))
get_score(".", "I", table)
print(max(2,2,2))
