# DNA_Protein_Alignment

# # Explanation of my approach

Input: DNA Sequence and protein Sequence

This program converts the DNA sequence into three forward reading frames and aligns them to the protein sequence using local alignment and blosum62 matrix. 
Affine gaps within frames are allowed, and eventually gaps between frames will also be allowed.

Right now the program is just attempting to align the first reading frame to the protein sequence, so a string is hardcoded for the DNA sequence
(But the logic for converting to the reading frames is finished)

protein.fasta contains the protein sequence, dna.fasta contains the DNA sequence (when not hardcoding)

affine.py has the Fill and Traceback code. Right now, one fill matrix and three pointer matrices are being used to try and align one reading frame to the protein sequence. A max score is calculated and a traceback is attempted. However, gaps are incorrect as it can end on a gap, and long sequences are not aligned properly.

The syntax of the affine alignment is smith_waterman_gotoh(seq1, seq2, seq3, prot_seq, scoring_table, go, ge, frame_pen)
