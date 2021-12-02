# DNA_Protein_Alignment

# Explanation of my approach
This program first converts the DNA sequence into three forward reading frames. I accomplished this by reading dna.fasta file and taking codons of three bases and converting them to their corresponding amino acids with a data structure created with the help of codon.txt. It has three starting points, hence the three forward reading frames.

I create the scoring matrix using a Blosum62 scoring table, allowing for stop codons and unknown purines and pyrimidines. I then pass the reading frames into my alignment algorithm.

My goal was to use the Smith Waterman algorithm with the Gotoh modification to allow for affine gaps. I created an S matrix (for match/mismatch), an X matrix, and a Y matrix (for horizontal and vertical affine gaps). Because this is localized, everything in the fill matrix is minimized at zero. For each letter in the reading frames and the protein sequences, I then use dynamic programming to produce the optimal alignment. I check to see if the current index has increased any of the sequence's lengths. If so, it is heavily penalized so that a frame shift can not be set to that sequence. If the current tracked sequence is diferent then the best scequence available, then the curr_seq is swapped and the frame penalty is assessed. The max score is calculated, and coordinates are updated for tracking purposes. To track affine gaps, an X and Y fill and pointer matrix are used. The pointer matrices try to determine where the highest score came from in order to go back later. The opening gap penalty is for the first time you inter the gap, and the extension is used for subsequent gaps. The traceback algorithm then iterates through the traceback matrices until the Fill matrix reaches a score of zero (for local alignment). The aligned sequences and the DNA bases are printed in blocks as you would see on an alignment program.

# Instructions

Input: DNA Sequence and protein Sequence in FASTA format
* Make sure the bases and amino acids are in capital letters.
The syntax of the affine alignment is smith_waterman_gotoh(seq1, seq2, seq3, prot_seq, dna_seq, scoring_table, go, ge, frame_pen)

1. Place your DNA sequence in dna.fasta
2. Place your Protein sequence in protein.fasta
3. Run `python3 align.py dna.fasta protein.fasta`

# What doesn't work
Unfortunately, there are some traceback issues. The correct score is calculated, but it swaps the gaps with the letters for some reason. For example, my program returned

YLILAL--V

YLIL-GALV

with the correct score, but the "AL" and the "--" should be swapped in the top sequence.

# Examples
You can run this with any DNA and Protein sequence you like, but here are some examples of the different results you could get

1. Perfect Match


Input DNA: ATACTTAATCTAGGCATTAGTGC Input Protein: ILNLGIS

Result:

Max score is 27

ALIGNMENT

Translated:  ILNLGIS

   Protein:  ILNLGIS
   
       DNA:  ATACTTAATCTAGGCATTAGT
       
2. Modifying the Protein to do a frame shift

Input DNA: ATACTTAATCTAGGCATTAGTGC  Input Protein: ILNLALV

Max score is 24

ALIGNMENT

Translated:  ILNLALV

   Protein:  ILNLALV
   
       DNA:  ATACTTAATCTAGCATTAGTG
       
 3. Introducing a Gap (Partially incorrect)

Max score is 21

ALIGNMENT

Translated:  ILN-L-GL

   Protein:  ILNLYRGI
   
       DNA:  ATACTTAATCTAGGCTTAAGT
       

Translated:  S

   Protein:  S
   
       DNA:  
       
 One of the gaps should swap with the L. 
 
 
# Thank you for this fun project, I hope you enjoy my work!




