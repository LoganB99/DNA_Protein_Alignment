def convert_to_amino(seq, codons):
    #read in dna sequences
    #set list comprehensions to variables
    #convert to codons (make the codon.txt into a structure
    #return the sequences
    # code the local alignment
    #figure out traceback
    #what to do with stop codons? use . in blosum

    seq1 = ""
    seq2 = ""
    seq3 = ""
    rf1 = ([seq[i:i+3] for i in range(0,len(seq)-2, 3)])
    for codon in rf1:
        seq1 += codons[codon]
    rf2 = ([seq[i:i+3] for i in range(1,len(seq)-2, 3)])
    for codon in rf2:
        seq2 += codons[codon]
    rf3 = ([seq[i:i+3] for i in range(2,len(seq)-2, 3)])
    for codon in rf3:
        seq3 += codons[codon]

    return [seq1, seq2, seq3]
