from get_score import *

#DP algorithm to fill matrices
def smith_waterman_gotoh(seq1, seq2, seq3, prot_seq, scoring_table, go, ge, frame_pen):
    #Main fill matrix for sequence 1
    #how do i do frame shift?
    seq1_fill = [[0 for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    #direction matrix for fill matrix
    seq1_ptr = [["" for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    #affine gap horizontally
    insert_x = [[float("-inf") for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    x_ptr = [["" for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    #affine gap vertically
    insert_y = [[float("-inf") for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    y_ptr = [["" for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    x_align = ""
    y_align = ""
    prot_align = ""
    max_score = 0
    max_i = 0
    max_j = 0
    for i in range(1, len(seq1)+1):
        for j in range(1, len(prot_seq)+1):
            seq1_fill_options = [ 
                    seq1_fill[i-1][j-1]+ get_score(seq1[i-1],prot_seq[j-1],scoring_table), #0
                    insert_x[i-1][j-1] + get_score(seq1[i-1],prot_seq[j-1],scoring_table), #1
                    insert_y[i-1][j-1] + get_score(seq1[i-1],prot_seq[j-1],scoring_table), #2
                    0]

            seq1_fill[i][j] = max(seq1_fill_options)
            if (seq1_fill[i][j] > max_score):
                max_score = seq1_fill[i][j]
                max_i = i
                max_j = j
            ptr = seq1_fill_options.index(max(seq1_fill_options))
            if (ptr == 0 or ptr == 3): # Comes from a match/mismatch
                seq1_ptr[i][j]="S"
            elif (ptr == 1): 
                seq1_ptr[i][j] = "X"
            elif (ptr == 2): 
                seq1_ptr[i][j] = "Y"
            #Gap open, gap extension
            insert_x_options = [seq1_fill[i-1][j] - go, insert_x[i-1][j] - ge]
            insert_x[i][j] = max(insert_x_options)
            ptr = insert_x_options.index(max(insert_x_options))
            if (ptr == 0):
                x_ptr[i][j] = "S"
            else:
                x_ptr[i][j] = "X"
            insert_y_options = [seq1_fill[i][j-1] - go, insert_y[i][j-1] - ge]
            insert_y[i][j] = max(insert_y_options)
            ptr = insert_y_options.index(max(insert_y_options))
            if (ptr == 0):
                y_ptr[i][j] = "S"
            else:
                y_ptr[i][j] = "Y"
    curr_i = max_i
    curr_j = max_j
    print(seq1_fill)
    print(max_score)
    goto = seq1_ptr[curr_i][curr_j]
    while(seq1_fill[curr_i][curr_j]!=0):
        if goto == "S":
            y_align = prot_seq[curr_j-1] + y_align
            x_align = seq1[curr_i-1] + x_align
            goto = seq1_ptr[curr_i-1][curr_j-1]
            curr_i -=1
            curr_j -=1
        elif goto == "X":
            x_align = seq1[curr_i-1] + x_align
            y_align = "-" + y_align
            goto = x_ptr[curr_i-1][curr_j]
            curr_i -=1
        elif goto == "Y":
            x_align = "-" + x_align
            y_align = prot_seq[curr_j-1] + y_align
            goto = y_ptr[curr_i][curr_j-1]
            curr_j -= 1

    print(y_align)
    print(x_align)
    return seq1_fill, insert_x, insert_y


