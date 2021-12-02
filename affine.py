from get_score import *

#DP algorithm to fill matrices
def smith_waterman_gotoh(seq1, seq2, seq3, prot_seq, dna_seq, scoring_table, go, ge, frame_pen):
    #Main fill matrix for sequence 1
    #how do i do frame shift?
    seq1_fill = [[0 for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    #direction matrix for fill matrix
    seq1_ptr = [["" for c in range(len(prot_seq)+1)] 
            for r in range(len(seq1)+1)]
    
    #affine gap horizontally
    insert_x = [[float("-inf") for c in range(len(prot_seq)+1)] 
            for r in range(max(len(seq1),len(seq2),len(seq3))+1)]
    x_ptr = [["" for c in range(len(prot_seq)+1)] 
            for r in range(max(len(seq1),len(seq2),len(seq3))+1)]
    #affine gap vertically
    insert_y = [[float("-inf") for c in range(len(prot_seq)+1)] 
            for r in range(max(len(seq1),len(seq2),len(seq3))+1)]
    y_ptr = [["" for c in range(len(prot_seq)+1)] 
            for r in range(max(len(seq1),len(seq2),len(seq3))+1)]
    x_align = ""
    y_align = ""
    prot_align = ""
    max_score = 0
    max_i = 0
    max_j = 0
    curr_seq = 0
    seq1 = list(seq1)
    seq2 = list(seq2)
    seq3 = list(seq3)
    for i in range(1, max(len(seq1),len(seq2),len(seq3))+1):
        for j in range(1, len(prot_seq)+1):
            curr_seq_start = curr_seq
            seq1_fill_options = [
                seq1_fill[i - 1][j - 1],  # 0
                insert_x[i - 1][j-1],  # 1
                insert_y[i-1][j - 1],  # 2
                0]
            seq1_score = 0
            seq2_score = 0
            seq3_score = 0
            if ((i-1) >= len(seq1)):
                seq1_score = float("-inf")
            else:
                seq1_score = get_score(seq1[i-1], prot_seq[j-1], scoring_table)
            if ((i-1) >= len(seq2)):
                seq2_score = float("-inf")
            else:
                seq2_score = get_score(seq2[i-1], prot_seq[j-1], scoring_table)
            if ((i-1) >= len(seq3)):
                seq3_score = float("-inf")
            else:
                seq3_score = get_score(seq3[i-1], prot_seq[j-1], scoring_table)

            frame_shift_options = [
                    seq1_score,
                    seq2_score,
                    seq3_score]
            fs_penalty = 0
            # figure out when to add frame_shift_options
            if(frame_shift_options[curr_seq] < (frame_shift_options[frame_shift_options.index(max(frame_shift_options))] - frame_pen)):
                #Frame shift is penalized
                fs_penalty = -1*frame_pen
                curr_seq = frame_shift_options.index(max(frame_shift_options))
            seq1_fill[i][j] = max(frame_shift_options) + fs_penalty + max(seq1_fill_options)
            #seq1_fill[i][j] += max(seq1_fill_options)
            seq1_fill[i][j] = max(seq1_fill[i][j], 0)
            if (seq1_fill[i][j] > max_score):
                max_score = seq1_fill[i][j]
                max_i = i
                max_j = j
            ptr = seq1_fill_options.index(max(seq1_fill_options))
            if (ptr == 0): # Comes from a match/mismatch
                seq1_ptr[i][j]=curr_seq
            elif (ptr == 1): 
                seq1_ptr[i][j] = "X"
            elif (ptr == 2): 
                seq1_ptr[i][j] = "Y"
            #Gap open, gap extension
            insert_x_options = [seq1_fill[i-1][j] - (go), insert_x[i-1][j] - ge]
            insert_x[i][j] = max(insert_x_options)
            ptr = insert_x_options.index(max(insert_x_options))
            if (ptr == 0):
                x_ptr[i][j] = curr_seq
            else:
                x_ptr[i][j] = "X"
            insert_y_options = [seq1_fill[i][j-1] - (go), insert_y[i][j-1] - ge]
            insert_y[i][j] = max(insert_y_options)
            ptr = insert_y_options.index(max(insert_y_options))
            if (ptr == 0):
                y_ptr[i][j] = curr_seq
            else:
                y_ptr[i][j] = "Y"
    curr_i = max_i
    curr_j = max_j
    #y_align = seq1[curr_i] #maybe switch
    #x_align = prot_seq[curr_j]
    goto = seq1_ptr[curr_i][curr_j]
    count=0
    dna_seq_final = ""
    while(seq1_fill[curr_i][curr_j]!=0):
        if goto == 0 or goto == 1 or goto == 2:
            curr_i -=1
            curr_j -=1
            x_align = prot_seq[curr_j] + x_align
            if goto == 0:
                curr_seq = 0
                y_align = seq1[curr_i] + y_align
                dna_seq_final = dna_seq[curr_i*3:3+curr_i*3] + dna_seq_final
            elif goto == 1:
                curr_seq = 1
                y_align = seq2[curr_i] + y_align
                dna_seq_final = dna_seq[1+curr_i*3:4+curr_i*3] + dna_seq_final
            else:
                curr_seq = 2
                y_align = seq3[curr_i] + y_align
                dna_seq_final = dna_seq[2+curr_i*3:5+curr_i*3] + dna_seq_final
            goto = seq1_ptr[curr_i][curr_j]
        elif goto == "X":
            #maybe switch
            curr_i -=1
            if curr_seq == 0:
                y_align = seq1[curr_i] + y_align
                dna_seq_final = dna_seq[curr_i*3:3+curr_i*3] + dna_seq_final
            elif curr_seq == 1:
                y_align = seq2[curr_i] + y_align
                dna_seq_final = dna_seq[1+curr_i*3:4+curr_i*3] + dna_seq_final
            else:
                y_align = seq3[curr_i] + y_align
                dna_seq_final = dna_seq[2+curr_i*3:5+curr_i*3] + dna_seq_final

            x_align = "-" + x_align
            goto = x_ptr[curr_i][curr_j]
        elif goto == "Y":
            curr_j -= 1
            y_align = "-" + y_align
            x_align = prot_seq[curr_j] + x_align
            goto = y_ptr[curr_i][curr_j]
        count += 1
            
    print("Max score is", max_score)

    print("ALIGNMENT")
    index = 0
    while (index < len(y_align)-8):
        print("Translated: ", y_align[index:index+8])
        print("   Protein: ", x_align[index:index+8])
        print("       DNA: ", dna_seq_final[3*index:3*(index+8)])
        index += 8
        print()
    print("Translated: ", y_align[index:])
    print("   Protein: ", x_align[index:])
    print("       DNA: ", dna_seq_final[3*index:])
    return seq1_fill, insert_x, insert_y


