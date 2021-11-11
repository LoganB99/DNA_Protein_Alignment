def get_score(amino1, amino2, table):
    index1 = table[0].index(amino1)
    index2 = table[0].index(amino2)
    return(int(table[index1][index2]))
