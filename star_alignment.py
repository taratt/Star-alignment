from copy import deepcopy
def global_align(x, y, s_match, s_mismatch, s_gap):

    A = []

    for i in range(len(y) + 1):

        A.append([0] * (len(x) + 1))

    for i in range(len(y) + 1):

        A[i][0] = s_gap * i

    for i in range(len(x) + 1):

        A[0][i] = s_gap * i

    for i in range(1, len(y) + 1):

        for j in range(1, len(x) + 1):

            A[i][j] = max(

                A[i][j - 1] + s_gap,

                A[i - 1][j] + s_gap,

                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (

                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (

                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)

            )

    align_X = ""

    align_Y = ""

    i = len(x)

    j = len(y)

    while i > 0 or j > 0:

        current_score = A[j][i]

        if i > 0 and j > 0 and (

                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or

                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][

                    i - 1] + s_mismatch) or

                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)

        ):

            align_X = x[i - 1] + align_X

            align_Y = y[j - 1] + align_Y

            i = i - 1

            j = j - 1

        elif i > 0 and (current_score == A[j][i - 1] + s_gap):

            align_X = x[i - 1] + align_X

            align_Y = "-" + align_Y

            i = i - 1

        else:

            align_X = "-" + align_X

            align_Y = y[j - 1] + align_Y

            j = j - 1

    return (align_X, align_Y, A[len(y)][len(x)] )

def get_input():
    sequences = []
    n = int(input())
    for i in range(n):
        sequences.append(input())

    return sequences

def get_center(sequences):
    score_matrix = {}
    for i in range(len(sequences)):
        score_matrix[i] = {}
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            if i != j:
                score_matrix[i][j] = global_align(sequences[i], sequences[j],3,-1,-2)

    center = 0
    center_score = float('-inf')
    for scores in score_matrix:
        sum = 0
        for i in score_matrix[scores]:
            sum += score_matrix[scores][i][2]
        if sum > center_score:
            center_score = sum
            center = scores
    alignments_needed = {}
    for i in sequences:
        if i != sequences[center]:
            s1,s2,sc = global_align(i, sequences[center],3,-1,-2)
            alignments_needed[sequences.index(i)] = (s2,s1,sc)
    #print(score_matrix)
    return center, alignments_needed

def align_gaps(seq1, seq2, aligneds, new):
    i = 0
    while i< max(len(seq1), len(seq2)):
        try:
            if i > len(seq1)-1:
                seq1 = seq1[:i] + "-" + seq1[i:]
                naligneds = []
                for seq in aligneds:
                    naligneds.append(seq[:i] + "-" + seq[i:])
                aligneds = naligneds
            elif i> len(seq2)-1:
                seq2 = seq2[:i] + "-" + seq2[i:]
                new = new[:i] + "-" + new[i:]
            elif (seq1[i] == "-" and i>= len(seq2)) or (seq1[i] == "-" and seq2[i] != "-") :
                seq2 = seq2[:i]+"-"+seq2[i:]
                new = new[:i]+"-"+new[i:]
            elif (seq2[i] == "-" and i>= len(seq1)) or (seq2[i] == "-" and seq1[i] != "-"):
                #print("me")
                seq1 = seq1[:i] + "-" + seq1[i:]
                naligneds = []
                for seq in aligneds:
                    naligneds.append(seq[:i] + "-" + seq[i:])
                aligneds = naligneds

        except:
            print("u oh")
        i+=1

    aligneds.append(new)
    #print(aligneds)
    return seq1, aligneds

def msa(alignments):
    aligned_center = alignments[list(alignments.keys())[0]][0]
    aligneds = []
    aligneds.append(alignments[list(alignments.keys())[0]][1])

    for seq in list(alignments.keys())[1:]:
        cent = alignments[seq][0]
        newseq = alignments[seq][1]
        aligned_center, aligneds = align_gaps(aligned_center, cent, aligneds, newseq)

    return aligneds, aligned_center

def calculate_scores(alignments):
    score = 0
    for i in range(len(alignments[0])):
        for j in range(len(alignments)-1):
            for k in range(j+1,len(alignments)):
                if alignments[j][i] == alignments[k][i]:
                    if alignments[j][i]!="-":
                        score +=3

                elif alignments[k][i]=="-" or alignments[j][i]=="-":
                    score-=2
                else:
                    score-=1
    return score

def order_results(aligneds, center_seq, center):
    i = 0
    j = 0
    results = []
    while i < len(sequences):
        if i == center:
            results.append(center_seq)
            i += 1
        else:
            results.append(aligneds[j])
            i += 1
            j += 1
    return results

def block_optimization(results, score):
    continue_bit = True

    while continue_bit:
        #print(results)
        block_columns = []
        blocks = []
        for i in range(len(results[0])):
            for seq in range(len(results) - 1):
                if (results[seq][i] != results[seq+1][i]) or results[seq][i]=="-":
                    block_columns.append(i)
                    break
        i = 0
        while i < len(block_columns):
            col = block_columns[i]
            counter = col
            while counter < len(results[0]):
                if counter+1 in block_columns:
                    counter+=1
                else:
                    break
            if counter != col:
                blocks.append((col,counter))
                i = counter
            else:
                i+=1
        #print(blocks)
        seq_blocks = []
        for b in blocks:
            block = []
            for seq in results:
                nseq = seq[b[0]:b[1]+1]
                block.append(nseq.replace("-",""))
            if "" not in block:
                seq_blocks.append(block)
        block = seq_blocks[0]
        counter = 0
        for block in seq_blocks:

            center, alignments = get_center(block)
            aligneds, center_seq = msa(alignments)
            res = order_results(aligneds, center_seq, center)

            index = seq_blocks.index(block)
            #print(res)
            news = []
            for sequence in results:
                news.append(sequence[:blocks[index][0]]+res[results.index(sequence)]+sequence[blocks[index][1]+1:])

            if calculate_scores(news) > score:
                results = news
                break
            counter+=1

        if counter == len(seq_blocks):
                continue_bit = False

    return results, calculate_scores(results)


if __name__ == '__main__':
    sequences = get_input()
    center, alignments = get_center(sequences)
    aligneds, center_seq = msa(alignments)
    #print(alignments)
    results = order_results(aligneds, center_seq, center)
    results, score = block_optimization(results,calculate_scores(results))
    print(calculate_scores(results))
    #print(score)
    for i in results:
        print(i)



