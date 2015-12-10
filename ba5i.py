match_score = 1
indel_penalty = -2
DIAG = 0
HORIZ = 1
VERT = 2

def calc_score_matrix(dna_s, dna_t):
    # dna_s = '_' + dna_s
    # dna_t = '_' + dna_t
    score_matrix = [[(0, -1) for x in range(len(dna_t))] for y in range(len(dna_s))]

    for i in range(1, len(dna_s)):
        for j in range(1, len(dna_t)):
            # get the costs
            if dna_s[i] == dna_t[j]:
                # match
                match = (score_matrix[i-1][j-1][0] + match_score, DIAG)
            else:
                # mismatch
                match = (score_matrix[i-1][j-1][0] + indel_penalty, DIAG)
            indel_h = (score_matrix[i][j-1][0] + indel_penalty, HORIZ)
            indel_v = (score_matrix[i-1][j][0] + indel_penalty, VERT)
            best_choice = max(match, indel_h, indel_v)
            score_matrix[i][j] = best_choice

    return score_matrix

def find_max_score(score_matrix):
    # look for max score in bottom row
    last_row = score_matrix[-1]
    max_score = max(last_row)
    i = len(score_matrix) - 1
    j = last_row.index(max_score)
    return max_score, i, j

def backtrack(dna_s, dna_t, score_matrix, i, j):
    # print 'start i, j: ', i, j
    aligned1 = ''
    aligned2 = ''
    while j > 0:
        if score_matrix[i][j][1] == DIAG:
            aligned1 += dna_s[i]
            aligned2 += dna_t[j]
            i -= 1
            j -= 1
        elif score_matrix[i][j][1] == HORIZ:
            aligned1 += '-'
            aligned2 += dna_t[j]
            j -= 1
        elif score_matrix[i][j][1] == VERT:
            aligned1 += dna_s[i]
            aligned2 += '-'
            i -= 1
        # print 'i, j: ', i, j
    return aligned1[::-1], aligned2[::-1]

if __name__ == '__main__':
    dna_s = raw_input().strip()
    dna_t = raw_input().strip()
    dna_s = '_' + dna_s
    dna_t = '_' + dna_t
    score_matrix = calc_score_matrix(dna_s, dna_t)
    max_score, i, j = find_max_score(score_matrix)
    print max_score[0]
    aligned1, aligned2 = backtrack(dna_s, dna_t, score_matrix, i, j)
    print aligned1
    print aligned2
