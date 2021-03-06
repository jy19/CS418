import sys
import collections

with open(sys.argv[1]) as genome_fasta:
        genome_fasta.next()
        dna = genome_fasta.next().strip().lower()
        dna += '$'

with open(sys.argv[2]) as read_fasta:
        pattern_names = []
        patterns = []
        for pattern in read_fasta:
                if pattern[0] == ">":
                        pattern_names.append(pattern.strip())
                else:
                        patterns.append(pattern.strip().lower())

# Make the lists

dna_and_pos = collections.OrderedDict()
for i in range(len(dna)):
        dna_and_pos.setdefault(dna[i:] + dna[:i], i)

first_and_SA = collections.OrderedDict()
sorted = sorted(dna_and_pos)
first = ""
for s in sorted:
        first_and_SA.setdefault(s, dna_and_pos[s])
        first += s[0]

BWT_and_LtoF = collections.OrderedDict()
BWT = ""
for sort in first_and_SA.keys():
        seq = sort[len(sort) - 1] + sort[:len(sort) - 1]
        BWT += seq[0]
        LtoF = first_and_SA.keys().index(seq)
        BWT_and_LtoF[seq] = LtoF

# Lists:

dna = list(dna)
pos = dna_and_pos.values()
first = list(first)
SA = first_and_SA.values()
BWT = list(BWT)
LtoF = BWT_and_LtoF.values()

# Find locations of patterns in genome

starting_positions = []
for pattern in patterns:
        list = SA
        for c in reversed(range(len(pattern))):
                if c == 0:
                        break
                temp_list = []
                for i in list:
                        if first[i] == pattern[c] and BWT[i] == pattern[c - 1]:

                           temp_list.append(LtoF[i])
                list = temp_list
        for i in list:
                starting_positions.append(SA[i])

print starting_positions
