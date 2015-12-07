from __future__ import print_function
import sys
from suffixtree import SuffixTree
from utils import generate_kmers, create_subscripts
import threading

# dirty hack for thread safe stdout, so can print as it goes in case program times out
print = lambda x: sys.stdout.write('%s\n' % x)

class MapThread(threading.Thread):
    def __init__(self, patterns, mapper, dna_name, kmer_size=None):
        self.dna_name = dna_name
        self.mapper = mapper
        self.patterns = patterns
        self.pattern_mappings = []
        self.kmer_size = kmer_size
        threading.Thread.__init__(self)

    def run(self):
        if self.kmer_size:
            for pattern in self.patterns:
                estimated_position = estimate_position(pattern[0], kmer_size)
                if estimated_position:
                    # entry = output_to_SAM(estimated_position, pattern[0], pattern[1], self.dna_name)
                    output_to_SAM(estimated_position, pattern[0], pattern[1], self.dna_name)
                    # self.pattern_mappings.append(entry)
        else:
            for pattern in self.patterns:
                positions = mapper.map(pattern[0])
                positions = [x + 1 for x in positions]
                if not positions:
                    continue
                # entry = output_to_SAM(positions[0], pattern[0], pattern[1], self.dna_name)
                output_to_SAM(positions[0], pattern[0], pattern[1], self.dna_name)
                # self.pattern_mappings.append(entry)

class Mapper:
    def __init__(self, dna):
        self.dna = dna
        self.suffix_tree = SuffixTree(len(dna))
        self.suffix_array = []
        self.first_col = []
        self.bwt = []
        self.ltof = []
        self.init_self()

    def init_self(self):
        """initializes lists needed for mapping"""
        for c in self.dna:
            self.suffix_tree.add_char(c)
        root = self.suffix_tree.nodes[self.suffix_tree.root]
        self.traverse_tree(root, root.start)
        self.first_col = create_subscripts([self.dna[x] for x in self.suffix_array])
        self.bwt = create_subscripts([self.dna[x - 1] if x > 0 else self.dna[-1] for x in self.suffix_array])
        self.ltof = [self.first_col.index(x) for x in self.bwt]

    def traverse_tree(self, node, char_depth):
        """traverse the tree recursively to generate suffix array"""
        if not node.edges:
            self.suffix_array.append(node.start - char_depth - 1)
            return
        chars = ['$', 'A', 'C', 'G', 'T']
        for char in chars:
            try:
                next_node = self.suffix_tree.nodes[node.edges[char]]
                new_depth = char_depth + (node.end - node.start)
                self.traverse_tree(next_node, new_depth)
            except KeyError:
                pass

    def get_position_range(self, char, list_type, start, end):
        """returns range of positions for a char in specified type of list"""
        positions = []
        for index in range(start, end + 1):
            if list_type[index].startswith(char):
                positions.append(index)
                break
        for index in range(end, start, -1):
            if list_type[index].startswith(char):
                positions.append(index)
                break
        return positions

    def map(self, pattern):
        """find a given pattern in the genome"""
        found_positions = []
        rev_pattern = pattern[::-1]
        current_positions = self.get_position_range(rev_pattern[0], self.first_col, 0, len(self.first_col) - 1)
        for i in range(1, len(rev_pattern)):
            try:
                bwt_positions = self.get_position_range(rev_pattern[i], self.bwt, current_positions[0],
                                                        current_positions[-1])
                ltof_positions = range(self.ltof[bwt_positions[0]], self.ltof[bwt_positions[-1]] + 1)
                current_positions = ltof_positions
            except IndexError:
                # may not be able to find, simply return so don't have to go through entire pattern
                return found_positions
        # after last char, push every position in SA to found_positions
        for i in current_positions:
            found_positions.append(self.suffix_array[i])
        return found_positions

def estimate_position(pattern, kmer_size):
    kmers = generate_kmers(pattern, kmer_size)
    possible_positions = {}
    kmer_offset = 0
    for kmer in kmers:
        kmer_positions = mapper.map(kmer)
        for pos in kmer_positions:
            try:
                possible_positions[pos - kmer_offset] += 1
            except KeyError:
                possible_positions[pos - kmer_offset] = 1
        kmer_offset += 1
    try:
        estimated_position = max(possible_positions, key=possible_positions.get) + 1
    except ValueError:
        estimated_position = None
    return estimated_position

def output_to_SAM(position, pattern, pattern_name, genome_name):
    entry = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}'.format(
        pattern_name, 0, genome_name, position, 255, str(len(pattern)) + 'M', '*', 0, 0, pattern.lower(), '*')
    print(entry)
    # return entry

# def print_all(pattern_mappings):
#     for mapping in pattern_mappings:
#         print mapping

if __name__ == '__main__':
    with open(sys.argv[1]) as genome_fasta:
        genome_name = genome_fasta.next().strip()[1:]
        dna = ""
        for line in genome_fasta:
            dna += line.strip().upper()
        dna += '$'

    with open(sys.argv[2]) as read_fasta:
        pattern_names = []
        patterns = []
        pattern_tuples = []
        for pattern in read_fasta:
            if pattern[0] == ">":
                pattern_names.append(pattern.strip()[1:])
            else:
                patterns.append(pattern.strip().upper())

        for i in range(len(patterns)):
            pattern_tuples.append((patterns[i], pattern_names[i]))

    thread_count = int(sys.argv[3])
    mapper = Mapper(dna)

    try:
        kmer_size = int(sys.argv[4])
    except IndexError:
        kmer_size = None

    threads = []
    increments = int(len(patterns)) / thread_count

    for i in range(thread_count):
        pattern_range = pattern_tuples[increments*i:increments*(i+1)]
        t = MapThread(pattern_range, mapper, genome_name, kmer_size)
        threads.append(t)

    print('Finished reading files and creating threads.')

    for thread in threads:
        thread.start()

    for thread in threads:
        thread.join()

    # for thread in threads:
    #     print_all(thread.pattern_mappings)
