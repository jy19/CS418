import sys
from suffixtree import SuffixTree

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
        return positions

    def map(self, pattern):
        found_positions = []
        rev_pattern = pattern[::-1]
        current_positions = self.get_position_range(rev_pattern[0], self.first_col, 0, len(self.first_col) - 1)
        for i in range(1, len(rev_pattern)):
            bwt_positions = self.get_position_range(rev_pattern[i], self.bwt, current_positions[0], current_positions[-1])
            ltof_positions = range(self.ltof[bwt_positions[0]], self.ltof[bwt_positions[-1]] + 1)
            current_positions = ltof_positions
        # after last char, push every position in SA to found_positions
        for i in current_positions:
            found_positions.append(self.suffix_array[i])
        return found_positions

def create_subscripts(text):
    """helper to append indices (subscripts) to a list of characters (text)"""
    subscripts = {}
    subscripted_text = []
    for i in range(len(text)):
        try:
            subscripts[text[i]]
        except KeyError:
            subscripts[text[i]] = 0
        subscripts[text[i]] += 1
        subscripted_text.append('{0}-{1}'.format(text[i], subscripts[text[i]]))
    return subscripted_text

if __name__ == '__main__':
    with open(sys.argv[1]) as genome_fasta:
        genome_fasta.next()
        dna = ""
        for line in genome_fasta:
            dna = line.strip().upper()
        dna += '$'

    with open(sys.argv[2]) as read_fasta:
        pattern_names = []
        patterns = []
        for pattern in read_fasta:
            if pattern[0] == ">":
                pattern_names.append(pattern.strip())
            else:
                patterns.append(pattern.strip().lower())

    # dna = 'CGTGATGCGCGGAC$'
    # patterns = ['GCG']

    mapper = Mapper(dna)
    positions = []
    for pattern in patterns:
        positions.extend(mapper.map(pattern))

    print positions
