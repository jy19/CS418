import sys
from suffixtree import SuffixTree

def traverse_tree(node, char_depth, suffix_tree, suffix_array):
    """traverse the tree recursively to generate suffix array"""
    # print node.start
    # print node.start, suffix_tree.text[node.start:node.end], char_depth
    if not node.edges:
        # print 'leaf node', node.start, char_depth
        suffix_array.append(node.start - char_depth - 1)
        return
    chars = ['$', 'A', 'C', 'G', 'T']
    for char in chars:
        try:
            next_node = suffix_tree.nodes[node.edges[char]]
            # new_depth = next_node.start - char_depth if char_depth > 0 else (node.end - node.start)
            new_depth = char_depth + (node.end - node.start)
            traverse_tree(next_node, new_depth, suffix_tree, suffix_array)
        except KeyError:
            pass

def create_suffix_array(suffix_tree):
    suffix_array = []
    root = suffix_tree.nodes[suffix_tree.root]
    traverse_tree(root, root.start, suffix_tree, suffix_array)
    return suffix_array

def create_matrices(suffix_array):
    bwt, first_col, ltof = [], [], []
    return bwt, first_col, ltof

def map():
    pass

if __name__ == '__main__':
    # with open(sys.argv[1]) as genome_fasta:
    #     genome_fasta.next()
    #     dna = genome_fasta.next().strip().lower()
    #     dna += '$'
    #
    # with open(sys.argv[2]) as read_fasta:
    #         pattern_names = []
    #         patterns = []
    #         for pattern in read_fasta:
    #                 if pattern[0] == ">":
    #                         pattern_names.append(pattern.strip())
    #                 else:
    #                         patterns.append(pattern.strip().lower())

    # dna = 'ACGGACT$'
    # dna = 'TCAGGTAGCTTACT$'
    dna = 'AACGCTTAGTC$'

    suffix_tree = SuffixTree(len(dna))
    for char in dna:
        suffix_tree.add_char(char)
    # suffix_tree.print_tree()
    suffix_array = create_suffix_array(suffix_tree)
    print suffix_array
