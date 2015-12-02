import sys
from suffixtree import SuffixTree

def traverse_tree(node, char_depth, suffix_tree, suffix_array):
    """traverse the tree recursively to generate suffix array"""
    if not node.edges:
        suffix_array.append(node.start - char_depth - 1)
        return
    chars = ['$', 'A', 'C', 'G', 'T']
    for char in chars:
        try:
            next_node = suffix_tree.nodes[node.edges[char]]
            new_depth = char_depth + (node.end - node.start)
            traverse_tree(next_node, new_depth, suffix_tree, suffix_array)
        except KeyError:
            pass

def create_suffix_array(suffix_tree):
    suffix_array = []
    root = suffix_tree.nodes[suffix_tree.root]
    traverse_tree(root, root.start, suffix_tree, suffix_array)
    return suffix_array

def create_matrices(dna, suffix_array):
    """create the other lists from SA needed to do mapping"""
    first_col = create_subscripts([dna[x] for x in suffix_array])
    bwt = create_subscripts([dna[x - 1] if x > 0 else dna[-1] for x in suffix_array])
    ltof = [first_col.index(x) for x in bwt]
    return bwt, first_col, ltof

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

def map():
    pass

if __name__ == '__main__':
    # with open(sys.argv[1]) as genome_fasta:
    #     genome_fasta.next()
    #     dna = genome_fasta.next().strip().lower()
    #     dna += '$'
    #
    # with open(sys.argv[2]) as read_fasta:
    #     pattern_names = []
    #     patterns = []
    #     for pattern in read_fasta:
    #         if pattern[0] == ">":
    #             pattern_names.append(pattern.strip())
    #         else:
    #             patterns.append(pattern.strip().lower())

    dna = 'CGTGATGCGCGGAC$'

    suffix_tree = SuffixTree(len(dna))
    for char in dna:
        suffix_tree.add_char(char)
    # suffix_tree.print_tree()
    suffix_array = create_suffix_array(suffix_tree)
    print create_matrices(dna, suffix_array)
