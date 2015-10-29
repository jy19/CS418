import sys
from utils import generate_kmers

def init_graph(kmers, edge_count, graph, counts):
    for kmer in kmers:
        # print kmer
        left = kmer[: -1]
        right = kmer[1:]
        edge_count += 1

        if left in graph:
            graph[left].append(right)
        else:
            graph[left] = [right]

        if left in counts:
            counts[left][1] += 1
        else:
            counts[left] = [0, 1]

        if right in counts:
            counts[right][0] += 1
        else:
            counts[right] = [1, 0]
    return graph, counts, edge_count

def build_graph(start, g, non_branching):
    global edge_count
    path = [start]
    cur_node = start

    while len(cur_node) > 0:
        next_node = g[cur_node][0]
        del g[cur_node][0]
        # if len( g[ cur_node ] ) == 0:
        #    del g[ cur_node ]

        edge_count -= 1

        # print cur_node
        if next_node in non_branching:
            # print "non branching:", cur_node
            path.append(next_node)
            cur_node = next_node
            continue
        else:
            path.append(next_node)
            break

    return path

def merge_nodes(nodes):
    contig = nodes[0]
    for i in range(1, len(nodes)):
        contig += nodes[i][-1]
    return contig


def has_outgoing(node):
    if len(graph[node]) > 0:
        return True
    else:
        return False

def generate_contigs(counts, graph):
    non_branching = set()
    contigs = []

    for key, item in counts.iteritems():
        if item[0] == 1 and item[1] == 1:
            non_branching.add(key)

    # print non_branching

    start = graph.keys()[0]
    while edge_count > 0:
        for i in graph.keys():
            if i in non_branching or len(graph[i]) == 0:
                continue
            start = i
            break
        # print "starting with:", start
        c = build_graph(start, graph, non_branching)
        contigs.append(c)
    return contigs

if __name__ == '__main__':
    k = int(sys.argv[2])
    reads = []
    # for line in sys.stdin:
    #     kmers.append(line.strip())
    with open(sys.argv[1]) as fh:
        for line in fh:
            if line[0] != '>':
                reads.append(line.strip())

    kmers = []
    for read in reads:
        kmers.extend(generate_kmers(read, k))
    kmers = list(set(kmers))
    graph = dict()
    counts = dict()
    edge_count = 0
    graph, counts, edge_count = init_graph(kmers, edge_count, graph, counts)

    contigs = generate_contigs(counts, graph)
    contigs = [merge_nodes(contig) for contig in contigs]

    # graph = dict()
    # counts = dict()
    # edge_count = 0
    # graph, counts, edge_count = init_graph(contigs, edge_count, graph, counts)
    # contigs2 = generate_contigs(counts, graph)
    # contigs2 = [merge_nodes(contig) for contig in contigs2]
    # print '\n'.join(contig for contig in sorted(contigs2))
    print '\n'.join(contig for contig in sorted(contigs))
