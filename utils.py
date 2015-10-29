from collections import Counter
import sys

def generate_kmers(read, k):
    kmers = []
    for i in range(len(read)):
        kmer = read[i:i+k]
        if len(kmer) < k:
            continue
        kmers.append(kmer)
    return kmers

def handle_errors(kmers, threshold):
    kmer_counts = dict(Counter(kmers))
    filtered_kmers = []
    for k, v in kmer_counts.iteritems():
        if v <= threshold:
            continue
        filtered_kmers.append(k)
    return filtered_kmers

def get_contig_info(contigs):
    longest, shortest = None, None
    avg_len, longest_len, shortest_len, total_len = 0, 0, sys.maxint, 0
    for contig in contigs:
        total_len += len(contig)
        if len(contig) > longest_len:
            longest_len = len(contig)
            longest = contig
        if len(contig) < shortest_len:
            shortest_len = len(contig)
            shortest = contig
    avg_len = total_len / len(contigs)
    contig_info = {'longest': longest, 'shortest': shortest, 'avg': avg_len}
    return contig_info
