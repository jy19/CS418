def generate_kmers(read, k):
    kmers = []
    for i in range(len(read)):
        kmer = read[i:i+k]
        if len(kmer) < k:
            continue
        kmers.append(kmer)
    return kmers

def handle_errors(kmers):
    pass
