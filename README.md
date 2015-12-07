To run the code for datasets without errors:

    python assembler.py input_file kmer_size

`kmer_size` is a range of values, so it can run multiple kmer size at once.
It should be formatted like `start-end`.

To run the code for datasets with errors:

    python assembler.py input_file kmer_size error_threshold

`error_threshold` is basically the number of time a kmer appears.
When generating kmers from reads, an error in a read will result in a kmer that does not appear as often as other kmers.
Try different `error_threshold` values to see if longest contig length changes.

To output to a file for easier reading, pipe stdout to file:

    python assembler.py input_file kmer_size [error_threshold] > output_file

To run the mapping project:

    python mapper.py genome_file reads_file thread_count

`thread_count` is the number of cores you'd like to run this program on (for speed-up purposes).
Note: if any of the arguments are missing the program will fail. 
To output to a SAM file to use with samtools:

    python mapper.py genome_file reads_file thread_count > file_name.sam

To run the mapping project with errors:

    python mapper.py genome_file reads_file thread_count kmer_size

`kmer_size` is the size of the kmer a read will be broken down into for error handling.