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