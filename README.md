# Demultiplexing

Illumina output FASTQ files contain multiplexed sequence data from different biological samples. This repo contains the code for an algorithm that will sort each record within a set of FASTQ files to be reassigned to files containing only the FASTQ records for that biological sample group. To do this, paired-end sequence reads must be sorted using their unique dual-matched index barcode, however, the code needs to account for barcodes having undetermined base calls (N) as well as index hopping, which results in some reads containing two different indexes on either side of the sequence.

Final script titled **Dplexer.py** can be found under assignment, the third.


