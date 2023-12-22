# Microbiome Analysis
Producing an example of microbiome analysis using the Dada2 pipeline. The codes from the Dada2 pipeline workflow tutorial (1.16) from https://benjjneb.github.io/dada2/tutorial.html have been adapted for use in this project. The taxonomy is assigned with reference to SILVA database v138.

Modifications from the Dada2 pipeline workflow tutorial (1.16):
1) The code from the tutorial has been modified for processing of big data, through incorporation of the optimised code from https://benjjneb.github.io/dada2/bigdata.html while keeping the rest of the workflow similar to that in the tutorial.
2) An additional line of code has been added to check that the number of forward and reverse fastq files matches before filtering the fastq files.
3) In the track reads table, an additional column has been added to show the percentage of reads retained for each sample
4) ASV IDs are assigned. The resultant sequence tables and taxonomy tables uses the shortened ASV IDs instead of displaying the full DNA sequence.
5) The ASV DNA sequences are stored in a FASTA file.
6) Code for creating phylogenic tree has been added in.
7) The microbiome dataset used is the same as the data for the mothur MiSeq SOP https://mothur.org/wiki/miseq_sop/, but the metadata used in subsequent data analysis has been modified.

Downstream Analysis:
