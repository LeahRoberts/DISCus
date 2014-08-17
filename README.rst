BWA_Bedops_aligner
==================

Script for aligning illumina paired end reads with BWA and counting read overlap using: (1) Bedops; and (2) read-pairs traversing the desired region.

How-To: Run Me
---------------

This bash scripts requires the reference sequences to have already been generated. Furthermore, the reads (illumina paired end - not interleaved) need to be in a directory below the reference files, and the bash script should be executed in the file containing the reads.

To run the script, simply type::

 bash <script> <REFERENCE> <BEDMAP_REFERENCE>

The REFERENCE will serve as the reference for the read mapping. 
The BEDMAP_REFERENCE will specify the location of the desired overlap regions.
 
*Note: this script has been modified to accept a particular format of data. See below for detailed formatting specifications.*

For the script to work, the fastq files should be named in a way similar to this::

 $name_1.fastq
 $name_2.fastq


The specifications for the reference files are:

1. A fasta file for mapping the reads to (REFERENCE).
2. A .bed file generated from a .gff file using gff2bed containing the exons you want to count reads overlapping (BEDMAP_REFERENCE).

You can generate this by using Bedops::

 $ gff2bed < genes.gff > genes.bed
 $ bam2bed < reads.bam > reads.bed

  
*It is very important that the nomenclature stays consistent between the reference sequences, particularly regarding the naming of the reference sequence. i.e. The fasta header for the reference should match that in the .bed file.*

For example, if the reference fasta header looks like this::

 >EC958_fimpromoter_off_on_extended
 ATAAACATTAAGTTAACCATATCCATACAAAATACAATGGTTTATGTTCTTCAAAATAAA
 TAAACAAAATCATTCATAAATTTACACATCACTTAAATTCTCCTGTTTCCGCACTTTTTT
 CTTTATTTTTTAAGCAACTGGAAGTTAATCCACTGCAATCTATTGTTATATTGAATCAAA

Then the bed file should look like this::

 EC958_fimpromoter_off_on_extended	1001	1011	.	.	+	artemis	exon	.	gene_id=exon:1002..1011
 EC958_fimpromoter_off_on_extended	1316	1326	.	.	+	artemis	exon	.	gene_id=exon:1317..1326
 EC958_fimpromoter_off_on_extended	3322	3332	.	.	+	artemis	exon	.	gene_id=exon:3323..3332
 EC958_fimpromoter_off_on_extended	3637	3647	.	.	+	artemis	exon	.	gene_id=exon:3638..3647


 
Errors with Script
--------------------

Strains with similar names (eg. HVM5 and HVM52) have a tendency to incorrectly align together. To solve this, I have been separating files with similar names and running the script in separate directories. 
