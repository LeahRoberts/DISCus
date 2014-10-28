DISCO - *D*NA *S*w*I*tch *CO*unter
========================================

Script for mapping illumina paired end reads with BWA and counting read overlap using: (1) Bedops; and (2) read-pairs traversing the desired region.

Installation Requirements
--------------------------

1. BWA (http://sourceforge.net/projects/bio-bwa/files/)
2. SAMtools (http://samtools.sourceforge.net/)
3. Bedops - bedmaps (http://bedops.readthedocs.org/en/latest/) 


File Requirements
------------------

1. REFERENCE in fasta format (see below - 'Construction of Reference')
2. BEDMAP_REFERENCE in bed format. (see below - 'Construction of Bedmap_reference')
3. illumina paired-end read files (fastq) NOT interleaved (see below - Fastq File Format)
4. Coordinates file (see below - 'Construction of Coordinates File')


Construction of Reference
--------------------------

The REFERENCE will serve as the reference for the read mapping and should be in fasta format. 

The REFERENCE should contain the invertible DNA sequence with 1000 bp flanking sequence for both orientations. For example:

------------------>>>>>>>>>>------------------/------------------<<<<<<<<<<------------------

* ------ = Flanking regions
* >>>>>> = Invertible DNA sequence, orientation 1
* <<<<<< = Invertible DNA sequence, orientation 2
 
The flanking sequences remain the same for each orientation, while the invertible DNA switch should be reverse complemented. In this way, both orientations are represented in the REFERENCE. 

Construction of Coordinates File:
-----------------------------------

In order to count reads traversing the bordering regions of the invertible DNA switch, it is necessary to assign reads to any of 6 regions in the REFERENCE, as in the above figure (you'll notice that there are six regions - 4 flanking regions and 2 invertible switch regions in opposing orientations). Assuming that these orientations are OFF and ON, the six regions are OFF-left-flank, OFF-switch-region, OFF-right-region, and similarly for the ON orientation. 

The start and end coordinates for each region is necessary for the assignation of reads to their correct region. Therefore, the script needs to read in a text file with these coordinate. The file needs to be in the below format, with the number changed accordingly::

	Region	Start	End
	A_left_flank	n/a	1000
	A_switch_region	1001	1313
	A_right_flank	1314	2313
	B_left_flank	2314	3313
	B_switch_region	3314	3626
	B_right_flank	3627	n/a
	
The 'n/a' regions are irrelevant as they represent the lower most and uppermost regions. The file needs to have a header, and should be **tab delimited**.

**The file needs to be called coordinates.txt** and should be in the directory above the reads.


Construction of Bedmap_reference
----------------------------------

The BEDMAP_REFERENCE will specify the location of the desired overlap regions and should be in .bed format, which can be generated from a .gff file using gff2bed.

You can generate this by using Bedops::

 $ gff2bed < genes.gff > genes.bed
 $ bam2bed < reads.bam > reads.bed


* It is very important that the nomenclature stays consistent between the reference sequences, particularly regarding the naming of the reference sequence. i.e. The fasta header for the reference should match that in the .bed file.*

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

The exons for the DNA switch analysis were created as 10 bp pseudo-exons overlapping each end of the DNA switch regions, totally 4 pseudo-exons regions. 

Fastq File Format
---------------------

For the script to work, the fastq files should be named in a way similar to this::

 $name_1.fastq
 $name_2.fastq
 $name_1.fastq.gz
 $name_2.fastq.gz

The read files can be zipped or unzipped. 


How-To: Run Me
---------------

This bash scripts requires the reference sequences to have already been generated. Furthermore, the reads (illumina paired end - not interleaved) need to be in a directory below the reference files, and the bash script should be executed in the file containing the reads.

To run the script, simply type::

 bash <script> <REFERENCE> <BEDMAP_REFERENCE>

Output
-------

The script will generate directories for each strain containing the BAM and BAI files, and the bedmaps results. 
**NOTE** that the script will delete the original fastq files and the SAM file.

Two other files will also be created:

1. Bedmap_results.csv - The concatenated results for the bedmaps counts of reads overlapping the provided exon locations
2. Paired_read_results.csv - The concatenated results for the paired-end read counts which traverse the region of interest

**NOTE**: The script works based on counting the overlapping reads for two orientations of an invertible DNA region. Thus, the input requires a REFERENCE with opposing orientations of an invertible DNA switch, arbitrarily named OFF and ON. The output assumes that the REFERENCE has been designed with OFF orientation first (i.e leftmost), and the ON orientation second. 
