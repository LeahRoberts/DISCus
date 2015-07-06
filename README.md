DISCus - *D*NA *I*nver*S*ion *C*ounter
==========================================

Simple Overview:
----------------
DISCus is a script for mapping illumina paired end reads with BWA to an invertible DNA pseudo-reference and counting read overlap using: (1) Bedops; and (2) read-pairs traversing the desired region.

Basic Commands:
---------------
To generate references:

	$ python DISCus_create_reference.py <file.fasta> <start_coordinate> <end_coordinate> <out_name>

To run DISCus (with full path to files):

	$ bash ~/PATH/TO/DISCus/DISCus_general.sh <ref.fasta> <ref.bed> <coordinates.txt>

Installation Requirements
--------------------------

1. BWA version: 0.7.12 (http://sourceforge.net/projects/bio-bwa/files/)
2. SAMtools version: 1.2 (using htslib 1.2) (http://samtools.sourceforge.net/)
3. Bedtools version: v2.23.0 (http://bedtools.readthedocs.org/en/latest/content/installation.html)
4. Bedops version: 2.4.14 (http://bedops.readthedocs.org/en/latest/content/installation.html) 
5. Python version: 2.7
6. Biopython version: 1.64 (http://biopython.org/wiki/Main_Page)

File Requirements
------------------

1. REFERENCE in fasta format (see below - 'Construction of Reference')
2. BEDMAP_REFERENCE in bed format. (see below - 'Construction of Reference' and 'Construction of Bedmap_reference')
3. illumina paired-end read files NOT interleaved, named: <strain>_1.fastq, _2.fastq (see below - Fastq File Format)
4. Coordinates file (see below - 'Construction of Reference' and 'Construction of Coordinates File')


Test Data
==========

**Step 1.** Make a directory called "reads".

**Step 2.** Download the paired fastq files from EMBL-EBI for the ST131 strains S37EC (ERR161302) and HVM1147 (ERR161318) into the "reads" directory:

	$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR161/ERR161302/ERR161302_*.fastq.gz
	$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR161/ERR161318/ERR161318_*.fastq.gz
	

**Step 3.** Above the reads directory, download the complete genome for *Escherichia coli* O25b:H4-ST131 str. EC958 from NCBI (NZ_HG941718.1, GI:802098267) (complete fasta sequence).

**Step 4.** Generate reference fasta, bed and coordinate files using the python script DISCus_create_reference.py:

	$ python ~/PATH/TO/DISCus_create_reference.py EC958_complete.fasta 5018228 5018540 EC958_OFF_ON

*If you get errors, make sure python and biopython are correctly installed.*
	

**Step 5.** Check that all the files are correctly placed (fastq files in reads directory, all other files in directory above reads). 

**Step 6.** Run DISCus from reads directory:

	$ bash ~/PATH/TO/DISCus_general.sh ../EC958_OFF_ON.fa ../EC958_OFF_ON.bed ../EC958_OFF_ON-coordinates.txt


Output results should be:

*Bedmaps results:*

	STRAIN,A_1,A_2,B_1,B_2
	ERR161302,61,70,3,0
	ERR161318,50,34,18,15

*Paired_read_results:*

	ERR161302,39,43,4,1
	ERR161318,17,17,9,0


Where A_1 and A_2 indicate reads overlapping the first orientation in the pseudo-reference (OFF), and B_1 and B_2 indicate reads overlapping the second orientation (ON). 

As a percentage:

(Bedmaps):
	ERR161302 = 98% OFF
	ERR161318 = 72% OFF

(Paired Reads):
	ERR161302 = 94% OFF
	ERR161318 = 79% OFF



DISCus Manual
=============

Introduction:
-------------

DISCus is a command-line tool created to interrogate inversion orientation of an invertible DNA region using Illumina paired-end sequencing data. First, a pseudo-reference is generated containing both orientations of a known invertible DNA region with 1 kb of flanking region for each orientation. Illumina paired-end reads are then mapped to this pseudo-reference, where reads will map primarily to whichever orientation was more represented in the bacterial population at the time of sequencing. Two methods are used to quantify the number of reads mapping to either orientation of the invertible region: (1) count the number of reads directly overlapping the unique bordering regions of both orientation, and (2) count the number of paired-reads that traverse from either flanking region to a switch region.

By quantifying the amount of reads mapping to either orientation of a known invertible DNA region, DISCus is able to sensitively quantify the proportion of 'forward' and 'reverse' orientations within a bacterial population. 


Construction of Reference
--------------------------

**Automated:**

The python script, *DISCus_create_reference.py*, can automatically generate a fasta pseudo-reference, bed file reference and coordinates file for analysis with DISCus using a fasta file of the genome of interest.
The script takes in four arguments and can be exected as shown below::

 	$ python DISCus_create_reference.py <genome_sequence.fasta> <start_coordinate> <end_coordinate> <filename>
 
Where **start_coordinate** is the start of the invertible DNA region of interest, and **end_coordinate** is the end of the invertible DNA region. The script also requires a filename, which will become the name of the output file as well as the fasta header. 

The fasta header generated will be::

 	> <filename> _ <start_coordinate> _ <end_coordinate>
	 Sequence...
 
The sequence will include 1000 bp of flanking region, as well as both orientations of the invertible DNA region.

**Example:**

Using the EC958_complete.fasta genome as the input, and wanting a 100 bp inversion region between 182100 and 182200, the command to execute the script would be::
 
  	$ python DISCus_create_reference.py EC958_complete.fasta 182100 182200 EC958_100bp

*Note: the script will not run unless all four arguments are given. The filename argument should be without spaces.*

The output of this will be:

1. a pseudo-reference with the header ">EC958_100bp_182100_182200"
2. a bedmaps reference file indicating 10 bp overlap regions on either end of the invertible DNA region of interest (further explained in "Construction of Bedmap_reference" section)
3. A coordinates file (txt) defining the regions of the pseudo-reference necessary for determining paired-end read traversal


**Manual:**

You can also construct the reference manually, as explained below. This section can also serve as further explanation regarding the pseudo-reference and theory of the script. 

The REFERENCE will serve as the reference for the read mapping and should be in fasta format. 

The REFERENCE should contain the invertible DNA sequence with 1000 bp flanking sequence for both orientations. For example:

------------------>>>>>>>>>>------------------/------------------<<<<<<<<<<------------------

* ------ = Flanking regions
* >>>>>> = Invertible DNA sequence, orientation 1
* <<<<<< = Invertible DNA sequence, orientation 2
 
The flanking sequences remain the same for each orientation, while the invertible DNA switch should be reverse complemented. In this way, both orientations are represented in the REFERENCE. 

Construction of Coordinates File:
-----------------------------------

This process is now automated by a python script, as described above.

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


Construction of Bedmap_reference
----------------------------------

This process is now automated by a python script, as described above.

The BEDMAP_REFERENCE will specify the location of the desired overlap regions and should be in .bed format.

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
Simple overview::

	$ bash PATH/TO/DISCus_general.sh <PATH_to_fasta_ref> <PATH_to_bedmaps_ref> <PATH_to_coordinates_file>
	

This bash scripts requires the reference sequences to have already been generated. Furthermore, the reads (illumina paired end - not interleaved) need to be in a directory below the reference files, and the bash script should be executed in the file containing the reads, according to this diagram:

![ScreenShot](https://github.com/LeahRoberts/DiSCus/blob/master/DISCus_how_to_run.png)


Output
-------

The script will generate directories for each strain containing the BAM and BAI files, and the bedmaps results. 


Two other files will also be created:

1. Bedmap_results.csv - The concatenated results for the bedmaps counts of reads overlapping the provided exon locations
2. Paired_read_results.csv - The concatenated results for the paired-end read counts which traverse the region of interest

**NOTE**: The script works based on counting the overlapping reads for two orientations of an invertible DNA region. Thus, the input requires a REFERENCE with opposing orientations of an invertible DNA switch, arbitrarily named OFF and ON. The output assumes that the REFERENCE has been designed with OFF orientation first (i.e leftmost), and the ON orientation second. 

Licence
--------

All scripts and files within this directory are under the MIT Licence (http://choosealicense.com/licenses/mit/)
