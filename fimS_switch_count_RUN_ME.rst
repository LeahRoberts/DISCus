8-5-14

How to run analysis of fimS switch 
===================================

bwa alignment + bedops counting of reads/counting of read-pair overlap from bam file coordinates 

Necessary Files
----------------

Reference Sequences:

	* EC958.OFF.fa *using bwa to map reads to this reference*
	* EC958_10bps_only.bed *using bedops to count the reads overlapping the exons defined in this file*

*Due to changes in the script, the OFF and ON alleles must now be separated into different references*

**NOTE**: The header for the reference fasta file must match that given in the .bed file (otherwise the counted read output will be zero)


Script:

	* <script to be renamed>  *does all of the read mapping and counting of reads*

**NOTE**: There are some glitches with this script, mainly that similarly named sequences (eg. HVM5 and HVM52) will be mapped together. Files like this should be separated into different directories until the script is fixed. 

How-To Run Me
--------------

1. Make sure the fastq files are paired-end and in a single directory.
2. Make sure both reference files are in the directory above the reads.
3. Run the script while in the reads directory:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
leah@binf-training:~/Raw_Data/FimS_Indian_strains_2$ bash <script> EC958.OFF.ON.alleles.fa EC958_10bps_only.bed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4. Wait for fim_OFF_bed_results.csv file to be generated. *Similarly, a 'fim_ON_bed_results.csv' file can also be generated via the sister script.*


**NOTE**: space is also a problem - May need to check that the script is still running, as it may halt/freeze/abort if it runs out of space. 
