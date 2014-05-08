8-5-14

How to run analysis of hyxR switch
===================================

bwa alignment + bedops counting of reads

Necessary Files
---------------

Reference Sequence:

	* hyxR_pseudo_reference.fa    *using bwa to maps reads to this reference*
	* hyxR_pseudo_reference_with_exons.bed    *using bedops to count the reads overlapping exons defined in this file*

**NOTE**: The header for the reference fasta file must match that given in the .bed file (otherwise the counted read output will be zero)


Script:

	* bwa_align_bedops_count.sh     *does all of the read mapping and counting of reads*

**NOTE**: There are some glitches with this script, mainly that similarly named sequences (eg. HVM5 and HVM52) will be mapped together. Files like this should be separated into different directories until the script is fixed. 


How-To Run Me
-------------

1. Make sure the fastq files are paired-end and in a single directory.
2. Make sure both reference files are in the directory above the reads.
3. Run the script while in the reads directory:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
leah@binf-training:~/Raw_Data/FimS_Indian_strains_2$ bash ~/bin/BWA_Bedops_aligner/bwa_align_bedops_count.sh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4. Wait for result.txt file to be generated.

**NOTE**: space is also a problem - May need to check that the script is still running, as it may halt/freeze/abort if it runs out of space.  
