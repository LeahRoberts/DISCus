BWA_Bedops_aligner
==================

Script for aligning reads with BWA and counting read overlap using Bedops

How-To: Run Me
---------------

This bash scripts requires the reference sequences to have already been generated. Furthermore, the reads need to be in a directory below the reference files, and the bash script should be executed in the file containing the reads.

The reference files required are:
  * A fasta file for mapping the reads to
  * A .bed file generated from a .gff file using gff2bed (see bedops) containing the exons you want to count reads overlapping
  
*It is very important that the nomenclature stays consistent between the reference sequences, particularly regarding the naming of the reference sequence. i.e. The fasta header for the reference should match that in the .bed file.*


Errors with Script
--------------------

The script iterates through each strain twice (generating .bam and .sam files twice, which overwrite). This doesn't create a problem for the results, it just means that the script takes twice as long to run.

Strains with similar names (eg. HVM5 and HVM52) have a tendency to incorrectly align together. To solve this, I have been separating files with similar names and running the script in separate directories. 
