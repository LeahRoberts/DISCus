*hyxR* Promoter Analysis
=========================

Analysis of the *hyxR* invertible DNA promoter (downstream of fimX in the PAI-X pathogenicity island) switching frequency in two datasets:

1. The ST131_100 Dataset
2. EC958 Fimbrial Switch Mutant Strains

Reference Files
----------------

1. hyxR_pseudo_reference.fa
2. hyxR_pseudo_reference_with_exons.bed

The reference files should be in the directory above the reads.

How-to: Run-Me
----------------

Execute script in the same directory as the fastq reads, formatted as below::

  name_1.fastq
  name_2.fastq
  
Or::

  name_1.fastq.gz
  name_2.fastq.gz

*NOTE*: Reads can be a combination of zipped or unzipped. The read pairs MUST have the same name. 

To run the script::
  
  $ bash DiSCO_hyxR.sh ../hyxR_pseudo_reference.fa ../hyxR_pseudo_reference_with_exons.bed
  
Output
-------

* Directories named for each strain in the analysis, containing BAM and BAI files, as well as the bedmaps results file.
* Two csv files, containing the concatenated results of the bedmaps analysis and the paired-end reads analysis:

1. hyxR_OFF_ON_bed_results.csv - results of the bedmaps analysis
2. hyxR_OFF_ON_paired-end.csv - results of the paired-read analysis

The csv files are formatted "Name,OFF_1,OFF_2,ON_1,ON_2".

* OFF_1/ON_2 refers to the left border of the OFF/ON orientation
* OFF_2/ON_2 refers to the right border
