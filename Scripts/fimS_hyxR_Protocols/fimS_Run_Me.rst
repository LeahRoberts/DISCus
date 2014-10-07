*fimS* Analysis 
=================

Analysis of the type 1 fimbriae invertible promoter (*fimS*) switching frequency in two datasets:

1. ST131 100 Dataset 
2. EC958 Fimbrial Switch Mutant Strains

Reference Files
----------------

1. EC958_OFF_ON_alleles.fa
2. EC958_10bps_only_off_on.bed

The reference files should be in the directory above the reads files. 

How-to: Run Me
---------------

Execute script in the directory with the fastq reads, formatted as below::

  name_1.fastq
  name_2.fastq
  
Or::

  name_1.fastq.gz
  name_2.fastq.gz
  
*NOTE*: reads can be a combination of zipped or unzipped. The read pairs MUST have the same name.

To run the script::

  $ bash FiSC_fimS.sh ../EC958_OFF_ON_alleles.fa ../EC958_10bps_only_off_on.bed

Output
-------

* Directories named for each strain in the analysis, containing the BAM and BAI files, as well as the bedmaps results file. 
* Two csv files, containing the concatenated results of the bedmaps analysis and the paired-end read analysis.

1. fimS_OFF_ON_bed_results.csv - results of the bedmaps analysis
2. fimS_OFF_ON_positions.csv - results of the paired-end read analysis

The csv files are formatted: "strain_name,OFF_1,OFF_2,ON_1,ON_2". 

* OFF_1/ON_1 refers to the left border of the OFF/ON orientation, where OFF_2/ON_2 refers to the right border. 


