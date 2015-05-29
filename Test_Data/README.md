Test Data for DISCus
====================

Simulated illumina paired-end read data was generated using GemSIM (McElroy et al., 2012) using the genome of the *Escherichia coli* ST131 strain EC958 (GI:802098267) with examples of the fimS region in either the OFF or ON orientation. Reads were then shuffled between orientations to give a known ratio of 60%:40% 'OFF' reads to 'ON' reads.

To generate reference files:
----------------------------

Using the script:

	$ python DISCus_create_reference.py EC958_complete.fasta 5018228 5018540 EC958_OFF60


To run DISCus:
---------------

The reference files should be in a directory above fastq files. Run the script in the folder with the fastq files:

	$ bash ~/PATH/TO/DISCUS/DISCus_general.sh ../EC958_OFF60.fa ../EC958_OFF60.bed ../EC958_OFF60-coordinates.txt	
