#!/bin/bash

# Bash script to map the reads from the 99 ST131 strains back to their velvet contigs

# Make all of the .sai files

for f in *
do
	if [[ $f == *.fastq ]]
	then

		name=$(ls $f | cut -f1 -d_)
		echo "processing " $name
	
		REFERENCE=../../ST131_99_ord/$name\_[[:digit:]]*\_Contigs.fas
		echo "Reference is " $REFERENCE
	
		bwa index $REFERENCE
	
		if [[ $f == *_1.fastq ]]
		then
			read1=$f
			name1=$(ls $f | cut -f1 -d.)
			bwa aln $REFERENCE $read1 > $name1.sai
			#echo $name1
			#echo $read1
		elif [[ $f == *_2.fastq ]]
		then
			read2=$f
			name2=$(ls $f | cut -f1 -d.)
			#echo "paired-end read " $name2
			bwa aln $REFERENCE $read2 > $name2.sai
			#echo $name2
			#echo $read2
		fi
	fi
done

for f in *
do
	if [[ $f == *fastq ]]
	then
#  Parse out just the name of the strain (again, in the format $strainname_1.fastq) 
        name=$(ls $f | cut -f1 -d_)
        echo $name       

		if [[ ! -e $name.bam ]]
        	then
         	echo "performing alignment on " $name

# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):

                	if [[ $(ls $f | cut -f1 -d.) == *_1* ]]
                	then
                        	name1=$(ls $f | cut -f1 -d.)
#                       	echo $name1

                	elif [[ $(ls $f | cut -f1 -d.) == *_2* ]]
                	then
                        	name2=$(ls $f | cut -f1 -d.)
#                       	echo $name2             

			 	if [[ $(ls $name1 | cut -f1 -d_) == $(ls $name2 | cut -f1 -d_) ]]
                        	then
                                	echo "creating sam file"
                                	bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
					echo "finished creating .sam file"
					samtools view -bS $name.sam > $name.bam
					samtools sort $name.bam $name.sorted
					samtools index $name.sorted.bam
					rm $name.bam $name.sam $name1.sai $name2.sai
					echo "finished read-mapping for " $name
				fi
			fi
		fi
	fi
done
