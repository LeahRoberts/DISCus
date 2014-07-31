#!/bin/bash

# Script to map reads to a reference using bwa, filter for reads mapping to the reference and count read overlap at particular regions 
# Reads are paired end but no interleaved
# reference is fasta format

# Need to have bwa and samtools installed
# need to execute in the directory which has the reads
# the reference needs to be one directory above
# The reads need to be in the format strainname_1.fastq for the script to work properly

# To run the script: bash bwa_aligner_read_overlap_count.sh $REFERENCE $IS

# The reference should be in the format ../reference.fa
# The IS should be in the form 'IS', and should match the name of the reference fasta file (to check this, simply cat the reference file and see what it's called). 

REFERENCE=$1
IS=$2

bwa index $REFERENCE

for f in *
do

# Used cut to parse out the strain name - THIS WILL ONLY WORK IF THE FASTQ FILE IS FORMATTED CORRECTLY (i.e. $strainname_1.fastq)
        echo "processing $(ls $f | cut -f1 -d.)"

                if [[ $f == *_1.fastq ]]
                then
			read1=$f
                        name1=$(ls $f | cut -f1 -d.)
                        bwa aln $REFERENCE $read1 > $name1.sai
                elif [[ $f == *_2.fastq ]]
		then
                        read2=$f
                        name2=$(ls $f | cut -f1 -d.)
                        bwa aln $REFERENCE $read2 > $name2.sai
        fi
done

# Write another loop that will take the .sai files from the previous section, and map them to the reference strain generating .sam and .bam files.

for f in *
do
	if [[ $f == *fastq ]]
	then
# Parse out just the name of the strain (again, in the format $strainname_1.fastq)
        	name=$(ls $f | cut -f1 -d_)

		if [[ ! -e $name.bam ]]
         	then
          		echo "performing alignment on " $name

# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):

                 	if [[ $(ls $f | cut -f1 -d.) == *_1* ]]
                 	then
                        	name1=$(ls $f | cut -f1 -d.)

	                 elif [[ $(ls $f | cut -f1 -d.) == *_2* ]]
        	         then
                		name2=$(ls $f | cut -f1 -d.)

				if [[ $(ls $name1 | cut -f1 -d_) == $(ls $name2 | cut -f1 -d_) ]]
                         	then
                                	echo "creating sam file"
                                 	bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
					echo "finished creating .sam file"
					samtools view -bS $name.sam > $name.bam
					samtools sort $name.bam $name.sorted
					samtools index $name.sorted.bam
					echo "Filtering reads that map to the reference"
					samtools view $name.sorted.bam $IS > $name.sorted.mapped.sam
				fi
			fi
		fi
	fi
done

# Cleanup the files you don't need (i.e. everything but the velvet.bam file):
mkdir tmp
for f in *
do
        if [[ $f == *.fastq ]]
        then
                mv $f tmp/

        elif [[ $f != *.sorted.mapped.sam ]]
        then
                rm $f
        fi
        mv tmp/* ./
        rm -r tmp/
done

# Create bam file to parse the read names that map to the desired reference

for f in *
do
	if [[ $f == *mapped.sam ]]
	then
		name=$(ls $f | cut -f1 -d.)
		echo "indexing fasta file to create .fai file"
		samtools faidx $REFERENCE
		echo "creating bam file from sam file"
		samtools view -bt $REFERENCE.fai $name.sorted.mapped.sam > $name.velvet.bam
		echo "finished creating bam file for velvet"
		rm $name.sorted.mapped.sam
	fi
done

# Take the velvet file and parse the read names from it using a python script. The python script must be in the ~/bin directory

for f in *
do
	if [[ $f == *.velvet.bam ]]
	then
		samtools view $f | cut -f1 -d$'\t' > list.txt
		cat list.txt > list_2.txt
		sed -i 's/$/\/1/g' list.txt
		sed -i 's/$/\/2/g' list_2.txt
		python ~/bin/fastq_parser.py $name\_1.fastq $name\_2.fastq
		echo "renaming output fastq files"
		mv $name\_reads_1.fastq $name\_1.fastq
		mv $name\_reads_2.fastq $name\_2.fastq

# The reads are parsed back in to new fastq files, which can then be remapped to another reference. This is done using another script, also located in the ~/bin directory:

		bash ~/bin/align_to_contigs.sh
	fi
done
