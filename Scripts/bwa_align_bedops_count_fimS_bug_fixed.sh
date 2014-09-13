#!/bin/bash

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name). Also check that the .bed reference contains exon features (not just CDS)

# Change the reference and bedmap reference files according to your study. These were located above the directory containing the fastq reads

REFERENCE=$1
BEDMAP_REFERENCE=$2

# Need to be in the reads directory
# If the reference has already been indexed once, there's no need to index it again and the next two commands can be commented (#) out

echo "indexing " $REFERENCE
bwa index $REFERENCE

echo "creating bam files"

for f in *
do
        if [[ $f == *fastq ]]
        then

# Parse out just the name of the strain (again, in the format $strainname_1.fastq) 
                name=$(echo $f | cut -f1 -d_)

                if [[ ! -e $name.bam ]]
                then
                	echo "checking for pairs - " $name
	
			if [[ $(echo $f | cut -f1 -d. | cut -f2 -d_) == "1" ]]
			then
				name1=$f
				echo "first paired read = " $name1
# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):
				for g in *
				do
					if [[ $(echo $g | cut -f1 -d. | cut -f2 -d_) == "2" ]] && [[ $(echo $g | cut -f1 -d_) == $name ]]
                        		then
						echo $f "and" $g "are a pair - performing alignment"
						/home/leah/bin/bwa/bwa mem $REFERENCE $f $g > $name.sam
						samtools view -bS -f 0x0002 -F 4 $name.sam > $name.bam
						samtools sort $name.bam $name.sorted
						samtools index $name.sorted.bam
                        		else
                        		        echo $f " and " $g " are not a pair"

                        		fi
				done
			
			elif [[ $(echo $f | cut -f1 -d. | cut -f2 -d_) == "2" ]]
			then
				name2=$f
				echo "second paired read = " $name2
				
				for g in *
				do
					if [[ $(echo $g | cut -f1 -d. | cut -f2 -d_) == "1" ]] && [[ $(echo $g | cut -f1 -d_) == $name ]]
					then
						echo $f "and" $g "are a pair"
						/home/leah/bin/bwa/bwa mem $REFERENCE $f $g > $name.sam
						samtools view -bS -f 0x0002 -F 4 $name.sam > $name.bam
						samtools sort $name.bam $name.bam.sorted
						samtools index $name.bam.sorted
						
					else
						echo $f "and" $g "are not a pair"
					fi
				done
               		fi
                fi
        fi
done
