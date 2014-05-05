#!/bin/bash

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name). Also check that the .bed reference contains exon features (not just CDS)

REFERENCE="../hyxR_pseudo_reference.fa"
BEDMAP_REFERENCE="../hyxR_pseudo_reference_with_exons.bed"
FILETYPE="fastq"
#DIRECTORY="./reads"

# Need to be in the reads directory
# If the reference has already been indexed once, there's no need to index it again and the next two commands can be commented (#) out

#echo "indexing " $REFERENCE
#bwa index $REFERENCE

for f in *
do
	echo "processing $(ls $f | cut -f1 -d.)"

		if [[ $f == *_1.fastq ]]
		then
			read1=$f
			name1=$(ls $f | cut -f1 -d.)
			bwa aln $REFERENCE $read1 > $name1.sai
#			echo $name1
#			echo $read1
		else
			read2=$f
			name2=$(ls $f | cut -f1 -d.)	
			echo "paired-end read " $name2
			bwa aln $REFERENCE $read2 > $name2.sai
#			echo $name2
#			echo $read2
	fi
done

for f in *
do
	name=$(ls $f | cut -f1 -d_)
	echo "performing alignment on " $name

		if [[ $(ls $f | cut -f1 -d.) == *_1* ]]
		then
			name1=$(ls $f | cut -f1 -d.)
#			echo $name1
	
		elif [[ $(ls $f | cut -f1 -d.) == *_2* ]]
		then
			name2=$(ls $f | cut -f1 -d.)
#			echo $name2		
		
			if [[ $(ls $name1 | cut -f1 -d_) == $(ls $name2 | cut -f1 -d_) ]]
	       	        then    
  	             	        bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
                        	samtools view -bS $name.sam > $name.bam
                    	    	samtools sort $name.bam $name.sorted
                     		samtools index $name.sorted.bam         

        	       	else    
    	                    	echo $name1 " and " $name2 " are not a pair"

			fi
		
		else
			echo "oh no something must be wrong"	

	fi
	
done

for f in *
do
	if [[ $f == *.sai ]]
	then
		echo "deleting " $f
		rm $f

	elif [[ $f == *.sorted.bam ]] 
	then
		echo "converting " $f " to bed file"
		bam2bed < $f > $f.bed
		bedmap --echo --count $BEDMAP_REFERENCE $f.bed > $f.result.bed
		echo "results for " $f " have finished compiling"
	fi
done

for f in *
do
        NAME=$(ls $f | cut -f1 -d.)

                if [[ $f == *.result.bed ]]
                then
                #       echo $NAME
                        EXONS=$(cut -f4 -d$'\t' $f)
                        COUNTS=$(cut -f2 -d\| $f)
                        
                        echo $NAME $'\n' $EXONS $'\n' $COUNTS >> result.txt
                #       echo $COUNTS 

        fi
done

echo "Finished!"
