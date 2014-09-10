#!/bin/bash

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name). Also check that the .bed reference contains exon features (not just CDS)

# Change the reference and bedmap reference files according to your study. These were located above the directory containing the fastq reads

REFERENCE=$1
BEDMAP_REFERENCE=$2

# Need to be in the reads directory
# If the reference has already been indexed once, there's no need to index it again and the next two commands can be commented (#) out

echo "indexing " $REFERENCE
bwa index $REFERENCE

# Write a loop that will generate .sai files by reading each of the paired strains to the reference strain using the tool BWA. These are outputted as $strain-name_1.sai or $strain-name_2.sai

for f in *
do

# Used cut to parse out the strain name - THIS WILL ONLY WORK IF THE FASTQ FILE IS FORMATTED CORRECTLY (i.e. $strainname.R1.fastq) 
	echo "processing $(echo $f | cut -f1 -d.)"

		if [[ $f == *.R1.fastq ]]
		then
			read1=$f
			name1=$(echo $f | cut -f1 -d.)
			bwa aln $REFERENCE $read1 > $name1.R1.sai
		elif [[ $f == *.R2.fastq ]]
		then
			read2=$f
			name2=$(echo $f | cut -f1 -d.)	
			echo "paired-end read " $name2
			bwa aln $REFERENCE $read2 > $name2.R2.sai
	fi
done

# Write another loop that will take the .sai files from the previous section, and map them to the reference strain generating .sam and .bam files.

for f in *
do
        if [[ $f == *fastq ]]
        then

# Parse out just the name of the strain (again, in the format $strainname_1.fastq) 
                name=$(echo $f | cut -f1 -d.)

                if [[ ! -e $name.bam ]]
                then
                echo "performing alignment on " $name
# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):

                        if [[ $f  == *R1* ]]
                        then
                                name1=$(echo $f | cut -f1-2 -d.)
#                               echo $name1
        
                        elif [[ $f == *R2* ]]
                        then
                                name2=$(echo $f | cut -f1-2 -d.)
#                               echo $name2             
                
                                if [[ $(echo $name1 | cut -f1 -d.) == $(ls $name2 | cut -f1 -d.) ]]
                                then
                                        bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
                                        samtools view -bS -F 4 $name.sam > $name.bam
                                        samtools sort $name.bam $name.sorted
                                        samtools index $name.sorted.bam
                                        rm $name1.sai $name2.sai $name.sam
                                        rm $name.bam

                                else
                                        echo $name1 " and " $name2 " are not a pair"

                                fi
                
                        else
                                echo "oh no something must be wrong"    
                        fi
                fi
        fi
done

# The above command should have generated the .sam and .bam alignment files for all the reads against the reference. 
# The next command then counts all the reads aligning the "exons" defined by the BEDMAP_REFERENCE file.
# This BEDMAP_REFERENCE file was generated using gff2bed (see bedops documentation). 

mkdir ../tmp/

for f in *
do
	
# Getting rid of the .sai files to clean up the directory

	if [[ $f == *.sorted.bam ]]
	then
		mv $f ../tmp
		
	elif [[ $f == *.bai ]]
	then
		mv $f ../tmp
	else
		echo "deleting " $f
		rm $f
	fi
	
	mv ../tmp/* ./
done	

rmdir ../tmp/

# Taking the sorted .bam files and converting them to .bed files, and then using bedmaps to count the reads overlapping the 'exon' regions (in this case, the borders of the invertible DNA switch for hyxR)
#The results for each strain are saved as $name.resut.bed

for f in *
do
	if [[ $f == *.sorted.bam ]] 
	then
		echo "converting " $f " to bed file"
		bam2bed < $f > $f.bed
		bedmap --echo --count $BEDMAP_REFERENCE $f.bed > $f.result.bed
		echo "results for " $f " have finished compiling"
	fi
done

# Loop to parse all of the $name.result.bed file contents into a single result.txt file with the strain identifier.

for f in *
do
       	NAME=$(echo $f | cut -f1 -d.)
        
	        if [[ $f == *.result.bed ]]
                then
                        POS_start=$(cut -f2 -d$'\t' $f)
			POS_end=$(cut -f3 -d$'\t' $f)
                        COUNTS=$(cut -f2 -d\| $f)
                        
                        echo $NAME','$POS_start'..'$POS_end','$COUNTS >> bedmap_results.csv

        fi
done

echo "Finished!"

