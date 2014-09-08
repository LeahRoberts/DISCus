#!/bin/bash

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name). Also check that the .bed reference contains exon features (not just CDS)

# Change the reference and bedmap reference files according to your study. These were located above the directory containing the fastq reads

REFERENCE=$1
BEDMAP_REFERENCE=$2

# Need to be in the reads directory
# If the reference has already been indexed once, there's no need to index it again and the next two commands can be commented (#) out

echo "indexing " $REFERENCE
/home/leah/bin/bwa/bwa index $REFERENCE

# Write a loop that will take the .sai files from the previous section, and map them to the reference strain generating .sam and .bam files.

for f in *
do
        if [[ $f == *fastq.gz ]]
        then

# Parse out just the name of the strain (again, in the format $strainname_1.fastq) 
                name=$(echo $f | cut -f1 -d.)

                if [[ ! -e $name.bam ]]
                then
                	echo "performing alignment on " $name
# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):

                	if [[ $(echo $name1 | cut -f1 -d.) == $(ls $name2 | cut -f1 -d.) ]]
                	then
                		/home/leah/bin/bwa/bwa mem $REFERENCE -p $f > $name.sam 
				samtools view -bS -F 4 $name.sam > $name.bam
                        	samtools sort $name.bam $name.sorted
                        	samtools index $name.sorted.bam
                        	rm $name1.sai $name2.sai $name.sam
                        	rm $name.bam
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

