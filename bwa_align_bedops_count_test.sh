#!/bin/bash

# The script should be run from inside the directory with the fq reads: 
# $ bash bwa_bedops_aligner_count.sh $REFERENCE $BEDMAP_REFERENCE

# the REFERENCE variable should be in the format ../reference.fa
# the BEDMAP_REFERENCE should be in the format ../bedmap_format.fa

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name). Also check that the .bed reference contains exon features (not just CDS)

REFERENCE=$1
BEDMAP_REFERENCE=$2

# First need to index the reference fasta file to be used in the mapping:
echo "indexing " $REFERENCE
bwa index $REFERENCE

# Write a loop that will generate .sai files by reading each of the paired strains to the reference strain using the tool BWA. These are outputted as $strain-name_1.sai or $strain-name_2.sai

for f in * 

do

# Used cut to parse out the strain name - THIS WILL ONLY WORK IF THE FASTQ FILE IS FORMATTED CORRECTLY (i.e. $strainname_1.fastq) 
	echo "processing $(ls $f | cut -f1 -d.)"

		if [[ $f == *_1.fastq ]]
		then
			read1=$f
			name1=$(ls $f | cut -f1 -d.)
			bwa aln $REFERENCE $read1 > $name1.sai
		else
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
		
# Want to only perform the alignment once - as there are two .sai files, this has the potential to interate through twice. This if statement prevents the script from iterating through more than once on the same strain:
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
  	             	        	bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
                        		samtools view -bS -F 4 $name.sam > $name.bam
                    	    		samtools sort $name.bam $name.sorted
                     			samtools index $name.sorted.bam         

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

# Getting rid of the .sai files to clean up the directory

for f in *
do
        if [[ $f == *.sam ]]
        then
                echo "deleting " $f
                rm $f

        elif [[ $f == *.sai ]]
        then
                echo "deleting " $f
                rm $f

# The next command then counts all the reads aligning the "exons" defined by the BEDMAP_REFERENCE file.
# This BEDMAP_REFERENCE file was generated using gff2bed (see bedops documentation). 

# Taking the sorted .bam files and converting them to .bed files, and then using bedmaps to count the reads overlapping the 'exon' regions (in this case, the borders of the invertible DNA switch for hyxR)
#The results for each strain are saved as $name.resut.bed

	elif [[ $f == *.sorted.bam ]] 
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
        NAME=$(ls $f | cut -f1 -d.)

                if [[ $f == *.result.bed ]]
                then
                        EXONS=$(cut -f4 -d$'\t' $f)
                        COUNTS=$(cut -f2 -d\| $f)
                        
                        echo $NAME $'\n'$EXONS $'\n'$COUNTS >> result.txt

        fi
done

echo "Finished!"
