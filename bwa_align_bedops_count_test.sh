#!/bin/bash

# This script is designed to align reads to a reference using BWA, and then count the 
# number of reads overlapping predefined regions. This script generates two outputs: one
# output describes the number of reads physically overlapping the predefined regions, and
# another output describing read-pairs that span across the desired regions.

# The script should be run from inside the directory with the fastq reads: 
# $ bash <script.sh> $REFERENCE $BEDMAP_REFERENCE

# It is very important that the fastq files are formatted correctly for the script to work:
# 	strain_1.fastq
# 	strain_2.fastq

# The REFERENCE is what the reads will be mapped to. It should be in the format ../reference.fa
# the BEDMAP_REFERENCE defines the desired regions. It should be in the format ../bedmap_format.fa

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name - refer to the README). 
# Also check that the .bed reference contains exon features (not just CDS).

REFERENCE=$1
BEDMAP_REFERENCE=$2

# First need to index the reference fasta file to be used in the mapping:
echo "indexing " $REFERENCE
bwa index $REFERENCE

# Write a loop that will generate .sai files by reading each of the paired strains to the 
# reference strain using the tool BWA. These are outputted as $strain-name_1.sai or $strain-name_2.sai.

for f in * 

do

# Used cut to parse out the strain name - THIS WILL ONLY WORK IF THE FASTQ FILE IS 
# FORMATTED CORRECTLY (i.e. $strainname_1.fastq).
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

# Write another loop that will take the .sai files from the previous section, and map them 
# to the reference strain generating .sam and .bam files.

for f in *
do
	if [[ $f == *fastq ]]
	then

# Parse out just the name of the strain (again, in the format $strainname_1.fastq). 
		name=$(ls $f | cut -f1 -d_)
		
# Want to only perform the alignment once - as there are two .sai files, this has the 
# potential to interate through twice. This if statement prevents the script from 
# iterating through more than once on the same strain:
		if [[ ! -e $name.bam ]]
		then
			echo "performing alignment on " $name

# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the 
# SAME strain files are being aligned to the reference (again using BWA):

			if [[ $(ls $f | cut -f1 -d.) == *_1* ]]
			then
				name1=$(ls $f | cut -f1 -d.)
	
			elif [[ $(ls $f | cut -f1 -d.) == *_2* ]]
			then
				name2=$(ls $f | cut -f1 -d.)
		
				if [[ $(ls $name1 | cut -f1 -d_) == $(ls $name2 | cut -f1 -d_) ]]
	       	        	then    
  	        
# The '-f 0x0002' flag will filter for only properly mapped reads, while the '-F 4' flag
# will get rid of reads that are not mapped to the reference.

		     	        	bwa sampe $REFERENCE $name1.sai $name2.sai $name1.fastq $name2.fastq > $name.sam
                        		samtools view -bS -f 0x0002 -F 4 $name.sam > $name.bam
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

# Getting rid of the .sai files to clean up the directory.

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

# Taking the sorted .bam files and converting them to .bed files, and then using bedmaps to count the 
# reads overlapping the 'exon' regions (in this case, the borders of the invertible DNA switch for hyxR).
# The results for each strain are saved as $name.resut.bed.

	elif [[ $f == *.sorted.bam ]] 
	then
		echo "converting " $f " to bed file"
		bam2bed < $f > $f.bed
		bedmap --echo --count $BEDMAP_REFERENCE $f.bed > $f.result.bed
		echo "results for " $f " have finished compiling"
	fi
done

# Loop to parse all of the $name.result.bed file contents into a single result.csv file with the strain identifier.

for f in *
do
        NAME=$(ls $f | cut -f1 -d.)

                if [[ $f == *.result.bed ]]
                then
                        EXONS=$(cut -f4 -d$'\t' $f)
                        COUNTS=$(cut -f2 -d\| $f)
                        
                        echo $NAME','$EXONS','$COUNTS >> fim_OFF_bed_results.csv

        fi
done

# This next section counts the number of read-pairs that overlap the 3 regions, namely Left Flank, Switch Region, 
# and Right Flank.

# Assigning of variables to tally the amount of reads in each region.
# Some of these variables also validates the correctness of the pre-generated bam file.
  
readcount=0

OFF_1=0
OFF_2=0

right_flank_no=0
left_flank_no=0
switch_region_no=0


# Reads the sorted bam file and obtains all of the read IDs.
# These are then sorted so that the paired reads are together.

echo "" > completed.reads

for f in *
do
	if [[ $f == *.sorted.bam ]]
	then
		samtools view $f | cut -f1 -d$'\t' > readnames
		sort readnames > readnames.sorted

		echo "assigning reads..."

# The next loop reads in the 'readnames.sorted' file which contains all of the read IDs. As there are two
# reads with the same ID, the script will generate twice as many overlaps as the actual amount.
# To account for this, the script checks the current read ID against a list of already analysed read IDs.
# If the read ID exists in the list, it doesn't proceed with the rest of the loop, and adds a '+1' to the readcount tally. 
# If it doesn't exist, it adds the read name to the list and carries out the appropriate analysis. 

		while read name
		do
			if grep -q $name completed.reads
			then
				readcount=$(($readcount + 1))
			else

# For each read that is fed through the while loop, its coordinate (as well as the coordinate of its pair) are 
# obtained from the sorted bam file using 'grep' and cutting the fourth field (which is the coordinate field). 
# This is temporarily saved in a 'positions.txt' file. Variable 'a' is then assigned one of these numbers,
# while variable 'b' is assigned the other. The script then determines what region the read lies in based on its
# starting coordinate, either 'left_flank' (0-999), 'switch_region' (1000-1313) or 'right_flank' (>=1314).
 
				echo $name >> completed.reads
				samtools view B36EC.sorted.bam | grep $name | cut -f4 -d$'\t' > position.txt
#				cat position.txt

				a=$(head -1 position.txt)
#				echo $a

				if [ $a -le 999 ]
				then
					read1='left_flank'
					left_flank_no=$(($left_flank_no + 1))
#					echo "read 1 is on the left"
		
				elif [ $a -ge 1000 -a $a -le 1313 ]
				then
					read1='switch_region'
					switch_region_no=$(($switch_region_no + 1))
#					echo "read 1 is in the switch region"
			
				elif [ $a -ge 1314 ]
				then
					read1='right_flank'
					right_flank_no=$(($right_flank_no + 1))
#					echo "read 1 is on the right"
				fi

				b=$(tail -1 position.txt)
#				echo $b

				if [ $b -le 999 ]
        			then
                			read2='left_flank'
					left_flank_no=$(($left_flank_no + 1))
#					echo "read 2 is on the left"        

        			elif [ $b -ge 1000 -a $b -le 1313 ]
        			then
                			read2='switch_region'
					switch_region_no=$(($switch_region_no + 1))
#					echo "read 2 is in the switch region"        

        			elif [ $b -ge 1314 ]
        			then
                			read2='right_flank'
        				right_flank_no=$(($right_flank_no + 1))
#					echo "read 2 is on the right"
				fi

# At this point the script counts the number of reads overlapping different region. 
# Since some of the reads can be in the same region, the first if statement filters only for 
# reads that are in different regions. 
# The next if statements encompass all of the possibilites for the read positions. Depending on
# where the reads lie, a '+1' is added to the tally for either OFF_1 (i.e. read pairs in the left flank and 
# in the switch region), or OFF_2 (i.e. read pairs in the right flank and in the switch region). 

				if [[ $read1 != $read2 ]]
				then
					if [[ $read1 == 'switch_region' ]] && [[ $read2 == 'right_flank' ]]
					then
						OFF_2=$(($OFF_2 + 1))
#						echo "adding to OFF_2"		
					elif [[ $read1 == 'switch_region' ]] && [[ $read2 == 'left_flank' ]]
					then
						OFF_1=$(($OFF_1 + 1))
#						echo "adding to OFF_1"
					elif [[ $read2 == 'switch_region' ]] && [[ $read1 == 'right_flank' ]]
					then
                        			OFF_2=$(($OFF_2 + 1))
#						echo "adding to OFF_2"
					elif [[ $read2 == 'switch_region' ]] && [[ $read1 == 'left_flank' ]]
					then
						OFF_1=$(($OFF_1 + 1))	
#						echo "adding to OFF_1"
					fi

# Creating a csv file for the results output

				NAME=$(ls $f | cut -f1 -d.)
				echo $NAME','$OFF_1','$OFF_2 >> fimS_OFF_positions.csv
# Printing out results

#				echo 'OFF_1 = ' $OFF_1
#				echo 'OFF_2 = ' $OFF_2
#				echo 'right_flank= ' $right_flank_no
#				echo 'left_flank= ' $left_flank_no
#				echo 'switch_region= ' $switch_region_no
#				echo 'read_count = ' $readcount

				fi

			fi
done <readnames.sorted

# Cleaning up files:

rm readnames
rm completed.reads
rm position.txt


echo "Finished!"
