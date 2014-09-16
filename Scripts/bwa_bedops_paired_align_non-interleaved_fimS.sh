#!/bin/bash

# This script is designed to align reads to a reference using BWA, and then count the 
# number of reads overlapping predefined regions. This script generates two outputs: one
# output describes the number of reads physically overlapping the predefined regions, and
# another output describing read-pairs that span across the desired regions.

# The script should be run from inside the directory with the fastq reads: 
# $ bash <script.sh> $REFERENCE $BEDMAP_REFERENCE

# It is very important that the fastq files are formatted correctly for the script to work:
# 	strain_R1.fastq
# 	strain_R2.fastq

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
	echo "processing $(echo $f | cut -f1 -d.)"

		if [[ $f == *_R1.fastq.gz ]]
		then
			gunzip $f
			read1=$(echo $f | cut -f1-2 -d.)
			name1=$(echo $f | cut -f1 -d.)
			bwa aln $REFERENCE $read1 > $name1.sai
		
		elif [[ $f == *_R2.fastq.gz ]]
		then
			gunzip $f
			read2=$(echo $f | cut -f1-2 -d.)
			name2=$(echo $f | cut -f1 -d.)	
			bwa aln $REFERENCE $read2 > $name2.sai
		
		elif [[ $f == *_R1.fastq ]]
		then
			name1=$(echo $f | cut -f1 -d.)
                        bwa aln $REFERENCE $f > $name1.sai
		
		elif [[ $f == *_R2.fastq ]]
		then	
			name2=$(echo $f | cut -f1 -d.)
                        bwa aln $REFERENCE $f > $name2.sai
	fi
done

# Write another loop that will take the .sai files from the previous section, and map them 
# to the reference strain generating .sam and .bam files.

for f in *
do
	if [[ $f == *fastq ]]
	then

# Parse out just the name of the strain (again, in the format $strainname_1.fastq). 
		name=$(echo $f | cut -f1 -d_)
		
# Want to only perform the alignment once - as there are two .sai files, this has the 
# potential to interate through twice. This if statement prevents the script from 
# iterating through more than once on the same strain:
		if [[ ! -e $name.bam ]]
		then
			echo "checking for pairs - " $name

            		if [[ $(echo $f | cut -f1 -d. | cut -f2 -d_) == "R1" ]]
            		then
            			name1=$f
            			echo "first paired read = " $name1
# Make sure that the names of paired files (.sai and .fastq) are the SAME and that the SAME strain files are being aligned to the reference (again using BWA):
                		for g in *
                		do
                			if [[ $(echo $g | cut -f2 -d_) == "R2.fastq" ]] && [[ $(echo $g | cut -f1 -d_) == $name ]]
                    			then
                    				echo $f "and" $g "are a pair - performing alignment"
                        			bwa sampe $REFERENCE $name\_R1.sai $name\_R2.sai $f $g > $name.sam
                        			samtools view -bS -f 0x0002 -F 4 $name.sam > $name.bam
                        			samtools sort $name.bam $name.sorted
                        			samtools index $name.sorted.bam
						rm $f
						rm $g
						rm $name.sam
                    			else
                        			echo $f " and " $g " are not a pair"

                    			fi
                		done

			elif [[ $(echo $f | cut -f1 -d. | cut -f2 -d_) == "R2" ]]
         		then
                    		name2=$f
                    		echo "second paired read = " $name2

                    		for g in *
                    		do
                    			if [[ $(echo $g | cut -f2 -d_) == "R1.fastq" ]] && [[ $(echo $g | cut -f1 -d_) == $name ]]
                        		then
                        			echo $f "and" $g "are a pair"
                            			bwa sampe $REFERENCE $name\_R1.sai $name\_R2.sai $g $f > $name.sam
                            			samtools view -bS -f 0x0002 -F 4 $name.sam > $name.bam
                            			samtools sort $name.bam $name.bam.sorted
                            			samtools index $name.bam.sorted
						rm $f
						rm $g
						rm $name.sam
                        		else
                        			echo $f "and" $g "are not a pair"
                        		fi
                		done
            		fi
		fi
	fi
	

done

# The above command should have generated the .sam and .bam alignment files for all the reads against the reference. 

# Getting rid of the .sai files to clean up the directory.

for f in *
do
        if [[ $f == *.sai ]]
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
		echo "results for "$f" have finished compiling"
	fi
done

# Loop to parse all of the $name.result.bed file contents into a single result.csv file with the strain identifier.

echo "STRAIN,OFF_1,OFF_2,ON_1,ON_2" > fimS_OFF_ON_bed_results.csv

for f in *
do
	if [[ $f == *.result.bed ]]
    	then
		echo $f
		NAME=$(echo $f | cut -f1 -d.) 
		OFF_1=$(head -1 $f | cut -f2 -d\|)
        	OFF_2=$(head -2 $f | tail -1 | cut -f2 -d\|)
        	ON_1=$(tail -2 $f | head -1 | cut -f2 -d\|)
        	ON_2=$(tail -1 $f | cut -f2 -d\|)
                        
    		echo $NAME','$OFF_1','$OFF_2','$ON_1','$ON_2 >> fimS_OFF_ON_bed_results.csv

    	fi
done

echo "finished creating csv file containing bedmaps results"

# This next section counts the number of read-pairs that overlap the 3 regions, namely Left Flank, Switch Region, 
# and Right Flank.

# Assigning of variables to tally the amount of reads in each region.
# Some of these variables also validates the correctness of the pre-generated bam file.
  
# Reads the sorted bam file and obtains all of the read IDs.
# These are then sorted so that the paired reads are together.


for f in *
do
	if [[ $f == *.sorted.bam ]]
	then
		samtools view $f | cut -f1 -d$'\t' > readnames
		sort readnames | uniq >> readnames.sorted
		echo "assigning reads to..."$f

		readcount=$(wc -l readnames.sorted)
	
		OFF_1=0
		OFF_2=0
		ON_1=0
		ON_2=0
		DIFF_REGION=0
		SAME_REGION=0
		
# The next loop reads in the 'readnames.sorted' file which contains all of the read IDs. As there are two
# reads with the same ID, the script will generate twice as many overlaps as the actual amount.
# To account for this, the script checks the current read ID against a list of already analysed read IDs.
# If the read ID exists in the list, it doesn't proceed with the rest of the loop, and adds a '+1' to the readcount tally. 
# If it doesn't exist, it adds the read name to the list and carries out the appropriate analysis. 

		while read name
		do
# For each read that is fed through the while loop, its coordinate (as well as the coordinate of its pair) are 
# obtained from the sorted bam file using 'grep' and cutting the fourth field (which is the coordinate field). 
# This is temporarily saved in a 'positions.txt' file. Variable 'a' is then assigned one of these numbers,
# while variable 'b' is assigned the other. The script then determines what region the read lies in based on its
# starting coordinate, either 'left_flank' (0-999), 'switch_region' (1000-1313) or 'right_flank' (>=1314).
 
			samtools view $f | grep $name | cut -f4 -d$'\t' > position.txt
				

			a=$(head -1 position.txt)

			if [ $a -le 1000 ]
			then
				read1='OFF_left_flank'
			
			elif [ $a -ge 1001 -a $a -le 1313 ]
			then
				read1='OFF_switch_region'
		
			elif [ $a -ge 1314 -a $a -le 2313 ]
			then
				read1='OFF_right_flank'
				
			elif [ $a -ge 2314 -a $a -le 3313 ]
			then
				read1='ON_left_flank'
				
			elif [ $a -ge 3314 -a $a -le 3626 ]
			then
				read1='ON_switch_region'

			elif [ $a -ge 3627 ]
			then
				read1='ON_right_flank'	
			fi

			b=$(tail -1 position.txt)

			if [ $b -le 999 ]
       			then
               			read2='OFF_left_flank'

       			elif [ $b -ge 1000 -a $b -le 1313 ]
       			then
               			read2='OFF_switch_region'

       			elif [ $b -ge 1314 -a $b -le 2313 ]
       			then
               			read2='OFF_right_flank'

			elif [ $b -ge 2314 -a $b -le 3313 ]
                        then
                                read2='ON_left_flank'
                               
			elif [ $b -ge 3314 -a $b -le 3626 ]
                        then
                                read2='ON_switch_region'
                                
                        elif [ $b -ge 3627 ]
                        then
                                read2='ON_right_flank'
			fi

# At this point the script counts the number of reads overlapping different region. 
# Since some of the reads can be in the same region, the first if statement filters only for 
# reads that are in different regions. 
# The next if statements encompass all of the possibilites for the read positions. Depending on
# where the reads lie, a '+1' is added to the tally for either OFF_1 (i.e. read pairs in the left flank and 
# in the switch region), or OFF_2 (i.e. read pairs in the right flank and in the switch region). 

			if [[ $read1 != $read2  ]]
			then
				DIFF_REGION=$(($DIFF_REGION + 1))
				if [[ $read1 == 'OFF_switch_region' ]] && [[ $read2 == 'OFF_right_flank' ]]
				then
					OFF_2=$(($OFF_2 + 1))
					echo $name >> mapped_reads_OFF.txt		
				
				elif [[ $read1 == 'OFF_switch_region' ]] && [[ $read2 == 'OFF_left_flank' ]]
				then
					OFF_1=$(($OFF_1 + 1))
					echo $name >> mapped_reads_OFF.txt
				
				elif [[ $read2 == 'OFF_switch_region' ]] && [[ $read1 == 'OFF_right_flank' ]]
				then
                       			OFF_2=$(($OFF_2 + 1))
					echo $name >> mapped_reads_OFF.txt
				
				elif [[ $read2 == 'OFF_switch_region' ]] && [[ $read1 == 'OFF_left_flank' ]]
				then
					OFF_1=$(($OFF_1 + 1))	
					echo $name >> mapped_reads_OFF.txt
				
				elif [[ $read1 == 'ON_switch_region' ]] && [[ $read2 == 'ON_right_flank' ]]
                                then    
                                        ON_2=$(($ON_2 + 1))
                                        echo $name >> mapped_reads_ON.txt                  
                                
				elif [[ $read1 == 'ON_switch_region' ]] && [[ $read2 == 'ON_left_flank' ]]
                                then    
                                        ON_1=$(($ON_1 + 1))
                                        echo $name >> mapped_reads_ON.txt
                                
				elif [[ $read2 == 'ON_switch_region' ]] && [[ $read1 == 'ON_right_flank' ]]
                                then    
                                        ON_2=$(($ON_2 + 1))
                                        echo $name >> mapped_reads_ON.txt
                                
				elif [[ $read2 == 'ON_switch_region' ]] && [[ $read1 == 'ON_left_flank' ]]
                                then    
                                        ON_1=$(($ON_1 + 1))   
                                        echo $name >> mapped_reads_ON.txt

				fi
			else
				SAME_REGION=$(($SAME_REGION + 1))
			fi
		
		done <readnames.sorted

### Printing out results (optional)

#	echo 'OFF_1 = ' $OFF_1
#       echo 'OFF_2 = ' $OFF_2
#	echo 'ON_1 = ' $ON_1
#	echo 'ON_2 = ' $ON_2
#       echo 'read_count = ' $readcount
#       echo "reads in the same region = " $SAME_REGION
#       echo "reads in different regions = " $DIFF_REGION   


# Creating a csv file for the results output

	rm readnames.sorted

	NAME=$(echo $f | cut -f1 -d.)
	echo $NAME','$OFF_1','$OFF_2','$ON_1','$ON_2 >> fimS_OFF_ON_positions.csv
	
	
	fi
done

# Cleaning up files:

rm readnames
rm readnames.sorted
rm position.txt

# Move all of the files into directories of their own

for f in *
do
	if [[ $f == *.result.bed ]]
	then
		NAME=$(echo $f | cut -f1 -d.)
		mkdir ana_$NAME
		mv $NAME* ana_$NAME
	fi
done


echo "Finished!"
