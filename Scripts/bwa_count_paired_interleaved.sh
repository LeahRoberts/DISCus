#!/bin/bash

# This script is designed to align reads to a reference using BWA, and then count the 
# number of reads overlapping predefined regions. This script generates one output: 
# counts the read-pairs that span across the desired regions.

# The script should be run from inside the directory with the fastq reads: 
# $ bash <script.sh> $REFERENCE 

# The REFERENCE is what the reads will be mapped to. It should be in the format ../reference.fa

# Make sure that the .fa and .bed references are compatible (i.e. same sequence name - refer to the README). 
# Also check that the .bed reference contains exon features (not just CDS).

REFERENCE=$1

# First need to index the reference fasta file to be used in the mapping:
echo "indexing " $REFERENCE
bwa index $REFERENCE

# Write another loop that will take the .sai files from the previous section, and map them 
# to the reference strain generating .sam and .bam files.

for f in *
do
	if [[ $f == *fastq.gz ]]
	then

# Parse out just the name of the strain (again, in the format $strainname_1.fastq). 
		name=$(echo $f | cut -f1 -d.)
		
# Want to only perform the alignment once - as there are two .sai files, this has the 
# potential to interate through twice. This if statement prevents the script from 
# iterating through more than once on the same strain:
		if [[ ! -e $name.bam ]]
		then
			echo "performing alignment on " $name

# The '-f 0x0002' flag will filter for only properly mapped reads, while the '-F 4' flag
# will get rid of reads that are not mapped to the reference.

    	        	/home/leah/bin/bwa/bwa mem $REFERENCE -p $f > $name.sam
                    	samtools view -bS -F 4 $name.sam > $name.bam
                    	samtools sort $name.bam $name.sorted
                     	samtools index $name.sorted.bam         
			rm $name.sam
			rm $name.bam
		fi
	fi
done

# The above command should have generated the .sam and .bam alignment files for all the reads against the reference. 

# This next section counts the number of read-pairs that overlap the 3 regions, namely Left Flank, Switch Region, 
# and Right Flank.

# Assigning of variables to tally the amount of reads in each region.
# Some of these variables also validates the correctness of the pre-generated bam file.
  
# Reads the sorted bam file and obtains all of the read IDs.
# These are then sorted so that the paired reads are together.

# Create a file for the output (with headers):

echo "STRAIN,OFF_1,OFF_2,ON_1,ON_2,Ambiguous" > fimS_OFF_ON_positions.csv

for f in *
do
	echo "" > completed.reads
	
	if [[ $f == *.sorted.bam ]]
	then
		samtools view $f | cut -f1 -d$'\t' > readnames
		sort readnames > readnames.sorted
		echo "assigning reads to..."$f

		readcount=0
				
		OFF_1=0
		OFF_2=0
		ON_1=0
		ON_2=0
		Ambiguous=0
		SAME_REGION=0	
#		right_flank=0
#		left_flank=0
#		switch_region=0

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
				samtools view $f | grep $name | cut -f4 -d$'\t' > position.txt
#				cat position.txt
				

				a=$(head -1 position.txt)
#				echo $a

				if [ $a -le 1000 ]
				then
					read1='OFF_left_flank'
				#	left_flank_no=$(($left_flank_no + 1))
#					echo "read 1 is on the left"
		
				elif [ $a -ge 1001 -a $a -le 1313 ]
				then
					read1='OFF_switch_region'
				#	switch_region_no=$(($switch_region_no + 1))
#					echo "read 1 is in the switch region"
			
				elif [ $a -ge 1314 -a $a -le 2313 ]
				then
					read1='OFF_right_flank'
				#	right_flank_no=$(($right_flank_no + 1))
#					echo "read 1 is on the right"
				
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
#				echo $b

				if [ $b -le 999 ]
        			then
                			read2='OFF_left_flank'
				#	left_flank_no=$(($left_flank_no + 1))
#					echo "read 2 is on the left"        

        			elif [ $b -ge 1000 -a $b -le 1313 ]
        			then
                			read2='OFF_switch_region'
				#	switch_region_no=$(($switch_region_no + 1))
#					echo "read 2 is in the switch region"        

        			elif [ $b -ge 1314 -a $b -le 2313 ]
        			then
                			read2='OFF_right_flank'
        			#	right_flank_no=$(($right_flank_no + 1))
#					echo "read 2 is on the right"

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
					if [[ $read1 == 'OFF_switch_region' ]] && [[ $read2 == 'OFF_right_flank' ]]
					then
						OFF_2=$(($OFF_2 + 1))
#						echo "adding to OFF_2"
						echo $name >> mapped_reads_OFF.txt		
					elif [[ $read1 == 'OFF_switch_region' ]] && [[ $read2 == 'OFF_left_flank' ]]
					then
						OFF_1=$(($OFF_1 + 1))
#						echo "adding to OFF_1"
						echo $name >> mapped_reads_OFF.txt
					elif [[ $read2 == 'OFF_switch_region' ]] && [[ $read1 == 'OFF_right_flank' ]]
					then
                        			OFF_2=$(($OFF_2 + 1))
#						echo "adding to OFF_2"
						echo $name >> mapped_reads_OFF.txt
					elif [[ $read2 == 'OFF_switch_region' ]] && [[ $read1 == 'OFF_left_flank' ]]
					then
						OFF_1=$(($OFF_1 + 1))	
#						echo "adding to OFF_1"
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
					elif [[ $read2 == 'ON*' ]] && [[ $read1 == 'OFF*' ]]
					then
						Ambiguous=$(($Ambiguous + 1))
					elif [[ $read2 == 'OFF*' ]] && [[ $read1 == 'ON*' ]]
					then
						Ambiguous=$(($Ambiguous + 1))
					
					fi
# Printing out results
				else
					SAME_REGION=$(($SAME_REGION + 1))
				fi

			fi

		done <readnames.sorted

echo 'OFF_1 = ' $OFF_1
echo 'OFF_2 = ' $OFF_2
echo 'ON_1 = ' $ON_1
echo 'ON_2 = ' $ON_2
echo 'Ambiguous= ' $Ambiguous
echo 'read_count = ' $readcount
echo 'Reads in same region = ' $SAME_REGION

# Creating a csv file for the results output

	NAME=$(ls $f | cut -f1 -d.)
	echo $NAME','$OFF_1','$OFF_2','$ON_1','$ON_2','$Ambiguous >> fimS_OFF_ON_positions.csv
	
	echo "completed.reads removed"	
#	rm completed.reads	
	fi
done

# Cleaning up files:

#rm readnames
#rm readnames.sorted
rm position.txt
#rm completed.reads
#rm *.fastq

# Move all of the files into directories of their own

for f in *
do
	if [[ $f == *.sorted.bam ]]
	then
		NAME=$(echo $f | cut -f1 -d.)
		mkdir ana_$NAME
		mv $NAME* ana_$NAME
	fi
done


echo "Finished!"

