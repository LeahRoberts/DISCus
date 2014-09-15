#!/bin/bash


# This script allows for the assigning of paired-end reads to different regions of a 
# reference sequence for the purpose of counting the number of paired reads that span multiple regions. 

# The script relies on the appropriate construction of a bam file (sorted) which filters out reads
# that are not properly paired or mapped. Samtools is also required to run this script. 

# The script can be run using the command "bash <script> <bam_file>"

# The script is currently set specifically for the Type 1 Fimbriae fimS switch in EC958. 
# However, with small manipulations it could be applied to other similar projects.  


# Assigning of variables to tally the amount of reads in each region.
# Some of these variables also validates the correctness of the pre-generated bam file.
  
readcount=0

OFF_1=0
OFF_2=0

right_flank_no=0
left_flank_no=0
switch_region_no=0

BAM_FILE=$1

# Reads the sorted bam file and obtains all of the read IDs.
# These are then sorted so that the paired reads are together.

samtools view $BAM_FILE | cut -f1 -d$'\t' > readnames
sort readnames > readnames.sorted

echo "assigning reads..."

echo "" > completed.reads

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
#		cat position.txt

		a=$(head -1 position.txt)
#		echo $a

		if [ $a -le 999 ]
		then
			read1='left_flank'
			left_flank_no=$(($left_flank_no + 1))
#			echo "read 1 is on the left"
		
		elif [ $a -ge 1000 -a $a -le 1313 ]
		then
			read1='switch_region'
			switch_region_no=$(($switch_region_no + 1))
#			echo "read 1 is in the switch region"
			
		elif [ $a -ge 1314 ]
		then
			read1='right_flank'
			right_flank_no=$(($right_flank_no + 1))
#			echo "read 1 is on the right"
		fi

		b=$(tail -1 position.txt)
#		echo $b

		if [ $b -le 999 ]
        	then
                	read2='left_flank'
			left_flank_no=$(($left_flank_no + 1))
#			echo "read 2 is on the left"        

        	elif [ $b -ge 1000 -a $b -le 1313 ]
        	then
                	read2='switch_region'
			switch_region_no=$(($switch_region_no + 1))
#			echo "read 2 is in the switch region"        

        	elif [ $b -ge 1314 ]
        	then
                	read2='right_flank'
        		right_flank_no=$(($right_flank_no + 1))
#			echo "read 2 is on the right"
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
#				echo "adding to OFF_2"		
			elif [[ $read1 == 'switch_region' ]] && [[ $read2 == 'left_flank' ]]
			then
				OFF_1=$(($OFF_1 + 1))
#				echo "adding to OFF_1"
			elif [[ $read2 == 'switch_region' ]] && [[ $read1 == 'right_flank' ]]
			then
                        	OFF_2=$(($OFF_2 + 1))
#				echo "adding to OFF_2"
			elif [[ $read2 == 'switch_region' ]] && [[ $read1 == 'left_flank' ]]
			then
				OFF_1=$(($OFF_1 + 1))	
#				echo "adding to OFF_1"
			fi
		fi
	fi
done <readnames.sorted

# Cleaning up files:

rm readnames
rm completed.reads
rm position.txt

# Printing out results

echo 'OFF_1 = ' $OFF_1
echo 'OFF_2 = ' $OFF_2
echo 'right_flank= ' $right_flank_no
echo 'left_flank= ' $left_flank_no
echo 'switch_region= ' $switch_region_no
echo 'read_count = ' $readcount

# Creating a csv file for result output

echo $OFF_1','$OFF_2 >> fimS_OFF_positions.csv
