#!/user/bin/env python

# CIRQUE_create_reference.py
# 
# Author: Leah Roberts
# Affiliations: Scott Beatson Lab Group, University of Queensland St Lucia
# Date: 29-04-2015

## Python Script to generate pseudo-reference and bed file for use with CIRQUE

#### Part 1: generate pseudo-reference file for mapping

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import os
import sys


# Takes in a fasta file of the genome from which the pseudo-reference is being constructed, 
# the coordinates of the desired invertible region, and also the name of the region (which 
# will determine the output file name and corresponding fasta header):

in_file = sys.argv[1]
start_pos = int(sys.argv[2])
end_pos = int(sys.argv[3])
name = sys.argv[4]

# Because python numbering starts at 0:
start_pos_2 = start_pos - 1
end_pos_2 = end_pos - 1

# Need to extract the invertible region and reverse complement it:

with open("out_file.fa", "w") as f:
	for seq_record in SeqIO.parse( in_file, "fasta"):
		f.write(str(seq_record.seq[start_pos_2:end_pos_2])) # prints sequence from start to end pos
		
with open("reverse_file.fa", "w'") as f:
	q = open("out_file.fa")
	r = Seq(q.read(), generic_dna)
	inv_region = r.reverse_complement() # Reverse complements the invertible region of interest
	f.write(str(inv_region))
			
# Now we should have two files - out_file.fa containing the invertible DNA region, and 
# reverse_file.fa containing the reverse complement of the former

# Extract flanking regions as well:

left_flank = start_pos - 1000
right_flank = end_pos + 1000

with open("left_flank.fa", "w") as f:
	for seq_record in SeqIO.parse( in_file, "fasta"):
		f.write(str(seq_record.seq[left_flank:start_pos_2]))

with open("right_flank.fa", "w") as f:
	for seq_record in SeqIO.parse(in_file, "fasta"):
		f.write(str(seq_record.seq[end_pos_2:right_flank]))

# Write a fasta header to the newfile

with open("new_file.fa", "w") as f:
	f.write(">" + name + "_" + str(start_pos) + "_" + str(end_pos) + "\n")
	
# Combine all of the files together:

os.system('cat left_flank.fa out_file.fa right_flank.fa left_flank.fa reverse_file.fa right_flank.fa >> new_file.fa')

# Get rid of newlines/change file name:

with open(name + ".fa", "w") as f:
	p = open("new_file.fa")
	final_form = p.read()
#	final_form = p.read().rstrip("\n")
	f.write(str(final_form))
	
# Print finishing statement and length of the final pseudo-reference sequence:

print "finished constructing pseudo-reference and bed file"

for seq_record in SeqIO.parse(name + ".fa", "fasta"):
	print "The length of the pseudo reference is: " + str(len(seq_record.seq)) + "bp"



#### Part 2: make a .bed reference file

# The bed reference file needs to contain 4 sets of coordinates for the bordering regions:

re = open("out_file.fa")
inv_length = len(re.read())


field_1 = name + "_" + str(start_pos) + "_" + str(end_pos)

field_2A = 995
field_2B = 1005
field_2C = field_2A + inv_length
field_2D = int(field_2C) + 10
field_2E = int(field_2C) + 2000
field_2F = int(field_2E) + 10
field_2G = field_2E + inv_length
field_2H = int(field_2G) + 10

with open(name + ".bed", "w") as f:
	f.write(str(field_1) + "\t" + str(field_2A) + "\t" + str(field_2B) + "\t" + "." + "\t" + "." + "\t" + "+" + "\t" + "artemis exon" + "\t" + "." + "\t" + "gene_id=exon:" + str(field_2A) + ".." + str(field_2B) + "\n")
	f.write(str(field_1) + "\t" + str(field_2C) + "\t" + str(field_2D) + "\t" + "." + "\t" + "." + "\t" + "+" + "\t" + "artemis exon" + "\t" + "." + "\t" + "gene_id=exon:" + str(field_2C) + ".." + str(field_2D) + "\n")
	f.write(str(field_1) + "\t" + str(field_2E) + "\t" + str(field_2F) + "\t" + "." + "\t" + "." + "\t" + "+" + "\t" + "artemis exon" + "\t" + "." + "\t" + "gene_id=exon:" + str(field_2E) + ".." + str(field_2F) + "\n")
	f.write(str(field_1) + "\t" + str(field_2G) + "\t" + str(field_2H) + "\t" + "." + "\t" + "." + "\t" + "+" + "\t" + "artemis exon" + "\t" + "." + "\t" + "gene_id=exon:" + str(field_2G) + ".." + str(field_2H) + "\n")

# Remove original files:

os.remove("out_file.fa")
os.remove("left_flank.fa")
os.remove("right_flank.fa")
os.remove("reverse_file.fa")
os.remove("new_file.fa")


#### Part 3: create coordinates file (if not using hyxR or fimS inversion regions)

A_left_flank = 1000
A_switch_A = 1001 
A_switch_B = 1000 + int(inv_length)
A_right_A = int(A_switch_B) + 1
A_right_B = int(A_switch_B) + 1000
B_left_A = int(A_right_B) + 1
B_left_B = int(A_right_B) + 1000
B_switch_A = int(B_left_B) + 1
B_switch_B = B_left_B + inv_length
B_right_A = int(B_switch_B) + 1

with open (name + "-coordinates.txt", "w") as f:
	f.write("Region" + "\t" + "Start" + "\t" + "End" + "\n")
	f.write("A_left_flank" + "\t" + "n/a" + "\t" + str(A_left_flank) + "\n")
	f.write("A_switch_region" + "\t" + str(A_switch_A) + "\t" + str(A_switch_B) + "\n")
	f.write("A_right_flank" + "\t" + str(A_right_A) + "\t" + str(A_right_B) + "\n")
	f.write("B_left_flank" + "\t" + str(B_left_A) + "\t" + str(B_left_B) + "\n")
	f.write("B_switch_region" + "\t" + str(B_switch_A) + "\t" + str(B_switch_B) + "\n")
	f.write("B_right_flank" + "\t" + str(B_right_A) + "\t" + "n/a" + "\n")


