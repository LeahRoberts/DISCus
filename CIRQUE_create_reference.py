#!/user/bin/env python

## Python Script to generate pseudo-reference and bed file for CIRQUE

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import os
import sys


# Needs to take in an fasta file with the coordinated of the invertible region:
# Try to make these imbedded:

in_file = sys.argv[1]
start_pos = int(sys.argv[2])
end_pos = int(sys.argv[3])
name = raw_input('please enter name of reference: ')

# Because python numbering starts at 0:
start_pos_2 = start_pos - 1
end_pos_2 = end_pos - 1

# Need to extract the invertible region and reverse complement it:

with open("out_file.fa", "w") as f:
	for seq_record in SeqIO.parse( in_file, "fasta"):
		f.write(str(seq_record.seq[start_pos_2:end_pos_2])) # prints sequencing from start to end pos
		
with open("reverse_file.fa", "w'") as f:
	q = open("out_file.fa")
	r = Seq(q.read(), generic_dna)
	inv_region = r.reverse_complement()
	f.write(str(inv_region))
			
# Now we should have two files - out_file.fa containing the invertible DNA region, and reverse_file.fa containing the reverse complement of the former

# Extract flanking regions as well:

left_flank = start_pos - 1000
right_flank = end_pos + 1000

with open("left_flank.fa", "w") as f:
	for seq_record in SeqIO.parse( in_file, "fasta"):
		f.write(str(seq_record.seq[left_flank:start_pos_2]))

with open("right_flank.fa", "w") as f:
	for seq_record in SeqIO.parse(in_file, "fasta"):
		f.write(str(seq_record.seq[end_pos_2:right_flank]))


with open("new_file.fa", "w") as f:
	f.write(">" + name + "," + str(start_pos) + ".." + str(end_pos) + "\n")
	
# Combine all of the files together:

os.system('cat left_flank.fa out_file.fa right_flank.fa left_flank.fa reverse_file.fa right_flank.fa >> new_file.fa')

# Remove original files:

os.remove("out_file.fa")
os.remove("left_flank.fa")
os.remove("right_flank.fa")
os.remove("reverse_file.fa")

# Get rid of newlines:

with open(name + ".fa", "w") as f:
	p = open("new_file.fa")
	final_form = p.read()
#	final_form = p.read().rstrip("\n")
	f.write(str(final_form))
	os.remove("new_file.fa")
	
# Print length of the final sequence:

print "finished constructing pseudo-reference"

for seq_record in SeqIO.parse(name + ".fa", "fasta"):
	print "The length of the pseudo reference is: " + str(len(seq_record.seq)) + "bp"

