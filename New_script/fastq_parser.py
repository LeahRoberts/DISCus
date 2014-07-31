#!/usr/bin/env python

# Fastq parser
# Parses out reads that match those in a txt file

import glob
import os
from Bio import SeqIO
import sys

# Make a list of the desired read names from a txt file:

with open('list.txt', 'r') as f:
        read_ids = [line.strip() for line in f]
        print "%d reads imported" % (len(read_ids))  
#       print len(read_ids)
#	print read_ids

# Takes the fastq read ID and compares it to the read names in the read_ids list
# If the read is found in the read_ids list, it is written out to a "new_reads.fastq" file

	in_file = sys.argv[1]
	name = in_file.split("_")[0]
#       print name
	print "Parsing " + in_file
        reads = SeqIO.parse(in_file, "fastq")
        output = open(name + "_reads_1.fastq", "w")
        for fastq_rec in reads:
#		print fastq_rec.id
		if fastq_rec.id in read_ids:
#			print fastq_rec.id + " exists!"
#			exit()
			SeqIO.write(fastq_rec, output, "fastq")
		
print "Finished parsing " + in_file
os.remove(in_file)

# Next loop does the same thing except with the second list and the second fastq file

with open('list_2.txt', 'r') as f:
        read2_ids = [line.strip() for line in f]
	print "%d reads imported" % (len(read_ids))
#       print read2_ids

        in_file = sys.argv[2]
        name2 = in_file.split("_")[0]
        print "Parsing " + in_file
#	print name2
        reads = SeqIO.parse(in_file, "fastq")
        output = open( name2 + "_reads_2.fastq", "w")
        for fastq_rec in reads:
#               print fastq_rec.id
                if fastq_rec.id in read2_ids:
#                       print fastq_rec.id + " exists!"
#                       exit()
                        SeqIO.write(fastq_rec, output, "fastq")

print "Finished parsing " + in_file
os.remove(in_file)
