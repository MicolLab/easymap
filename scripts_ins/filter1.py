#This module filters the reads from a SAM file and extracts the unaligned reads with aligned pairs in a fastq file.

import argparse

#We create the input and output objects
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
args = parser.parse_args()

#Now we select in the .sam file the unpaired reads whose mates are paired 
with open(str(args.input), 'r') as f1: 					#We open the input file and create the output as objects (f1, f2)
	with open(str(args.output), 'w') as f2:
		for line in f1:								#To read through the lines of the file
			if not line.startswith('@'):			#We create a contition to eliminate the headers in the sam file
				sp = line.split() 					#Now we split each line into an array
				if sp[2] != '*' and sp[5] == '*': 	#Sp[2] is not an asterisk because the read takes the contig name of its mate when they are not aligned. The CIGAR (sp[5]) reveals if the read hasnt been aligned with an asterisk.
						f2.write('@'+sp[0] + '\n' + sp[9] + '\n' + '+' + '\n' + sp[10] + '\n' ) #The selected reads are written as a fastq file