#This module filters the reads from a SAM file extracting the localy aligned reads to a fastq file.

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
args = parser.parse_args()

#We will select in the .sam file the locally aligned reads
with open(str(args.input), 'r') as f1: 																#We open the input file and create the output, both as objects (f1, f2)
	with open(str(args.output), 'w') as f2: 
		for line in f1:						
			if not line.startswith('@'): 															#We create a contition to discard the headers in the sam file
				sp = line.split() 																	#Now we split each line into an array 
				if 'S' in sp[5]: 																	#If the read is locally aligned the cigar will contain an "S" 
						f2.write('@'+sp[0] + '\n' + sp[9] + '\n' + '+' + '\n' + sp[10] + '\n' ) 	#The selected reads are written as a fastq file
