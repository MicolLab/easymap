
'''
This program reads a fasta file with one or more contigs and calculates
the total length of the squence(s) 

'''

import argparse

# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gnm', action='store', dest='genome', required=True)
args = parser.parse_args()
genome = args.genome

# Function to parse fasta file
def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith('>'):
			if name: yield (name, ''.join(seq))
			name, seq = line, []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

genome_length = 0

# Read fasta file
with open(genome) as fp:
	for name, seq in read_fasta(fp):
		contig_length = len(seq)
		genome_length += contig_length

print genome_length