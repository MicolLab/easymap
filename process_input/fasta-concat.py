
#
# This script receives a list of fasta files, reads the content of each one
# and writes it to a single fasta-formatted file.
#

import argparse, os, shutil, fnmatch

# Process command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gnm', action="store", dest='basename')
parser.add_argument('-out_dir', action="store", dest='output_dir')
args = parser.parse_args()

basename = args.basename
output_dir = args.output_dir

# Create a subdirectory to place the reads. If it already exists, remove it first
if os.path.exists(output_dir):
	shutil.rmtree(output_dir)
os.makedirs(output_dir)

# Function to parse fasta file (based on one of the Biopython IOs)
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

# Function to divide a long string ('data') into chunks of defined length ('batch_size')
def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]

# Create list with all the files in user_data folder
input_files = sorted(os.listdir('./user_data'))

# Create a list with only the files that match the basename provided by the user and end in '.fa'
ref_files = fnmatch.filter(input_files, basename + '*.[Ff][Aa]') # fnmatch filters a list using a string that accepts wildcards

# Create and open output file
output = open(output_dir + '/genome.fa', 'w')

# Get the content of each file and append it to output fasta file
for ref_file in ref_files:
	with open('user_data/' + ref_file) as fp:
		for name_contig, seq_contig in read_fasta(fp):
			split_name_contig = name_contig.split(' ')
			output.write(split_name_contig[0] + '\n')
			
			# Write to file a small chunk of the contig sequence in each line
			for chunk in batch_gen(seq_contig, 80):
				output.write(chunk + '\n')

# Close output file
output.close()
