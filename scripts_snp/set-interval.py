#SDL 01/20
#Snippet to roughly determine the size of the candidate interval accodring to the size of the whole genome
import argparse
#We create the input and output objects
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
args = parser.parse_args()

#A function to parse the genome fasta
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

#Creates a dictionary with name_contig: seq_contig
ch = {}
with open(str(args.input), 'r') as fp: 
	for name_contig, seq_contig in read_fasta(fp):
		ch[name_contig[1:len(name_contig)]] = len(seq_contig)

#Calculates average length of contigs > 4mb
tot_len=0
n_chr=0
for chromosome in ch: 
	if int(ch[chromosome]) > 4000000:
		n_chr += 1
		tot_len += int(ch[chromosome])
try:
	av_chr = tot_len/n_chr
except:
	av_chr = 40000000

#Sets default interval with as 4mb and switches it in three categories accodding to average chromosome size
interval_width = 4000000
if av_chr <= 40000000 and av_chr > 1 : interval_width = 4000000
if av_chr > 40000000 and av_chr < 100000000 : interval_width = 10000000
if av_chr >= 100000000 : interval_width = 20000000

print interval_width
