import subprocess
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-genome', action="store", dest='genome', required=True)
parser.add_argument('-bam', action="store", dest='bam', required=True)
parser.add_argument('-out',action = "store", dest ='out', required=True)
args = parser.parse_args()

contig_source = args.genome
bam = args.bam

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

with open(contig_source) as fp:
	dic = {}
	for name_contig, seq_contig in read_fasta(fp):
		if name_contig not in dic:
			dic[name_contig] = ""
		l = len(seq_contig)	
		dic[name_contig] = l 


t = open(args.out,"w")
subprocess.call(['./samtools1/samtools', 'index', '-b', bam, bam[:-3]+'bai'])		
for contig in dic:
	if dic[contig] > 200000:
		depth_output = subprocess.check_output(['./samtools1/samtools', 'depth', '-a', '-r',contig[1:]+":100000-200000", bam])
		t.write(depth_output)
	
t.close()