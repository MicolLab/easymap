#Aim:  Using the variants file generated from the parental control (in case you are dealing with an outcross in a reference background and you have sequenced the polimorfic parental as a control.) and
#      using the gnm_ref file, the variants located in the control will be replaced in the reference genome.

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-var', '-v', action = 'store', dest = 'v_file')
parser.add_argument('-gnm_ref', '-r', action = 'store', dest = 'gnm_ref')
parser.add_argument('-out',action='store',dest = 'output')

args = parser.parse_args()

v_file = args.v_file
gnm_ref= args.gnm_ref
out = args.output

#data in var file should be translated into chromosome: [position,change] format. A list of list will be used.
variants_list = []
with open(v_file,"r") as variants:
	for line in variants:
		if line.startswith("#"):
			pass
		line.rstrip()
		l = []
		spl = line.split("\t")
		l.extend([spl[0],spl[1],spl[3]]) # 0 = chromosme, 1 = position, 3= Alt Base
		variants_list.append(l)

# Function to divide a long string ('data') into chunks of defined length ('batch_size')
def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]

#data of the genome reference will be stored in a list of contigs using the Biopython function

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

#Iteration through contigs. Sequences will be separed in a list, the elements in the list will be replaced
with open(gnm_ref,"r") as fp:
	for contig in read_fasta(fp): #contig[0] = name, contig[1] = sequence
		contig_mutated = list(contig[1]) #String to a list in order to change positions 		
		for variant in variants_list:
			if variant[0] == contig[0][1:]: #If the contig is the same (deleting >)
				contig_mutated[int(variant[1])-1] = variant[2].rstrip() #Change the position for the variant one
		contig_mutated = "".join(contig_mutated) #recreate the string
		#Write the brand new fasta
		with open(out,"a") as output_file:
			output_file.write(contig[0]+"\n")
			for chunk in batch_gen(contig_mutated,80):
				output_file.write(chunk+"\n")