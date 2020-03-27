import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-in', action="store", dest = 'input')
parser.add_argument('-out', action="store", dest = 'out')
args = parser.parse_args()


#Input 
input = args.input
f1 = open(input, 'r')

#Output
output = args.out
f3 = open(output, 'w')

#Acceptable characters
nts = ['A', 'C', 'T', 'G', 'N']

#Fasta cleaning
fasta = ''
for line in f1: 

	
	if line.startswith('>'):
		header = line.strip().lower()
		f3.write( header + '\n')

	else:
		content = (''.join(i for i in line if not i.isdigit())).strip()
		seq = 'yes'
		for i in content: 
			if i not in nts:
				seq = 'no'
				break
		if seq == 'yes' and len(content.strip()) >= 1:
			f3.write(content)
			f3.write('\n')

f1.close()
f3.close()
