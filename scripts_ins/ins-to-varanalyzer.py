# This script formats the information gathered from the mapping annalysis so that it can be processed by varanalyzer. For each insertion, the script picks the most probable insertion site by checking the RD values.

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest='input')
parser.add_argument('-b', action="store", dest='output')

args = parser.parse_args()

#input
f1 = open(args.input, 'r')
lines = f1.readlines()

#output
f2 = open(args.output, 'w')
f2.write('#data' + '\t' + 'contig' + '\t' + 'pos' + '\t' + 'ref' + '\t' + 'alt' + '\n')

#First I create a list of the insertions
insertion_list = list()
for i, line in enumerate(lines):
	if not line.startswith('@'):
		sp = line.split()
		if sp[0] == 'LOCAL':
			insertion = sp[2].strip()
			if insertion not in insertion_list:
				insertion_list.append(insertion)
				

#Then I loop through the insertion list, and take the highest RD value:
for ins in insertion_list:
	max_rd = 0
	for line in lines: 
		if not line.startswith('@'):
			sp = line.split()
			if sp[0] == 'LOCAL' and sp[2].strip() == ins and max_rd < int(sp[4]) and 'TOTAL' not in str(sp[5]):
				max_rd = int(sp[4])
				chrom = sp[1].strip()
				pos = sp[3].strip()
				direction = str(sp[5]).strip()
	
	if direction == 'LEFT':
		f2.write('lim' + '\t' + chrom + '\t' + pos + '\t' + '-' + '\t' + '-' + '\n')
	if direction == 'RIGHT':
		pos = str(int(pos)-1)
		f2.write('lim' + '\t' + chrom + '\t' + pos + '\t' + '-' + '\t' + '-' + '\n')
