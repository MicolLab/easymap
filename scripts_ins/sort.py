# This script sorts the information data outputed by the mapping analysis into insertion clusters (assigns each nucleotide to an insertion number). Then, the clusters are filtered to eliminate false positives. Finally, the script determines a candidate region that contains each insertion. 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')#Input 'output_analysis.txt'
parser.add_argument('-b', action="store", dest = 'finput')#Input '34k_genome_2c.fa'
parser.add_argument('-c', action="store", dest = 'output1')#Output1 'output_ordered.csv'
parser.add_argument('-d', action="store", dest = 'output2')#Output2 'sorted_insertions.txt'
parser.add_argument('-m', action="store", dest = 'mode', default = 'pe')
args = parser.parse_args()

#Input 'output_analysis.txt
input = str(args.input)
f1 = open(input, 'r')
lines = f1.readlines()	

#Input '34k_genome_2c.fa'
finput = str(args.finput)
f2 = open(finput, 'r') 
fasta_lines = f2.readlines()

#Output 'output_ordered.csv'
output1 = str(args.output1)
f3 = open(output1, 'w')
f3.write('@' + 'DATA\t' + 'Contig' + '\tins'+ '\t' + 'NT' + '\t' + '  RD' + '\t' + 'Position' + '\n')

###################################################################################################################################################################
#																																								  #
#														Sort input file (output analysis.txt > output ordered.csv)												  #
#																																								  #
###################################################################################################################################################################
#Lists
contigs = []
data = []

#Create a list with all the genome contigs
for i, line in enumerate(fasta_lines):
	if line.startswith('>'): #fasta sequences start with '>'
		sp = line.split(' ')  #because some names have whitespaces and extra info that is not written to sam file
		cont = sp[0].strip()  #strip() is to remove the '\r\n' hidden chars
		cont = cont[1:]       #to remove the first char of string (>)
		if cont not in contigs:
			contigs.append(cont)

#Create a list from the input file
for i, line in enumerate(lines):
	if not line.startswith('@'):
		data.append(line.split())
		
#Sort list and write to file
import operator
sorted_data = sorted(data, key=lambda e: (e[1], int(e[2]))) 

import csv
with open(output1, 'wb') as f3:
    writer = csv.writer(f3)
    writer.writerows(sorted_data)
   
f1.close()
f2.close()
f3.close()


###################################################################################################################################################################
#																																								  #
#														Sort data into insertions (output_ordered.csv > sorted_insertions.txt)									  #
#																																								  #
###################################################################################################################################################################

insertion_id = 1

#Input file 
input = str(args.output1)
f1 = open(input, 'r')
lines = f1.readlines()	

#Output file
output2 = str(args.output2)
f2 = open(output2, 'w')
f2.write('@' + 'DATA\t' + 'Contig' + '\tins'+ '\t' + 'NT' + '\t' + '  RD' + '\t' + 'Direction' + '\n')

insertion_id = 1

#Insertion sorting
if args.mode == 'pe': 
	for i, line in enumerate(lines):
		if not line.startswith('@'):
			sp = line.split(',')
			#if str(sp[0]).strip() == 'PAIRED' and  str(sp[4]).strip() == 'TOTAL':
			p = int(sp[2])
			contig = sp[1].strip('\t')
			try:
				d = abs(int(p2) - int(p))
				if d > 1000 or contig != contig2:
					insertion_id = insertion_id + 1
					f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4])
					p2 = p
					contig2 = contig

				else:
					f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
					p2 = p
					contig2 = contig
			except:
				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
				p2 = p
				contig2 = contig
		
			#else:
			#	f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )

elif args.mode == 'se':
	for i, line in enumerate(lines):
		if not line.startswith('@'): 
			sp = line.split(',')
			#if str(sp[0]).strip() == 'LOCAL_RD' and  str(sp[4]).strip() == 'TOTAL_RD':
			if 'LOCAL' in str(sp[0]).strip(): # and  'TOTAL' in str(sp[4]).strip():
				p = int(sp[2])
				contig = sp[1].strip('\t')
				try:
					d = abs(int(p2) - int(p))
					if d > 500 or contig != contig2:

						insertion_id = insertion_id + 1
						f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4])
						p2 = p
						contig2 = contig

					else:
						f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
						p2 = p
						contig2 = contig
				except:
					f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )
					p2 = p
					contig2 = contig
			else:
				f2.write(sp[0] + '\t' + sp[1] + '\t' + str(insertion_id) + '\t' + sp[2] + '\t' + sp[3] + '\t' + sp[4] )

f1.close()
f2.close()

###################################################################################################################################################################
#																																								  #
#															Filter insertions and rewrite sorted_insertions.txt													  #
#																																								  #
###################################################################################################################################################################

#Input file 
input = str(args.output2)
f1 = open(input, 'r')
lines = f1.readlines()	

insertions_final = list()
insertions_raw = list()

for e in range(1, (insertion_id + 1)):
	insertions_raw.append(e)

if args.mode == 'pe': 
	for insertion in insertions_raw:
		max_RD = 0
		directions = list()
		max_pos = 0
		min_pos = float('inf')
		local = "false"
		span = None

		for l, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split()
				if int(sp[2]) == insertion: 														# and sp[0].strip() == "PAIRED":	
					#1st criterion: insertion must have forward and reverse supporting reads
					read_direction = sp[5].strip("_RD")
					if read_direction not in directions and "TOTAL" not in read_direction:
						directions.append(read_direction)

					#2nd criterion: we calculate the maximum read depth in the data corresponding to the insertion
					if int(sp[4]) > max_RD:
						max_RD = int(sp[4])

					#3rd criterion: insertion data span 
					if int(sp[3]) > max_pos:
						max_pos = int(sp[3])
					if float(sp[3]) < min_pos:
						min_pos = int(sp[3])
					span = max_pos - min_pos

				#4th criterion: there must be at least one local alignment supporting the insertion
				if int(sp[2]) == insertion and "LOCAL" in str(sp[0].strip()):	
					local = "true"

		threshold = 0
		if len(directions) >= 2:
			threshold = threshold + 1
		if max_RD >= 3:
			threshold = threshold + 1
		if span > 250:
			threshold = threshold + 1

		if threshold >= 2:											#threshold >= 2 for default filtering 
			if local == "true":
				insertions_final.append(insertion)

elif args.mode == 'se': 
	for insertion in insertions_raw:
		max_RD = 0
		directions = list()
		max_pos = 0
		min_pos = float('inf')
		span = None

		for l, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split()
				if int(sp[2]) == insertion and sp[0].strip() == "LOCAL_RD":
					#1st criterion: insertion must have forward and reverse supporting reads
					read_direction = sp[5].strip("_RD")
					if read_direction not in directions and "TOTAL" not in read_direction:
						directions.append(read_direction)

					#2nd criterion: we calculate the maximum read depth in the data corresponding to the insertion
					if int(sp[4]) > max_RD:
						max_RD = int(sp[4])

					#3rd criterion: insertion data span 
					if int(sp[3]) > max_pos:
						max_pos = int(sp[3])
					if float(sp[3]) < min_pos:
						min_pos = int(sp[3])
					span = max_pos - min_pos

		threshold = 0
		if len(directions) >= 2:
			threshold = threshold + 1
		if max_RD >= 3:
			threshold = threshold + 1
		if span > 200:
			threshold = threshold + 1

		if threshold >= 2:												#threshold >= 2 for default filtering 
			insertions_final.append(insertion)

f1.close()
f1 = open(input, 'w')
new_id = 1

for fin_ins in insertions_final:
	for i, line in enumerate(lines):
		if not line.startswith('@'):
			sp = line.split()
			if int(sp[2]) == int(fin_ins):
				newline = str( sp[0].strip() + '\t' + sp[1].strip() + '\t' + str(new_id) + '\t' + sp[3].strip() + '\t' + sp[4].strip() + '\t' + sp[5].strip() + '\n' )
				f1.write(newline)

	new_id = new_id + 1
f1.close()


###################################################################################################################################################################
#																																								  #
#															Create candidate region for each insertions 														  #
#																																								  #
###################################################################################################################################################################

if args.mode == 'pe': 
	f2 = open(args.output2, 'r')
	lines = f2.readlines()
	candidate_regions = list() #This list will have the format: list(list(d1, d2, e))

	for ins in range(1, len(insertions_final) + 1):
		d1 = float('inf')
		d2 = 0
		for i, line in enumerate(lines):
			if not line.startswith('@'):
				sp = line.split()
				if sp[0].strip() == 'PAIRED' and int(sp[2].strip()) == int(ins): 

					#The following module creates a candidate region for the insertion
					#Delimiters:
					if sp[5].strip() == 'R': #reverse
						p = int(sp[3])
						if p < d1: 
							d1 = p

					if sp[5].strip() == 'F': #forward
						p2 = int(sp[3])
						if p2 > d2:
							d2 = p2
	
		cr = list()
		cr.append(d1)
		cr.append(d2)
		cr.append(ins)
		
		candidate_regions.append(cr)

	f2 = open(args.output2, 'a')
	for i in candidate_regions: 
		f2.write('@#')
		f2.write(((str(i).strip('[')).strip(']')) + '\n')
