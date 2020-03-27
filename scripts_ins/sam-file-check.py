#This module checks SAM files to see if they are correctly written, it does three comprobations:
#	- Searches for nucleotides in the sequence column of the SAM file
#	- Length of sequence and quality strings comparation
#	- Checks that the SAM file is not empty

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
args = parser.parse_args()
input = str(args.input)

error = 0

#______________________________________FUNCTIONS________________________________________________
#Nucleotide detection:
def nts():
	line.strip()
	sp = line.split('\t')
	if 'A' in sp[9] or 'T' in sp[9] or 'G' in sp[9] or 'C' in sp[9]: #The function looks for either one of the nucleotides in sp[9], which corresponds to the sequence of the read
		pass
	else:
		error = 1
		quit()
	return

#To compare the lenght of the strings sequence and quality:	
def lencheck():
	sp = line.split('\t')
	nts = len(sp[9])
	qs = len(sp[10])
	if qs == nts:
		pass
	else:
		error = 1
		quit()
	return
#_________________________________________________________________________________________________

#Line count
with open(input) as f1:
	linesum = sum(1 for _ in f1) #Counts the lines in the sam file

f1 = open(input, 'r')
lines = f1.readlines()	

#We apply the functions to check the sam file for a maximum of 100 lines 
if linesum == 0:
	error = 1
	quit()
	
if linesum >= 100:
	for i, line in enumerate(lines):
		for i in range(0,100):
			if not line.startswith('@'):
				nts()
				lencheck()
		
if 0<linesum<100:
	for i, line in enumerate(lines):
		if not line.startswith('@'):
			nts()
			lencheck()

print error
f1.close()

