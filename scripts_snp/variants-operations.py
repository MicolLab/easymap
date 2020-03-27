# This script performs different operations with polymorphism data contained in VA files
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input1')
parser.add_argument('-b', action="store", dest = 'input2')
parser.add_argument('-c', action="store", dest = 'output')
parser.add_argument('-mode', action="store", dest = 'mode', default='N', choices=['I','A','B', 'U', 'N']) #['I','A','B', 'U', 'N']) (Intersection, vcf1 - vcf2, vcf2 - vcf1, union, None)
parser.add_argument('-primary', action="store", dest = 'primary', default=1, choices=['1','2']) #1: info from input -a ; 2: info from input -b 
args = parser.parse_args()

mode = args.mode
primary = int(args.primary)

#Datasets A (VCF1) and B (VCF2)
#Input 1  
input1 = args.input1
f1 = open(input1, 'r')
lines1 = f1.readlines()	

#Input 2 
input2 = args.input2
f2 = open(input2, 'r')
lines2 = f2.readlines()	

#Output
output = args.output
f3 = open(output, 'w')

c = '-'

if mode == 'N':
	if primary == 1:
		for i, line in enumerate(lines1): 
			f3.write(line)

	elif primary == 2: 
		for i, line in enumerate(lines2):
			f3.write(line)

if mode == 'I':
	if primary == 1:
		vcf1 = list() #List of vcf1 positions
		vcf2 = list() #List of vcf2 positions

		for i, line in enumerate(lines1): #We add the positions of f1 to a list 
			if not line.startswith('#'): 
				sp = line.split('\t')
				vcf1.append(sp[0].strip() + '-' + sp[1].strip())

		for i, line in enumerate(lines2): #We add the positions of f2 to a list
			if not line.startswith('#'):
				sp = line.split('\t')
				vcf2.append(sp[0].strip() + '-' + sp[1].strip())
	
		intersection =list(set(vcf1).intersection(vcf2))

		for i, line in enumerate(lines1):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (sp[0].strip()+'-'+sp[1].strip()) in intersection:
					f3.write(line)

	elif primary == 2:
		vcf1 = list() #List of vcf1 positions
		vcf2 = list() #List of vcf2 positions

		for i, line in enumerate(lines2): #We add the positions of f1 to a list 
			if not line.startswith('#'): 
				sp = line.split('\t')
				vcf1.append(sp[0].strip() + '-' + sp[1].strip())

		for i, line in enumerate(lines1): #We add the positions of f2 to a list
			if not line.startswith('#'):
				sp = line.split('\t')
				vcf2.append(sp[0].strip() + '-' + sp[1].strip())
	
		intersection =list(set(vcf1).intersection(vcf2))

		for i, line in enumerate(lines2):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (sp[0].strip()+'-'+sp[1].strip()) in intersection:
					f3.write(line)

elif mode == 'U':
	if primary == 1:
		vcf1 = list() #List of vcf1 positions
		vcf2 = list() #List of vcf2 positions

		for i, line in enumerate(lines1): #We add the positions of f1 to a list 
			if not line.startswith('#'): 
				sp = line.split('\t')
				vcf1.append(sp[0].strip() + '-' + sp[1].strip())
	
		for i, line in enumerate(lines2): #We add the positions of f2 to a list
			if not line.startswith('#'):
				sp = line.split('\t')
				vcf2.append(sp[0].strip() + '-' + sp[1].strip())
			
		union = list(set(vcf1) | set(vcf2))

		for i, line in enumerate(lines1):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (sp[0].strip()+'-'+sp[1].strip()) in union:
					f3.write(line)

		input3 = 'op_output.vcf'
		f3 = open(input3, 'r')
		lines3 = f3.readlines()
		vcf3 = list()

		for i, line in enumerate(lines3):
			if not line.startswith('#'):
				sp = line.split('\t')
				vcf3.append(sp[0].strip() + '-' + sp[1].strip())	

		R = list(set(union).difference(vcf3)) 

		input3 = 'op_output.vcf'
		f3 = open(input3, 'a')	

		for i, line in enumerate(lines2):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (sp[0].strip()+'-'+sp[1].strip()) in R:
					f3.write(line)

	elif primary == 2:			
		vcf1 = list() #List of vcf1 positions
		vcf2 = list() #List of vcf2 positions

		for i, line in enumerate(lines2): #We add the positions of f1 to a list 
			if not line.startswith('#'): 
				sp = line.split('\t')
				vcf1.append(sp[0].strip() + '-' + sp[1].strip())
	
		for i, line in enumerate(lines1): #We add the positions of f2 to a list
			if not line.startswith('#'):
				sp = line.split('\t')
				vcf2.append(sp[0].strip() + '-' + sp[1].strip())
			
		union = list(set(vcf1) | set(vcf2))

		for i, line in enumerate(lines2):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (sp[0].strip()+'-'+sp[1].strip()) in union:
					f3.write(line)

		input3 = 'op_output.vcf'
		f3 = open(input3, 'r')
		lines3 = f3.readlines()
		vcf3 = list()

		for i, line in enumerate(lines3):
			if not line.startswith('#'):
				sp = line.split('\t')
				vcf3.append(sp[0].strip() + '-' + sp[1].strip())	

		R = list(set(union).difference(vcf3))

		input3 = 'op_output.vcf'
		f3 = open(input3, 'a')	

		for i, line in enumerate(lines1):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (sp[0].strip()+'-'+sp[1].strip()) in R:
					f3.write(line)

elif mode == 'A':
	vcf1 = list() 
	vcf2 = list() 
	#We add the positions of f1 to a list 
	for i, line in enumerate(lines1):
		if not line.startswith('#'): 
			sp = line.split('\t')
			vcf1.append(sp[0].strip() + '-' + sp[1].strip())
			
	#We add the positions of f2 to a list
	for i, line in enumerate(lines2):
		if not line.startswith('#'):
			sp = line.split('\t')
			vcf2.append(sp[0].strip() + '-' + sp[1].strip())
		
	A = list(set(vcf1).difference(vcf2))

	for i, line in enumerate(lines1):
		if not line.startswith('#'):
			sp = line.split('\t')
			if (sp[0].strip()+'-'+sp[1].strip()) in A:
				f3.write(line)
			
elif mode == 'B':			
	vcf1 = list() 
	vcf2 = list() 
	#We add the positions of f1 to a list 
	for i, line in enumerate(lines1):
		if not line.startswith('#'): 
			sp = line.split('\t')
			vcf1.append(sp[0].strip() + '-' + sp[1].strip())
			
	#We add the positions of f2 to a list
	for i, line in enumerate(lines2):
		if not line.startswith('#'):
			sp = line.split('\t')
			vcf2.append(sp[0].strip() + '-' + sp[1].strip())
		
	B = list(set(vcf2).difference(vcf1))

	for i, line in enumerate(lines2):
		if not line.startswith('#'):
			sp = line.split('\t')
			if (sp[0].strip()+'-'+sp[1].strip()) in B:
				f3.write(line)