#This module will process the information in the .sam file to obtain the read depth per nucleotide aligned. 

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest = 'input')
parser.add_argument('-b', action="store", dest = 'output')
parser.add_argument('-c', action="store", dest = 'finput')
args = parser.parse_args()

#Input file 
input = str(args.input)
f1 = open(input, 'r')
lines = f1.readlines()	

#fasta input
fasta_input = str(args.finput)
fasta_f1 = open(fasta_input, 'r')
fasta_lines = fasta_f1.readlines()

#Output file
output = str(args.output)
f2 = open(output, 'w')
f2.write('@' + 'ANALYSIS\t' + 'Contig' + '\t' + 'NT' + '\t' + 'RD' + '\t' + 'Direction/Position' + '\n') #We write the header to the file

#A list of the values that make up the flag for each read
flags = (2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1, 0)

#We create a list with all the contigs in the refference genome
contigs = []
for i, line in enumerate(fasta_lines):
	if line.startswith('>'): #fasta sequences start with '>'
		sp = line.split(' ')  #because some names have whitespaces and extra info that is not written to sam file
		cont = sp[0].strip()  #strip() is to remove the '\r\n' hidden chars
		cont = cont[1:]       #to remove the first char of string (>)
		if cont not in contigs:
			contigs.append(cont)

#analyze SAM file
for c in contigs:  
	di_f = dict() 										#Dictionary for the forward reads
	di_r = dict() 										#Dictionary for the reverse reads
	di_total = dict() 									#Dictionary for all reads
	for i, line in enumerate(lines): 
		if not line.startswith('@'):
			sp = line.split('\t')
			contig = sp[2] 								#This variable reads the name of the contig in each line..
			if c == contig: 							#.. then compares it to the "c" to group the reads accorting to which contig they belong
														#We define the CIGAR parameter, from which we will extract the aligned nucleotides
				cigar = sp[5].strip()
				if cigar != '*': 						#Discards unaligned reads
				
					#Total dictionary______________________________________________________________________________________________
										
					x = '' 								#The varialbe X will generate a list of each part of the CIGAR code
					l = 0 								#The variable l will contain the sum of the nucleotides aligned in the genome for each read
			
					#The following lines read through the CIGAR separating by \t each letter with its corresponding digits
					for i in cigar: 
						if i == 'M' or i == 'D' or i == 'I' or i == 'N' or i == 'S' or i == 'H' or i == 'P' or i == 'X' : 
							x += str(i) + '\t'
						else:
							x += str(i)
					sp2 = x.split()
			
					#We read through the split X variable and perform the operations necesary to determine l:
					for i in sp2:
						if 'M' in i: 					#For aligned nucleotides we add the number
							num = i.replace('M', '')
							l = int(l) + int(num)
						if 'D' in i: 					#For deletions in the read we add the number of deleted nucleotides
							num = i.replace('D', '')
							l = int(l) + int(num)
						if 'I' in i: 					#For insertions in the read we substract the number of inserted nucleotides
							num = i.replace('I', '')
							l = int(l) - int(num)
					
					#We determine the start and ending of an aligment and write in the dictionary the accumulated read depth of each nucleotide 
					if sp[3] != '0':
						p = int(sp[3]) 					#Initial position
						pf = p + int(l) - 1 			#Final position
						for i in range(p,(pf+1)):  		#i will take the value of each nucleotide in the refference sequence aligned to complete a dictionary
							try: 						#with the read depths of each nucleotide read
								di_total[i] =  1 + di_total[i]
							except KeyError:
								di_total[i] = 1
				

				
					#Determining alignment direction
					f = int(sp[1]) 						#We add the flag to a variable
					l=[] 								#This list will contain the deconstructed flag
					
					while f not in flags: 				#We substract gradually the values from the flags list from the flag of the read adding each substracted value to the list
						for i in flags:
							if i <= f: 
								l.append(i)
								f = f - i 
								break
					else: 
						l.append(f)

					if 16 in l: 						#The flag 16 marks reverse reads, using this information we will determine the direction of each read
						direction = 'R'
					else: #
						direction = 'F'
					
					#Forward dictionary______________________________________________________________________________________________
					if direction == 'F':
						x = '' 
						l = 0
					
						for i in cigar: 
							if i == 'M' or i == 'D' or i == 'I' or i == 'N' or i == 'S' or i == 'H' or i == 'P' or i == 'X' : 
								x += str(i) + '\t'
							else:
								x += str(i)
						sp2 = x.split()
					
						for i in sp2:
							if 'M' in i:
								num = i.replace('M', '')
								l = int(l) + int(num)
							if 'D' in i:
								num = i.replace('D', '')
								l = int(l) + int(num)
							if 'I' in i:
								num = i.replace('I', '')
								l = int(l) - int(num)
							
						if sp[3] != '0':
							p = int(sp[3]) 
							pf = p + int(l) - 1 
							for i in range(p,(pf+1)):
								try: 
									di_f[i] =  1 + di_f[i]
								except KeyError:
									di_f[i] = 1
						
					#Reverse dictionary______________________________________________________________________________________________
					elif direction == 'R':
						x = '' 
						l = 0
					
						for i in cigar: 
							if i == 'M' or i == 'D' or i == 'I' or i == 'N' or i == 'S' or i == 'H' or i == 'P' or i == 'X' : 
								x += str(i) + '\t'
							else:
								x += str(i)
						sp2 = x.split()
					
						for i in sp2:
							if 'M' in i:
								num = i.replace('M', '')
								l = int(l) + int(num)
							if 'D' in i:
								num = i.replace('D', '')
								l = int(l) + int(num)
							if 'I' in i:
								num = i.replace('I', '')
								l = int(l) - int(num)
							
						if sp[3] != '0':
							p = int(sp[3]) 
							pf = p + int(l) - 1 
							for i in range(p,(pf+1)):
								try: 
									di_r[i] =  1 + di_r[i]
								except KeyError:
									di_r[i] = 1
						

	#Finally we sort the dictionaries and write it to an output file	
	for key,value in sorted(di_total.items(), key=lambda i: int(i[0])):
		f2.write('PAIRED\t' +  c + '\t' + str(key) + '\t'+ str(value) + '\tTOTAL\n')

	for key,value in sorted(di_f.items(), key=lambda i: int(i[0])):
		f2.write('PAIRED\t' + c + '\t' + str(key) + '\t'+ str(value) + '\tF\n')
	
	for key,value in sorted(di_r.items(), key=lambda i: int(i[0])):
		f2.write('PAIRED\t' +  c + '\t' + str(key) + '\t'+ str(value) + '\tR\n')