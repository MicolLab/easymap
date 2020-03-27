#
#
# This program simulates a mutagenesis. It takes a DNA sequence and introduces mutations.
# 
# PARAMETERS EXPLANATION:
# 
# -nbr (integer): Total number of mutations. It must be an integer. Absolute value. The ratio between
#       the number of mutations asked to create and the length of the fasta sequence is limited to 0.02.
# 
# -mod [d, e, li]: Mutation mode. 'd': Any base can be mutated and any of the three alternative bases
#       can be the variant, like in genetic drift; 'e': Only CG > AT transitions are created, as caused by
#       ethyl methanesulfonate (EMS); 'li': The sequence of a large insertion provided by the user is inserted
#       at random positions.
# 
# -con: Source of the contig sequence to mutate. It must be in fasta format. It only supports fasta files
#       with a single contig. It does not support sequences that contain bases not stored in
#       'options_basecall_errors'.
# 
# -ins: Source of the insertion sequence. It must be in fasta format and contain a single contig.
#       Only required in mode 'large insertions'. 
# 
# -out: Name of the file (or relative path) where mutated fasta data will be written
# 
# 
# 2016 - David Wilson - dws1985@hotmail.com 2017 - Sergio Andreu (revision) 


import argparse, os, shutil
from random import randint


# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-nbr', action="store", dest='total_nbr_mutations', type=int, required=True)
parser.add_argument('-mod', action="store", dest='mutator_mode',
choices=set(('d','e','li')), required=True) #Choose between d, e, li (d = genetic drift SNP mutations, e = EMS SNP mutations GC>AT, li = long insertion (t-dna, transposon))
parser.add_argument('-con', action="store", dest='contig_source', required=True)
parser.add_argument('-ins', action="store", dest='insertion_source')
parser.add_argument('-out', action="store", dest='output')
parser.add_argument("-causal_mut", action = "store", dest="causal_mut")

args = parser.parse_args()

total_nbr_mutations = int(args.total_nbr_mutations)
mutator_mode = args.mutator_mode   
contig_source = args.contig_source
insertion_source = args.insertion_source
output_folder = args.output

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


# Read insertion fasta file
with open(contig_source) as fp:
	genome_length = []
	for name, seq in read_fasta(fp):
		contig_length = len(seq)
		genome_length.append(contig_length)

if args.causal_mut != None:
	causal_mut_split = args.causal_mut.split("-")

	if len(causal_mut_split) == 1:
		causal_mut_split = args.causal_mut.split(",")
		causal_mut_position = int(causal_mut_split[1]) -1 
		causal_mut_chromosome = int(causal_mut_split[0])
		try:
		 	if genome_length[causal_mut_chromosome-1]<int(causal_mut_split[1]) -1: quit('Quit. Selected mutation is outside the size of the chosen contig.')
		except: quit('Quit. Selected mutation is in a contig not found in fasta.')
	elif len(causal_mut_split) == 2: 
		causal_mut_s = args.causal_mut.split("-") 
		n = 0
		for mutation in causal_mut_s:
			causal_mut_split = mutation.split(",")
			if n == 0:
				causal_mut_position = int(causal_mut_split[1]) -1 
				causal_mut_chromosome = int(causal_mut_split[0])
				try: 
					if genome_length[causal_mut_chromosome-1]<int(causal_mut_split[1]) -1: quit('Quit. Selected mutation one is outside the size of the chosen contig.')
				except: quit('Quit. Selected mutation one is in a contig not found in fasta.')		
			else:
				causal_mut2_position = int(causal_mut_split[1]) -1 
				causal_mut2_chromosome = int(causal_mut_split[0])
				try: 
					if genome_length[causal_mut2_chromosome-1]<int(causal_mut_split[1]) -1: quit('Quit. Selected mutation two is outside the size of the chosen contig.')
				except: quit('Quit. Selected mutation two is in a contig not found in fasta.')

			n +=1
	else:
		quit("Quit. More than two causal motations were given.")


if mutator_mode == 'li' and insertion_source is None:
	quit('Quit. Selected mode is "large insertions" but no insertion sequence source was provided. See program description.')		

# DEPRECATED. IS ONLY USED WHEN PROGRAM WAS USED STANDALONE
# Process 'args.output'
#if '/' in output:
#	slash_pos = output.rfind('/')
#	output_folder = output[:slash_pos]
#	output_file_1 = output[slash_pos+1:]
#	if os.path.exists(output_folder): # If the output directory already exists, remove it first
#		shutil.rmtree(output_folder)
#	os.makedirs(output_folder)
#else:
#	output_folder = './'
#	output_file_1 = output


# Create subdirectories to place the output files. If they already exist, remove them first
output_folder_seq = output_folder + '/mutated_genome/'
output_folder_info = output_folder + '/info/'

if os.path.exists(output_folder):
	shutil.rmtree(output_folder)

os.makedirs(output_folder_seq)
os.makedirs(output_folder_info)





# Function to divide a long string ('data') into chunks of defined length ('batch_size')
def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]


# Read contig fasta file
contigs = []
all_lenght = []
total_lenght = 0
t = 0
with open(contig_source) as fp:
	for name_contig, seq_contig in read_fasta(fp):
		listed = []
		lenght = []
		listed.append(name_contig)
		listed.append(seq_contig)
		if t == 0:
			lenght = [name_contig, 0,len(seq_contig)]
		else:
			lenght = [name_contig,all_lenght[-1][-1], len(seq_contig)+ all_lenght[-1][-1]]
		all_lenght.append(lenght)
		contigs.append(listed)
		total_lenght += len(seq_contig)
		t+= 1

muta_contig= {}
number_chromosome =len(contigs)
for i in range(0,total_nbr_mutations):
	election =randint(0, total_lenght-1)
	for window in all_lenght:
		if election >= window[1] and election < window[2]:
			if window[0] not in muta_contig:
				muta_contig[window[0]]= 0
			muta_contig[window[0]] += 1 

# Check that ratio (nbr of mutations / contig length) is below threshold and stop program if necessary
ratio_ml = float(total_nbr_mutations) / float(total_lenght)
if ratio_ml > 0.02:
	quit('Quit. The ratio (nbr of mutations /  total contig length) is > 0.02')


no_mutated= []
n_chrom= 0
#Calculation of the total contig size
for chromosome in contigs:
	seq_contig = chromosome[1]
	chromosome = chromosome[0]
	try: 
		total_nbr_mutations = muta_contig[chromosome]
	except:
		lis= [chromosome,seq_contig] 
		output2 = open(output_folder_seq + '/mutated_genome.fa', 'a')
		if n_chrom ==0:
			output2.write(lis[0])
		else:
			output2.write("\n"+lis[0])
		for chunk in batch_gen(lis[1].upper(), 80):
			output2.write('\n' + chunk)
		output2.close()	
		n_chrom += 1		
		continue
	# Calculate contig length
	contig_length = len(seq_contig)
	# Convert contig sequence to uppercase
	seq_contig_uppercase = seq_contig.upper()

	
	# Randomly choose mutation positions, add them to a list, and finally sort the list
	all_mut_pos = list() # Will store the list of mutated positions

	iter = 0
	while iter < total_nbr_mutations:
		mut_pos = randint(1, contig_length - 2) # I deliberately exclude first and last nucleotides
		if mut_pos not in all_mut_pos:
			if mutator_mode == 'd' or mutator_mode == 'li': # If mode is drift, any base can be mutated
				iter +=1
				all_mut_pos.append(mut_pos)
			if mutator_mode == 'e': # If mode is EMS, ony G and A can be mutated
				wt_base = seq_contig_uppercase[mut_pos:mut_pos+1]
				if wt_base == 'G' or wt_base == 'C':
					iter +=1
					all_mut_pos.append(mut_pos)
	if args.causal_mut != None:				
		if n_chrom+1 == causal_mut_chromosome and causal_mut_position not in all_mut_pos:
			if mutator_mode == "d" or mutator_mode == "li":
				all_mut_pos.append(causal_mut_position)
			if mutator_mode =="e":
				iter = 0
				while True:
					wt_base = seq_contig_uppercase[causal_mut_position+iter:causal_mut_position+iter+1]
					iter +=1
					if wt_base == 'G' or wt_base == 'C':
						all_mut_pos.append(causal_mut_position+iter-1)
						break
		if len(args.causal_mut.split(",")) > 2:
			if n_chrom+1 == causal_mut2_chromosome and causal_mut2_position not in all_mut_pos:
				if mutator_mode == "d" or mutator_mode == "li":
					all_mut_pos.append(causal_mut2_position)
				if mutator_mode =="e":
					iter = 0
					while True:
						wt_base = seq_contig_uppercase[causal_mut2_position+iter:causal_mut2_position+iter+1]
						iter +=1
						if wt_base == 'G' or wt_base == 'C':
							all_mut_pos.append(causal_mut2_position+iter-1)
							break


	all_mut_pos.sort()


	# If mode is SNPs (drift or EMS), use list 'all_mut_pos' and determine the mutant allele at each position
	if mutator_mode == 'd' or mutator_mode == 'e':
		all_mut_info = list() # Will store mutated positions, wt base, and mut base

		# Iterate over list, determine the wt base for each position in list, and randomly choose a mt base
		for i in all_mut_pos:
			wt_base = seq_contig_uppercase[i:i+1]

			if mutator_mode == 'd': # If mode is drift, the mutation can be any base
				options_drift = {'A': ['T', 'C', 'G'], 'T': ['C', 'G', 'A'], 'C': ['G', 'A', 'T'], 'G': ['A', 'T', 'C']}
				try:
					possible_mt_bases = options_drift[wt_base] # Lookup the alterantive bases of "wt_base" in "options_drift"
					mt_base = possible_mt_bases[randint(0, 2)] # Randomly pick between the possible mut bases
				except KeyError:
					mt_base = 'N'

			if mutator_mode == 'e': # If mode is ems, only GC>AT possible
				options_ems = {'G': 'A', 'C': 'T'}
				try:
					mt_base = options_ems[wt_base]
				except KeyError:
					mt_base = 'N'

			mut_info = (i, wt_base, mt_base) 
			all_mut_info.append(mut_info)

		# Iterate over all_mut_info list, mutate the wt sequence, and write the mutant version to file	
		mutant_seq_contig = list(seq_contig_uppercase) # I use a list because python strings do not support item assignment
			
		for mut_info in all_mut_info:
			# Write to info file
			mutant_seq_contig[mut_info[0]] = mut_info[2] # For each mutation in list "all_mut_info", appropriately substitute in list "mutant_seq"


			

		output2 = open(output_folder_seq + '/mutated_genome.fa', 'a')
		if n_chrom == 0:
			output2.write(chromosome + '__mutated')
		else:
			output2.write("\n"+chromosome + '__mutated')

		# Write to file a small chunk of the contig sequence in each line
		mutant_seq_contig = ''.join(mutant_seq_contig)
		for chunk in batch_gen(mutant_seq_contig, 80):
			output2.write('\n' + chunk)

		output2.close()
		
		
	if mutator_mode == 'li':
	# How insertion positions are handled:
	#  ATCG     .ATCG     A.TCG     AT.CG     ATC.G     ATCG.
	#  01234    0          1          2          3          4

		mutant_seq_contig = ''

		# Read insertion fasta file
		with open(insertion_source) as fp:
			for name_ins, seq_ins in read_fasta(fp):
				#print (name_ins, seq_ins)
				pass

		# Create a new list that contains the insertion points and the beggining and end coordinates 
		all_mut_pos_plus_ends = list(all_mut_pos)

		all_mut_pos_plus_ends.insert(0, 0) # Add the value "0" at position 0 of the list
		all_mut_pos_plus_ends.append(len(seq_contig)) # Add the last base of the contig to the end of the list 

		# Iterate over list "all_mut_pos_plus_ends", take the coord values and substring the contig sequence, and build the mutant sequence
		iter_ins = 0
		while iter_ins < len(all_mut_pos_plus_ends):
			try:
				fragment = seq_contig_uppercase[all_mut_pos_plus_ends[iter_ins]:all_mut_pos_plus_ends[iter_ins + 1]]
				#print fragment_len
				mutant_seq_contig += fragment
				if iter_ins < len(all_mut_pos_plus_ends) - 2:
					mutant_seq_contig += seq_ins
			except IndexError:
				pass
			
			iter_ins +=1

		output2 = open(output_folder_seq + '/mutated_genome.fa', 'a')
		if n_chrom ==0:
			output2.write(chromosome +'__mutated')
		else:
			output2.write("\n"+chromosome + '__mutated')
		# Write to file a small chunk of the contig sequence in each line
		for chunk in batch_gen(mutant_seq_contig, 80):
				output2.write('\n' + chunk)
				
		output2.close()	


	# Iterate over list "all_mut_info" and create file with info of mutations
	output1 = open(output_folder_info + '/info_all_mutations.txt', 'a')

	# If mode is SNPs, write header with 3 columns and use list "all_mut_info".
	if mutator_mode == 'd' or mutator_mode == 'e':
		if n_chrom == 0:
			output1.write("chromosome"+"\t"+'position\twt_base\tmt_base\n')
		for mut_info in all_mut_info:
			pos = int(mut_info[0]) + 1
			output1.write(chromosome[1:len(chromosome)]+"\t"+str(pos) + '\t' + str(mut_info[1]) + '\t' + str(mut_info[2]) + '\n')

	# If mode is large insertions, write name of inserted sequence, single column header, and just use list "all_mut_pos".
	if mutator_mode == 'li':
		if n_chrom == 0:
			output1.write('name of inserted sequence: ' + name_ins[1:] + '\n\npositions of insertions\n')
		for mut_pos in all_mut_pos:
			output1.write(chromosome[1:len(chromosome)]+ "\t" + str(mut_pos) + '\n')
	output1.close()

	n_chrom +=1

if len(no_mutated)!= 0:
	for chrom in no_mutated: 
			output2 = open(output_folder_seq + '/mutated_genome.fa', 'a')
			if n_chrom ==0:
				output2.write(chrom[0])
			else:
				output2.write("\n"+chrom[0])
			for chunk in batch_gen(chrom[1].upper(), 80):
				output2.write('\n' + chunk)
	output2.close()			

