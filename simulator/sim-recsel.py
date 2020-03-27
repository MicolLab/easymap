# 
# 
# This program simulates the F2 offspring of a cross between two different isogenic lines and the
# selection of individuals with a particular trait. The original purpose was to create mapping
# populations. It takes two chromosomes from Arabidopsis thaliana and creates recombinant versions.
# The recombination frequency distribution is based in Salome et al 2012, Heredity (2012) 108,
# 447-455; doi:10.1038/hdy.2011.95. (The program can be modified to be used with other species
# mainly by modifying the list(s) 'chr_xo_freq'). The positions of the crossovers are random.
# The program then selects the recombinant chromosomes based on whether they carry or not a given
# mutation selected by the user (the phenotype-causing mutation).
# 
# PARAMETERS EXPLANATION:
# 
# -outdir: Name of directory relative to current location where to place the recombinant chromosomes
# 		created.
# 
# -rec_freq_distr: List of ';'-separated pairs (Nbr XO events, freq (%)). Example: "0,15; 1,20; 2,18; 3,10; ..."  ; --> -   sin comillas
#			Sum of freq values must be 100%. This distribution is determined experimentally (e.g. in
#			Arabidopsis: Salome et al 2012, Heredity (2012) 108, 447-455; doi:10.1038/hdy.2011.95.) 
# 
# -parmut: Parental A. In modes 'r', 'd', and 'di', this parental must be the mutant parental. It
# 		must be a fasta-formatted file with a single contig. The length of the contig must be equal to the
# 		length of -parpol.
# 
# -parpol: Parental B. Polymorphic to parental A. In modes 'r', 'd', and 'di', this parental must
# 		be the non-mutant parental. It must be a fasta-formatted file with a single contig. The length of
# 		the contig must be equal to the length of -parmut.
# 
# -mutapos (integer): Chromosome and position of the causal mutation in -parmut.
# 
# -mutbpos (integer): Position of the causal mutation in -parpol. Only required if mode is 'dr'.
# 
# -smod [r, d, di, dr]: Selection mode. 'r': Recessive mutation and selection of the mutant phenotype
# 		(recessive phenotypic class); 'd': Dominant mutation and selection of the mutant phenotype (dominant
# 		phenotypic class); 'di': Dominant mutation and selection of the wild type phenotype (recessive
# 		phenotypic class); 'wt': Recessive mutation and selection of the wild type phenotype (dominant
#		phenotypic class); 'dr': Two recessive mutations and selection of the double mutant phenotype.
#						'mc': mutant cross (mutation 1 and two come from different parentals)
# 		In the first three modes, the causal mutation is provided in parental -parmut. In the last mode,
# 		each causal mutation is provided in a different parental.
# 
# -nrec (integer): Number of recombinant chromosomes to create. This number corresponds to the number
# 		of chromosomes after the selection has been performed.



import argparse, os, shutil
from random import randint


# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-outdir', action="store", dest='out_dir', required=True)
parser.add_argument('-rec_freq_distr', action ="store", dest = "rec_freq_distr", required = True)
parser.add_argument('-parmut', action="store", dest='parental_a_mutant', required=True) #mutated genome
parser.add_argument('-parpol', action="store", dest='parental_b_polymorphic', required=True) #polymorphic genome
parser.add_argument('-mutpos', action="store", dest='mut_pos',type=str ,required=True)
#parser.add_argument('-mutbpos', action="store", dest='mut_b_pos', type=str)
parser.add_argument('-smod', action="store", dest='selection_mode',
required=True, choices=set(('r','d','di','dr','wt','mc'))) #Choose between... These 4 options are avaiable for the standalone program, 'r' and 'd' for now. 
parser.add_argument('-nrec', action="store", dest='nbr_rec_chrs', type=int, required=True)
args = parser.parse_args()

out_dir = args.out_dir
rfd = args.rec_freq_distr
parental_a_mutant = args.parental_a_mutant
parental_b_polymorphic = args.parental_b_polymorphic
mut_a = args.mut_pos.split("-")[0].split(",")
mut_a_chrom = int(mut_a[0])
mut_a_pos = int(mut_a[1])
selection_mode = args.selection_mode # "r" = recessive, "d" = dominant mt-phe, "di" = dominant wt-phe, "dr" = double recessive
if selection_mode == "dr":
	mut_b =  args.mut_pos.split("-")[1].split(",")
	mut_b_chrom = int(mut_b[0])
	mut_b_pos = int(mut_b[1])
if selection_mode == "mc":
	mut_b =  args.mut_pos.split("-")[1].split(",")
	mut_b_chrom = int(mut_b[0])
	mut_b_pos = int(mut_b[1])
nbr_haploid_recombinants = args.nbr_rec_chrs

if selection_mode == 'dr' and mut_b_pos is None:
	quit('Quit. Selected mode is "double recessive" but no position for mutation b was provided. See program description.')



# Functions

# Function to transform 'rfd' input into 'chr_xo_freq'
def chr_freq_generator(rfd):
	chr_xo_freq = []
	for contig_fr in rfd.split("/"):
		i = 0
		f = 0
		former_value = 0
		chr_xo= []
		for items in contig_fr.split("-"):
			items = items.split(",")
			position = []
			if f != 0:
				i = i + int(former_value)
			position.append(i)
			former_value = items[1]
			f += int(items[1])
			position.append(f)
			position.append(int(items[0]))
			chr_xo.append(position)
		chr_xo_freq.append(chr_xo)
	return chr_xo_freq

# Function to parse fasta files

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




# Function to create a recombinant chromosome from two parental chromosomes and
# the list 'crossover_positions' with XO positions.
# I decide to group this code in function to be used in the different selection modes
def create_rec_seq(crossover_positions, starting_parental):
	rec_chr = ""
	for key,val in enumerate(crossover_positions):	
		if starting_parental == 0:
			if key % 2 == 0:
				try: rec_chr += seq_parental_a[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass
			if key % 2 != 0:
				try: rec_chr += seq_parental_b[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass
		if starting_parental == 1:
			if key % 2 == 0:
				try: rec_chr += seq_parental_b[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass
			if key % 2 != 0:
				try: rec_chr += seq_parental_a[crossover_positions[key]:crossover_positions[key+1]]
				except IndexError: pass					
	return rec_chr

# Function to divide a long string ('data') into chunks of defined length ('batch_size')
def batch_gen(data, batch_size):
	for i in range(0, len(data), batch_size):
		yield data[i:i+batch_size]

# Function to create Fasta file and write the sequence of a recombinant chromosome to it.
def create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1):
	

	with open(out_dir + '/rec_gnm_' + str(iter1 + 1) + '.fa', 'a') as output_file:
		if n_contig == 0:
			output_file.write('>'+ contig_parental_a[n_contig][0][1:])


		else:

			output_file.write("\n"+'>'+ contig_parental_a[n_contig][0][1:])

		for chunk in batch_gen(rec_chr, 80): # Write to file a small chunk of the contig sequence in each line
			output_file.write('\n' + chunk)


# Create a subdirectory to place the recombinant chromosomes. If it already exists, remove it first
if os.path.exists(out_dir): # In easymap, out_dir points to './project/0_input/sim_data/sim_recsel_output'
	shutil.rmtree(out_dir)
os.makedirs(out_dir)

# Read fasta files of parentals
contig_parental_a = []
total_a= 0
with open(parental_a_mutant) as fp:
	for name_parental_a, seq_parental_a in read_fasta(fp):
		lists= [name_parental_a, seq_parental_a]
		contig_parental_a.append(lists)
		total_a += len(seq_parental_a)
contig_parental_b = []
total_b= 0
with open(parental_b_polymorphic) as fp:
	for name_parental_b, seq_parental_b in read_fasta(fp):
		lists = [name_parental_b, seq_parental_b]
		total_b += len(seq_parental_b)
		contig_parental_b.append(lists)


# Check that both parentals have identical length. If not, quit.
if total_a != total_b:
	quit('Quit. The lengths of the two parental sequences provided are not identical')


chr_xo_freq_total = chr_freq_generator(rfd)
# Create a defined number of recombinant chromosomes 
iter1 = 0

nbr_contigs = len(contig_parental_a)

while iter1 < nbr_haploid_recombinants:
	iter2 = 0
	is_in = []
	while iter2 < nbr_contigs:
		crossover_positions = []
		list_seq_parental_a=contig_parental_a[iter2][1]
		list_seq_parental_b=contig_parental_b[iter2][1]

		chr_len = len(list_seq_parental_a)
		chr_xo_freq = chr_xo_freq_total[iter2]
		
		# Reset variables to 0 ('does not contain the mutation') at the beggining of each loop
		chr_carries_mutation_a = False
		chr_carries_mutation_b = False

		# Randomly choose which of the two parentals the program will choose as "seed".
		# The probabilities are 0.5/0.5. Use later.
		starting_parental = randint(0, 1) # 0: start with mutated (parmut) genome, 1: start with polymorphic (parpol) genome
		# Create a radom number and compare it to the XO frequency table of the user-chosen chr 
		# The XO frequency table represents the frequency of each number of XOs in a chromosome
		# By comparing a series of random numbers with this table, a series of XOs can be created 
		# numbers that follow exactly the real frequencies
		
		rand_nbr = randint(1, 100) # 100 different options (values 1 and 100 are included)
		for i in chr_xo_freq:
			if rand_nbr > i[0] and rand_nbr <= i[1]:
				nbr_crossovers = i[2] # This is the number of XOs the current rec chr will have 
		# Randomly create the genomic position of each XO event
		iter3 = 0
		crossover_positions = []
		while iter3 < nbr_crossovers:
			crossover_positions.append(randint(1, chr_len))
			iter3 +=1
		
		# Add the beggining and end coordinates of the chromosome to the list 'crossover_positions' 
		crossover_positions.append(0)
		crossover_positions.append(chr_len)
		crossover_positions.sort()

		
		# Determine if the resulting recombinant chr carries the primary (a) causal mutation

		if mut_a_chrom == iter2 +1:
			for key,val in enumerate(crossover_positions):
				if starting_parental == 0 and key % 2 == 0:
					if mut_a_pos > crossover_positions[key] and mut_a_pos <= crossover_positions[key+1]:
						chr_carries_mutation_a = True
				
				if starting_parental == 1 and key % 2 != 0:
					if mut_a_pos > crossover_positions[key] and mut_a_pos <= crossover_positions[key+1]:
						chr_carries_mutation_a = True
		
		
		# Determine if the resulting recombinant chr carries the secondary (b) causal mutation
		
		if selection_mode == "dr":
		 	if mut_b_chrom == iter2+1:
				for key,val in enumerate(crossover_positions):
					
					if starting_parental == 0 and key % 2 == 0:
						if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
							chr_carries_mutation_b = True

					if starting_parental == 1 and key % 2 != 0:
						if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
							chr_carries_mutation_b = True


		if selection_mode == "cm": 
			if mut_b_chrom == iter2+1:
				for key,val in enumerate(crossover_positions):
					
					if starting_parental == 0 and key % 2 != 0:
						if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
							chr_carries_mutation_b = True

					if starting_parental == 1 and key % 2 == 0:
						if mut_b_pos > crossover_positions[key] and mut_b_pos <= crossover_positions[key+1]:
							chr_carries_mutation_b = True

		# Chromosome selection
		# If recombinant chr contains the desired mutation(s), execute 'create_rec_seq()'
		# and 'create_rec_chr_file()' functions.
		# The first function creates a Fasta file using the info in 'crossover_positions', and
		# the second writes the output to a file.
		
		# Select all chromosomes that carry mutation A (all phenotypically mutant plants)
		n_contig = iter2
		if selection_mode == 'r':
			if chr_carries_mutation_a == True:
				seq_parental_a =  list_seq_parental_a
				seq_parental_b =  list_seq_parental_b
				rec_chr = create_rec_seq(crossover_positions, starting_parental)
				create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
				is_in.append(n_contig)
				iter2 += 1

			elif mut_a_chrom != iter2+1:
				seq_parental_a =  list_seq_parental_a
				seq_parental_b =  list_seq_parental_b
				rec_chr = create_rec_seq(crossover_positions, starting_parental)
				create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
				is_in.append(n_contig)
				iter2 +=1

		
		# Select all chromosomes that carry mutation A and also some that do not, so the final
		# proportion of chrs that carry the mutation is 0.67 (all phenotypically mutant plants).
		# This happens when selecting mutants with dominant mutation based on phenotype.
		if selection_mode == 'd':
			if chr_carries_mutation_a == True:
				seq_parental_a =  list_seq_parental_a
				seq_parental_b =  list_seq_parental_b
				rec_chr = create_rec_seq(crossover_positions, starting_parental)
				create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
				is_in.append(n_contig)
				iter2 += 1

			else:
				chr_filtering_threshold = randint(1,100)
				if chr_filtering_threshold > 50:
					seq_parental_a =  list_seq_parental_a
					seq_parental_b =  list_seq_parental_b
					rec_chr = create_rec_seq(crossover_positions, starting_parental)
					create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
					is_in.append(n_contig)
					iter2 += 1
		

		# Select all chromosomes that do not carry mutation A (all phenotypically wild type plants)
		if selection_mode == 'di':
			if chr_carries_mutation_a == False:
				seq_parental_a =  list_seq_parental_a
				seq_parental_b =  list_seq_parental_b
				rec_chr = create_rec_seq(crossover_positions, starting_parental)
				create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
				is_in.append(n_contig)
				iter2 += 1
		
		# Select all chromosomes that do not carry mutation A and also some that do, so the final
		# proportion of chrs that carry the mutation is 0.33 (all phenotypically wt plants).
		# This happens when selecting wild type plants in from a recessive cross.	
		if selection_mode == 'wt':
			if chr_carries_mutation_a == False:
				seq_parental_a =  list_seq_parental_a
				seq_parental_b =  list_seq_parental_b
				rec_chr = create_rec_seq(crossover_positions, starting_parental)
				create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
				is_in.append(n_contig)
				iter2 += 1

			else:
				chr_filtering_threshold = randint(1,100)
				if chr_filtering_threshold > 50:
					seq_parental_a =  list_seq_parental_a
					seq_parental_b =  list_seq_parental_b
					rec_chr = create_rec_seq(crossover_positions, starting_parental)
					create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
					is_in.append(n_contig)
					iter2 += 1

		# Select all chromosomes that carry mutations A and B (all phenotipically double recessive mutants)
		if selection_mode == 'dr' or selection_mode == "cm":

			if mut_a_chrom == mut_b_chrom and mut_a_chrom == iter2+1:
				if chr_carries_mutation_a == True and chr_carries_mutation_b == True:

					seq_parental_a =  list_seq_parental_a
					seq_parental_b =  list_seq_parental_b
					rec_chr = create_rec_seq(crossover_positions, starting_parental)
					create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
					is_in.append(n_contig)
					iter2 += 1
					

			elif mut_a_chrom == iter2+1 and mut_b_chrom != iter2+1:

				if chr_carries_mutation_a == True:
					seq_parental_a =  list_seq_parental_a
					seq_parental_b =  list_seq_parental_b
					rec_chr = create_rec_seq(crossover_positions, starting_parental)
					create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
					is_in.append(n_contig)
					iter2 += 1


			elif mut_a_chrom != iter2+1 and mut_b_chrom == iter2+1:
				if chr_carries_mutation_b == True:
					seq_parental_a =  list_seq_parental_a
					seq_parental_b =  list_seq_parental_b
					rec_chr = create_rec_seq(crossover_positions, starting_parental)
					create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
					is_in.append(n_contig)
					iter2 += 1

			else:
					seq_parental_a =  list_seq_parental_a
					seq_parental_b =  list_seq_parental_b
					rec_chr = create_rec_seq(crossover_positions, starting_parental)
					create_rec_chr_file(n_contig,rec_chr, contig_parental_a,iter1)
					is_in.append(n_contig)
					iter2 += 1




		if len(is_in) == nbr_contigs:
			iter1+=1