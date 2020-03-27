# 
# 
# This program simulates a massively parallel sequencer. It accepts either a single DNA sequence
# or a population, as in mapping by sequencing experiments. The input format is fasta, while the
# output format is fastq. The user can choose to simulate either single- or paired-end libraries.
# The program can simulate fixed- (e.g. HiSeq2000) and variable-length (e.g. IonProton) reads.
# Basecalling errors can be introduced. High-throughput sequencing is mostly random across a template,
# but several sources of bias exist. A well-documented source is the preference for neutral GC content
# fragments during library PCR amplification. This program uses GC content as the input for a bias
# mechanism. The strength of this bias is adjustable by the user. Quality scores (phred scaled) are
# all set to 40 (the maximum value in Phred scale). In Sanger encoding, this value corresponds to 'I'.
# 
# PARAMETERS EXPLANATION:
# 
# -if: Input folder. The name of folder relative to the script location that contains the template
# 		genome(s). Only template genome(s) can be present in the folder for the program to work properly.
# 		In addition, the program assumes that all template genomes are equal in length. If not, a warning
# 		message is show but the program continues execution. The template genomes must be in fasta format.
# 		To simulate reads from a single sample, simply place in this folder one fasta template. Each file
# 		must contain only one contig. It only accepts fasta files with the bases stored in
# 		'options_basecall_errors'.
# 
# -out: Name of directory relative to current location where to place the reads created.
# 
# -mod [se, pe]: Mode [single-end or paired end]. If paired-end mode chosen, the user must provide
# 		a value for '-flm' and '-fls' parameters.
# 
# -rd (integer): Read depth. The average number of reads by which a given locus in the genome is
# 		represented.
# 
# -rlm (integer): Read length mean.
# 
# -rls (integer): Read length standard deviation. To simulate reads with fixed length, set this
# 		parameter to 0; to simulate variable length reads with a normal distribution, set it to >0.
# 
# -flm (integer): Fragment length mean. flm >= 2*rlm. Required if -mod=pe.
# 
# -fls (integer): Fragment length standard deviation. Fragment length follows a normal distribution.
# 		Required if -mod=pe.
# 
# -ber (integer) [0-10%]: Basecalling error rate. Percentage of bases incorrectly called. If base
# 		calling occurs and template base is A, output can be C, G or T, each of them with 1/3 probability.
# 		Basecalling errors consist only in base substitutions, not indels.
# 
# -gbs (integer) [0 - 100]: GC/AT bias strength. The amount of effect caused by non-neutral GC/AT
# 		content during library amplification (see above). 
#
# 



import argparse, os, re, shutil
from random import randint, gauss
from string import maketrans

# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-input_folder', action="store", dest='input_folder', required=True)
parser.add_argument('-mod', action="store", dest='mode',
choices=set(('se','pe')), required=True) # Choose between se (= single end) and pe (= paired end)
parser.add_argument('-rd', action="store", dest='read_depth', type=int, required=True)
parser.add_argument('-rlm', action="store", dest='read_length_mean', type=int, required=True)
parser.add_argument('-rls', action="store", dest='read_length_sd', type=int, required=True)
parser.add_argument('-flm', action="store", dest='fragment_length_mean', type=int)
parser.add_argument('-fls', action="store", dest='fragment_length_sd', type=int)
parser.add_argument('-ber', action="store", dest='basecalling_error_rate', type=int, required=True)
parser.add_argument('-gbs', action="store", dest='gc_bias_strength', type=int, required=True)
parser.add_argument('-out', action="store", dest='out', required=True)

args = parser.parse_args()

# Assing argument parameters to variables
input_folder = args.input_folder
mode = args.mode
read_depth = args.read_depth
read_length_mean = args.read_length_mean # Must be > 0
read_length_sd = args.read_length_sd
basecalling_error_rate = args.basecalling_error_rate # As percentage [0-10%]. It is constrained because >10-20% increases dramatically processing required
gc_bias_strength = args.gc_bias_strength # As percentage [0-100%]

if mode == 'pe':
	if args.fragment_length_mean == None:
		quit('Quit. Paired end mode selected but no fragment length mean value was provided. See help.')
	if args.fragment_length_sd == None:
		quit('Quit. Paired end mode selected but no fragment length standard deviation value was provided. See help.')
	else:
		fragment_length_mean = int(args.fragment_length_mean) # >= 2 * read_length_mean. only required in pe mode
		fragment_length_sd = int(args.fragment_length_sd) # Only required in pe mode

#if it doesn't exist, create folder named 'output' inside the folder that contains the program
#output_folder = args.out # Set to '.' to output to same folder as .py file.
#if not os.path.exists(output_folder):
#   os.makedirs(output_folder)

# Create a subdirectory to place the reads. If it already exists, remove it first
output_folder = args.out
if os.path.exists(output_folder):
	shutil.rmtree(output_folder)
os.makedirs(output_folder)


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


#Function to calcule GC content of a read and filter it or not accordingly the the result and
#to the stringency chosen by the user
def depth_gc_bias(gc_bias_input_seq, gc_bias_input_length):
	read_gc_content = 100 * (gc_bias_input_seq.count('G') + gc_bias_input_seq.count('C')) / gc_bias_input_length
	distance_to_neutral_gc = abs(50 - read_gc_content)
	
	gc_random1 = randint(0,50)
	gc_random2 = randint(0,100)
	
	#print read_gc_content, 'dist:', distance_to_neutral_gc, 'random:', gc_random1
	
	#Compare read deviation from neutral GC content (possible values are 0-50) to a random
	#number between 0 and 50. The more neutral the GC content, the more probability that
	#the random number is the biggest. If so, the read passes the filter.
	if gc_random1 >= distance_to_neutral_gc:
		read_passes_gcbias_filter = True         
	
	#The following code allows the user to tune the stregth of the GC bias influence.
	#It basically allows to soften the effect of the filter above. 
	#If the read did not pass the GC filter in the first opportunity, it still can do it
	#but in a user-controlled way, by regulating "gc_bias_strength".
	#If gc_bias_strength is 0, all reads that did not pass int he first opportunity, will pass
	#this time (second opportunity). If 100, non of them will. If 50, only half of them 
	#(randomly determined by randint function) will.
	elif gc_random1 < distance_to_neutral_gc:
		if gc_random2 >= gc_bias_strength:
			read_passes_gcbias_filter = True
		else:
			read_passes_gcbias_filter = False
	else:
		read_passes_gcbias_filter = False
			
	#If the read passed the filter in first or second opportunity, the function returns True.
	#If not, the function returns False. This result is handled downstream in this script to
	#include or skip the read in the final reads dataset.
	
	return read_passes_gcbias_filter 


#Function to introduce basecalling errors in a read
#Important: Behaviour only checked for 1-5% error rate and 100-200 pb reads
def create_basecalling_errors(bce_seq, bce_length):
		
	#Define a different error rate for each execution based on user paramenter "basecalling_error_rate"
	error_rate_current_read = round(gauss(basecalling_error_rate, basecalling_error_rate))

	#Convert negative values to zero.
	if error_rate_current_read < 0:
		error_rate_current_read = 0
		
	#Calculate the absolute number of errors in current read	 
	number_of_errors_current_read = int(bce_length * (error_rate_current_read / 100))	
	
	#Create x="number_of_errors_current_read" error positions in read and store them in list "error_positions"
	error_positions = []
	error_iterator = 0
	while error_iterator < number_of_errors_current_read:
		new_error_position = randint(0, bce_length - 1) #0-based
		if new_error_position not in error_positions:
			error_positions.append(new_error_position)
			error_iterator += 1
	
	#Convert the string "bce_seq" in a list to enable introduction of basecalling errors
	#by item assignment. Downstream I convert it back to string.
	bce_seq_as_list = list(bce_seq)
	
	#For each error position in "error_positions" read real base and randomly determine error base.
	for error_position in error_positions:
		base = bce_seq[error_position:error_position+1] #0-based		
		#Define what are the alternative bases depending on the reference base. Ambiguous bases are ignored deliberately.
		options_basecall_errors = {
											'A': ['T', 'C', 'G'],'T': ['C', 'G', 'A'],'C': ['G', 'A', 'T'],'G': ['A', 'T', 'C'],
											'R': ['R', 'R', 'R'],'Y': ['Y', 'Y', 'Y'],'S': ['S', 'S', 'S'],'W': ['W', 'W', 'W'],
											'K': ['K', 'K', 'K'],'M': ['M', 'M', 'M'],'B': ['B', 'B', 'B'],'D': ['D', 'D', 'D'],
											'H': ['H', 'H', 'H'],'V': ['V', 'V', 'V'],'N': ['N', 'N', 'N'],'.': ['.', '.', '.'],
											'-': ['-', '-', '-']
										  }
		
		possible_error_bases = options_basecall_errors[base] #Lookup the alternative bases of "base" in "options_basecall_errors"
		error_base = possible_error_bases[randint(0, 2)] #Randomly pick between the possible mut bases
		bce_seq_as_list[error_position] = error_base #Substitute real base for error base

	#Convert "bce_seq" from list back to string	
	processed_bce_seq = ''.join(bce_seq_as_list) #I overwrite "bce_seq" variable with the new version with basecalling errors
	
	return processed_bce_seq


#Function to obtain the reverse complementary of a DNA sequence
def reverse_complementary(seq):
	revcomp = seq.translate(maketrans('ACGT', 'TGCA'))[::-1]
	return revcomp


#Miscellaneous processing steps before the actual sequencing:
#Parse input information, check input, create output folder, calculate the number of reads

#Create list with all tamplate genomes
template_genomes = sorted(os.listdir('./' + input_folder))

#Create list with templates lengths. This list is downstream used for checking format and length of data
template_lengths = []
for template_genome in template_genomes:
	seq_templates= []
	with open(input_folder + '/' + template_genome) as fp:
		for name, seq in read_fasta(fp):
			seq_templates.append(seq)
		seq_templates = "".join(seq_templates) 
		template_lengths.append(len(seq_templates))

#Check if input file is fasta formatted. If file is not in fasta format, function read_fasta()
#fails to create array with name(s) and sequence(s) of contig(s) in file, and list with template
#lengths is empty. If so, warn user that the input has not fasta format.
if len(template_lengths) == 0:
	quit('Warning: The program quit because input folder does not contain fasta files')

#Check length of all template genomes. If not all have same length, warn user
unique_template_lenghts = list(set(template_lengths))
num = len(unique_template_lenghts)
if num > 1 and num != 1:
	print 'Warning: Not all the templates provided have the same length. The program assumes all the templates have the same length to distribute reads across them. Therefore, read depth may not be homogeneous across contigs. If this is a concern, please provide contigs with equal lengths.\nContinuing anyway.'

#Calculate the number of reads needed to reach the chosen read depth
#Calculate the number of input template genomes provided by user and the number of reads needed
#per individual genome
genome_size = sum(template_lengths)/len(template_lengths)
number_of_genomes = len(template_lengths)
number_of_reads = (read_depth * genome_size) / read_length_mean
number_of_reads_per_genome = number_of_reads / number_of_genomes

#print 'number of genomes: ' + str(number_of_genomes)
#print 'number of reads per genome: ' + str(number_of_reads_per_genome)
#print 'number of total reads: ' + str(number_of_reads)
#print 'template lengths: ' + str(template_lengths)
#print 'mean genome length: ' + str(genome_size)
#print 'creating reads...'


#####################
# Do the sequencing #
#####################

#Split program in SE and PE modes because they have significant differences
if mode == 'se':

	#Create output fastQ file
	output = open(output_folder + '/se_reads.fq', 'w')

	#For each input genome to sequence, get its sequence and store the amount of reads to create in
	#'number_of_reads_per_genome'
	for genome_index, template_genome in enumerate(template_genomes):
		seq_template = []		
		with open(input_folder + '/' + template_genome) as fp:
			for name, seq in read_fasta(fp):
				seq_template.append(seq.upper())
		seq_template = "".join(seq_template)
		#Create reads
		read_count = 1
		while read_count <= number_of_reads_per_genome:
			genome_length = len(seq_template)  
	        #Calculate genome length for each genome, in case it is not always the same 
			#Determine length, start and end positions, and sequence of read
			read_length = int(round(gauss(read_length_mean, read_length_sd)))
			#Function "gauss()" might return not real-world values. I compensate for that.	
			if read_length < 5: read_length = 5
			#Define read start and end positions and take substring from template genome.
			read_start = randint(1, genome_length - read_length + 1)
			read_end = read_start + read_length
			read_seq = seq_template[read_start-1:read_end-1]                               #My values are 1-based, but python is 0-based. I compensate subtracting one position.
			
			#Subject read to GC bias filter. If it does not pass the filter, discard it and move on
			#to the next read without further operations.
			if gc_bias_strength > 0:
				gc_bias_result = depth_gc_bias(read_seq, read_length)
				if gc_bias_result == False:
					continue			
			
			#If user has chosen to introduce some amount of base calling errors, do so
			if basecalling_error_rate > 0:
				read_seq = create_basecalling_errors(read_seq, read_length)

			#Since reads can come from each strand with 0.5:0.5 probability,
			#Randomly reverse complement 50% of them.
			strand = randint(0,1)
			if strand == 1:
				read_seq = reverse_complementary(read_seq)
			
			#Create a string that contains the phred scores of the current read
			qual_seq = re.sub('\w', 'I', read_seq)
			
			#Write read to FastQ file
			output.write('@gn' + str(genome_index + 1) + '-rd' + str(read_count) + '\n' + read_seq + '\n+\n' + qual_seq + '\n')
			
			read_count += 1
			
	#Close FastQ file	
	output.close()

if mode == 'pe':

	#Create output forward and reverse fastQ files
	output_for = open(output_folder + '/pe-for_reads.fq', 'w')
	output_rev = open(output_folder + '/pe-rev_reads.fq', 'w')

	#For each input genome to sequence, get its sequence and store the amount of reads to create
	#in 'number_of_reads_per_genome'
	for genome_index, template_genome in enumerate(template_genomes):
		seq_template = []		
		with open(input_folder + '/' + template_genome) as fp:
			for name, seq in read_fasta(fp):
				seq_template.append(seq.upper())
		seq_template = "".join(seq_template)
		
		#Create reads
		fragment_count = 1
		while fragment_count <= number_of_reads_per_genome / 2:
			genome_length = len(seq_template)                                             #Calculate genome length for each genome, in case it is not always the same 
			
			#Determine length, start and end positions, and sequence of read pair
			fragment_length = int(round(gauss(fragment_length_mean, fragment_length_sd)))
			#Function "gauss()" might return not real-world values. I compensate for that.	
			if fragment_length < 2 * read_length_mean:
				fragment_length = 2 * read_length_mean
			#Define read pair start and end positions
			fragment_start = randint(1, genome_length - fragment_length + 1)
			fragment_end = fragment_start + fragment_length
			
			#Subject fragment to GC bias filter. If it does not pass the filter, discard it and move on
			#to the next fragment without further operations.
			#First, obtain sequence of current fragment
			if gc_bias_strength > 0:
				fragment_seq = seq_template[fragment_start-1 : fragment_end-1]
				gc_bias_result = depth_gc_bias(fragment_seq, fragment_length)
				if gc_bias_result == False:
					continue
			
			#Define lenghts of for and rev reads
			read_length_for = int(round(gauss(read_length_mean, read_length_sd)))
			read_length_rev = int(round(gauss(read_length_mean, read_length_sd)))
			#Function "gauss()" might return not real-world values. I compensate for that.
			if read_length_for < 5: read_length_for = 5
			if read_length_rev < 5: read_length_rev = 5
			
			#Define read start and end positions and read sequences in a single step			
			read_for_seq = seq_template[fragment_start-1 : fragment_start+read_length_for-1]    #My values are 1-based, but python is 0-based. I compensate subtracting one position.
			read_rev_seq = seq_template[fragment_end-read_length_rev-1 : fragment_end-1]
			
			#If user has chosen to introduce some amount of base calling errors, do so
			if basecalling_error_rate > 0:
				read_for_seq = create_basecalling_errors(read_for_seq, read_length_for)
				read_rev_seq = create_basecalling_errors(read_rev_seq, read_length_rev)

			#Obtain reverse complemetary of reverse read in pair
			read_rev_seq = reverse_complementary(read_rev_seq)

			#Determine which read will be first and which second in the pair
			#If necesssary, swap them
			pos_in_pair = randint(0,1) #randomly choose between two options
			if pos_in_pair == 1:
				read_for_seq, read_rev_seq = read_rev_seq, read_for_seq
			
			#Create a string that contains the phred scores of the current read
			qual_for_seq = re.sub('\w', 'I', read_for_seq)
			qual_rev_seq = re.sub('\w', 'I', read_rev_seq)
			
			#Write read pair to FastQ files
			output_for.write('@gn' + str(genome_index + 1) + '-rd' + str(fragment_count) + '/1\n' + read_for_seq + '\n+\n' + qual_for_seq + '\n')
			output_rev.write('@gn' + str(genome_index + 1) + '-rd' + str(fragment_count) + '/2\n' + read_rev_seq + '\n+\n' + qual_rev_seq + '\n')
			
			fragment_count += 1

	#Close FastQ files	
	output_for.close()
	output_rev.close()
