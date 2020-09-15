#
# This program accepts a list of variants (two modes available: SNPs [snp] / large insertions [lim])
# and returns their putative effect on genes and proteins. For that it needs a gff file and a fasta file
# of the genome.
# 
# It is optimized for samples where a significant part of the mutation positions do not
# lie within transcribed regions. It does a first screening to detect all mutation positions that
# do not lie within transcribed regions (here only 'mrna' features from file gff are loaded, saving work)
# and then does a second (much more computationally demanding) analysis only focused on mutations 
# that do lie within transcribed regions. First, mutation positions are tested for untranscribed exonic
# regions and for introns (here all cds and UTR features from gff are loaded into memory). Finally, only
# if mutation positions are in coding sequences, it is analyzed their effect on a protein (computationally
# intensive and requires use of fasta chromosomes). This way of conditional operation avoids doing unnecessary
# work. The goal is to reduce the long execution time of previous version of this program, so it can be used
# with hundreds of thousands of mutations.
#
# The input variants must be in the following format (... = other fields can be present here):
# snp mode: CHR\tPOS\tREF\tALT\t...
# lim mode: CHR\tPOS\t-\t-\t...
#
# All arguments are required except '-ann'.
#
# Adding functional info to the gene hits is optional and requires a 2-column file with the
# following format: gene\tdescription. To turn it on use argument '-ann <source of functional descriptions>'.
#
# TO DO: maybe check
# a few ref letters in the input to see if they match the reference, consider the possibility of analyzing
# indels as if they were lims (just tell where the protein is broken).
#
# David Wilson - dws1985@hotmail.com
#
#

# UPDATE - Fixed splicing detection, SDL




import argparse
from string import maketrans


# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-pname', action="store", dest='project_name', required=True)
parser.add_argument('-out', action="store", dest='output', required=True)
parser.add_argument('-itp', action="store", dest='input_type', choices=set(('snp','lim')), required=True)
parser.add_argument('-con', action="store", dest='contigs_source', required=True)
parser.add_argument('-gff', action="store", dest='gff_source', required=True)
parser.add_argument('-var', action="store", dest='variants_source', required=True)
parser.add_argument('-rrl', action="store", dest='regulatory_region_length', required=True) # To turn off, set to 0
parser.add_argument('-ann', action="store", dest='gene_ann_source')

args = parser.parse_args()


project = args.project_name
input_type = args.input_type
contigs_source = args.contigs_source
gff_source = args.gff_source
variants_source = args.variants_source
regulatory_region_length = int(args.regulatory_region_length)
gene_ann_source = args.gene_ann_source
output = args.output

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


# Function to translate DNA into protein
#CHECK THAT ALL INFO IS CORRECT
def dna_to_prot(dna_seq):
	genetic_code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
                   "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                   "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
                   "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
                   "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                   "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                   "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                   "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                   "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
                   "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                   "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                   "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                   "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
                   "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                   "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                   "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
   
	prot_seq = []
	start = 0
	while start + 2 < len(dna_seq):
		codon = dna_seq[start:start+3]
		try:
			if genetic_code[codon] == "STOP":
				prot_seq.append('*')
				break
			prot_seq.append(genetic_code[codon])
		except KeyError: # To prevent non ATGC nts cause the script to throw an error and stop
			prot_seq.append('_UNKNOWCODON')
			break
		start += 3
	prot_seq = ''.join(prot_seq)
	return prot_seq


# Function to obtain the reverse complementary of a DNA sequence
def reverse_complementary(seq):
	revcomp = seq.translate(maketrans('ACGT', 'TGCA'))[::-1]
	return revcomp


# If input is a SNPs, extract and store in array only selected info (type, chr, pos, ref, alt)
# If input is a custom file comming from the analysis of large insertion mutations (LIM),
# simply store it in mut_array. This is currently under development
mut_array = []
if input_type == 'snp':
	# Extract needed info from pseudovcf input file and store it in the array 'mut_array'
	# muts_array.append('@type\tchr\tpos\tref\talt')
	with open(variants_source) as input_mut:
		for line_mut in input_mut:
			if not line_mut.startswith('#'):
				fields_mut = line_mut.split('\t')
				useful_mut_info = fields_mut[0].lower(), fields_mut[1].lower(), int(fields_mut[2]), fields_mut[3].upper(), fields_mut[4].upper().strip()
				mut_array.append(useful_mut_info)

if input_type == 'lim':
	with open(variants_source) as input_mut:
		for line_mut in input_mut:
			if not line_mut.startswith('#'):
				fields_mut = line_mut.split('\t')
				useful_mut_info = fields_mut[0].lower(), fields_mut[1].lower(), int(fields_mut[2]), fields_mut[3], fields_mut[4].strip()
				mut_array.append(useful_mut_info)
				

# Extract needed info from template gff file: only mRNAs and certain columns
# The following code only checks whether a mutation position lies within a mRNA sequence
gff_array1 = []
with open(gff_source) as input_gff:
	for line_gff in input_gff:
		if not line_gff.startswith('#'):
			fields_gff = line_gff.split('\t')
			if fields_gff[2].lower() == 'mrna':
				useful_gff_info = fields_gff[2].lower(), fields_gff[0].lower(), int(fields_gff[3]), int(fields_gff[4]), fields_gff[6], fields_gff[8].split(';')[0][3:]
				gff_array1.append(useful_gff_info)

# Check whether each mutation position lies within a mRNA sequence or a putative regulatory region of the template gff file
variants_info = []
for variant in mut_array:

	# Reset 'hit' variable
	is_hit = False
	
	# For each input mutation/variant, go through the gff info and detect hits. A single mutation can
	# affect simultaneously more than one transcription unit or putatve regulatory sequence. This is 
	# handled in the code. The output is written to 'variants_info' array.
	for mrna in gff_array1:
	
		if variant[1] == mrna[1] and variant[2] >= mrna[2] and variant[2] <= mrna[3]: # mrna[1]: mRNA chrom, mrna[2]: mRNA left coord, mrna[3]: mRNA right coord
			is_hit = True
			hit = 'tu'
			dumped_info = variant[0], variant[1], variant[2], variant[3], variant[4], hit, mrna[2], mrna[3], mrna[4], mrna[5]
			variants_info.append(dumped_info)
		
		if regulatory_region_length > 0:
		
			if mrna[4] == '+':
				if variant[1] == mrna[1] and variant[2] < mrna[2] and variant[2] >= (mrna[2] - regulatory_region_length):
					is_hit = True
					hit = 'rr'
					dumped_info = variant[0], variant[1], variant[2], variant[3], variant[4], hit, mrna[2], mrna[3], mrna[4], mrna[5], 'promoter', '-', '-', '-'
					variants_info.append(dumped_info)

			if mrna[4] == '-':
				if variant[1] == mrna[1] and variant[2] > mrna[3] and variant[2] <= (mrna[3] + regulatory_region_length):
					is_hit = True
					hit = 'rr'
					dumped_info = variant[0], variant[1], variant[2], variant[3], variant[4], hit, mrna[2], mrna[3], mrna[4], mrna[5], 'promoter', '-', '-', '-'
					variants_info.append(dumped_info)

	# If no mRNA or putative regulatory region is hit, parse the mutation and add its info to 'variants_info'
	if is_hit == False:
		hit = 'nh'
		dumped_info = variant[0], variant[1], variant[2], variant[3], variant[4], hit, '-', '-', '-', '-', '-', '-', '-', '-'
		variants_info.append(dumped_info)

del mut_array, gff_array1


# From input gff file, load (cds, exon, UTRs x chrom, feat, start, end, model) in array
# I add exons so I can later use this array to analyze mutation within introns, although this
# makes more complicated the detection on the functional element affected by the mutation just
# after "if variant_info[5] == 'tu':".
gff_array2 = []
with open(gff_source) as input_gff:
	for line_gff in input_gff:
		if not line_gff.startswith("#"):
			fields_gff = line_gff.split('\t')
			if fields_gff[2].lower() == 'cds' or fields_gff[2].lower() == 'exon' or fields_gff[2].lower() == 'five_prime_utr' or fields_gff[2].lower() == 'three_prime_utr':
				useful_gff_info = fields_gff[0].lower(), int(fields_gff[3]), int(fields_gff[4]), fields_gff[2].lower(), fields_gff[8]
				gff_array2.append(useful_gff_info)

# Analyze variants that are marked as interrupting a mRNA or putative regulatory region
variants_info2 = []
for variant_info in variants_info:

	if variant_info[5] == 'tu':
		
		# Check if mutation position lies in 'UTRs' (untranslated regions) or introns
		feature_hit = 'intron' #Since gff files do not contain introns info, I set 'feature_hit' default value to 'intron'
		
		for feature in gff_array2:
			
			if feature[3] != 'exon': # This if statements is necessary because array also coatains exon records (besides cds and utrs)
				if variant_info[9] in feature[4] and variant_info[2] >= feature[1] and variant_info[2] <= feature[2]:
					feature_hit = feature[3]
		
		# Obtain and store the exons coordinates of the current gene (this is needed for both exon and intron muations)
		exon_coords_list = []
		for feature in gff_array2:
			if variant_info[9] in feature[4] and feature[3] == 'exon':
				exon_coords = feature[1],feature[2]
				exon_coords_list.append(exon_coords)
		
		number_of_exons = len(exon_coords_list)

		# Reverse the list with exon coordinates to more easily calculate intron boundaries downstream
		if variant_info[8] == '-':
			exon_coords_list.reverse()

		# If feature hit is an exon evaluate if splicing signals are affected
		if feature_hit != 'intron':

			# Only take into account the first and last bases of each exon (Brent & Guigo 2004 - Recent advances in gene structure prediction. Current Opinion in Structural Biology)
			numberOfExonBasesConsidered = 3

			exon_counter = 1
			exon_left_end_hit = False
			exon_right_end_hit = False
			
			# Loop through the exon coordinates and compare them with the mutation position to find possible splicing sites affected
			exon_splicing_hit = False
			while exon_counter <= number_of_exons:

				exon_left_coord = exon_coords_list[exon_counter-1][0]
				exon_right_coord = exon_coords_list[exon_counter-1][1]
				
				if exon_counter == 1:
					if abs(int(variant_info[2]) - exon_right_coord) < numberOfExonBasesConsidered:
						exon_affected = exon_counter
						exon_end_affected = 'right'
						exon_splicing_hit = True
						break
				elif exon_counter == number_of_exons:
					if abs(int(variant_info[2]) - exon_left_coord) < numberOfExonBasesConsidered:
						exon_affected = exon_counter
						exon_end_affected = 'left'
						exon_splicing_hit = True
						break
				else:
					if abs(int(variant_info[2]) - exon_left_coord) < numberOfExonBasesConsidered:
						exon_affected = exon_counter
						exon_end_affected = 'left'
						exon_splicing_hit = True
						break
					if abs(int(variant_info[2]) - exon_right_coord) < numberOfExonBasesConsidered:
						exon_affected = exon_counter
						exon_end_affected = 'right'
						exon_splicing_hit = True
						break

				exon_counter += 1
			
			if exon_splicing_hit == True:
				if variant_info[8] == '-':
					exon_affected = number_of_exons - exon_affected + 1
					if exon_end_affected == 'left': exon_end_affected = '3\''
					if exon_end_affected == 'right': exon_end_affected = '5\''
				else:
					if exon_end_affected == 'left': exon_end_affected = '5\''
					if exon_end_affected == 'right': exon_end_affected = '3\''

				exonSplicingSignal = ', putative splicing signal in the ' + exon_end_affected + ' end of exon ' + str(exon_affected) + ' affected'
			
			else:
				exonSplicingSignal = ''


			# If feature hit is CDS, evaluate if mutation has an efect on the aminoacid sequence
			if feature_hit == 'cds':

				# Create a list with the start and end coordinates of each cds stretch of the gene
				cds_list = []
				for feature in gff_array2:
					if variant_info[9] in feature[4] and feature[3] == 'cds':
						cds_coords = feature[1], feature[2]
						cds_list.append(cds_coords)

				# Reconstruct the coding sequence of the wild type gene.
				target_fasta_header = '>' + variant_info[1]
						
				with open(contigs_source) as fp:
					for name_contig, seq_contig in read_fasta(fp): # Create an array with the names and sequences of the contigs in the fasta input using the function 'read_fasta(fp)'
						if name_contig.lower() == target_fasta_header.lower(): # Only work with the contig where the mutation lies
							cds_seq_list_wt = []
							cds_seq_list_mt = []
							for cds in cds_list:
								cds_seq = seq_contig[cds[0]-1 : cds[1]]
								
								if cds[0]-1 <= variant_info[2] and cds[1] >= variant_info[2]:
									# Calculate the position of the mutation in the CDS sequence
									relative_mut_pos = variant_info[2] - cds[0]
									
									# Replace wt-base for mut-base at mut-pos (if input is large insertions, 
									# substitute wt-base for the symbol'-').
									cds_seq_as_list = list(cds_seq)									
									cds_seq_as_list[relative_mut_pos] = variant_info[4]
									cds_seq_mut = ''.join(cds_seq_as_list)
									
									# Append cds wt seq to wt list and cds mut seq to mt list
									# If mRNA is in the reverse strand, reverse complement the cds sequences
									# when adding them to the lists 'cds_seq_list_wt' and 'cds_seq_list_mt'
									if variant_info[8] == '+':
										cds_seq_list_wt.append(cds_seq)
										cds_seq_list_mt.append(cds_seq_mut)
									if variant_info[8] == '-':
										cds_seq_list_wt.append(reverse_complementary(cds_seq))
										cds_seq_list_mt.append(reverse_complementary(cds_seq_mut))
								
								else:
									if variant_info[8] == '+':
										cds_seq_list_wt.append(cds_seq)
										cds_seq_list_mt.append(cds_seq)
									if variant_info[8] == '-':
										cds_seq_list_wt.append(reverse_complementary(cds_seq))
										cds_seq_list_mt.append(reverse_complementary(cds_seq))					
				
				# Reconstruct the coding sequence of the mutant gene
				full_cds_seq_wt = (''.join(cds_seq_list_wt)).upper()
				full_cds_seq_mt = (''.join(cds_seq_list_mt)).upper()
					
				if input_type == 'snp':	
					# Translate the coding sequences of the wild type and mutant genes
					prot_wt = dna_to_prot(full_cds_seq_wt)
					prot_mt = dna_to_prot(full_cds_seq_mt)
					
					# Determine if protein has an amino acid change. If so, store its position and the wt and mut aas
					aa_change = False
					aa_position = 1
					for aa_wt, aa_mt in zip(prot_wt, prot_mt):
						if aa_wt != aa_mt:
							aa_change = True
							if aa_mt == "*": aa_mt = "STOP"
							result_aa_wt, result_aa_mt = aa_wt, aa_mt
							result_aa_position = aa_position								
						aa_position += 1
					
					if aa_change == False:
						result_aa_wt, result_aa_mt, result_aa_position = '-', '-', 'no aa change'
					
					# Write info as a comma-separated list to the list 'variants_info2'
					condensed_info = variant_info[0], variant_info[1], variant_info[2], variant_info[3], variant_info[4], variant_info[5], variant_info[6], variant_info[7], variant_info[8], variant_info[9], feature_hit + exonSplicingSignal, result_aa_position, result_aa_wt, result_aa_mt
					variants_info2.append(condensed_info)
				
				if input_type == 'lim':
					#Determine position of insertion in protein sequence
					result_nt_position = int(float(full_cds_seq_mt.find('-') + 1)/3)
					
					# Write info as a comma-separated list to the list 'variants_info2'
					condensed_info = variant_info[0], variant_info[1], variant_info[2], variant_info[3], variant_info[4], variant_info[5], variant_info[6], variant_info[7], variant_info[8], variant_info[9], feature_hit, result_nt_position, '-', '-'
					variants_info2.append(condensed_info)
			
			else:
				condensed_info = variant_info[0], variant_info[1], variant_info[2], variant_info[3], variant_info[4], variant_info[5], variant_info[6], variant_info[7], variant_info[8], variant_info[9], feature_hit + exonSplicingSignal, '-', '-', '-'
				variants_info2.append(condensed_info)

		# If intron is hit, evaluate if splicing is affected
		elif feature_hit == 'intron' and input_type == 'snp':

			# Obtain a list with the index, start and end coordinates of the introns of the hit gene
			# First, set some variables to starting values
			number_of_introns = len(exon_coords_list) - 1
			intron_counter = 1
			intron_left_end_hit = False
			intron_right_end_hit = False
			
			# Loop through the intron coordinates (determined on the fly based on exon coordinates
			# in 'exon_coords_list') and compare them with mutation position to find possible splicing
			# sites affected
			while intron_counter <= number_of_introns:
				intron_left_coord = exon_coords_list[intron_counter-1][1] + 1 # Because the first intron starts 1 position after the first exon ends
				distance_left = variant_info[2] - intron_left_coord           # To calculate the distance between the intron beginning and the position of the mutation
				
				if 0 <= distance_left < 8:
					intron_left_end_hit = True
					result_intron_number = intron_counter
					break
				
				intron_right_coord = exon_coords_list[intron_counter][0] - 1
				distance_right = intron_right_coord - variant_info[2]
				
				if 0 <= distance_right < 8:
					intron_right_end_hit = True
					result_intron_number = intron_counter
					break
					
				intron_counter += 1
			
			if intron_left_end_hit == True or intron_right_end_hit == True:
				if variant_info[8] == '+':
					intron_left_end, intron_right_end = 'donor', 'acceptor'
				else:
					intron_left_end, intron_right_end = 'acceptor', 'donor'
					result_intron_number = number_of_introns - result_intron_number +1
			
			if intron_left_end_hit == True:
				intron_result = 'intron, putative splicing ' + str(intron_left_end) + ' sequence of intron ' + str(result_intron_number) + ' affected'
			elif intron_right_end_hit == True:
				intron_result = 'intron, putative splicing ' + str(intron_right_end) + ' sequence of intron ' + str(result_intron_number) + ' affected'
			else:
				intron_result = 'intron'
			
			# Write info as a comma-separated list to the list 'variants_info2'
			condensed_info = variant_info[0], variant_info[1], variant_info[2], variant_info[3], variant_info[4], variant_info[5], variant_info[6], variant_info[7], variant_info[8], variant_info[9], intron_result, '-', '-', '-'
			variants_info2.append(condensed_info)

		# Introns in and lim input type
		else:
			# If transcriptional unit is hit but is neither in cds or in intron
			condensed_info = variant_info[0], variant_info[1], variant_info[2], variant_info[3], variant_info[4], variant_info[5], variant_info[6], variant_info[7], variant_info[8], variant_info[9], feature_hit, '-', '-', '-'
			variants_info2.append(condensed_info)

	else:
		# If no transcriptional unit has been hit, simply copy 'variants_info' to the new array 'variants_info2'
		variants_info2.append(variant_info)

del input_mut, input_gff, variants_info
	

# Retrieve gene functional annotation and create final output

# Create output file
output = open(output, 'w')

if input_type == 'snp': header_pos = 'aa_pos'
if input_type == 'lim': header_pos = 'nt_pos'

# If no gene annotation file provided, simply print 'variants_info2' to output file.
if gene_ann_source == 'user_data/n/p':
	output.write('@type\tcontig\tposition\tref_base\talt_base\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\t' + header_pos + '\taa_ref\taa_alt\tgene_annotation_info\n')
		
	for variant in variants_info2:	
		for index, field in enumerate(variant):
			if index == 0:
				output.write(str(field))
			else:
				output.write('\t' + str(field))
		output.write('\t-\n')
	output.close()	

# If genome annotation file provided, merge that info with variants_info2 i a new list called 'variants_info3'
else:
	# Open gene annotation file and load contents in array
	ann_array = []
	with open(gene_ann_source) as gene_ann_input:
		for line_ann in gene_ann_input:
			ann_array.append(line_ann)
	
	# Create an array that will contain all the info
	variants_info3 = []
	
	# Iterate over 'variants_info2', get gene,
	for variant_info2 in variants_info2:
		
		# First, determine if the variant has interrupted a gene
		hit_gene = ''
		try:
			hit_gene = variant_info2[9][:-2]
		except IndexError:
			pass
		
		# If the variant interrupts a gene, search for the gene in 'ann_array', and, if found,
		# get any functional info available. Then combine info in 'variants_info3'
		
		if hit_gene != '':
			# Default value
			ann_result = 'Gene not found in the annotation file'
			
			for ann_gene in ann_array:
				ann_fields = ann_gene.split('\t')
				if hit_gene == ann_fields[0].strip('\n'):
					ann_info_fields = ann_fields[1:]
					ann_info_string = '; '.join(ann_info_fields)
					if ann_info_string == '': ann_info_string = 'Information not found in the annotation file'
					ann_result = ann_info_string.strip('\n')
					ann_result = ann_result.replace('\t', "; ")

			# Combine info and append it to 'variants_info3'
			condensed_info2 = variant_info2[0], variant_info2[1], variant_info2[2], variant_info2[3], variant_info2[4], variant_info2[5], variant_info2[6], variant_info2[7], variant_info2[8], variant_info2[9], variant_info2[10], variant_info2[11], variant_info2[12], variant_info2[13], ann_result
			variants_info3.append(condensed_info2)
		
		# If the variant does not interrupt any gene, just copy its info from variants_info2
		else:
			tmp_variant_info2_list = list(variant_info2); tmp_variant_info2_list.append('-'); variant_info2 = tuple(tmp_variant_info2_list)
			variants_info3.append(variant_info2)

	output.write('@type\tcontig\tposition\tref_base\talt_base\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\t' + header_pos + '\taa_ref\taa_alt\tgene_annotation_info\n')	
	
	for variant in variants_info3:		
		for index, field in enumerate(variant):
			if index == 0:
				output.write(str(field))
			else:
				output.write('\t' + str(field))
		
		output.write('\n')
	output.close()

