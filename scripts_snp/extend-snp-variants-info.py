
'''

This script was written to add more info to the list of variants that is presented to the user
when snp mode is used.

The new info is qual, ref_count, alt_count, af, dist_to_selected pos

The process from the first list of candidate variants to the final list that is given to user
is a bit complicated. The following scripts and files are involved:

Modes snp and ins: transform a .va file into a format suitable for varanalyzer.py

file.va --> ins-to-varanalyzer.py --> ins_to_varanalyzer.txt --> varanalyzer.py
file.va --> snp-to-varanalyzer.py --> snp_to_varanalyzer.txt --> varanalyzer.py

Example 1_intermediate_files/ins_to_varanalyzer.txt:
1	data
2	contig
3	pos
4	ref
5	alt

Example 1_intermediate_files/snp_to_varanalyzer.txt:
1	data
2	contig
3	pos
4	ref
5	alt
6	qual
7	ref_count
8	alt_count
9	af

Only the first 5 columns are taken by varanalyzer.py, no matter what the mode is. varanalyzer.py adds
another 9 columns to the output it produces

snp mode {
	varanalyzer.py: 1_intermediate_files/***_to_varanalyzer.txt --> 1_intermediate_files/varanalyzer_output.txt
	6	mrna_start
	7	mrna_end
	8	strand
	9	gene_model
	10	gene_element
	11	aa_pos
	12	aa_ref
	13	aa_alt
	14	gene_funct_annotation

	primers.py adds more columns
	primers.py: 1_intermediate_files/varanalyzer_output.txt --> 1_intermediate_files/primers_output.txt
	16	forward primer                                                                                      
	17	Tm forward                                                                                         
	18	reverse primer
	19	Tm reverse
	
	And here is where this script acts:
	extend-snp-variants-info.py: 1_intermediate_files/primers_output.txt --> 3_workflow_output/snp_variants_report.txt
	20	qual								(from_snp_to_varanalyzer.txt)
	21	ref_count						(from_snp_to_varanalyzer.txt)
	22	alt_count						(from_snp_to_varanalyzer.txt)
	23	af									(from_snp_to_varanalyzer.txt)
	24 dist_to_selectes_pos			(calculated from map_info.txt)
}

ins mode {
	varanalyzer.py: 1_intermediate_files/***_to_varanalyzer.txt --> 1_intermediate_files/varanalyzer_output.txt
	6	mrna_start
	7	mrna_end
	8	strand
	9	gene_model
	10	gene_element
	11	aa_pos
	12	aa_ref
	13	aa_alt
	14	gene_funct_annotation

	In ins mode primers.py add more columns that in snp mode
	primers.py: 1_intermediate_files/varanalyzer_output.txt --> 3_workflow_output/ins_variants_report.txt
	15	forward primer
	16	Tm forward
	17	insertion primer 5'
	18	Tm insertion primer 5'
	19	insertion primer 3'
	20	Tm insertion primer 3'
	21	reverse primer
	22	Tm reverse
}

'''


import argparse
from string import maketrans

# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('--project-name', action="store", dest='project_name', required=True)
parser.add_argument('--variants', action="store", dest='variants_input', required=True)
parser.add_argument('--snp-info', action="store", dest='snp_info_input', required=True)
parser.add_argument('--map-info', action="store", dest='map_info', required=True)
parser.add_argument('--output-file', action="store", dest='output_file', required=True)
parser.add_argument('--region', action="store", dest='region', required=True)

args = parser.parse_args()

project = args.project_name
variants_input = args.variants_input
snp_info_input = args.snp_info_input
map_info_input = args.map_info
output_file = args.output_file

# Create a list that will contain the extended information of every snp
list_extended_info = []

# Create two 2-level list to store the content of the two inut files
in1_array = []
in2_array = []

with open(variants_input, 'r') as in1:
	for line1 in in1:
		if not line1.startswith('@'):
			fields1 = line1.split('\t')
			in1_array.append(fields1)

with open(snp_info_input, 'r') as in2:
	for line2 in in2:
		if not line2.startswith('@'):
			fields2 = line2.split('\t')
			in2_array.append(fields2)

# Calculate distance to predicted position of mutation according to map_info.txt (created with map-mutation.py)
with open(map_info_input) as map_info:
	for map_info_line in map_info:
		map_fields = map_info_line.split('\t')
		if map_fields[0] == '?': # This line stores the information I am retrieving
			# 'selected_position' is the center of the candidate interval calculated by map-mutation.py
			selected_position = (int(map_fields[2]) + int(map_fields[3])) / 2 # map_fields[2] and map_fields[3] are the left and right boundaries of the candidate interval

# Iterate over in1_array and fo reach record:
# 	-iterate over in2_array to grab quality, counts and af (from snp-to-varanalyzer.txt)
#	-calculte distance to predicted position (from map-info.txt)
for line1 in in1_array:
	id1 = line1[1].strip() + ':' + line1[2].strip()
	
	# Get quality, ref_count, alt_count, and allele_frequency
	for line2 in in2_array:
		id2 = line2[1].strip() + ':' + line2[2].strip()
		
		quality, ref_count, alt_count, allele_frequency = '-', '-', '-', '-'
		
		if id1.lower() == id2.lower():
			quality = line2[5].strip()
			ref_count = line2[6].strip()
			alt_count = line2[7].strip()
			allele_frequency = line2[8].strip()
			break
		
	# If, for some reason, the script fails to assign its original quality, counts and freq to a variant,
	# something wrong is happening (every variant should always be in both input files).
	# In this event, stop the workflow by sending 'error' to the workflow
	if quality == '-' or ref_count == '-' or alt_count == '-' or allele_frequency == '-':
			print 'error'
			quit()
	
	# Calculate the distance from the variant position to the selected position
	if str(args.region) == "CR":
		dist_to_selected_position = int(line1[2]) - selected_position
	elif str(args.region) == "total":
		dist_to_selected_position = '-'
	
	extended_info = quality, ref_count, alt_count, allele_frequency, str(dist_to_selected_position)
	line1.append(extended_info)
	
	# Append line to the final list of SNPs	
	list_extended_info.append(line1)

# Create output file
output = open(output_file, 'w')
output.write('@type\tcontig\tposition\tref_base\talt_base\tquality\tref_count\talt_count\talt_allele_freq\tdist_to_selected_pos\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\taa_pos\taa_ref\taa_alt\tgene_funct_annot\tf_primer\ttm_f_primer\tr_primer\ttm_r_primer\n')    

for line in list_extended_info:
	final_fields = line[0], line[1], line[2], line[3], line[4], line[19][0], line[19][1], line[19][2], line[19][3], line[19][4], line[5], line[6], line[7], line[8], line[9], line[10], line[11], line[12], line[13], line[14], line[15], line[16], line[17], line[18]
	for index, field in enumerate(final_fields):
			if index == 0:
				output.write(str(field))
			else:
				output.write('\t' + str(field))	

output.close()


# Add flanking sequences to the candidate SNPs 

fp = open(project + '/1_intermediate_files/gnm_ref_merged/genome.fa', 'r')
input_file = open(output_file, 'r')

# This function parses the information from the fasta files
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

# We create a list of all the contigs with the format [[contig_name, sequence], [contig_name, sequence] ...]
fastalist = list()
for name_contig, seq_contig in read_fasta(fp):
	fastalist.append([name_contig.lower(), seq_contig])

# We retrieve the upstream and downstream sequences of each polymorphism from the fastalist. We also create a new list with the complete lines
final_lines = list()
for line in input_file:
	if not line.startswith('@'):
		sp = line.split()
		chromosome = str(sp[1])
		position = str(sp[2])
		for chrm in fastalist:
			if chrm[0].strip() == '>'+chromosome.strip().lower():
				upstream =  chrm[1][int(position)-51:int(position)-1]
				downstream =  chrm[1][int(position):int(position)+50]
				line = line.strip('\n')
				line = line + '\t' + upstream + '\t' + downstream
				final_lines.append(line)

# We re-write the file with the extended information to the output file from the final_lines list
output_file = open(output_file, 'w')

output_file.write('@type\tcontig\tposition\tref_base\talt_base\tquality\tref_count\talt_count\talt_allele_freq\tdist_to_selected_pos\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\taa_pos\taa_ref\taa_alt\tgene_funct_annot\tf_primer\ttm_f_primer\tr_primer\ttm_r_primer\tupstream\tdownstream\n')    
for line in final_lines:
	output_file.write(line + '\n')

# If the program reaches the end, emit 'success' to let the workflow know that the qual, counts and af
# info was found fo all variants 
print 'success'




