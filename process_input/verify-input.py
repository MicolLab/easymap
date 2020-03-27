
# Multiple fasta files can be provided (one per contig). The program asks for the folder that
# contains the fasta inputs. The same happens with fastq files (f and r / sample 1 and sample2)
#
#
#


import argparse, os, fnmatch


parser = argparse.ArgumentParser()
parser.add_argument('-gnm', action="store", dest='gnm')
parser.add_argument('-ins', action="store", dest='ins_source')
parser.add_argument('-fq', action="store", dest='fq_source')
parser.add_argument('-gff', action="store", dest='gff_source')
parser.add_argument('-ann', action="store", dest='ann_source')
parser.add_argument('-fa_match', action="store", dest='fa_match')
parser.add_argument('-gff_match', action="store", dest='gff_match')
args = parser.parse_args()

gnm_basename = args.gnm
ins_source = args.ins_source
fq_source = args.fq_source
gff_source = args.gff_source
ann_source = args.ann_source
fa_match = args.fa_match
gff_match = args.gff_match


# If gnm argument provided, check fasta file
if gnm_basename != None:
	gnm_result = 0 # 0:pass ; 1:error

	# Create list with all the files in user_data folder
	input_files = sorted(os.listdir('./user_data'))

	# Create a list with only the files that match the basename provided by the user and end in '.fa'
	gnm_files = fnmatch.filter(input_files, '*' + gnm_basename + '.fa') # fnmatch filters a list using a string that accepts wildcards 

	for gnm_file in gnm_files:
		if os.stat('user_data/' + gnm_file).st_size == 0: # Checks whether the file is empty
			gnm_result = 1
		else:
			gnm_contents = open('user_data/' + gnm_file, 'r')
			for index, line in enumerate(gnm_contents):
				if index == 0:
					if not line.startswith('>'):
						gnm_result = 1
				if index == 1:
					break	
			gnm_contents.close()
	
	print gnm_result
	

# If ins argument provided, check fasta file
if ins_source != None:
	ins_result = 0 # 0:pass ; 1:error

	if os.stat(ins_source).st_size == 0: # Checks whether the file is empty
		ins_result = 1
	else:
		ins_contents = open(ins_source, 'r')
		for index, line in enumerate(ins_contents):
			if index == 0:
				if not line.startswith('>'):
					ins_result = 1
			if index == 1:
				break	
		ins_contents.close()
	
	print ins_result


# If fq argument provided, check fastq file(s)
if fq_source != None:
	fq_result = 0
	
	if os.stat(fq_source).st_size == 0:
		fq_result = 1
	else:
		fq_contents = open(fq_source, 'r')
		for index, line in enumerate(fq_contents):
			if index == 0:
				if not line.startswith('@'):
					fq_result = 1
			if index == 2:
				if not line.startswith('+'):
					fq_result = 1
			if index == 3:
				break
		fq_contents.close()
	
	print fq_result
	
	
# If gff argument provided, check gff3 file
if gff_source != None:
	gff_result = 0 # 0:pass ; 1:error
	if os.stat(gff_source).st_size == 0:
		gff_result = 1
	else:
		i = 0 # Counter to move to the tenth line in the file that does not start with "#"
		gff_contents = open(gff_source, 'r')
		for index, line in enumerate(gff_contents):
			if not line.startswith("#"):
				if i == 10:
					fields = line.split('\t')
					try:
						y = int(fields[3]) # This raises an exception if the variable inside int() is not an integer
					except:
						gff_result = 1
						break
					try:
						y = int(fields[4])
					except:
						gff_result = 1
						break
					if fields[6].strip() != '+' and fields[6].strip() != '-' and fields[6].strip() != '0':
						gff_result = 1
				if index == 11:
					break
				i += 1
		gff_contents.close()
	print gff_result


# The format of the annotation file is flexible, so this block is nowadays not called by process-input.sh
# If ann argument provided, check functional annotation file file
if ann_source != None:
	ann_result = 0 # 0:pass ; 1:error
	ann_contents = open(ann_source, 'r')
	
	if os.stat(ann_source).st_size == 0:
		ann_result = 1
	else:		
		for line in ann_contents:
			line_content = line.strip()
		
			if line_content:
				fields = line.split('\t')
			
				if len(fields) < 2:
					ann_result = 1
			
				if len(fields) >= 2:
					if not fields[1].strip():
						ann_result = 1
						print line
	
		ann_contents.close()

	print ann_result
	
	
# If fa_match and gff_match arguments provided, check contigs match
if fa_match != None and gff_match != None:
	
	match_result = 0 
	
	# First, verify that the files are not empty	
	if os.stat(fa_match).st_size == 0 or os.stat(gff_match).st_size == 0:
		match_result = 1
	
	# If files are not empty...
	
	else:
		# Retrieve the name of the contigs in the fasta file and store them in a list
		fa_match_contents = open(fa_match, 'r')
		fa_contigs = []
		
		for fa_line in fa_match_contents:
			if fa_line.startswith('>'):
				fa_fields = fa_line.split('\t')
				fa_contigs.append(fa_fields[0][1:].lower().strip())
		
		fa_match_contents.close()
		
		# Retrieve the name of the contigs (unique) in the gff file and store them in a list
		gff_match_contents = open(gff_match, 'r')
		gff_contigs = []
		
		for gff_line in gff_match_contents:
			gff_fields = gff_line.split('\t')
			contig_name = gff_fields[0].lower().strip()
			if contig_name not in gff_contigs:
				gff_contigs.append(contig_name)
		
		gff_match_contents.close()
		
		# Check if all the contigs in fasta file are also in gff file
		for fa_contig in fa_contigs:
			if fa_contig not in gff_contigs:
				match_result = 2

		# Check if all the contigs in fasta file are also in gff file
		if match_result != 2:
			for gff_contig in gff_contigs:
				if gff_contig not in fa_contigs:
					match_result = 3
	
	print match_result
	
	# 0: pass
	# 1: fa, gff, or both files are empty
	# 2: FASTA file has contigs not present in GFF3 file. Since this is potentially very
	#    dangerous, if exit_code=2, process-input.sh stops.
	# 3: GFF3 file has contigs not present in FASTA file. In this case, just warn the user
	#	  about the presence of extra contigs in the GFF3 file. The analysis will keep going.		 	
	

# Every 'if block' returns only one variable with a numeric value. This numeric value
# is read by the bash wrapper, which decides what to do accordingly
