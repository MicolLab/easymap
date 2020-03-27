
#python ./graphic_output/report.py -variants ./user_projects/project/3_workflow_output/candidate_variants.txt -log ./user_projects/project/2_logs/log.log -output_html ./user_projects/project/3_workflow_output/report.html -project user_projects/project  -mut_type snp -files_dir ./user_projects/project/3_workflow_output/
#python ./graphic_output/report.py -variants ./user_projects/project/3_workflow_output/insertions_output.txt -log ./user_projects/project/2_logs/log.log -output_html ./user_projects/project/3_workflow_output/report.html -project user_projects/project  -mut_type lin -files_dir ./user_projects/project/3_workflow_output/

import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-mut_type', action="store", dest = 'mut_type')
parser.add_argument('-variants', action="store", dest = 'variants')
parser.add_argument('-log', action="store", dest = 'log')
parser.add_argument('-output_html', action="store", dest = 'output_html')
parser.add_argument('-project', action="store", dest = 'project')
parser.add_argument('-files_dir', action="store", dest = 'files_dir')

args = parser.parse_args()

#Log input
input_log = args.log

#VAriants input
input_var = args.variants


#HTML output
output_html = args.output_html
output = open(output_html, 'w')

#Others
mut_type = args.mut_type

#______________________________________________List output '+files_dir+'__________________________________________________
from os import listdir
from os.path import isfile, join

files_dir = args.files_dir

files = [f for f in listdir(files_dir) if isfile(join(files_dir, f))]

#______________________________________________Log info___________________________________________________________
with open(input_log, 'r') as f1:
	for line in f1:
		if line.startswith('Project name'):
			sp = line.split()
			project_name = str(sp[-1])

		if line.startswith('Reference sequence:'):
			sp = line.split()
			ref_files = str(sp[-1])

		if line.startswith('Insertion sequence:'):
			sp = line.split()
			ins_file = str(sp[-1])

		if line.startswith('Library type (test sample):'):
			sp = line.split()
			reads_type = str(sp[-1])

		if line.startswith('Single-end reads (test sample):'):
			sp = line.split()
			reads_s = str(sp[-1])

		if line.startswith('Forward reads (test sample):'):
			sp = line.split()
			reads_f = str(sp[-1])

		if line.startswith('Reverse reads (test sample):'):
			sp = line.split()
			reads_r = str(sp[-1])

		if line.startswith('Library type (control sample):'):
			sp = line.split()
			reads_type_control = str(sp[-1])

		if line.startswith('Single-end reads (control sample):'):
			sp = line.split()
			reads_s_control = str(sp[-1])

		if line.startswith('Forward reads (control sample):'):
			sp = line.split()
			reads_f_control = str(sp[-1])

		if line.startswith('Reverse reads (control sample):'):
			sp = line.split()
			reads_r_control = str(sp[-1])

		if line.startswith('GFF file:'):
			sp = line.split()
			gff_file = str(sp[-1])

		if line.startswith('Annotation file:'):
			sp = line.split()
			if str(sp[-1]) == "n/p": ann_file = "Not provided"
			else: ann_file = str(sp[-1])

		if line.startswith('Data source:'):
			sp = line.split()
			data_source = str(sp[-1])

		if line.startswith('Simulator (sim-recsel.py) command:'):
			try:
				sp = line.split('+')
				selected_position = str(sp[-3])
			except:
				selected_position = "n/p"

		if line.startswith('Type of cross [bc/oc]:'):
			sp = line.split()
			if str(sp[-1]) == 'bc': cross_type = 'backcross'
			if str(sp[-1]) == 'oc': cross_type = 'outcross'

		if line.startswith('Mutant strain [ref/noref]:'):
			sp = line.split()
			if str(sp[-1]) == 'ref': mut_background = 'reference'
			if str(sp[-1]) == 'noref': mut_background = 'non-reference'

		if line.startswith('SNP analysis type [par/f2wt]:'):
			sp = line.split()
			if str(sp[-1]) == 'par': snp_analysis_type = 'parental'
			if str(sp[-1]) == 'f2wt': snp_analysis_type = 'wild type F2'

		if line.startswith('Parental used as control [mut/nomut/np]:'):
			sp = line.split()
			if str(sp[-1]) == 'mut': parental_used_as_control = 'mutant'
			if str(sp[-1]) == 'nomut': parental_used_as_control = 'wild type'

if mut_type == 'snp':
	#SNP mappint control samples
	if snp_analysis_type == 'parental' and parental_used_as_control == 'mutant' : control = ' parental of the mutant strain.'
	if snp_analysis_type == 'parental' and parental_used_as_control == 'wild type' : control = ' wild type parental of the mapping population.'
	if snp_analysis_type == 'wild type F2' : control = ' wild type F2 of the mapping population.'

if data_source == 'sim':
	with open(input_log, 'r') as f1:
		for line in f1:
			if line.startswith('Simulator (sim-mut.py) command:'):
				sp = line.split()
				number_mutations = str(sp[-1].split('+')[0])

			if line.startswith('Simulator (sim-recsel.py) command:'):
				sp = line.split()
				rec_chromosomes = str(sp[-1].split('+')[-1])

			if line.startswith('Simulator (sim-seq.py) command:'):
				sp = line.split()
				read_depth = str(sp[-1].split('+')[0])



if mut_type == 'lin':
	#______________________________________________Get list of insertions_____________________________________________
	insertions_list = list()

	infile = './'+args.project+'/3_workflow_output/sorted_insertions.txt'
	sorted_insertions = open(infile, 'r')

	for line in sorted_insertions:
		if not line.startswith('@'):
			sp = line.split()
			if sp[2].strip() not in insertions_list:
				insertions_list.append(sp[2].strip())

	insertions_pos_list = list()
	with open(input_var) as f:
		for line in f:
			if not line.startswith('@'):
				sp = line.split()



				ins_localizer = str(sp[1].strip().lower() + '-' + sp[2].strip())
				
				ins_localizer_1_1 = str((sp[1]).strip().lower() + '-' + str((int(sp[2]))).strip())
				ins_localizer_1_2 = str((sp[1]).strip().lower() + '-' + str((int(sp[2])+1)).strip())
				ins_localizer_1_3 = str((sp[1]).strip().lower() + '-' + str((int(sp[2])-1)).strip())


				with open(infile) as f2:
					for line in f2:
						if not line.startswith('@'):
							sp2 = line.split()


							ins_localizer_2 = str((sp2[1]).strip().lower() + '-' + str((int(sp2[3]))).strip())

							if ins_localizer_2 == ins_localizer_1_1:
								for insertion in insertions_list:
									if insertion == sp2[2]:
										sublist = [insertion, str(sp2[3]), sp[1]]
										if sublist not in insertions_pos_list:
											insertions_pos_list.append(sublist)
								break
							
							elif ins_localizer_2 == ins_localizer_1_2 or ins_localizer_2 == ins_localizer_1_3 :
								for insertion in insertions_list:
									if insertion == sp2[2]:
										sublist = [insertion, str(sp2[3]), sp[1]]
										if sublist not in insertions_pos_list:
											insertions_pos_list.append(sublist)
								break
#______________________________________________Writting header____________________________________________________

#HTML/CSS stuff
output.write(
'<!DOCTYPE html>' + '\n'
'<html>' + '\n'
'<head>' + '\n'

'	<meta charset="utf-8" />' + '\n'
'	<title>Easymap - report</title>' + '\n'
'	<style>' + '\n'
'		#wrapper {' + '\n'
'			width: 100%;' + '\n'
'			max-width: 1050px;' + '\n'
'			margin: 0 auto;' + '\n'
'			font-family: arial, helvetica;' + '\n'
'		}	' + '\n'
'		.easymap {' + '\n'
'			border: 0;' + '\n'
'			color: rgb(139, 167, 214);' + '\n'
'			background-color: rgb(139, 167, 214);' + '\n'
'			height: 5px;' + '\n'
		
'		}	' + '\n'
'		.img { width: auto; max-width:100% }' + '\n'
'		.result1 { max-width: 864px; }' + '\n'
'		.result2 { max-width: 1000px; }' + '\n'


'		table {border-collapse:collapse; table-layout:fixed;}' + '\n'
'		table td {border:solid 0px #fab;  word-wrap:break-word; vertical-align:top;}' + '\n'
'       tr:hover { background-color: #ededed; }' + '\n'
'		#t {  border: 0px solid red; word-wrap:break-word; table-layout:fixed; line-height: 24px;}' + '\n'
'		#candidates { vertical-align:top; line-height: normal; text-align:left; word-wrap:break-word; line-height: 14px;; }' + '\n'

'	</style>' + '\n'

'</head>' + '\n'
)

#Header and run summary
output.write(
'<body>' + '\n'
'	<div id="wrapper">' + '\n'
'		<hr class="easymap">' + '\n'
'		<h1>' +  project_name + '</h1>' + '\n'
'		<hr class="easymap">' + '\n'

)

#Exp/sim and read files
output.write(
'		<h2>Run summary</h2>' + '\n'
'		<table id="t" >' + '\n'
'		<col width="300">' + '\n'
'		<col width="700">' + '\n'
	)

output.write(

'		<tr>' + '\n'
'			<td> <b>Input genome files ID:</b></td>' + '\n'
'			<td>' + ref_files + '</td>' + '\n'
'		</tr>' + '\n'

	)


#Insertion sequence
if mut_type == 'lin' : 
	output.write(

	'		<tr>' + '\n'
	'			<td> <b>Input insertion file:</b></td>' + '\n'
	'			<td>' + ins_file + '</td>' + '\n'
	'		</tr>' + '\n'

		)

#Reads files
if mut_type == 'lin' and data_source == 'exp': 
	if reads_type == 'pe':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input read files:</b></td>' + '\n'
'			<td>' + reads_f + ', &nbsp;&nbsp;&nbsp;&nbsp;' + reads_r + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type == 'se':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input read files:</b></td>' + '\n'
'			<td>' + reads_s + '</td>' + '\n'
'		</tr>' + '\n'
			)


if mut_type == 'snp' and data_source == 'exp': 
	if reads_type == 'pe':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input test read files:</b></td>' + '\n'
'			<td>' + reads_f + ', &nbsp;&nbsp;&nbsp;&nbsp;' + reads_r + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type == 'se':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input test read files:</b></td>' + '\n'
'			<td>' + reads_s + '</td>' + '\n'
'		</tr>' + '\n'
			)

	if reads_type_control == 'pe':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input control read files:</b></td>' + '\n'
'			<td>' + reads_f_control + ', &nbsp;&nbsp;&nbsp;&nbsp;' + reads_r_control + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif reads_type_control == 'se':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Input control read files:</b></td>' + '\n'
'			<td>' + reads_s_control + '</td>' + '\n'
'		</tr>' + '\n'
			)


if data_source == 'sim':
	if mut_type == 'lin':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Number of insertions (simulation):</b></td>' + '\n'
'			<td>' + number_mutations + '</td>' + '\n'
'		</tr>' + '\n'
'		<tr>' + '\n'
'			<td> <b>Read depth (simulation):</b></td>' + '\n'
'			<td>' + read_depth + '</td>' + '\n'
'		</tr>' + '\n'
		)

	elif mut_type == 'snp':
		output.write(
'		<tr>' + '\n'
'			<td> <b>Selected causal mutation (simulation):</b></td>' + '\n'
'			<td>' + selected_position.split(',')[1] + ', Contig ' + selected_position.split(',')[0] + '</td>' + '\n'
'		</tr>' + '\n'
'		<tr>' + '\n'
'			<td> <b>Number of mutations (simulation):</b></td>' + '\n'
'			<td>' + number_mutations + '</td>' + '\n'
'		</tr>' + '\n'
'		<tr>' + '\n'
'			<td> <b>Number of recombinant chromosomes (simulation):</b></td>' + '\n'
'			<td>' + rec_chromosomes + '</td>' + '\n'
'		</tr>' + '\n'
'		<tr>' + '\n'
'			<td> <b>Read depth (simulation):</b></td>' + '\n'
'			<td>' + read_depth + 'x</td>' + '\n'
'		</tr>' + '\n'
		)
		

if mut_type == 'snp': 
	output.write(

'		<tr>' + '\n'
'			<td> <b>Experimental design: </b></td>' + '\n'
'			<td> Mutation in ' + mut_background + ' genetic background. A ' + cross_type + ' was performed to obtain the mapping population. The control sample is the ' + control +'</td>' + '\n'
'		</tr>' + '\n'

		)

#Gff and ann '+files_dir+'
output.write(

'		<tr>' + '\n'
'			<td> <b>Structural annotation file:</b></td>' + '\n'
'			<td>' + gff_file + '</td>' + '\n'
'		</tr>' + '\n'

'		<tr>' + '\n'
'			<td> <b>Functional annotation file:</b></td>' + '\n'
'			<td>' + ann_file + '</td>' + '\n'
'		</tr>' + '\n'
'		</table>' + '\n'

)

#Link to log file
output.write(
'		<br>' + '\n'
'		<a href="../2_logs/log.log" target="_blank">Click to see log file</a>' + '\n'
'		<br><br>' + '\n'
'		<hr class="easymap">' + '\n'
)


#Input data quality assessment
output.write(
'		<h2>Input data quality assessment</h2>' + '\n'
	)
if data_source == 'exp':

	if mut_type == 'lin': 
		if reads_type == 'pe':
			output.write(
	'		<b>Paired end reads quality assessment<br></b>' + '\n'
	'		<p>Forward reads<br></p>' + '\n'
	'		<center> <img class="img" src="paired-end-problem-forward-reads-qual-stats.png"></center> ' + '\n'
	'		<p>Reverse reads<br></p>' + '\n'
	'		<center><img class="img" src="paired-end-problem-reverse-reads-qual-stats.png"></center> ' + '\n'
			)

		if reads_type == 'se':
			output.write(
	'		<b>Single end reads quality assessment<br></b>' + '\n'
	'		<center><img class="img" src="single-end-problem-reads-qual-stats.png" width="500"></center> ' + '\n'
			)


	if mut_type == 'snp': 
		if reads_type == 'pe':
			output.write(
	'		<b>Test reads quality assessment<br></b>' + '\n'
	'		<p>Forward reads<br></p>' + '\n'
	'		<center> <img class="img" src="'+'paired-end-problem-forward-reads-qual-stats.png" > </center> ' + '\n'
	'		<p>Reverse reads<br></p>' + '\n'
	'		<center> <img class="img" src="paired-end-problem-reverse-reads-qual-stats.png"  > </center> ' + '\n'
			)

		if reads_type == 'se':
			output.write(
	'		<b>Test reads quality assessment<br></b>' + '\n'
	'		<center> <img class="img" src="single-end-problem-reads-qual-stats.png" width="500" > </center> ' + '\n'
			)

		if reads_type_control == 'pe':
			output.write(
	'		<b>Control reads quality assessment<br></b>' + '\n'
	'		<p>Forward reads<br></p>' + '\n'
	'		<center> <img class="img" src="'+'paired-end-control-forward-reads-qual-stats.png" > </center> ' + '\n'
	'		<p>Reverse reads<br></p>' + '\n'
	'		<center> <img class="img" src="paired-end-control-reverse-reads-qual-stats.png"  > </center> ' + '\n'
			)

		if reads_type_control == 'se':
			output.write(
	'		<b>Control reads quality assessment<br></b>' + '\n'
	'		<center> <img class="img" src="single-end-control-reads-qual-stats.png" width="500" > </center> ' + '\n'
			)

#Read depth distribution graphics
if mut_type == 'snp': 
	output.write(
	'		<b>Test and control samples read depth distribution<br></b>' + '\n'
	'		<p>Test sample<br></p>' + '\n'
	'		<center> <img class="img" src="frequence_depth_alignment_distribution_sample.png" > </center> ' + '\n'
	'		<p>Control sample<br></p>' + '\n'
	'		<center> <img class="img" src="frequence_depth_alignment_distribution_control.png" > </center> ' + '\n'
	'		<hr class="easymap">' + '\n'
		)

if mut_type == 'lin': 
	output.write(
	'		<b>Test sample read depth distribution<br></b>' + '\n'
	'		<center> <img class="img" src="frequence_depth_alignment_distribution_sample.png" > </center> ' + '\n'
	'		<hr class="easymap">' + '\n'
		)


#__________________________________LIN cartographic report________________________________________________________________
if mut_type == 'lin': 
	#Genome overview
	for f in sorted(files): 
			if 'insertions_overview' in str(f):
				output.write(
				'		<h2>Mapping analysis overview</h2>' + '\n'					
				'		<h3>Genomic overview</h3>' + '\n'
				'		<center> <img class="img" src="'  +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" >  </center>' + '\n'
				)

	output.write(
		#Insertions summary
		'		<h3>Insertions summary</h3>' + '\n'
		)


	#first we check that there are insertions 
	n_candidates = 0
	with open(input_var) as candidates:
		for line in candidates:
			if not line.startswith('@'):
				n_candidates = n_candidates + 1

	#then if the number of insertions is > 0, we write the candidates in a table

	if n_candidates == 0:
		output.write(
		'		<center> <p>No insertions were found<br></p> <center/>' + '\n'
		)


	if n_candidates > 0: 
		output.write(
			#Table
			'		<table id="candidates" border="0" align="center" cellpadding="10">' + '\n'
			'		  <tr>' + '\n'
			'		    <th>Ins</th>' + '\n'
			'		    <th>Contig</th>' + '\n'
			'		    <th>Position</th>' + '\n'
			'		    <th>Gene (gene element)</th>' + '\n'
			'		    <th>Wt amino acids</th>' + '\n'
			'		  </tr>' + '\n'
			)

		with open(input_var) as candidates:
			variants_list=list()
			for line in candidates:
				if not line.startswith('@'):
					sp = line.split('\t')
					contig = str(sp[1]).strip()
					position = str(sp[2]).strip()
					aminoacid = str(sp[11]).strip() 
					primer_f = str(sp[15]).strip()
					primer_r = str(sp[21]).strip()
					primer_5 = str(sp[17]).strip()
					primer_3 = str(sp[19]).strip()
					upstream = str(sp[23]).strip()
					downstream = str(sp[24]).strip()
					if str(sp[9]).strip() !='-':
						gene = str(sp[9]).strip() + ' (' + str(sp[10]).strip() + ')'
					else:
						gene = '-'
					if str(sp[14]).strip() != '-':
						annotation = str(sp[14]).strip()
					else:
						annotation = 'Functional annotation not available'

					for i in insertions_pos_list:
						if str(i[2]).lower().strip() == contig.lower():
							if int(i[1]) == int(position) or (int(i[1])+1) == int(position) or (int(i[1])-1) == int(position):
								ins = str(i[0])

					variants_list.append([ins, contig, position, gene, aminoacid, primer_f, primer_r, primer_5, primer_3, upstream, downstream, annotation])

					output.write(
					'		  <tr>' + '\n'
					'		    <td>'+ins+'</th>' + '\n'
					'		    <td>'+contig+'</th>' + '\n'
					'		    <td>'+position+'</th>' + '\n'
					'		    <td>'+gene+'</th>' + '\n'
					'		    <td>'+aminoacid+'</th>' + '\n'
					'		  </tr>' + '\n'
						)

			output.write(
				'	</table>' + '\n'
				)

		#Link to variants file
		output.write(
		'		<br>' + '\n'
		'		<br><a href="insertions_output.txt" target="_blank">Click to see extended information</a>' + '\n'
		'		<br><br>' + '\n'

		)

		past_ins = None
		#Insertions
		for ins in variants_list:
			if past_ins == ins[0]:
				past_ins = ins[0]

			if past_ins != ins[0]:
				for f in sorted(files): 
					if ('_ins_' + ins[0] + '.png') in str(f):

						output.write(
						'		<hr class="easymap">' + '\n'
						'		<h2> Insertion   ' +  str(ins[0]) +'</h2>' + '\n'
						'		<center> <img class="img" src="'  +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" >  </center>' + '\n'
						)

				for f in sorted(files):
					if '_lin_' + ins[0] + '_gene' in str(f):
						fname = str(f).split('_')[5]
						spf=fname.split('.')
						gene1=""
						for b in spf[:-1] : gene1=gene1+b+"."
						gene=gene1[:-1]
						#DEPRECATED: gene = str(f).split('_')[5].split('.')[0]+'.'+str(f).split('_')[5].split('.')[1]
						for i in variants_list:
							if gene in str(i[3]):
								output.write(
								'		<h3>' + gene + '</h3>' + '\n'
								'		<center> <img class="img" src="'  +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" >  </center>' + '\n'
								'		<table id="t">' + '\n'
								'		<col width="300">' + '\n'
								'		<col width="700">' + '\n'

								'		<tr>' + '\n'
								'			<td> <b>Insertion position:</b></td>' + '\n'
								'			<td style="font-family:Lucida Console, monospace">' + str(ins[2]) + ', ' + str(ins[1]) + '</td>' + '\n'
								'		</tr>' + '\n'

								)

								if ann_file != "Not provided": 
									output.write(
									'		<tr>' + '\n'
									'			<td> <b>Functional annotation:</b></td>' + '\n'
									'			<td>' + str(i[11]) + '</td>' + '\n'
									'		</tr>' + '\n'
									)

								output.write(

								'		<tr>' + '\n'
								'			<td> <b>Forward primer:</b></td>' + '\n'
								'			<td style="font-family:Lucida Console, monospace">' + str(i[5]) + '</td>' + '\n'
								'		</tr>' + '\n'

								'		<tr>' + '\n'
								'			<td> <b>Reverse primer:</b></td>' + '\n'
								'			<td style="font-family:Lucida Console, monospace">' + str(i[6]) + '</td>' + '\n'
								'		</tr>' + '\n'

								'		<tr>' + '\n'
								'			<td> <b>Insertion 5 primer:</b></td>' + '\n'
								'			<td style="font-family:Lucida Console, monospace">' + str(i[7]) + '</td>' + '\n'
								'		</tr>' + '\n'

								'		<tr>' + '\n'
								'			<td> <b>Insertion 3 primer:</b></td>' + '\n'
								'			<td style="font-family:Lucida Console, monospace">' + str(i[8]) + '</td>' + '\n'
								'		</tr>' + '\n'


								'		<tr>' + '\n'
								'			<td> <b>Flanking sequences:</b></td>' + '\n'
								'			<td style="font-family:Lucida Console, monospace">' + str(i[9])[15:] + '<font style="font-family:Lucida Console, monospace" color="red">[Ins' + str(ins[0]) + ']</font>' + str(i[10])[0:35] + '</td>' + '\n'
								'		</tr>' + '\n'


								'		</table>' + '\n'
								'		<br>' + '\n'

								)

					switch = "no"
					for f in sorted(files):
						if '_lin_' + ins[0] + '_gene' in str(f):
							switch = "yes"

					if switch == "no": 
						output.write(
						'		<table id="t">' + '\n'
						'		<col width="300">' + '\n'
						'		<col width="700">' + '\n'
						'		<tr>' + '\n'
						'			<td> <b>Insertion position:</b></td>' + '\n'
						'			<td style="font-family:Lucida Console, monospace">' + str(ins[2]) + ', ' + str(ins[1]) + '</td>' + '\n'
						'		</tr>' + '\n'
						'		</table>' + '\n'
						'		<br>' + '\n'
						)
						break

				past_ins = ins[0]


		for f in sorted(files):
			if 'gene_plot_lin' in str(f):
				output.write(
				'		<left> <img class="img" src="gene_legend_ins.png" align="middle" >  </left>' + '\n'
				'		<br>' + '\n'
				)
				break

		#Link to images 
		output.write(
		'		<br><a href="report_images.zip" target="_blank">Click to download all image files</a>' + '\n'
		'		<hr class="easymap">' + '\n'
		)


#__________________________________SNP cartographic report________________________________________________________________
if mut_type == 'snp':

	#Chromosomes FA vs POS
	#Mapping
	output.write(
	'		<h2>Mapping analysis overview</h2>' + '\n'
	'		<p>All input contigs are displayed, with all the polymorphisms used for mapping the causal mutation and their linear description (SMA, boost values or AF difference between mapping populations). The candidate region determined by the analysis is highlighted.</p>' + '\n'
		) 
	for f in sorted(files):
		if 'mapping' in str(f):
			output.write(
			'		<left> <img class="img" src="' +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" > </left>' + '\n'
			)

	#Candidates
	#selected chromosome:

	for f in sorted(files):
		if 'candidates' in str(f) and 'zoom' in str(f):
			sp = str(f).split('_')
			selected_chromosome = str(sp[1]).strip()

	output.write(
	'		<hr class="easymap">' + '\n'
	'		<h2>Candidate polymorphisms overview</h2>' + '\n'
	'		<p>The chromosome containing the candidate region is displayed along with the total polymorphisms in the test sample, highlighting the window that contains the candidate variants and representing the selected position as a dashed line. </p>' + '\n'
		) 

	#Candidates 
	for f in sorted(files):
		if 'candidates_' in str(f) and selected_chromosome.lower().strip() == ((str(f).split('candidates_')[1]).split('.png')[0]).lower() and 'zoom' not in str(f):
			output.write(
			'		<left> <img class="img" src="'  +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" > </left>' + '\n'
			)
	#Candidates zoom
	for f in sorted(files):
		if 'candidates' in str(f) and 'zoom' in str(f):
			output.write(
			'		<p>Zooming-in the window containing the candidate polymorphisms, including only typical EMS changes. </p>' + '\n'
			'		<left> <img class="img" src="'  +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" > </left>' + '\n'
			)
			

	#Legend
	##PathInterfaceToImages = '../user_projects/' + project_name + '/3_workflow_output/'
	#ReplacementImage = 'onerror="this.src=\'./legend.png'"'
	output.write(
	'		<h3>Legend</h3>' + '\n'
	'		<left> <img class="img" src="legend.png" align="middle" >  </left>' + '\n'
	)


	#Candidates table:
	output.write(
	'		<hr class="easymap">' + '\n'
	'		<h2>Candidate region analysis</h2>' + '\n'
	)

	#first we check that there are candidates 
	n_candidates = 0
	with open(input_var) as candidates:
		for line in candidates:
			if not line.startswith('@'):
				sp = line.split()
				if sp[10].strip() != "nh":
					n_candidates = n_candidates + 1

	#then if the number of candidates is > 0, we write the candidates in a table

	if n_candidates == 0:
		output.write(
		'		<center> <p>No candidate variants found</p> <center/>' + '\n'
		'		<p>Click to see a list of <a href="candidate_variants_total.txt" target="_blank">all variants</a>. </p>' + '\n'
		
		)


	if n_candidates > 0: 
		output.write(
			#Table
			'		<table id="candidates" border="0" align="center" cellpadding="10">' + '\n'
			'		  <tr>' + '\n'
			'		    <th>ID</th>' + '\n'
			'		    <th>Contig</th>' + '\n'
			'		    <th>Position</th>' + '\n'
			'		    <th>AF</th>' + '\n'
			'		    <th>DTP</th>' + '\n'
			'		    <th>Nucleotide (Ref/Alt)</th>' + '\n'
			'		    <th>Gene (gene element)</th>' + '\n'
			'		    <th>Amino acid (Ref/Alt)</th>' + '\n'
			'		  </tr>' + '\n'
			)

		with open(input_var) as candidates:
			i=1
			variants_list=list()
			for line in candidates:
				if not line.startswith('@'):
					sp = line.split('\t')
					contig = str(sp[1]).strip()
					position = str(sp[2]).strip()
					AF = str(sp[8]).strip()
					DTP = str(sp[9]).strip()
					nucleotide = str(sp[3]).strip() + ' &rarr; ' + str(sp[4]).strip()
					alt_nt = str(sp[4]).strip()
					aminoacid = str(sp[17]).strip() + ' &rarr; ' + str(sp[18]).strip()
					if aminoacid == '- &rarr; -': aminoacid = '-'
					primer_f = str(sp[20]).strip()
					primer_r = str(sp[22]).strip()
					upstream = str(sp[24]).strip()
					downstream = str(sp[25]).strip()
					if str(sp[14]).strip() != '-':
						gene = str(sp[14]).strip() + ' (' + str(sp[15]).strip() + ')'
					else:
						gene = '-'
					if str(sp[19]).strip() != '-':
						annotation = str(sp[19]).strip() 
					else:
						annotation = ' Functional annotation not available'

					variants_list.append([str(i), contig, position, AF, DTP, nucleotide, gene, aminoacid, primer_f, primer_r, upstream, downstream, annotation, alt_nt])

					output.write(
					'		  <tr>' + '\n'
					'		    <td>'+str(i)+'</th>' + '\n'
					'		    <td>'+contig+'</th>' + '\n'
					'		    <td>'+position+'</th>' + '\n'
					'		    <td>'+AF+'</th>' + '\n'
					'		    <td>'+DTP+'</th>' + '\n'
					'		    <td>'+nucleotide+'</th>' + '\n'
					'		    <td>'+gene+'</th>' + '\n'
					'		    <td>'+aminoacid+'</th>' + '\n'
					'		  </tr>' + '\n'
						)

					i=i+1

			output.write(
				'	</table>' + '\n'
				'	<br><p>Click to see <a href="candidate_variants.txt" target="_blank">extended information</a>  or  <a href="candidate_variants_total.txt" target="_blank"> all variants</a>. </p>' + '\n'
				)

	output.write(
	'		<hr class="easymap">' + '\n'
	)


	#Candidate SNPs
	if n_candidates > 0:
		output.write(
		'		<h2>Candidate variants</h2>' + '\n'
		'		<p>This section contains a list of the candidate mutations affecting gene open reading frames.</p>' + '\n'
			) 
		for var in variants_list:
			gene_name = var[6].split(' (')[0]
			for f in sorted(files):
				if 'gene_plot_snp' in str(f) and gene_name in str(f) and var[2] in str(f):
					output.write(
					'		<h3>ID ' + str(var[0]) + ': ' + gene_name + '</h3>' + '\n'
					'		<center> <img class="img" src="'  +  str(f).split('3_workflow_output/')[-1]  + ' " align="middle" >  </center>' + '\n'
					'		<table id="t">' + '\n'
					'		<col width="300">' + '\n'
					'		<col width="700">' + '\n'
					)

					output.write(
					'		<tr>' + '\n'
					'			<td> <b>Position:</b></td>' + '\n'
					'			<td style="font-family:Lucida Console, monospace">' + var[2] + ', '+ var[1] + '</td>' + '\n'
					'		</tr>' + '\n'

					'		<tr>' + '\n'
					'			<td> <b>Allelic frequence:</b></td>' + '\n'
					'			<td style="font-family:Lucida Console, monospace">' + var[3] + '</td>' + '\n'
					'		</tr>' + '\n'

					'		<tr>' + '\n'
					'			<td> <b>Distance to peak (DTP):</b></td>' + '\n'
					'			<td style="font-family:Lucida Console, monospace">' + var[4] + '</td>' + '\n'
					'		</tr>' + '\n'
					)

					if ann_file != "Not provided":
						output.write(
						'		<tr>' + '\n'
						'			<td> <b>Functional annotation:</b></td>' + '\n'
						'			<td>' + var[12] + '</td>' + '\n'
						'		</tr>' + '\n'
						)

					output.write(
					'		<tr>' + '\n'
					'			<td> <b>Forward primer:</b></td>' + '\n'
					'			<td style="font-family:Lucida Console, monospace">' + var[8] + '</td>' + '\n'
					'		</tr>' + '\n'

					'		<tr>' + '\n'
					'			<td> <b>Reverse primer:</b></td>' + '\n'
					'			<td style="font-family:Lucida Console, monospace">' + var[9] + '</td>' + '\n'
					'		</tr>' + '\n'

					'		<tr>' + '\n'
					'			<td> <b>Flanking sequences:</b></td>' + '\n'
					'			<td style="font-family:Lucida Console, monospace">' + var[10][10:] + '<font style="font-family:Lucida Console, monospace" color="red">' + var[13] + '</font>' + var[11][0:40] + '</td>' + '\n'
					'		</tr>' + '\n'

					'		</table>' + '\n'
					'		<br>' + '\n'
					'		<br>' + '\n'

					)

		output.write(
		'		<left> <img class="img" src="gene_legend_snp.png" align="middle" >  </left>' + '\n'
		'		<br>' + '\n'
		)

	#Link to images 
	output.write(
	'		<br><a href="report_images.zip" target="_blank">Click to download all image files</a>' + '\n'
	)

output.write('</div>')
output.close()

'''


'		<tr>' + '\n'
'			<td> <b>Upstream sequence:</b></td>' + '\n'
'			<td style="font-family:Lucida Console, monospace">' + var[10] + '</td>' + '\n'
'		</tr>' + '\n'

'		<tr>' + '\n'
'			<td> <b>Downstream sequence:</b></td>' + '\n'
'			<td style="font-family:Lucida Console, monospace">' + var[11] + '</td>' + '\n'
'		</tr>' + '\n'

'		<tr>' + '\n'
'			<td> <b>Upstream sequence:</b></td>' + '\n'
'			<td style="font-family:Lucida Console, monospace">' +  str(i[9]) + '</td>' + '\n'
'		</tr>' + '\n'

'		<tr>' + '\n'
'			<td> <b>Downstream sequence:</b></td>' + '\n'
'			<td style="font-family:Lucida Console, monospace">' +  str(i[10]) + '</td>' + '\n'
'		</tr>' + '\n'
'''





# Create a copy (rfeed.html) of the report.html file where the paths to all images and links to server locations
# have been replaced. report.html is intended to be viewed by someone using the command-line interface. rfeed.html
# is used to feed view-report.php in the web-interface with the correct path to resources (images, files, links) 

htmlFile2 = open('user_projects/' + project_name + '/3_workflow_output/rfeed.html', 'w')
#htmlFile2 = open('user_projects/project/3_workflow_output/rfeed.html', 'w')

with open(output_html, 'r') as htmlFile1:
		for line in htmlFile1:
			 line = line.replace('<img class="img" src="', '<img class="img" src="../user_projects/' + project_name + '/3_workflow_output/')
			 line = line.replace('a href="', 'a href="../user_projects/' + project_name + '/3_workflow_output/')
			 htmlFile2.write(line)

htmlFile2.close()


# TO DO:
# -Show RD graphs in report
# -Validate HTML
