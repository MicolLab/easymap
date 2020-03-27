# This script contains the functions used for drawing the programs output images. The functions are called from graphic-output.py when they are needed. 

import argparse, math
from PIL import Image, ImageDraw, ImageFont, ImageOps
from StringIO import StringIO

#Common arguments
parser = argparse.ArgumentParser()
parser.add_argument('-my_mut', action="store", dest = 'my_mut')			#snp or lin
parser.add_argument('-m', action="store", dest = 'mode', default = 'P')
parser.add_argument('-pname', action="store", dest='project_name')
parser.add_argument('-cross', action="store", dest='my_cross')
parser.add_argument('-snp_analysis_type', action="store", dest='my_snp_analysis_type')
parser.add_argument('-gff', action="store", dest = 'gff')				#Genome feature file
parser.add_argument('-iva', action="store", dest = 'input_va')	 		#Output de varanalyzer
parser.add_argument('-rrl', action="store", dest = 'rrl') 				#Regulatory region lenght
parser.add_argument('-f', action="store", dest = 'output_html')

#Arguments for point mutation mapping graphic output
parser.add_argument('-asnp', action="store", dest = 'input_snp')		
parser.add_argument('-bsnp', action="store", dest = 'input_f_snp')		#Fasta genome input
parser.add_argument('-interval_width', action="store", dest = 'interval_width')

#Arguments for large insertions mapping graphic output
parser.add_argument('-a', action="store", dest = 'input')		
parser.add_argument('-b', action="store", dest = 'input_f')				#Fasta genome input
parser.add_argument('-ins_pos', action="store", dest = 'ins_pos')

args = parser.parse_args()
project = args.project_name

def red(p):
	if len(str(p)) <= 3:
		r = p

	elif len(str(p)) == 4:
		r = str(p)[:1]  + ' kb'

	elif len(str(p)) == 5:
		r =  str(p)[:2]  + ' kb'

	elif len(str(p)) == 6:
		r =  '0' + '.' + str(p)[:2] + ' Mb'

	elif len(str(p)) == 7:
		r = str(p)[:1] + '.' + str(p)[1:3] + ' Mb'

	elif len(str(p)) == 8:
		r = str(p)[:2] + ' Mb'

	elif len(str(p)) == 9:
		r = str(p)[:3] + ' Mb'
	return r; 

args = parser.parse_args()


#############################################################################################################
#																											#
# 							SNP - Alelic frequence VS Chromosome position									#
#																											#
#############################################################################################################
def fa_vs_pos():
	#Input 1
	input1 = args.input_snp
	f1 = open(input1, 'r')
	lines = f1.readlines()	

	#Input 2
	input2 = args.input_f_snp
	f2 = open(input2, 'r')
	lines_f = f2.readlines()	
	contig_source = args.input_f_snp


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


	# Read contig fasta file
	contig_lengths = list()

	with open(contig_source) as fp:
		fastalist = list()
		for name_contig, seq_contig in read_fasta(fp):
			innerlist = list()
			innerlist.append(name_contig.strip('>'))
			innerlist.append(len(seq_contig))
			fastalist.append(innerlist)
			contig_lengths.append(len(seq_contig))

	max_contig_len = 0
	for i in contig_lengths:
		if int(i) > max_contig_len:
			max_contig_len = int(i)

	#FA vs POS graphs 
	for i in fastalist:
		if int(i[1]) > 2000000: 
			wide=int(880*float(i[1])/max_contig_len) + 120					#<-------------------------------------------------------------------------------- SET IMAGE SIZE
			height=500
			im = Image.new("RGB", (1000, int(height)), (255,255,255))		# Esta puesto para da imagenes de tamano fijo (manteniendo la escala). Para tener ancho variable: im = Image.new("RGB", (wide, int(height)), (255,255,255))
			draw = ImageDraw.Draw(im)
			
			#get fonts from foler 'fonts'
			fnt2 = ImageFont.truetype('fonts/VeraMono.ttf', 14)

			r = red(int(i[1]))

			if 'Mb' in r:
				max_graph_x = int(i[1]) + 10000

			elif 'kb' in r: 
				max_graph_x = int(i[1])

			#Scaling factors
			scaling_factor_x = (max_graph_x)/(wide - 120) 								#nts/pixel        
			scaling_factor_y = (1.001/(63/100.0*height))  								#af/pixels

			#Candidate region 
			if args.my_mut == 'snp':
				binput = open(project + '/1_intermediate_files/map_info.txt', 'r')
				
				#Retrieving candidate region coordinates
				for line in binput:
					if line.startswith('?'):
						sp = line.split()						
						chromosome = sp[1].strip().lower()
						chromosome_candidate = chromosome
						if chromosome == i[0].lower():
							if int(sp[2]) > 0 :
								cr_start = int(sp[2])  
							else:
								cr_start = 0
							if  int(sp[3]) < int(i[1]) :
								cr_end = int(sp[3]) 
							else:
								cr_end = int(i[1])

				#Drawing candidate region:
				if chromosome_candidate == i[0].lower():
					cr_start_im = int(cr_start/scaling_factor_x) + 70
					cr_end_im = int(cr_end/scaling_factor_x) + 70
					draw.rectangle( [cr_start_im, int(15/100.0*height), cr_end_im, int(80/100.0*height)], fill=(249, 222, 252) )

			# af_candidates: framing the candidate region
			if args.my_mut == 'af_candidates':
				binput = open(project + '/1_intermediate_files/map_info.txt', 'r')
				#Retrieving candidate region coordinates
				for line in binput:
					if line.startswith('?'):
						sp = line.split()						
						chromosome = sp[1].strip().lower()
						chromosome_candidate = chromosome
						if chromosome == i[0].lower():
							cr_start_raw = int(sp[2])
							cr_end_raw = int(sp[3])
							if int(sp[2]) > 0 :
								cr_start = int(sp[2])  
							else:
								cr_start = 0
							if  int(sp[3]) < int(i[1]) :
								cr_end = int(sp[3]) 
							else:
								cr_end = int(i[1])

				if chromosome_candidate == i[0].lower():
					#Drawing a frame for the candidate region:
					cr_start_im = int(cr_start/scaling_factor_x) + 70
					cr_end_im = int(cr_end/scaling_factor_x) + 70
					fa_img_08 = int(80/100.0*height) - int(0.8/scaling_factor_y) - 1
					draw.rectangle( [cr_start_im, int(15/100.0*height)+1, cr_end_im+1, fa_img_08], fill=(249, 222, 252), outline=(112, 112, 112) )

					#Drawing a dotted line in the frame
					cr_middle = ((cr_start_raw + cr_end_raw)/2)/scaling_factor_x + 70
					h = int(16/100.0*height)
					while h in range(int(15/100.0*height), fa_img_08):
						draw.line((cr_middle, h) + (cr_middle, h+5), fill=(255, 0, 0, 0), width=1)
						h = h + 10
					
			#snps
			r, g, b = 31, 120, 180
			if args.my_mut == 'af_control':
				r, g, b = 245, 120, 44

			for l, line in enumerate(lines):
				sp = line.split()
				if i[0].lower() == sp[0].lower() :
					fa = float(sp[6])/(float(sp[6])+float(sp[5]))
					fa_img = int(80/100.0*height) - int(fa/scaling_factor_y) - 1
					pos_img = int(int(sp[1])/scaling_factor_x) + 70
					draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(r, g, b))

			if args.my_snp_analysis_type == 'f2wt' and args.my_mut == 'snp':
				

				#Filler variants 
				problem_var = open(project + '/1_intermediate_files/filler_variants.va', 'r')
				for line in problem_var:
					sp = line.split()
					if i[0].lower() == sp[0].lower():
						#f2 mut
						fa = float(sp[8])/(float(sp[8])+float(sp[7]))
						fa_img = int(80/100.0*height) - int(fa/scaling_factor_y)
						pos_img = int(int(sp[1])/scaling_factor_x) + int(70)
						draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(237, 194, 168)) 
						#f2wt snps
						fa = float(sp[6])/(float(sp[6])+float(sp[5]))
						fa_img = int(80/100.0*height) - int(fa/scaling_factor_y) - 1
						pos_img = int(int(sp[1])/scaling_factor_x) + 70
						draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(167, 190, 206))


				'''
				# Problem variants
				problem_var = open(project + '/1_intermediate_files/F2_filtered.va', 'r')
				for line in problem_var:
					sp = line.split()
					if i[0].lower() == sp[0].lower():
						fa = float(sp[6])/(float(sp[6])+float(sp[5]))
						fa_img = int(80/100.0*height) - int(fa/scaling_factor_y)
						pos_img = int(int(sp[1])/scaling_factor_x) + int(70)
						draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(167, 190, 206)) 

				#Control variants
				problem_var = open(project + '/1_intermediate_files/control_filtered.va', 'r')
				for line in problem_var:
					sp = line.split()
					if i[0].lower() == sp[0].lower():
						fa = float(sp[6])/(float(sp[6])+float(sp[5]))
						fa_img = int(80/100.0*height) - int(fa/scaling_factor_y)
						pos_img = int(int(sp[1])/scaling_factor_x) + int(70)
						draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(237, 194, 168)) 
				'''


				#Mapping variants
				for l, line in enumerate(lines):
					sp = line.split()
					if i[0].lower() == sp[0].lower():
						#f2 mut
						fa = float(sp[8])/(float(sp[8])+float(sp[7]))
						fa_img = int(80/100.0*height) - int(fa/scaling_factor_y)
						pos_img = int(int(sp[1])/scaling_factor_x) + int(70)
						draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(245, 120, 44)) 
						#f2wt snps
						fa = float(sp[6])/(float(sp[6])+float(sp[5]))
						fa_img = int(80/100.0*height) - int(fa/scaling_factor_y) - 1
						pos_img = int(int(sp[1])/scaling_factor_x) + 70
						draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(31, 120, 180))


			my_cross = str(args.my_cross)
			#Boost / mm 																						
			if args.my_mut == 'snp':
				binput = open(project + '/1_intermediate_files/map_info.txt', 'r')
				blines = binput.readlines()

				#Boost line
				if my_cross == 'oc' :
					for b, bline in enumerate(blines):
						sp = bline.split()
						if bline.startswith('!'):
							boost_max = float(sp[3])
					for b, bline in enumerate(blines):
						sp = bline.split()					
						if bline.startswith('@') and sp[4].lower().strip('>') == i[0].lower():
							boost_value = float(sp[3].strip())/boost_max
							boost_value_img = int(80/100.0*height) - int(boost_value/scaling_factor_y )

							window_position = int(sp[1])
							window_position_img = int(window_position/scaling_factor_x) + 70

							try:
								draw.line(((window_position_img, boost_value_img) + (window_position_img_2, boost_value_img_2)), fill=(255, 0, 0, 0), width=1)	
								window_position_img_2 = window_position_img
								boost_value_img_2 = boost_value_img

							except:
								window_position_img_2 = window_position_img
								boost_value_img_2 = boost_value_img

					window_position_img = None 
					boost_value_img = None
					window_position_img_2 = None 
					boost_value_img_2 = None

				#MM line
				if my_cross == 'oc' :
					for b, bline in enumerate(blines):
						sp = bline.split()					
						if bline.startswith('@') and sp[4].lower().strip('>') == i[0].lower():
							mm_value = float(sp[2].strip())
							mm_value_img = int(80/100.0*height) - int(mm_value/scaling_factor_y )
							window_position = int(sp[1])
							window_position_img = int(window_position/scaling_factor_x) + 70
							try:
								draw.line(((window_position_img, mm_value_img) + (window_position_img_2, mm_value_img_2)), fill=(46, 255, 0), width=1)	
								window_position_img_2 = window_position_img
								mm_value_img_2 = mm_value_img
							except:
								window_position_img_2 = window_position_img
								mm_value_img_2 = mm_value_img

				if my_cross == 'bc' :
					r, g, bl = 46, 255, 0
					if args.my_snp_analysis_type == 'f2wt':
						r, g, bl = 255, 0, 255
					for b, bline in enumerate(blines):
						sp = bline.split()					
						if bline.startswith('@') and sp[3].lower().strip('>') == i[0].lower():
							mm_value = float(sp[2].strip())
							mm_value_img = int(80/100.0*height) - int(mm_value/scaling_factor_y )
							window_position = int(sp[1])
							window_position_img = int(window_position/scaling_factor_x) + 70
							try:
								draw.line(((window_position_img, mm_value_img) + (window_position_img_2, mm_value_img_2)), fill=(r, g, bl), width=1)	
								window_position_img_2 = window_position_img
								mm_value_img_2 = mm_value_img
							except:
								window_position_img_2 = window_position_img
								mm_value_img_2 = mm_value_img

				window_position_img = None 
				mm_value_img = None
				window_position_img_2 = None 
				mm_value_img_2 = None


			#Axes
			draw.line(((wide - 49), int(15/100.0*height)) + ((wide - 49), int(80/100.0*height)), fill=(255, 255, 255, 0), width=2)  #cleanup
			draw.line((68, int(15/100.0*height)) + (68, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#Y axis
			draw.line((68, int(80/100.0*height)) + ((wide - 50), int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#X axis
			draw.line(((wide - 50), int(15/100.0*height)) + ((wide - 50), int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#-Y axis
			draw.line((68, int(15/100.0*height)) + ((wide - 50), int(15/100.0*height)), fill=(0, 0, 0, 0), width=1)	#-X axis

			
			#Axis rulers_____________________
			#X Axis

			if max_contig_len > 1000000  and max_contig_len <= 40000000 :
				mbs = int(0/scaling_factor_x) + 68
				x_tag = 0
				while mbs in range(68, wide-50):
					draw.line((mbs, int(81/100.0*height) ) + (mbs, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	
					if len(str(x_tag)) == 1:
						draw.text(((mbs - 4), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					elif len(str(x_tag)) == 2: 
						draw.text(((mbs - 8), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					
					mbs = mbs + 1000000/scaling_factor_x +1
					x_tag = x_tag + 1

			# fix 8/19
			if max_contig_len > 40000000 and max_contig_len <= 90000000:
				mbs = int(0/scaling_factor_x) + 68
				x_tag = 0
				while mbs in range(68, wide-50):
					draw.line((mbs, int(81/100.0*height) ) + (mbs, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	
					if len(str(x_tag)) == 1:
						draw.text(((mbs - 4), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					elif len(str(x_tag)) == 2: 
						draw.text(((mbs - 8), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					
					mbs = mbs + 5000000/scaling_factor_x +1
					x_tag = x_tag + 5

			# fix 8/19
			if max_contig_len > 90000000 :
				mbs = int(0/scaling_factor_x) + 68
				x_tag = 0
				while mbs in range(68, wide-50):
					draw.line((mbs, int(81/100.0*height) ) + (mbs, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	
					if len(str(x_tag)) == 1:
						draw.text(((mbs - 4), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					if len(str(x_tag)) == 2:
						draw.text(((mbs - 8), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					elif len(str(x_tag)) == 3: 
						draw.text(((mbs - 12), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					
					mbs = mbs + 10000000/scaling_factor_x +1
					x_tag = x_tag + 10


			elif max_contig_len <= 1000000:
				mbs = int(0/scaling_factor_x) + 68
				x_tag = 0
				while mbs in range(68, wide-50):
					draw.line((mbs, int(81/100.0*height) ) + (mbs, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	
					draw.text(((mbs - 4*len(str(x_tag))), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
					mbs = mbs + 100000/scaling_factor_x +1
					x_tag = x_tag + 100000


			#Y axis
			fa_img_0 = int(80/100.0*height) - int(0/scaling_factor_y) - 1		
			fa_img_1 = int(80/100.0*height) - int(1/scaling_factor_y) - 1
			fa_img_05 = int(80/100.0*height) - int(0.5/scaling_factor_y) - 1
			fa_img_025 = int(80/100.0*height) - int(0.25/scaling_factor_y) - 1
			fa_img_075 = int(80/100.0*height) - int(0.75/scaling_factor_y) - 1

			draw.line(( 68 , fa_img_0 +1) + ( 63 , fa_img_0 +1 ), fill=(0, 0, 0, 0), width=1)	
			draw.line(( 68 , fa_img_1 ) + ( 63 , fa_img_1 ), fill=(0, 0, 0, 0), width=1)	
			draw.line(( 68 , fa_img_05 ) + ( 63 , fa_img_05 ), fill=(0, 0, 0, 0), width=1)	
			draw.line(( 68 , fa_img_025 ) + ( 65 , fa_img_025 ), fill=(0, 0, 0, 0), width=1)	
			draw.line(( 68 , fa_img_075 ) + ( 65 , fa_img_075 ), fill=(0, 0, 0, 0), width=1)	

			draw.text(((48), fa_img_0-6), ( '0' ), font=fnt2, fill=(0,0,0,255))
			draw.text(((32), fa_img_1-8), ( '1.0' ), font=fnt2, fill=(0,0,0,255))
			draw.text(((32), fa_img_05-8), ( '0.5' ), font=fnt2, fill=(0,0,0,255))


			#Y axis label
			txt=Image.new('L', (140, 20))
			d = ImageDraw.Draw(txt)
			d.text( (0, 0), "Allele frequency",  font=fnt2, fill=255)
			w=txt.rotate(90,  expand=1)
			im.paste( ImageOps.colorize(w, (0,0,0), (0,0,0)), (2,150),  w)

			#X axis label
			if int(i[1]) > 1000000: x_title = str(i[0]) + ' (Mb)'
			if int(i[1]) <= 1000000: x_title = str(i[0]) + ' (bp)'
			w, h = draw.textsize(str(x_title))
			draw.text((( (wide-120)/2- w/2 +70), (int(87/100.0*height))), (x_title), font=fnt2, fill=(0,0,0,255))


			#Crop and save image, specifying the format with the extension
			w, h = im.size
			if args.my_mut == 'snp':
				im.crop((0, 60, w-0, h-40)).save(project + '/3_workflow_output/mapping_' + str(i[0]) + '.png')

			if args.my_mut == 'af_control':
				im.crop((0, 60, w-0, h-40)).save(project + '/3_workflow_output/control_' + str(i[0]) + '.png')

			if args.my_mut == 'af_sample':
				im.crop((0, 60, w-0, h-40)).save(project + '/3_workflow_output/problem_' + str(i[0]) + '.png')

			if args.my_mut == 'af_candidates':
				im.crop((0, 60, w-0, h-40)).save(project + '/3_workflow_output/candidates_' + str(i[0]) + '.png')

#############################################################################################################
#																											#
# 				SNP - Alelic frequence VS Chromosome position - Zoom in candidate region					#
#																											#
#############################################################################################################
def candidates_zoom():
	#Input 1
	input1 = open(project + '/3_workflow_output/candidate_variants.txt', 'r')
	lines = input1.readlines()	

	#Input 2
	input2 = open(project + '/1_intermediate_files/map_info.txt', 'r')
	lines_map = input2.readlines()	

	#Input fasta
	contig_source = args.input_f_snp


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


	# Read contig fasta file
	contig_lengths = list()

	with open(contig_source) as fp:
		fastalist = list()
		for name_contig, seq_contig in read_fasta(fp):
			innerlist = list()
			innerlist.append(name_contig.strip('>'))
			innerlist.append(len(seq_contig))
			fastalist.append(innerlist)
			contig_lengths.append(len(seq_contig))

	for line in lines_map:
		if line.startswith('?'):
			sp = line.split()
			chromosome = str(sp[1]).lower().strip()
			for i in fastalist:
				if i[0].strip().lower() == chromosome:
					chromosome_length = int(i[1])

			reg_min = sp[2]
			reg_min_real = reg_min
			if int(reg_min) < 0 : reg_min = '0'
			
			reg_max = sp[3]
			reg_max_real = reg_max

			if int(reg_max) > chromosome_length : reg_max = chromosome_length


	#Create image
	wide=1000								
	height=500
	im = Image.new("RGB", (wide, int(height)), (255,255,255))
	draw = ImageDraw.Draw(im)
	
	#get fonts from foler 'fonts'
	fnt2 = ImageFont.truetype('fonts/VeraMono.ttf', 14)

	max_graph_x = int(int(reg_max) - int(reg_min))

	#Scaling factors
	scaling_factor_x = (max_graph_x)/(wide - 120)									#nts/pixel        
	scaling_factor_y = float(0.201/(63/100.0*height))								#fa/pixels
	scaling_factor_y_1 = float(1.001/(63/100.0*height))								#fa/pixels

	#shading
	draw.rectangle( [70, int(15/100.0*height), wide-50, int(80/100.0*height)], fill=(249, 222, 252) )

	#snps
	r, g, b = 31, 120, 180
	for l, line in enumerate(lines):
		if not line.startswith('@'):
			sp = line.split()
			fa = float(sp[8]) - 0.8
			fa_img = int(80/100.0*height) - int(fa/scaling_factor_y)
			pos_img = int((int(sp[2])-int(reg_min))/scaling_factor_x) + 70
			draw.ellipse((pos_img-2, fa_img-2, pos_img+2, fa_img+2), fill=(r, g, b))

	
	#Peak line
	peak = (((int(reg_max_real) + int(reg_min_real))/2) - int(reg_min))/scaling_factor_x + 70

	if 	(((int(reg_max) + int(reg_min_real))/2) - int(reg_min)) < 0: peak = 70

	h = int(16/100.0*height)
	while h in range(int(15/100.0*height), int(80/100.0*height)):
		draw.line((peak, h) + (peak, h+7), fill=(255, 0, 0, 0), width=1)
		h = h + 14
	

	#Axes
	draw.line(((wide - 49), int(15/100.0*height)) + ((wide - 49), int(80/100.0*height)), fill=(255, 255, 255, 0), width=2)  #cleanup
	draw.line((68, int(15/100.0*height)) + (68, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#Y axis
	draw.line((68, int(80/100.0*height)) + ((wide - 50), int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#X axis
	draw.line(((wide - 50), int(15/100.0*height)) + ((wide - 50), int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	#-Y axis
	draw.line((68, int(15/100.0*height)) + ((wide - 50), int(15/100.0*height)), fill=(0, 0, 0, 0), width=1)	#-X axis
	
	#Axis rulers_____________________
	#X Axis
	interval_width = int(args.interval_width)
	if interval_width == 4000000:
		mark = 68
		mark_2 = 68 + 500000/scaling_factor_x/5
		x_tag = int(reg_min)
		while mark in range(68, wide-50):
			draw.line((mark, int(81/100.0*height) ) + (mark, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)	
			w, h = draw.textsize(str(x_tag))
			draw.text(((mark - w/2 -4), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
			mark = mark + 500000/scaling_factor_x
			x_tag = x_tag + 500000

	if interval_width == 20000000:
		mark = 68
		mark_2 = 68 + 4000000/scaling_factor_x/5
		x_tag = int(reg_min)
		while mark in range(68, wide-50):
			draw.line((mark, int(81/100.0*height) ) + (mark, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)
			w, h = draw.textsize(str(x_tag))
			draw.text(((mark - w/2 -4), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
			mark += 4000000/scaling_factor_x
			x_tag += 4000000

	if interval_width == 10000000:
		mark = 68
		mark_2 = 68 + 2000000/scaling_factor_x/5 
		x_tag = int(reg_min)
		while mark in range(68, wide-50):
			draw.line((mark, int(81/100.0*height) ) + (mark, int(80/100.0*height)), fill=(0, 0, 0, 0), width=1)
			w, h = draw.textsize(str(x_tag))
			draw.text(((mark - w/2 -4), (int(81.8/100.0*height))), (str(x_tag).strip()), font=fnt2, fill=(0,0,0,255))
			mark += 2000000/scaling_factor_x
			x_tag += 2000000


	#Y axis
	fa_img_0 = int(80/100.0*height) - int(0/scaling_factor_y) -1
	fa_img_09 = int(80/100.0*height) - int(0.1/scaling_factor_y) -1
	fa_img_1 = int(80/100.0*height) - int(0.2/scaling_factor_y) -1
	fa_img_95 = int(80/100.0*height) - int(0.15/scaling_factor_y) -1
	fa_img_85 = int(80/100.0*height) - int(0.05/scaling_factor_y) -1

	draw.line(( 68 , fa_img_0 +1) + ( 63 , fa_img_0 +1 ), fill=(0, 0, 0, 0), width=1)	
	draw.line(( 68 , fa_img_09 +1) + ( 63 , fa_img_09 +1 ), fill=(0, 0, 0, 0), width=1)	
	draw.line(( 68 , fa_img_1 +1) + ( 63 , fa_img_1 +1 ), fill=(0, 0, 0, 0), width=1)	
	draw.line(( 68 , fa_img_95 +1) + ( 65 , fa_img_95 +1 ), fill=(0, 0, 0, 0), width=1)	
	draw.line(( 68 , fa_img_85 +1) + ( 65 , fa_img_85 +1 ), fill=(0, 0, 0, 0), width=1)	

	draw.text(((32), fa_img_0-8), ( '0.8' ), font=fnt2, fill=(0,0,0,255))
	draw.text(((32), fa_img_1-8), ( '1.0' ), font=fnt2, fill=(0,0,0,255))
	draw.text(((32), fa_img_09-8), ( '0.9' ), font=fnt2, fill=(0,0,0,255))

	#Y axis label
	txt=Image.new('L', (140, 20))
	d = ImageDraw.Draw(txt)
	d.text( (0, 0), "Allele frequency",  font=fnt2, fill=255)
	w=txt.rotate(90,  expand=1)
	im.paste( ImageOps.colorize(w, (0,0,0), (0,0,0)), (2,150),  w)

	#X axis label
	x_title = chromosome + ' (bp)'
	w, h = draw.textsize(str(x_title))
	draw.text((( (wide-120)/2- w/2 +70), (int(87/100.0*height))), (x_title), font=fnt2, fill=(0,0,0,255))

	#Crop and save image, specifying the format with the extension
	w, h = im.size
	if args.my_mut == 'snp':
		im.crop((0, 60, w-0, h-40)).save(project + '/3_workflow_output/candidates_' + chromosome + '_zoom.png')


#############################################################################################################
#																											#
# 							SNP - Alelic frequence VS Chromosome position - LEGEND							#
#																											#
#############################################################################################################


def legend():

	wide = 250
	high = 310

	im = Image.new("RGB", (wide, high), (255,255,255))
	draw = ImageDraw.Draw(im)
	fnt2 = ImageFont.truetype('fonts/VeraMono.ttf', 14)

	w = 10
	h = 10
	length = high - 10
	width = wide - 10

	#legend box
	draw.line((w, h) + (width, h), fill=256, width=1)
	draw.line((w, h) + (w, length), fill=256, width=1)
	draw.line((width, length) + (width, h), fill=256, width=1)
	draw.line((w, length) + (width, length), fill=256, width=1)

	#legend items
	draw.text((w+20, h+20), 'Legend:', font=fnt2, fill=(0,0,0,255))
	
	draw.ellipse((w+40-2, h+60-2, w+40+2, h+60+2), fill=(31, 120, 180))
	draw.text((w+60, h+52), 'F2 problem SNPs', font=fnt2, fill=(0,0,0,255))

	draw.ellipse((w+40-2, h+90-2, w+40+2, h+90+2), fill=(245, 120, 44))
	draw.text((w+60, h+82), 'Control SNPs', font=fnt2, fill=(0,0,0,255))

	draw.line((w+38, h+120) + (w+44, h+120), fill=(46, 255, 0), width=2)
	draw.text((w+60, h+112), 'SMA', font=fnt2, fill=(0,0,0,255))

	draw.line((w+38, h+150) + (w+44, h+150), fill=(255, 0, 0), width=2)
	draw.text((w+60, h+142), 'Boost', font=fnt2, fill=(0,0,0,255))

	draw.line((w+38, h+180) + (w+44, h+180), fill=(255, 0, 255), width=2)
	draw.text((w+60, h+172), 'AF difference', font=fnt2, fill=(0,0,0,255))

	draw.line((w+37, h+210) + (w+45, h+210), fill=(255, 0, 0), width=2)
	draw.line((w+40, h+210) + (w+42, h+210), fill=(255, 255, 255), width=2)
	draw.text((w+60, h+202), 'Selected position', font=fnt2, fill=(0,0,0,255))

	draw.line((w+38, h+240) + (w+44, h+240), fill=(249, 222, 252), width=8)
	draw.text((w+60, h+232), 'Candidate region', font=fnt2, fill=(0,0,0,255))
	draw.rectangle( [w+38, h+244, w+44, h+236], fill=None, outline=(0,0,0) )

	draw.line((w+38, h+270) + (w+44, h+270), fill=(255, 252, 232), width=8)
	draw.text((w+60, h+262), 'Selected chromosome', font=fnt2, fill=(0,0,0,255))
	draw.rectangle( [w+38, h+274, w+44, h+266], fill=None, outline=(0,0,0) )

	im.save(project + '/3_workflow_output/legend.png')


#############################################################################################################
#																											#
# 											Gene plot - LEGEND												#
#																											#
#############################################################################################################

def gene_legend():

	wide = 250
	high = 230

	im = Image.new("RGB", (wide, high), (255,255,255))
	draw = ImageDraw.Draw(im)
	fnt2 = ImageFont.truetype('fonts/VeraMono.ttf', 14)

	w = 10
	h = 10
	length = high - 10
	width = wide - 10

	#legend box
	draw.line((w, h) + (width, h), fill=256, width=1)
	draw.line((w, h) + (w, length), fill=256, width=1)
	draw.line((width, length) + (width, h), fill=256, width=1)
	draw.line((w, length) + (width, length), fill=256, width=1)

	#legend items
	draw.text((w+20, h+20), 'Legend:', font=fnt2, fill=(0,0,0,255))
	
	draw.line((w+30, h+60) + (w+44, h+60), fill=(59, 119, 214), width=9)
	draw.text((w+60, h+52), 'Exon', font=fnt2, fill=(0,0,0,255))

	draw.line((w+30, h+90) + (w+44, h+90), fill=(188, 209, 242), width=9)
	draw.text((w+60, h+82), 'UTR', font=fnt2, fill=(0,0,0,255))

	draw.line((w+30, h+120) + (w+44, h+120), fill=(14, 54, 119), width=2)
	draw.text((w+60, h+112), 'Intron', font=fnt2, fill=(0,0,0,255))

	draw.line((w+30, h+150) + (w+34, h+150), fill=(14, 54, 119), width=2)
	draw.line((w+40, h+150) + (w+44, h+150), fill=(14, 54, 119), width=2)
	draw.text((w+60, h+142), 'Putative promoter', font=fnt2, fill=(0,0,0,255))

	# SNP arrow
	#draw.line((w+37, h+180+4+1) + (w+37, h+180-4+1), fill=(255, 0, 0), width=3)
	#draw.polygon( (w+37-4, h+180-5+2) + (w+37+4, h+180-5+2) +  (w+37, h+180-5-5+2), fill=(255,0,0))
	#draw.text((w+60, h+172), 'Variant position', font=fnt2, fill=(0,0,0,255))

	# Ins triangle
	draw.polygon((w+37-7, h+180+4) + (w+37+7, h+180+4) +  (w+37, h+180-5-5+2), fill=(255,0,0))
	draw.text((w+60, h+172), 'Insertion site', font=fnt2, fill=(0,0,0,255))

	im.save(project + '/3_workflow_output/gene_legend_ins.png')


#############################################################################################################
#																											#
# 									LIN - GENOME OVERVIEW & HISTOGRAMS										#
#																											#
#############################################################################################################

def insertions_overview_and_histograms():

	#Input 1
	input = args.input
	f1 = open(input, 'r')
	lines = f1.readlines()	

	#__________________________________________Insertions overview image_____________________________________________________________
	#________________________________________________________________________________________________________________________________

	#Input 2
	finput = args.input_f
	f2 = open(finput, 'r')
	flines = f2.readlines()	
	#define a superlist with innerlists, each of them containing all the info of each contig 
	superlist = list()


	#Input 3
	input_pos = args.ins_pos
	f3 = open(input_pos, 'r')
	lines_pos = f3.readlines()	


	#Create a list with all the genome contigs
	contigs = []
	length = 0
	n = 0
	#dict_contigs = dict()
	lengthlist = list()


	contig_source = args.input_f

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

	long_contigs=list()
	# Read contig fasta file
	with open(contig_source) as fp:
		fastalist = list()
		for name_contig, seq_contig in read_fasta(fp):
			innerlist = list()
			innerlist.append(name_contig.strip('>'))
			innerlist.append(len(seq_contig))
			if int(len(seq_contig)) > 700000:
				fastalist.append(innerlist)
				long_contigs.append(name_contig.strip('>').lower())
	try:
		max_list = list()
		for c in fastalist:
			max_list.append(int(c[1]))
		max_length = max(max_list)

	except:
		max_length = fastalist[0][1]

	mb_max =  int(math.ceil(float(max_length)/1000000))


	#Calculate the length of the image acording to the number of contigs
	num_contigs = int(len(fastalist))
	
	contigs_image_length = 65 * num_contigs + 60
	im = Image.new("RGB", (1000, contigs_image_length+120), (255,255,255))

	#read fasta file to determine number of contigs
	for i, line in enumerate(flines):
		if line.startswith('>'): 		#fasta sequences start with '>'
			sp = line.split(' ')  		#because some names have whitespaces and extra info that is not written to sam file
			cont = sp[0].strip()  		#strip() is to remove the '\r\n' hidden chars
			cont = cont[1:]       		#to remove the first char of string (>)
			if cont not in contigs:
				contigs.append(cont)
				innerlist = list()
				innerlist.append(cont)
				if cont.lower() in long_contigs:
					superlist.append(innerlist)

	for c in superlist:
		for l in fastalist:
			if c[0] == l[0]:
				c.append(l[1])

	contigs_scaling_factor = mb_max*1000000 / 850.0

	#translate real chromosome lengths into image coordinates
	contig_counter = 1
	contig_xi_coord = 100

	for c in superlist:
		contig_xf_coord = int(contig_xi_coord + c[1]/contigs_scaling_factor)
		contig_y_coord = ((contigs_image_length / num_contigs) * contig_counter)
		contig_counter +=1
		c.append(contig_y_coord)
		c.append(contig_xi_coord)
		c.append(contig_xf_coord)

	#add insertions aproximate position to superlist
	for c in superlist: 
		positions_list = list()
		for i, line in enumerate(lines_pos):
			if not line.startswith('#'): 
				sp = line.split('\t')
				contig = str(sp[1].strip()).lower()
				insertion_pos = int(sp[2])
				if contig == c[0].strip().lower() and insertion_pos not in positions_list:
					positions_list.append(insertion_pos)
		c.append(positions_list)

	#initialize draw
	draw = ImageDraw.Draw(im)

	#get fonts from foler 'fonts'
	fnt3 = ImageFont.truetype('fonts/VeraMono.ttf', 14)
	tab = 50
	number = 1

	#sort superlist (8-7-19 fix)
	superlist.sort(key = lambda x: x[0]) 

	#Drawing the chromosomes:
	for c in superlist:
		previous_pos_x = 0
		previous_chr = 'none'
		draw.line((c[3], c[2]) + (c[4], c[2]), fill=(31, 120, 180), width=13)
		draw.text(((c[3]), (c[2] -55)), c[0], font=fnt3, fill=(0,0,0,255))

		for i in c[5]:
			d = abs(i/contigs_scaling_factor+contig_xi_coord - previous_pos_x)
			if d > 21 or previous_chr != c[0]: 
				draw.polygon([((i/contigs_scaling_factor+contig_xi_coord), (c[2]-8)), (((i/contigs_scaling_factor+contig_xi_coord)+10), (c[2]-18)), ((i/contigs_scaling_factor + contig_xi_coord)-10, c[2]-18)], fill = (200, 0, 0, 200))
				draw.line((((i/contigs_scaling_factor+contig_xi_coord), (c[2]-8)), ((i/contigs_scaling_factor+contig_xi_coord)+10), (c[2]-18)), fill=256, width=1)
				draw.line((((i/contigs_scaling_factor+contig_xi_coord)+10), (c[2]-18)) + ((i/contigs_scaling_factor + contig_xi_coord)-10, (c[2]-18)), fill=256, width=1) 
				draw.line(((i/contigs_scaling_factor+contig_xi_coord), (c[2]-8)) + ((i/contigs_scaling_factor + contig_xi_coord)-10, (c[2]-18)), fill=256, width=1) 
				draw.text(((i/contigs_scaling_factor+contig_xi_coord -4), (c[2]-35)), ( str(number)), font=fnt3, fill=(0,0,0,255))
			else:
				draw.polygon([((i/contigs_scaling_factor+contig_xi_coord), (c[2]+8)), (((i/contigs_scaling_factor+contig_xi_coord)-10), (c[2]+18)), ((i/contigs_scaling_factor + contig_xi_coord)+10, c[2]+18)], fill = (200, 0, 0, 200))
				draw.line((((i/contigs_scaling_factor+contig_xi_coord), (c[2]+8)), ((i/contigs_scaling_factor+contig_xi_coord)-10), (c[2]+18)), fill=256, width=1)
				draw.line((((i/contigs_scaling_factor+contig_xi_coord)-10), (c[2]+18)) + ((i/contigs_scaling_factor + contig_xi_coord)+10, (c[2]+18)), fill=256, width=1) 
				draw.line(((i/contigs_scaling_factor+contig_xi_coord), (c[2]+8)) + ((i/contigs_scaling_factor + contig_xi_coord)+10, (c[2]+18)), fill=256, width=1) 
				draw.text(((i/contigs_scaling_factor+contig_xi_coord -4), (c[2]+20)), ( str(number)), font=fnt3, fill=(0,0,0,255))

			number = number + 1
			previous_pos_x = i/contigs_scaling_factor+contig_xi_coord
			previous_chr = c[0]

		#Chromosome caps
		image_file_left = StringIO(open("./fonts/left_cap.png",'rb').read())
		image_file_right = StringIO(open("./fonts/right_cap.png",'rb').read())

		cap_left = Image.open(image_file_left)
		cap_right = Image.open(image_file_right)

		im.paste(cap_left, (c[3], c[2]-6))
		im.paste(cap_right, (c[4], c[2]-6))

	#Axis
	draw.line((100, contigs_image_length + 50) + (950, contigs_image_length + 50), fill=256, width=1)

	mb = 850.0/mb_max
	m = 100.0
	tag = 'Position (Mb)'
	num = 0
	draw.text((480, contigs_image_length + 88), str(tag), font=fnt3, fill=256)

	if mb_max < 40:
		while int(m) in range(100, 951):
			draw.line(( int(m), contigs_image_length + 48) + (int(m), contigs_image_length + 52), fill=256, width=1)
			w, h = draw.textsize(str(num))
			draw.text((int(m) - w/2 -2, contigs_image_length + 57), str(num), font=fnt3, fill=256)
			num = num + 1
			m = m + mb

	if mb_max >= 40:
		while int(m) in range(100, 951):
			draw.line(( int(m), contigs_image_length + 48) + (int(m), contigs_image_length + 52), fill=256, width=1)
			w, h = draw.textsize(str(num))
			draw.text((int(m) - w/2 -2, contigs_image_length + 57), str(num), font=fnt3, fill=256)
			num = num + 5
			m = m + mb*5

	im.save(project + "/3_workflow_output/insertions_overview.png")

	#_________________________________________________________Local and paired analysis graphs________________________________________________________________	
	#_________________________________________________________________________________________________________________________________________________________
	if args.mode == 'pe':
		insertions = list()
		for i, line in enumerate(lines):
			if not line.startswith('@'):	
				sp = line.split('\t')
				insertion = str(sp[2]).strip()
				if insertion not in insertions and insertion != '-':
					insertions.append(insertion)

		for e in insertions:
			try:
				del region_min
			except:
				pass
			region_max = 0
			rd_max_paired = 0
			rd_max_local = 0
			for i, line in enumerate(lines):
				if not line.startswith('@'):
					sp = line.split('\t')
					#Max and min for genome region in graphic
					if sp[2] == e:		# and sp[0] == 'PAIRED' : 
						if int(sp[3]) > region_max:
							region_max = int(sp[3])
						else:
							try:
								if sp[3] < region_min: 
									region_min = int(sp[3])
							except:
								region_min = int(sp[3])
			
					#Max and min read depth 
					if sp[2] == e and sp[0] == 'PAIRED' and sp[5].strip() != 'TOTAL': 
						if int(sp[4]) > rd_max_paired:
							rd_max_paired = int(sp[4])		

					if sp[2] == e and sp[0] == 'LOCAL_RD' and sp[5].strip() != 'TOTAL_RD': 
						if int(sp[4]) > rd_max_local:
							rd_max_local = int(sp[4])		

			rd_max_paired = rd_max_paired + 1
			rd_max_local = rd_max_local + 1
			region_max = region_max + 100
			if region_min > 200:
				region_min = region_min - 100

			#Images, axes and title 
			im = Image.new("RGB", (1000, 1000), (255,255,255))
			draw = ImageDraw.Draw(im)
			draw.line((120, 449) + (900, 449), fill=256, width=1) 								#x axis paired
			draw.line((120, 150) + (120, 449), fill=256, width=1)   							#y axis paired
			draw.line((120, 754) + (900, 754), fill=256, width=1) 								#x axis local
			draw.line((120, 455) + (120, 754), fill=256, width=1)   							#y axis local
			
			draw.line((120, 150) + (900, 150), fill=256, width=1) 								#-x axis paired
			draw.line((900, 150) + (900, 449), fill=256, width=1)   							#-y axis paired
			draw.line((120, 455) + (900, 455), fill=256, width=1) 								#-x axis local
			draw.line((900, 455) + (900, 754), fill=256, width=1)   							#-y axis local

			draw.text(((450), (795)), ('Nucleotide'), font=fnt3, fill=(0,0,0,255))
			draw.text(((140), (155)), ('Flanking unpaired alignments'), font=fnt3, fill=(0,0,0,255))
			draw.text(((140), (460)), ('Flanking local alignments'), font=fnt3, fill=(0,0,0,255))

			#Y axis label
			txt=Image.new('L', (500, 50))
			d = ImageDraw.Draw(txt)
			d.text( (0, 0), "Read depth (x)",  font=fnt3, fill=255)
			w=txt.rotate(90,  expand=1)
			im.paste( ImageOps.colorize(w, (0,0,0), (0,0,0)), (35,0),  w)

			#Scaling factors 
			nucleotides = region_max - region_min
			scaling_factor_x = nucleotides/780.0
			scaling_factor_y_paired = rd_max_paired/280.0
			scaling_factor_y_local = rd_max_local/280.0

			#LOCAL/PAIRED GRAPHICS
			for i, line in enumerate(lines):
				if not line.startswith('@'):
					sp = line.split('\t')
					ins_contig = sp[1]
					if sp[2] == e and sp[0].strip() == 'PAIRED' and sp[5].strip() != 'TOTAL':
						raw_x_position = int(sp[3])
						img_x_position = int(raw_x_position/scaling_factor_x)
						img_relative_x_position = img_x_position - int(region_min/scaling_factor_x) + 121

						raw_y_position = int(sp[4])
						img_y_position_p = 450 - int(raw_y_position/scaling_factor_y_paired)
				
						#draw
						if sp[5].strip() == 'R':
							draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_y_position_p), fill=(64, 159, 65, 100), width=2)	
						
						elif sp[5].strip() == 'F':
							draw.line((img_relative_x_position, 448) + (img_relative_x_position, img_y_position_p), fill=(31, 120, 180, 100), width=2)

			img_relative_y_position_2_r = float('inf')
			img_relative_y_position_2_l = float('inf')
			cand_pos_l = 'none'
			cand_pos_r = 'none'
			for i, line in enumerate(lines):
				if not line.startswith('@'):
					sp = line.split('\t')
					ins_contig = sp[1]
					if sp[2] == e and sp[0].strip() == 'LOCAL_RD' and sp[5].strip() != 'TOTAL_RD':
						raw_x_position = int(sp[3])
						img_x_position = int(raw_x_position/scaling_factor_x)
						img_relative_x_position = img_x_position - int(region_min/scaling_factor_x) + 121

						raw_y_position = int(sp[4])
						img_y_position_l = int(raw_y_position/scaling_factor_y_local)
						img_relative_y_position = 755 - img_y_position_l 

						#draw
						if sp[5].strip() == 'RIGHT_RD':
							draw.line((img_relative_x_position, 753) + (img_relative_x_position, img_relative_y_position), fill=(64, 159, 65, 100), width=2)
							if img_relative_y_position < img_relative_y_position_2_r:
								cand_pos_r = img_relative_x_position
								img_relative_y_position_2_r = img_relative_y_position
						if sp[5].strip() == 'LEFT_RD':
							draw.line((img_relative_x_position, 753) + (img_relative_x_position, img_relative_y_position), fill=(31, 120, 180, 100), width=2)
							if img_relative_y_position <= img_relative_y_position_2_l:
								cand_pos_l = img_relative_x_position
								img_relative_y_position_2_l = img_relative_y_position
			'''
			#Candidate regions
			for i, line in enumerate(lines):
				if line.startswith('@#'):
					if 'inf' not in line:
						sp = line.split(',')
						if int(sp[2].strip()) == int(e):
								cr = [int(sp[0].strip('@#')), int(sp[1].strip())]
								cr_min = min(cr)
								cr_max = max(cr)
								draw.line((((120 +int(sp[0].strip('@#'))/scaling_factor_x - int(region_min/scaling_factor_x)) , 448) + ((120 +int(sp[0].strip('@#'))/scaling_factor_x - int(region_min/scaling_factor_x)) , 151)), fill=(147, 147, 147, 0), width=1)
								draw.line((((120 +int(sp[1].strip())/scaling_factor_x - int(region_min/scaling_factor_x)) , 448) + ((120 +int(sp[1].strip())/scaling_factor_x - int(region_min/scaling_factor_x)) , 151)), fill=(147, 147, 147, 0), width=1)
			'''

			#Candidate position:
			if cand_pos_r != 'none' and cand_pos_l == 'none':
				draw.line((cand_pos_r-1 , 456) + (cand_pos_r-1 , 753), fill=(147, 147, 147, 0), width=1)		
			if cand_pos_l != 'none' and cand_pos_r == 'none':
				draw.line((cand_pos_l-1 , 456) + (cand_pos_l-1 , 753), fill=(147, 147, 147, 0), width=1)
			if cand_pos_r != 'none' and cand_pos_l != 'none':
				draw.line((cand_pos_r-1 , 456) + (cand_pos_r-1 , 753), fill=(147, 147, 147, 0), width=1)

			#Axis anotations
			#x Axis
			x_p = 120 + int(100/scaling_factor_x)
			x_p_2 = 120 + int(200/scaling_factor_x)
			ruler = region_min + 100

			while x_p in range(120, 900):
				draw.line((x_p, 755) + (x_p, 762), fill=256, width=1)
				w, h = draw.textsize(str(ruler))
				draw.text((x_p - w/2 - 5, 766), (str(ruler)), font=fnt3, fill=(0,0,0,255))  
				ruler = ruler + 200
				x_p = int(x_p + (200/scaling_factor_x))			#Ruler with 200 nts separations
			
			while x_p_2 in range(120, 900):
				draw.line((x_p_2, 755) + (x_p_2, 758), fill=256, width=1)
				x_p_2 = int(x_p_2 + (200/scaling_factor_x))		#Ruler with 100 nts separations

			#y Axis - paired
			if rd_max_paired > 20:
				y_p = 450 - int(5/scaling_factor_y_paired)
				counter = 10
				while y_p in range(150, 451): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 10
					y_p = int(y_p - (10/scaling_factor_y_paired))

			if 20 >= rd_max_paired > 10:
				y_p = 450 - int(5/scaling_factor_y_paired)
				counter = 5
				while y_p in range(150, 451): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 5
					y_p = int(y_p - (5/scaling_factor_y_paired))

			if rd_max_paired <= 10:
				y_p = 450 - int(1/scaling_factor_y_paired)
				counter = 1
				while y_p in range(150, 451): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 1
					y_p = int(y_p - (1/scaling_factor_y_paired))
			
			#y Axis - local
			if rd_max_local > 20:
				y_p = 755 - int(5/scaling_factor_y_local)
				counter = 5
				while y_p in range(455, 751): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 5
					y_p = int(y_p - (5/scaling_factor_y_local))

			if 20 >= rd_max_local > 10:
				y_p = 755 - int(5/scaling_factor_y_local)
				counter = 5
				while y_p in range(455, 751): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 5
					y_p = int(y_p - (5/scaling_factor_y_local))

			if rd_max_local <= 10:
				y_p = 755 - int(1/scaling_factor_y_local)
				counter = 1
				while y_p in range(455, 751): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 1
					y_p = int(y_p - (1/scaling_factor_y_local))


			#Legend paired________________________________________________________________________________________
			#w and h can be used to re-position the legend in the figure
			w = 690
			h = 160

			#legend box
			draw.polygon([(w,h), (w,h+100), (w+200,h+100), (w+200,h) ], fill = (255, 255, 255, 0))
			draw.line((w, h) + (w+200, h), fill=256, width=1)
			draw.line((w, h) + (w, h+100), fill=256, width=1)
			draw.line((w+200, h) + (w+200, h+100), fill=256, width=1)
			draw.line((w, h+100) + (w+200, h+100), fill=256, width=1)
			draw.text((w+10, h+10), 'Legend:', font=fnt3, fill=(0,0,0,255))
			#legend items
			draw.text((w+45, h+32), 'Forward reads', font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+32+7) + (w+35, h+32+7), fill=(31, 120, 180), width=10)
			draw.text((w+45, h+52), 'Reverse reads ', font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+52+7) + (w+35, h+52+7), fill=(64, 159, 65), width=10)
			draw.text((w+45, h+72), 'Candidate region', font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+72+8) + (w+35, h+72+8), fill=(147, 147, 147), width=1)

			#Legend local________________________________________________________________________________________
			#w and h can be used to re-position the legend in the figure
			w = 690
			h = 465

			#legend box
			draw.polygon([(w,h), (w,h+100), (w+200,h+100), (w+200,h) ], fill = (255, 255, 255, 0))
			draw.line((w, h) + (w+200, h), fill=256, width=1)
			draw.line((w, h) + (w, h+100), fill=256, width=1)
			draw.line((w+200, h) + (w+200, h+100), fill=256, width=1)
			draw.line((w, h+100) + (w+200, h+100), fill=256, width=1)
			draw.text((w+10, h+10), 'Legend:', font=fnt3, fill=(0,0,0,255))
			#legend items
			draw.text((w+45, h+32), """3' clipped reads """, font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+32+7) + (w+35, h+32+7), fill=(31, 120, 180), width=10)
			draw.text((w+45, h+52), """5' clipped reads """, font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+52+7) + (w+35, h+52+7), fill=(64, 159, 65), width=10)
			draw.text((w+45, h+72), 'Predicted position', font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+72+8) + (w+35, h+72+8), fill=(147, 147, 147), width=1)

			#save image, specifying the format with the extension
			w, h = im.size
			im.crop((0, 100, w, h-100)).save(project + '/3_workflow_output/img_1_ins_' + str(e) + '.png')

	#_________________________________________________________________________________________________________________________________________________________
	if args.mode == 'se':
		insertions = list()
		for i, line in enumerate(lines):
			if not line.startswith('@'):	
				sp = line.split('\t')
				insertion = str(sp[2]).strip()
				if insertion not in insertions and insertion != '-':
					insertions.append(insertion)
						
		for e in insertions:
			try:
				del region_min
			except:
				pass
			region_max = 0
			rd_max_local = 0
			for i, line in enumerate(lines):
				if not line.startswith('@'):
					sp = line.split('\t')
					#Max and min for genome region in graphic
					if sp[2] == e: 
						if int(sp[3]) > region_max:
							region_max = int(sp[3])
						else:
							try:
								if sp[3] < region_min: 
									region_min = int(sp[3])
							except:
								region_min = int(sp[3])
			
					#Max read depth 
					if sp[2] == e and sp[0] == 'LOCAL_RD' and sp[5].strip() != 'TOTAL_RD': 
						if int(sp[4]) > rd_max_local:
							rd_max_local = int(sp[4]) + 1

			region_max = region_max + 100
			if region_min > 200:
				region_min = region_min - 100

			#Images, axes and title 
			im = Image.new("RGB", (1000, 600), (255,255,255))
			draw = ImageDraw.Draw(im)
			draw.line((120, 450) + (900, 450), fill=256, width=1) 								#x axis 
			draw.line((120, 150) + (120, 450), fill=256, width=1)   							#y axis 
			draw.line((120, 150) + (900, 150), fill=256, width=1) 								#-x axis 
			draw.line((900, 150) + (900, 450), fill=256, width=1)   							#-y axis 

			draw.text(((450), (500)), ('Nucleotide'), font=fnt3, fill=(0,0,0,255))

			#Y axis label
			txt=Image.new('L', (150, 30))
			d = ImageDraw.Draw(txt)
			d.text( (0, 0), "Read depth (x)",  font=fnt3, fill=255)
			w=txt.rotate(90,  expand=1)
			im.paste( ImageOps.colorize(w, (0,0,0), (0,0,0)), (35,200),  w)

		
			#Scaling factors 
			nucleotides = region_max - region_min
			scaling_factor_x = nucleotides/780.0
			scaling_factor_y_local = rd_max_local/280.0

			img_relative_y_position_2_r = float('inf')
			img_relative_y_position_2_l = float('inf')
			cand_pos_l = 'none'
			cand_pos_r = 'none'
			#LOCAL GRAPHICS
			for i, line in enumerate(lines):
				if not line.startswith('@'):
					sp = line.split('\t')
					ins_contig = sp[1]
					if sp[2] == e and sp[0].strip() == 'LOCAL_RD' and sp[5].strip() != 'TOTAL_RD':
						raw_x_position = int(sp[3])
						img_x_position = int(raw_x_position/scaling_factor_x)
						img_relative_x_position = img_x_position - int(region_min/scaling_factor_x) + 121

						raw_y_position = int(sp[4])
						img_y_position_l = int(raw_y_position/scaling_factor_y_local)
						img_relative_y_position = 450 - img_y_position_l 

						#draw
						if sp[5].strip() == 'RIGHT_RD':
							draw.line((img_relative_x_position, 449) + (img_relative_x_position, img_relative_y_position), fill=(64, 159, 65, 100), width=3)
							if img_relative_y_position < img_relative_y_position_2_r:
								cand_pos_r = img_relative_x_position
								img_relative_y_position_2_r = img_relative_y_position

						if sp[5].strip() == 'LEFT_RD':
							draw.line((img_relative_x_position, 449) + (img_relative_x_position, img_relative_y_position), fill=(31, 120, 180, 100), width=3)
							if img_relative_y_position <= img_relative_y_position_2_l:
								cand_pos_l = img_relative_x_position
								img_relative_y_position_2_l = img_relative_y_position

			#Candidate position
			if cand_pos_r != 'none' and cand_pos_l == 'none':
				draw.line((cand_pos_r-1 , 449) + (cand_pos_r-1 , 151), fill=(147, 147, 147, 0), width=1)		
			if cand_pos_l != 'none' and cand_pos_r == 'none':
				draw.line((cand_pos_l-1 , 449) + (cand_pos_l-1 , 151), fill=(147, 147, 147, 0), width=1)
			if cand_pos_r != 'none' and cand_pos_l != 'none':
				draw.line((cand_pos_r-1 , 449) + (cand_pos_r-1 , 151), fill=(147, 147, 147, 0), width=1)

			#Axis anotations
			#x Axis
			x_p = 120 + int(25/scaling_factor_x)
			x_p_2 = 120 + int(50/scaling_factor_x)

			ruler = region_min + 25
		
			while x_p in range(120, 900):
				draw.line((x_p, 450) + (x_p, 457), fill=256, width=1)
				w, h = draw.textsize(str(ruler))
				draw.text((x_p - w/2 - 5, 460), (str(ruler)), font=fnt3, fill=(0,0,0,255))  
				ruler = ruler + 50
				x_p = int(x_p + (50/scaling_factor_x)) #Ruler with 50 nts separations

			while x_p_2 in range(120, 900):
				draw.line((x_p_2, 450) + (x_p_2, 455), fill=256, width=1)
				x_p_2 = int(x_p_2 + (50/scaling_factor_x)) #Ruler with 25 nts separations

			#y Axis 
			if rd_max_local > 8:
				y_p = 450 - int(5/scaling_factor_y_local)
				counter = 5
				while y_p in range(150, 451): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 5
					y_p = int(y_p - (5/scaling_factor_y_local))

			if rd_max_local < 8:
				y_p = 450 - int(1/scaling_factor_y_local)
				counter = 1
				while y_p in range(150, 451): 
					draw.line((120, y_p) + (115, y_p), fill=256, width=1)
					draw.text((90, y_p-8), ( str(counter)), font=fnt3, fill=(0,0,0,255))
					counter = counter + 1
					y_p = int(y_p - (1/scaling_factor_y_local))

			#Legend local________________________________________________________________________________________				<---------------------------------------------
			#w and h can be used to re-position the legend in the figure
			w = 690
			h = 160

			#legend box
			draw.polygon([(w,h), (w,h+100), (w+200,h+100), (w+200,h) ], fill = (255, 255, 255, 0))
			draw.line((w, h) + (w+200, h), fill=256, width=1)
			draw.line((w, h) + (w, h+100), fill=256, width=1)
			draw.line((w+200, h) + (w+200, h+100), fill=256, width=1)
			draw.line((w, h+100) + (w+200, h+100), fill=256, width=1)
			draw.text((w+10, h+10), 'Legend:', font=fnt3, fill=(0,0,0,255))
			#legend items
			draw.text((w+45, h+32), """3' clipped reads """, font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+32+7) + (w+35, h+32+7), fill=(31, 120, 180), width=10)
			draw.text((w+45, h+52), """5' clipped reads """, font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+52+7) + (w+35, h+52+7), fill=(64, 159, 65), width=10)
			draw.text((w+45, h+72), 'Predicted position', font=fnt3, fill=(0,0,0,255))
			draw.line((w+10, h+72+8) + (w+35, h+72+8), fill=(147, 147, 147), width=1)

			#save image, specifying the format with the extension
			w, h = im.size
			im.crop((0, 100, w, h-50)).save(project + '/3_workflow_output/img_1_ins_' + str(e) + '.png')

#############################################################################################################
#																											#
# 											GENE PLOT														#
#																											#
#############################################################################################################

def gene_plot(): 

	if args.my_mut == 'lin':
		#Input 1
		input = args.input
		f1 = open(input, 'r')
		lines = f1.readlines()	

	#Input varanalyzer
	input = args.input_va
	f3 = open(input, 'r')
	lines_va = f3.readlines()	

	#Input gff
	input = args.gff
	f4 = open(input, 'r')
	lines_gff = f4.readlines()	

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



	#We create an 'intermediate list', which will contain the information necesary for the gene plot gathered from the output file of the varanalyzer module, the 
	#sorted_insertions.txt file and the genome feature file. Intermediate list format:
	#	['Chr', 'insertion/snp position', 'insertion number / -', 'gene name', [ref-aa, alt-aa, pos-aa, ref-base, alt-base, strand], [list of gene features: 'type', 'start', 'end'], [list of positions(required for calculations)]]
	intermediate_list = list()
	for i, line in enumerate(lines_va):
		if not line.startswith('@'):
			sp = line.split('\t')
			if sp[5].strip() != 'nh':
				temp_list = list()
				temp_list.append(sp[1])
				temp_list.append(sp[2])
				temp_list.append('-')
				temp_list.append(sp[9])
				refalt = [sp[12].strip(), sp[13].strip(), sp[11].strip(), sp[3].strip(), sp[4].strip(), sp[8].strip()]
				temp_list.append(refalt)
				intermediate_list.append(temp_list)

	if args.my_mut == 'lin':
		for p in intermediate_list: 
			for i, line in enumerate(lines): 
				sp = line.split()
				if p[0].lower().strip() == sp[1].lower().strip() and sp[2] not in p[2]:	#and 'TOTAL' in sp[5]
					if int(p[1]) == int(sp[3]) or int(p[1]) == (int(sp[3]) - 1):
						p[2] = sp[2]

	for p in intermediate_list:
		features = list()
		positions = list()
		refalt = list()
		for g, line in enumerate(lines_gff):
			if not line.startswith('#'):
				sp = line.split('\t')
				if (p[3] + ';')  in (sp[8].strip() + ';'):
					feature = [sp[2], sp[3], sp[4]]
					positions.append(int(sp[3]))
					positions.append(int(sp[4]))
					features.append(feature)
		p.append(features)
		p.append(positions)

	for p in intermediate_list:
		p[5].append(['rr', int(args.rrl)])
		if p[4][5].strip() == '+':
			p[6].append(min(p[6]) - int(args.rrl))
		if p[4][5].strip() == '-':
			p[6].append(max(p[6]) +  int(args.rrl))

	#Drawing the genes:
	for p in intermediate_list:
		wide=1000 																							#<-------------------------------------------------------------------------------- SET IMAGE SIZE
		height=(35/100.0)*wide
		im = Image.new("RGB", (wide, int(height)), (255,255,255))
		draw = ImageDraw.Draw(im)
		gene_max_raw = max(p[6])
		gene_min_raw = min(p[6])
		gene_length = gene_max_raw - gene_min_raw
		gene_px_length = float((0.7)*wide)
		gene_scaling_factor = gene_length/gene_px_length  #bp/pixel

		#Fonts
		fnt2 = ImageFont.truetype('fonts/arial.ttf', int(0.016*wide))
		fnt3 = ImageFont.truetype('fonts/arial.ttf', int(0.024*wide))
		fnt4 = ImageFont.truetype('fonts/arial.ttf', int(0.02*wide))

		#Gene name
		draw.text((int(0.05*wide), int(0.03*wide)), (str(p[3])), font=fnt3, fill=(0,0,0,255))

		#Gene baseline
		if p[4][5] == '+':
			draw.line((int(0.15*wide) + int(int(args.rrl)/gene_scaling_factor), int(180/350.0*height)) + (int(0.15*wide) + gene_px_length, int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
		if p[4][5] == '-':
			draw.line((int(0.15*wide), int(180/350.0*height)) + (int(0.15*wide) + gene_px_length - int(int(args.rrl)/gene_scaling_factor), int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
		
		#Gene features
		atg_list = list()
		if p[4][5] == '+':
			for e in p[5]:		
				if e[0].strip() == 'rr':
					inicio = int(0.15*wide)
					fin = int(0.15*wide) + int(e[1])/gene_scaling_factor
					step = int(0.005*wide)
					s = inicio
					while s in range(inicio, int(fin)):
						draw.line((s, int(180/350.0*height)) + (s + step, int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
						s = s + step*2

		if p[4][5] == '-':
			for e in p[5]:		
				if e[0].strip() == 'rr':
					fin = int(0.85*wide)
					inicio = int(0.85*wide) - int(e[1])/gene_scaling_factor
					step = int(0.005*wide)
					s = int(inicio)
					while s in range(int(inicio), int(fin)):
						draw.line((s, int(180/350.0*height)) + (s + step, int(180/350.0*height)), fill=(14, 54, 119), width=int(0.004*wide))
						s = s + step*2

		for e in p[5]:
			if e[0].strip() == 'exon':
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide)
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(59, 119, 214), width=int(0.02*wide))
				#draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)	
				#draw.line((inicio, int(190/350.0*height)) + (fin+1, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#if p[4][5] == '+':
				#	draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#if p[4][5] == '-':
				#	draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
		
		for e in p[5]:
			if 'utr' in (e[0].strip()).lower() and 'five' in (e[0].strip()).lower():																		# Backup UTR drawing
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(188, 209, 242), width=int(0.02*wide))
				#draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
			if 'utr' in (e[0].strip()).lower() and 'three' in (e[0].strip()).lower():																		# Backup UTR drawing
				inicio = int((int(e[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				fin = int((int(e[2]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
				draw.line((inicio, int(180/350.0*height)) + (fin, int(180/350.0*height)), fill=(188, 209, 242), width=int(0.02*wide))
				#draw.line((inicio, int(170/350.0*height)) + (fin, int(170/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#draw.line((inicio, int(190/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#if p[4][5] == '+':
				#	draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#if p[4][5] == '-':
				#	draw.line((fin, int(170/350.0*height)) + (fin, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)
				#	draw.line((inicio, int(170/350.0*height)) + (inicio, int(190/350.0*height)), fill=(0, 4, 71, 0), width=2)

		#Gene direction
		if p[4][5] == '+':
			draw.polygon([(int(0.84*wide), int(168/350.0*height)), (int(0.851*wide) , int(168/350.0*height)), (int(0.851*wide), int(181/350.0*height))], fill = (255, 255, 255, 0))
			draw.polygon([(int(0.84*wide), int(192/350.0*height)), (int(0.851*wide) , int(192/350.0*height)), (int(0.851*wide), int(180/350.0*height))], fill = (255, 255, 255, 0))
			#draw.line((int(0.841*wide), int(170/350.0*height)) + (int(0.85*wide), int(180/350.0*height)), fill=(0, 4, 71, 0), width=2)
			#draw.line((int(0.841*wide), int(190/350.0*height)) + (int(0.851*wide), int(180/350.0*height)), fill=(0, 4, 71, 0), width=2)

				
		if p[4][5] == '-':
			draw.polygon([(int(0.16*wide), int(168/350.0*height)), (int(0.149*wide) , int(168/350.0*height)), (int(0.149*wide), int(180/350.0*height))], fill = (255, 255, 255, 0))
			draw.polygon([(int(0.16*wide), int(192/350.0*height)), (int(0.149*wide) , int(192/350.0*height)), (int(0.149*wide), int(180/350.0*height))], fill = (255, 255, 255, 0))
			#draw.line((int(0.158*wide), int(170.5/350.0*height)) + (int(0.148*wide), int(180.5/350.0*height)), fill=(0, 4, 71, 0), width=2)
			#draw.line((int(0.159*wide), int(190/350.0*height)) + (int(0.149*wide), int(180/350.0*height)), fill=(0, 4, 71, 0), width=2)

			draw.line((int(0.15*wide), int(169/350.0*height)) + (int(0.16*wide), int(169/350.0*height)), fill=(255, 255, 255, 0), width=1)

		#Scale bar
		px_scale_test = float(100/gene_scaling_factor)
		if px_scale_test >=4: 
			scale = 100
			scale_tag = '100 bp'
		if px_scale_test < 4: 
			scale = 500
			scale_tag = '500 bp'

		if args.my_mut == 'snp':
			#scale = 100
			#scale_tag = '100 bp'
			w, h = draw.textsize(str(scale_tag))
			px_scale = float(scale/gene_scaling_factor)
			draw.line((int(0.91*wide) - int(px_scale) - w/2 + px_scale/2, int(110/350.0*height)) + (int(0.91*wide) - w/2 + px_scale/2, int(110/350.0*height)), fill=(0, 0, 0, 0), width=int(0.002*wide))
			draw.text((int(0.87*wide), int(117.8/350.0*height)), (scale_tag), font=fnt2, fill=(0,0,0,255))

		if args.my_mut == 'lin':
			#scale = 100
			#scale_tag = '100 bp'
			w, h = draw.textsize(str(scale_tag))
			px_scale = float(scale/gene_scaling_factor)
			draw.line((int(0.91*wide) - int(px_scale) - w/2 + px_scale/2, int(250/350.0*height)) + (int(0.91*wide) - w/2 + px_scale/2, int(250/350.0*height)), fill=(0, 0, 0, 0), width=int(0.002*wide))
			draw.text((int(0.87*wide), int(257.8/350.0*height)), (scale_tag), font=fnt2, fill=(0,0,0,255))

		#Insertion triangle and info
		if args.my_mut == 'lin':
			ins_pos = int((int(p[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			draw.polygon([(ins_pos, int(170/350.0*height)), (ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide))], fill = (200, 0, 0, 200))
			draw.line((ins_pos, int(170/350.0*height)) + (ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
			draw.line((ins_pos, int(170/350.0*height)) + (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
			draw.line((ins_pos - int(0.02*wide), int(170/350.0*height) - int(0.025*wide)) + (ins_pos + int(0.02*wide), int(170/350.0*height) - int(0.025*wide)), fill=(0, 0, 0, 0), width=1)
			draw.text((ins_pos - int(0.04*wide), int(0.115*wide)), ('Insertion ' + str(p[2])), font=fnt4, fill=(0,0,0,255))


		#SNP arrow and info
		if args.my_mut == 'snp':
			snp_pos = int((int(p[1]) - gene_min_raw)/gene_scaling_factor)  + int(0.15*wide) 
			draw.line((snp_pos, int(194/350.0*height)) + (snp_pos , int(194/350.0*height) + int(0.03*wide)), fill=(180, 0, 0, 0), width=int(0.005*wide))						
			draw.polygon([(snp_pos, int(191/350.0*height)), (snp_pos - int(0.01*wide), int(191/350.0*height) + int(0.01*wide)), (snp_pos + int(0.01*wide), int(191/350.0*height) + int(0.01*wide))], fill = (200, 0, 0, 200))

			#Aa change
			if p[4][0].strip() != '-' : 												
				aach = 'yes'
				draw.text((int(snp_pos - int(0.092*wide)), int(0.75*height)), (
					str(p[4][0])+ ' (' + str(p[4][2]) +')' +  '        '  +
					str(p[4][1])), font=fnt4, fill=(0,0,0,255))   				
				str_len = len(str(p[4][2]))
			else:
				aach = 'no'

			#Base change # SDL V2 update - fix for indels
			msg_1 = str(p[4][3])
			msg_2 = str(p[4][4])
			w1, h1 = draw.textsize(msg_1, font=fnt4)
			w2, h2 = draw.textsize(msg_2, font=fnt4)

			draw.text(((  snp_pos - (w1 + 22)   ),(int(0.67*height))), msg_1, font=fnt4, fill=(0,0,0,255))
                        draw.text(((  snp_pos +       31    ),(int(0.67*height))), msg_2, font=fnt4, fill=(0,0,0,255))

			# Deprecated
			#draw.text((int(snp_pos - (int(0.036*wide) + 16*(len(str(p[4][3]).strip())-1) )), int(0.67*height)), (
			#	str(p[4][3]) +   '        '  +
			#	str(p[4][4])), font=fnt4, fill=(0,0,0,255))   


			#Arrows
			#Base
			image_file = StringIO(open("./fonts/arrowright.png",'rb').read())
			arrow = Image.open(image_file)
			arrow = arrow.resize((47, 47), Image.ANTIALIAS)
			arrow = arrow.crop((12, 17, 47, 34))
			im.paste(arrow, (int(snp_pos - int(0.013*wide)), int(0.685*height)))

			#Aa
			if aach == "yes":					
				image_file = StringIO(open("./fonts/arrowright.png",'rb').read())
				arrow = Image.open(image_file)
				arrow = arrow.resize((47, 47), Image.ANTIALIAS)
				arrow = arrow.crop((12, 17, 47, 34))

				#We paste the arrow in a different position depending of the number of characters that the aa position has:
				if str_len == 1: im.paste(arrow, (int(snp_pos - int(0.041*wide)), int(0.765*height)))
				if str_len == 2: im.paste(arrow, (int(snp_pos - int(0.030*wide)), int(0.765*height)))
				if str_len == 3: im.paste(arrow, (int(snp_pos - int(0.019*wide)), int(0.765*height)))
				if str_len == 4: im.paste(arrow, (int(snp_pos - int(0.008*wide)), int(0.765*height)))
				if str_len == 5: im.paste(arrow, (int(snp_pos + int(0.003*wide)), int(0.765*height)))

		#save image, specifying the format with the extension. For SNP images we save them with diferent sizes depending on if theres an aminoacid change or not
		w, h = im.size
		if args.my_mut == 'lin':
			im.crop((70, 100, w-20, h-60)).save(project + '/3_workflow_output/gene_plot_' + str(args.my_mut) + '_' + str(p[2]) + '_gene_' + str(p[3])+ '.png')

		if args.my_mut == 'snp' and aach == 'no':
			im.crop((70, 100, w-20, h-70)).save(project + '/3_workflow_output/gene_plot_' + str(args.my_mut) + '_' + str(p[1]) + '_gene_' + str(p[3])+ '.png')
		
		if args.my_mut == 'snp' and aach == 'yes':
			im.crop((70, 100, w-20, h-40)).save(project + '/3_workflow_output/gene_plot_' + str(args.my_mut) + '_' + str(p[1]) + '_gene_' + str(p[3])+ '.png')
