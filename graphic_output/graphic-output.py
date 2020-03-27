# This script calls on functions from draw.py to draw the images required for each annalysis

#from __future__ import division
import argparse
from draw import fa_vs_pos, insertions_overview_and_histograms, gene_plot, legend, candidates_zoom, gene_legend
parser = argparse.ArgumentParser()

#INPUT VARIABLES FOR SNP
parser.add_argument('-my_mut', action="store", dest = 'my_mut')			#snp or lin
parser.add_argument('-asnp', action="store", dest = 'input_snp')		
parser.add_argument('-bsnp', action="store", dest = 'input_f_snp')		#Fasta genome input
parser.add_argument('-snp_analysis_type', action="store", dest='my_snp_analysis_type')
parser.add_argument('-cross', action="store", dest='my_cross')
parser.add_argument('-interval_width', action="store", dest='interval_width')

#INPUT VARIABLES FOR LIN
parser.add_argument('-a', action="store", dest = 'input')		
parser.add_argument('-b', action="store", dest = 'input_f')				#Fasta genome input
parser.add_argument('-gff', action="store", dest = 'gff')				#Genome feature file
parser.add_argument('-m', action="store", dest = 'mode', default = 'pe')
parser.add_argument('-ins_pos', action="store", dest = 'ins_pos')

#SHARED VARIABLES
parser.add_argument('-iva', action="store", dest = 'input_va')	 		#Output de varanalyzer
parser.add_argument('-rrl', action="store", dest = 'rrl') 				#Regulatory region lenght
parser.add_argument('-pname', action="store", dest='project_name')

args = parser.parse_args()
project = args.project_name


if args.my_mut == 'af_control' or args.my_mut == 'af_sample' or args.my_mut == 'af_candidates' :
	#FA vs POS image
	fa_vs_pos()

if args.my_mut == 'snp':
	#FA vs POS image
	fa_vs_pos()

	#Candidates
	candidates_zoom()	

	#Gene gene_plot
	gene_plot()	
	#Legend
	#legend()


elif args.my_mut == 'lin':
	#Insertions overview and histogram
	insertions_overview_and_histograms()

	#Gene gene_plot
	gene_plot()

#gene_legend()
