# This script formats the SNP information for varanalyzer
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a', action="store", dest='input')
parser.add_argument('-b', action="store", dest='output')

args = parser.parse_args()

#input
f1 = open(args.input, 'r')
lines = f1.readlines()

#output
f2 = open(args.output, 'w')
f2.write('#data\tcontig\tpos\tref\talt\tqual\tref_count\talt_count\taf\n')

for line in lines:
	sp = line.split()

	allele_frequency = float(sp[6].strip()) / ( float(sp[5].strip()) + float(sp[6].strip()) )
	f2.write('snp' + '\t' + sp[0].strip() + '\t' + sp[1].strip() + '\t' + sp[2].strip() + '\t' + sp[3].strip() + '\t' + sp[4].strip() + '\t' + sp[5].strip() + '\t' + sp[6].strip() + '\t' + "{0:.2f}".format(allele_frequency) + '\n')
