import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-input', action="store", dest='input')
parser.add_argument('-out', action="store", dest='out')

args = parser.parse_args()
out = open(args.out, "w")
for line in open(args.input, "r"):
	new_line = line.strip("\n") +"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\n"
	out.write(new_line)
