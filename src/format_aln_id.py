#!/usr/bin/python

'''
This script finds "." and "|" is sequence IDs and converts these characters to "_". 
This scripts works strictly with alignment files in FASTA format.  

Author: Dariya K. Sydykova
'''

import argparse

def format(in_file, out_file):

	f=open(in_file,"r")
	out=open(out_file,"w")
	for line in f:
		if line.startswith(">"):
			if "." in line:
				line=line.replace(".","_")
			if "|" in line:
				line=line.replace("|","_")
			out.write(line)
		else:
			out.write(line)
			continue
			
def main():
	'''
	format sequence IDs in a FASTA file
	'''
	#creating a parser
	parser = argparse.ArgumentParser(description='Reformat sequence IDs in FASTA files.')
	#adding arguments 
	parser.add_argument('-i', metavar='<raw_aln.fasta>', type=str,
                    help='input alignment file')
	parser.add_argument('-o', metavar='<reformatted_aln.fasta>', type=str,
                    help='output alignment file with reformatted sequence IDs')

	args = parser.parse_args()

	#set up output file name if none is given
	if args.o is None:
		reformatted_fasta = 'reformatted_'+args.i
	else:
		reformatted_fasta = args.o
		
	raw_fasta=args.i
	
	format(raw_fasta,reformatted_fasta)

if __name__ == "__main__":
    main()