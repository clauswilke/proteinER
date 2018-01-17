#!/usr/bin/python

'''
This script reformats sequence IDs in the FASTA file. Specifically, the script replaces "." and "|" with a "_". The sequences are left entirely untouched. The output keeps the FASTA format and the order of the sequences from the input file. 

Author: Dariya K. Sydykova
'''

import argparse
import textwrap

def format_aln(in_aln_file, out_aln_file):
	
	f=open(in_aln_file,"r")
	out=open(out_aln_file,"w")
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
	Reformat sequence IDs in the FASTA file
	'''
	#creating a parser
	parser = argparse.ArgumentParser(
	        formatter_class=argparse.RawDescriptionHelpFormatter,
			description='Reformat sequence IDs in the FASTA file',
	        epilog=textwrap.dedent('''\
            This script keeps the FASTA format and the order of the sequences 
            '''))

	#adding arguments 
	parser.add_argument('-a', metavar='<aln.fasta>', type=str, help='input alignment file')
	parser.add_argument('-o', metavar='<reformated_aln.fasta>', type=str, help='output alignment file')

	args = parser.parse_args()

	#set up output file name if none is given
	if args.o is None:
		out_aln_file = "reformatted_aln.fasta"
	else:
		out_aln_file = args.o
		
	in_aln_file=args.a

	format_aln(in_aln_file, out_aln_file)

if __name__ == "__main__":
	main()