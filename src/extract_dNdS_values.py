#!/usr/bin/python

'''
This script extracts dN/dS values from HyPhy output file for the first sequence in the alignment

Author: Dariya K. Sydykova
'''
import argparse, sys, os
from Bio import SeqIO

def extract_dNdS(aln_file, rates_file, out_rates):	
	records = list(SeqIO.parse(aln_file, "fasta")) #read fasta file
	first_seq=records[0].seq #get the first sequence in the alignment and turn it into a Seq object
	
	r=open(rates_file,"r") #read the rates file
	out=open(out_rates,"w") #open output file
	
	for line in r:
		print(line)
		print(len(first_seq))
		token=line.split(",")
		if line.startswith("Site"):
			out.write('fasta_position,fasta_aa,'+",".join(token[1:])) #write a new heading
			continue
		
		site=int(token[0])
		if first_seq[site-1]!='-': #if the site is not a gap, write the fasta position, amino acid, and dN/dS value to the output file
			aa=first_seq[site-1]
			
			new_line=str(site)+','+aa+','+",".join(token[1:])
			out.write(new_line)
	
	total_sites=len(first_seq)
	if total_sites!=site:
		print('The length of the first sequence does not match the total output positions')
		os.remove(out_rates)
		sys.exit()

def main():
	'''
	Extract dN/dS values from HyPhy output for non-gap sites
	'''
	#creating a parser
	parser = argparse.ArgumentParser(description='Assign dN/dS of 0 to conserved sites.')
	#adding arguments 
	parser.add_argument('-a', metavar='<aa_aln.fasta>', type=str, help='input amino acid alignment file')
	parser.add_argument('-r', metavar='<rates.csv>', type=str, help='HyPhy FEL1 file')
	parser.add_argument('-o', metavar='<processed_rates.csv>', type=str, help='output processed rates file')

	args = parser.parse_args()

	#set up output file name if none is given
	if args.o is None:
		out_rates = "processed_"+args.r
	else:
		out_rates = args.o
		
	aln_file=args.a
	rates_file=args.r

	extract_dNdS(aln_file, rates_file, out_rates)
	

if __name__ == "__main__":
	main()