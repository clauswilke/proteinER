#!/usr/bin/python

'''
This script fixes HyPhy's FEL1 site-wise dN/dS values. The script finds convserved sites in a sequence alignment and assigns these sites dN/dS of 0.  

Author: Dariya K. Sydykova
'''

import argparse
from Bio import AlignIO

def find_unchanged_sites(aln, rates_file, out_file):

	r=open(rates_file,"r")
	out=open(out_file,"w")
	
	total_col=len(aln[0]) #total number of sites
	
	conserved_sites_lst=[]
	for i in range(total_col):
		#a list of all amino acids at site i
		col = aln[:,i]
		
		#checks if the list of amino acids at site i is identical to 
		#the first amino acid in a list repeated the number of times the list is
		if col == len(col) * col[0]:
			conserved_sites_lst.append(True)
		else:
			conserved_sites_lst.append(False)
		
	site=0
	for line in r:
		if line.startswith("dN/dS"):
			out.write(line)
			continue
		
		conserved = conserved_sites_lst[site]
		token=line.split(",")
		
		if conserved:
			new_line = str(0)+","+",".join(token[1:])
		else:
			new_line = line
		
		out.write(new_line)
		site+=1

def main():
	'''
	find conserved sites in FASTA file and assign dN/dS of 0 to those sites.
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

	aln = AlignIO.read(aln_file, "fasta") 
	find_unchanged_sites(aln, rates_file, out_rates)

if __name__ == "__main__":
	main()