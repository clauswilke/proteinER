#!/usr/bin/python

'''
This script changes HyPhy's FEL one-rate dN/dS values. The script finds convserved sites in a sequence alignment and assigns these sites dN/dS of 0.  

Author: Dariya K. Sydykova
'''

import argparse
import textwrap
from Bio import AlignIO

def fix_dNdS(aln, rates_file, out_file):

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
			conserved_sites_lst.append(i+1)
		else:
			continue
	
	for line in r:
		if line.startswith("Site"):
			token=line.split(",")
			new_header = token[0]+",dN/dS,"+",".join(token[4:])
			out.write(new_header)
			continue
		
		token=line.split(",")
		site=token[0]
		dS=float(token[1])
		dN=float(token[2])
		
		if dS==0:
			dN_dS=0
		elif dS==1:
			dN_dS = dN
		else:
			print("dS is not equal 1")
			sys.exit()
		
		if site in conserved_sites_lst:
			new_line = site+",0,"+",".join(token[4:])
		else:
			new_line = site+","+str(dN_dS)+","+",".join(token[4:])
		
		out.write(new_line)

def main():
	'''
	Find conserved sites and assign them dN/dS=0
	'''
	#creating a parser
	parser = argparse.ArgumentParser(
	        formatter_class=argparse.RawDescriptionHelpFormatter,
			description='Find conserved sites and assign them dN/dS=0',
	        epilog=textwrap.dedent('''\
            This script produces a CSV with the following columns:
            
            Column name           Description
            ===================================================================
            Site                  Site number, extracted from the alignment 
                                  FASTA file.

            dN/dS                 Site-wise dN/dS, calculated from HyPhy output.
                                  dN='beta'
                                  dS='alpha'

            LRT                   Likelihood ration test statistic for 
                                  beta = alpha, versus beta does not equal alpha
                       
            p-value               Likelihood ration test statistic for 
                                  beta = alpha, versus beta does not equal alpha 
            
            Total_branch_length   The total length of branches contributing 
                                  to inference at this site, and used to 
                                  scale dN-dS
            '''))

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
	fix_dNdS(aln, rates_file, out_rates)

if __name__ == "__main__":
	main()