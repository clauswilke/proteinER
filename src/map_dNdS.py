#!/usr/bin/python

'''
This script extracts dN/dS values from HyPhy output file for the first sequence in the alignment

Author: Dariya K. Sydykova
'''
import argparse
import sys
import os
import textwrap
from Bio import SeqIO

def extract_dNdS(aln_file, rates_file, out_rates):	
	records = list(SeqIO.parse(aln_file, "fasta")) #read fasta file
	first_seq=records[0].seq #get the first sequence in the alignment and turn it into a Seq object
	
	r=open(rates_file,"r") #read the rates file
	out=open(out_rates,"w") #open output file
	
	fasta_pos=0
	for line in r:
		token=line.split(",")
		if line.startswith("Site"):
			out.write('fasta_position,fasta_aa,'+",".join(token[1:])) #write a new heading
			continue
		
		site=int(token[0])
		if first_seq[site-1]!='-': #if the site is not a gap, write the fasta position, amino acid, and dN/dS value to the output file
			aa=first_seq[site-1]
			fasta_pos+=1

			new_line=str(fasta_pos)+','+aa+','+",".join(token[1:])
			out.write(new_line)
	
	total_sites=len(first_seq.ungap("-"))
	if total_sites!=fasta_pos:
		print('The length of the first sequence does not match the total output positions')
		os.remove(out_rates)
		sys.exit()

def main():
	'''
	Extract dN/dS values from HyPhy output for non-gap sites
	'''
	#creating a parser
	parser = argparse.ArgumentParser(
		    formatter_class=argparse.RawDescriptionHelpFormatter,
			description='Extract dN/dS values from HyPhy output for non-gap sites',
	        epilog=textwrap.dedent('''\
            This script produces a CSV with the following columns:
            
            Column name           Description
            ===================================================================
            fasta_position        Site number of the reference sequence
                                  
            fasta_aa              Amino acid of the reference sequence
            
            dN/dS                 Site-wise dN/dS, calculated from HyPhy output.
                                  dN='beta'
                                  dS='alpha'

            LRT                   Likelihood ration test statistic for 
                                  beta = alpha, versus beta does not equal alpha
                       
            p-value               p-value for the LRT 
            
            Total_branch_length   The total length of branches contributing 
                                  to inference at this site, and used to 
                                  scale dN-dS
            '''))

	#adding arguments 
	parser.add_argument('-a', metavar='<aa_aln.fasta>', type=str, help='input amino acid alignment file')
	parser.add_argument('-r', metavar='<rates.csv>', type=str, help='HyPhy FEL file')
	parser.add_argument('-o', metavar='<trimmed_dNdS.csv>', type=str, help='output trimmed rates file')

	args = parser.parse_args()

	#set up output file name if none is given
	if args.o is None:
		out_rates = "trimmed_dNdS.csv"
	else:
		out_rates = args.o
		
	aln_file=args.a
	rates_file=args.r

	extract_dNdS(aln_file, rates_file, out_rates)
	

if __name__ == "__main__":
	main()