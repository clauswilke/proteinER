#!/usr/bin/python

'''
This script takes in nucleotide msa (in fasta format) and converts it to amino acid msa (also in fasta format)

Author: Dariya K. Sydykova
'''

import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

def translate(nuc_aln_file,aa_aln_file):	
	nuc_aln = SeqIO.parse(nuc_aln_file, "fasta") 
	
	aa_seq_records = list()
	for nuc_seq in nuc_aln:
		aa_seq = nuc_seq.seq.translate(to_stop=True) #translate nucleotide seq to amino acids, returns a Seq object
		aa_seq_r = SeqRecord(aa_seq) #create a SeqRecord object from amino acid Seq object to start a list 
		aa_seq_r.id = nuc_seq.id #set SeqRecord.id to be printed in the new fasta file
		aa_seq_r.description = '' #set SeqRecord.description to an empty string to avoid <unknown description> being printed in the new fasta file 
		aa_seq_records.append(aa_seq_r) 

	#write sequences
	SeqIO.write(aa_seq_records, aa_aln_file, "fasta") 

def main():
	'''
	Translate nucleotide sequences to amino acid sequences.
	'''
	#creating a parser
	parser = argparse.ArgumentParser(description='Translate nucleotide sequences into amino acid sequences')
	#adding arguments 
	parser.add_argument('-i', metavar='<nuc_aln.fasta>', type=str,
                    help='input nucleotide alignment file')
	parser.add_argument('-o', metavar='<aa_aln.fasta>', type=str,
                    help='output amino acid alignment file')

	args = parser.parse_args()

	#set up output file name if none is given
	if args.o is None:
		aa_aln_file = "aa_aln.fasta"
	else:
		aa_aln_file = args.o
		
	nuc_aln_file = args.i

	translate(nuc_aln_file,aa_aln_file)

if __name__ == "__main__":
	main()