#!/usr/bin/python

'''
This script takes in aligned amino acid msa (in fasta format) and converts it to nucleotide msa (also in fasta format) while preserving the alignment

Author: Dariya K. Sydykova
'''

import argparse
import sys
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
    
def back_translate(aa_aln,codon_dict,out_aln_file,gencode):
	new_codon_aln=[] #create an empty list to store new sequences
	for aa_seq in aa_aln: #iterate over sequences in the alignment
		new_codon_seq="" #initialize an empty string to store the new back translated sequence
		
		if aa_seq.id in codon_dict: #check if sequence ID from amino acid alignment is present in the codon sequence file
			old_codon_rec=codon_dict[aa_seq.id] #get the codon sequence 
		else:
			print("Amino acid alignment sequence IDs do not match codon sequence IDs")
			sys.exit()
			
		old_codon_seq=str(old_codon_rec.seq) #turn the codon sequence into a string
		if "-" in old_codon_seq: #remove gaps if there are any
			old_codon_seq = old_codon_seq.replace("-","")
					
		i=0
		for aa in aa_seq.seq: 
			if aa=="-": #change the amino acid gaps to match codon gaps
				new_codon_seq += "---"
			else: 
				codon = old_codon_seq[i:i+3] #get the codon that corresponds to the amino acid
				if gencode[codon]!=aa: #check if the codon belongs to the amino acid
					print('codon does not match the aa')
					sys.exit()
				else: #add the codon to the back translated sequence if it belongs to the correct amino acid
					new_codon_seq += codon
				i+=3 
		
		codon_seq=Seq(new_codon_seq) #turn back translate sequence string into a Seq object
		codon_seq_rec = SeqRecord(codon_seq) #create a SeqRecord object from Seq
		codon_seq_rec.id = aa_seq.id #set SeqRecord.id to be printed in the new fasta file
		codon_seq_rec.description = '' #set SeqRecord.description to an empty string to avoid <unknown description> being printed in the new fasta file 
		new_codon_aln.append(codon_seq_rec) #add a new codon sequence to the list
	
	new_codon_msa = MultipleSeqAlignment(new_codon_aln) # create an MSA object
	AlignIO.write(new_codon_msa, out_aln_file, "fasta") #write a new alignment

def main():
	'''
	Back translate aligned amino acid sequences to nucleotide sequences while maintaining the alignment
	'''
	#creating a parser
	parser = argparse.ArgumentParser(
	        formatter_class=argparse.RawDescriptionHelpFormatter,
			description='Translate amino acid sequences into nucleotide sequences while maintaining the amino acid alignment and the order of sequences.')
	#adding arguments 
	parser.add_argument('-a', metavar='<aa_aln.fasta>', type=str, 
	help='input amino acid alignment file')
	parser.add_argument('-n', metavar='<nuc_aln.fasta>', type=str, 
	help='input unaligned/aligned nucleotide/codon sequence file')
	parser.add_argument('-o', metavar='<codon_aln.fasta>', type=str, 
    help='output codon alignment file') 
	args = parser.parse_args()
    
    #set up output file name if none is given
	if args.o is None:
		nuc_aln_file = "codon_aln.fasta"
	else:
		nuc_aln_file = args.o
    
	aa_aln_file = args.a
	nuc_seq_file = args.n
	
	aa_aln = AlignIO.read(aa_aln_file, "fasta") 

	nuc_seq = open(nuc_seq_file) #read in codon file
	codon_dict = SeqIO.to_dict(SeqIO.parse(nuc_seq, "fasta")) #read in codon sequence file as a dictionary, key=sequence ID, value=sequence itself

	back_translate(aa_aln,codon_dict,nuc_aln_file,gencode)

if __name__ == "__main__":
	main()