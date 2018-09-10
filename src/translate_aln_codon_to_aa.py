#!/usr/bin/python

'''
This script takes in a nucleotide multiple sequence alignment (MSA) in FASTA format and converts it to an amino acid MSA in FASTA format.

Author: Dariya K. Sydykova
'''

# load packages required to run this script
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

# this function translates nucleotide sequence to amino acid sequence
# the inputted sequence can be aligned or not aligned
def translate(nuc_aln_file, aa_aln_file):
    # read in a nucleotide MSA as a SeqIO object
    nuc_aln = SeqIO.parse(nuc_aln_file, "fasta") 
    
    # create an empty list to store translated amino acid sequences
    aa_seq_records = list()
    
    # loop over nucleotide sequence present in nucleotide MSA
    for nuc_seq in nuc_aln:
        # translate nucleotide sequence to amino acids
        aa_seq = nuc_seq.seq.translate(to_stop=True) # returns a Seq object
        # create a SeqRecord object from amino acid Seq object to start a list 
        aa_seq_r = SeqRecord(aa_seq)
        # set SeqRecord.id to be printed in the new fasta file
        aa_seq_r.id = nuc_seq.id
        # set SeqRecord.description to an empty string to avoid <unknown description> being printed in the new fasta file 
        aa_seq_r.description = ''
        # add SeqRecrod object to the list
        aa_seq_records.append(aa_seq_r) 

    # write an amino acid MSA in FASTA format
    SeqIO.write(aa_seq_records, aa_aln_file, "fasta") 

def main():
    '''
    Translate nucleotide sequences to amino acid sequences.
    '''
    # creating a parser
    parser = argparse.ArgumentParser(description = 'Translate nucleotide sequences (either aligned or not aligned) into amino acid sequences')
    
    # adding arguments 
    parser.add_argument('-n', metavar = '<nuc_aln.fasta>', type = str, help = 'input nucleotide/codon alignment file')
    parser.add_argument('-o', metavar = '<aa_aln.fasta>', type = str, help = 'output amino acid alignment file')

    args = parser.parse_args()

    # set up output file name if none is given
    if args.o is None:
        aa_aln_file = "aa_aln.fasta"
    else:
        aa_aln_file = args.o

    # store inputted arguments as variables
    nuc_aln_file = args.n

    # translate inputted alignment
    translate(nuc_aln_file, aa_aln_file)

if __name__ == "__main__":
    main()