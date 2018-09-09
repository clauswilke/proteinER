#!/usr/bin/python

'''
This script reformats sequence IDs in the FASTA file. Specifically, the script replaces "." and "|" with a "_". The sequences are left entirely untouched. The output keeps the FASTA format and the order of the sequences from the input file. 

Author: Dariya K. Sydykova
'''

# load packages required to run this script
import argparse
import textwrap

# this function replaces "." and "|" with "_" in sequence IDs of an alignment
def format_aln(in_aln_file, out_aln_file):
    
    # open input and output files
    f = open(in_aln_file,"r")
    out = open(out_aln_file,"w")

    # loop over the lines in the input file
    for line in f:

        # check if a line starts with ">"
        # files in FASTA format contain ">" that are followed by a sequence ID
        # if a line starts with ">", reformat the sequence ID contained in that line
        if line.startswith(">"):
            # if the line contains ".", replace every instance of "." with "_"
            if "." in line:
                line=line.replace(".","_")
            # if the line contains "|", replace every instance of "|" with "_"
            if "|" in line:
                line=line.replace("|","_")
            # if the line doesn't contain "." or "|", move on to the next line
            out.write(line)
        # if a line doesn't start with ">", write the line to the output file
        else:
            out.write(line)
            continue

def main():
    '''
    Reformat sequence IDs in the FASTA file
    '''
    #creating a parser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description = 'Reformat sequence IDs in the FASTA file',
            epilog = textwrap.dedent('''\
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

    # store inputted arguments as variables
    in_aln_file = args.a

    # reformat alignment's sequence IDs
    format_aln(in_aln_file, out_aln_file)

if __name__ == "__main__":
    main()