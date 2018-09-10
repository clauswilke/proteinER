#!/usr/bin/python

'''
This script calculates site-specific dN/dS using HyPhy's FEL one-rate output. 

HyPhy's FEL one-rate method assigns dS=1 to sites with substitutions. This script simply takes a site's dN and dS values and divides them dN/dS. This script also finds convserved sites in a sequence alignment and assigns these sites dN/dS=0.  

Author: Dariya K. Sydykova
'''

# load packages required to run this script
import argparse
import textwrap
import sys
from Bio import AlignIO

# this function calculates site-wise dN/dS from site's dN and site's dS
def calc_dNdS(aln, rates_file, out_file):
    
    # open input and output files
    r = open(rates_file, "r")
    out = open(out_file, "w")
    
    # get a total number of sites
    total_col = len(aln[0]) 
    
    # create a list to store position of sites that are conserved
    conserved_sites_lst = []

    # loop over the alignment to find conserved sites
    for i in range(total_col):
        # get a list of all amino acids at site position i
        col = aln[:,i]
        
        # check if a site is conserved
        # the check is done by comparing the first amino acid at a site to the rest of the amino acids at that site
        if col == len(col) * col[0]:
            conserved_sites_lst.append(i+1) # if a site is conserved, add its position to the list
        else:
            continue # else, move on to the next site
    
    # loop over a file that contains HyPhy rates
    for line in r:
        
        # if the line is the header, add an extra column that will store dN/dS values
        if line.startswith("Site"):
            token = line.split(",")
            new_header = token[0] + ",dN/dS," + ",".join(token[4:])
            out.write(new_header)
            continue
        
        # make a list of individual values stored in each line
        token = line.split(",")
        
        # index site's position, dN, and dS values
        site = token[0]
        dS = float(token[1])
        dN = float(token[2])
        
        # calculate a site's dN/dS
        # if site dS = 0, set dN/dS = 0
        if abs(dS - 0.0) <= 0.000000000001:
            dN_dS = 0
        # else, calculate dN/dS by dividing site's dN by the site's dS
        else:
            # check if a site's dS was set to 1 by the HyPhy output
            # if it wasn't, this script will stop running because this script is intended to work only with HyPhy's one-rate inference method
            if abs(dS - 1.0) > 0.000000000001: 
                print("dS is not equal to 1")
                print("Check HyPhy output")
                sys.exit()
            else:
                dN_dS = dN/dS
            
        # if a sites position is stored in the list with conserved sites, set dN/dS = 0 and write a new line that contains all of the HyPhy output and dN/dS of 0
        if site in conserved_sites_lst:
            new_line = site+",0,"+",".join(token[4:])
        # else, write a new line that contains all of the HyPhy output and dN/dS calculated earlier
        else:
            new_line = site+","+str(dN_dS)+","+",".join(token[4:])
        
        # write the new line to the output file
        out.write(new_line)

def main():
    '''
    Calculate site-specific dN/dS from HyPhy's FEL one-rate model output
    '''
    # creating a parser
    parser = argparse.ArgumentParser(
            formatter_class = argparse.RawDescriptionHelpFormatter,
            description = "Calculate site-specific dN/dS from HyPhy's FEL one-rate model output",
            epilog = textwrap.dedent('''\
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

    # adding arguments 
    parser.add_argument('-a', metavar = '<aa_aln.fasta>', type = str, help = 'input amino acid alignment file')
    parser.add_argument('-r', metavar = '<rates.csv>', type = str, help = 'HyPhy FEL file')
    parser.add_argument('-o', metavar = '<processed_rates.csv>', type = str, help = 'output processed rates file')

    args = parser.parse_args()

    # set up output file name if none is given
    if args.o is None:
        out_rates = "processed_dNdS.csv"
    else:
        out_rates = args.o
        
    # store inputted arguments as variables
    aln_file = args.a
    rates_file = args.r

    # read an alignment file in FASTA format
    aln = AlignIO.read(aln_file, "fasta")
    
    # calculate site-wise dN/dS 
    calc_dNdS(aln, rates_file, out_rates)

if __name__ == "__main__":
    main()