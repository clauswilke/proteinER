#!/usr/bin/python

'''
This script parses the raw rate4site output into a CSV.

Author: Benjamin Jack
'''

import argparse
import os
import textwrap
import pandas as pd

def main():
    '''
    Parse a fixed-width rate4site output file into a CSV.
    '''
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Extract rate values from raw rate4site output and write rates to a CSV.',
        epilog=textwrap.dedent('''\
            This script produces a CSV with the following columns: 
			
            Column name     Description
            ===================================================================
            fasta_position  (Defined in Rate4Site file as POS column)
                            Site number, extracted from the alignment 
                            FASTA file
            
            fasta_aa        (Defined in Rate4Site file as SCORE column)
                            The amino acid in the reference sequence in one 
                            letter code.
            
            r4s_rate        (Defined in Rate4Site file as SCORE column)
                            The conservation scores. lower value = higher 
                            conservation.
            '''))    
    
    parser.add_argument('rates', metavar='<r4s_rates>', type=str,
                        help='rate file output from rate4site')
    parser.add_argument('-o', metavar='<output file>', type=str,
                        help='name of output file')
    args = parser.parse_args()

    if args.o is None:
        outfile = 'extracted_' + \
            os.path.splitext(os.path.basename(args.rates))[0] + '.csv'
    else:
        outfile = args.o
    
    # Import r4s output as dataframe
    rates = pd.read_fwf(args.rates, 
                        skiprows=13, # Skip r4s header junk
                        skipfooter=2, # Skip mean and std dev footer
                        widths=[5, 5, 9], # Specifiy column widths
                        usecols=[0,1,2], # Grab the first 4 columns
                        names=['fasta_position', 'fasta_aa', 'r4s_rate'])
    # Write dataframe to file
    rates.to_csv(outfile, index=False)

if __name__ == "__main__":
    main()
    
    
