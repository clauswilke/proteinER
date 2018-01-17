#!/usr/bin/python

'''
Parse HyPhy output file = protein rate JSON into a CSV

Author: Stephanie J. Spielman, Dariya K. Sydykova
'''

import json
import argparse
import textwrap

def parse_json(relprot_json_file, relprot_csv_file):

    with open (relprot_json_file , "rU") as f:
        parsed_json = json.load(f)


    site_block =  parsed_json[ "MLE" ][ "content" ][ "0" ]
    final_header = "Site,Rate,lower_bound_95CI,upper_bound_95CI\n"
    with open(relprot_csv_file, "w") as f:
        f.write(final_header)
        for site in range(len(site_block)):
            site_row = site_block[site]
            f.write( ",".join( [ str(site+1), str(site_row[0]), str(site_row[1]), str(site_row[2]) ] ) + "\n")

def main():
    '''
    Parse relative protein rates json file to write csv file with site-wise rates.
    '''
    
    #creating a parser
    parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='Parse relative protein rates json file to write csv file with site-wise rates.',
    epilog=textwrap.dedent('''\
            This script produces a CSV with the following columns: 
            
            Column name       Description
            ===================================================================
            Site              Site number, extracted from the alignment FASTA 
                              file
                                  
            Rate              Maximum likelihood estimate of the rate
            
            lower_bound_95CI  Lower bound of the 95% confidence interval of the
                              rate estimate
            
            upper_bound_95CI  Upper bound of the 95% confidence interval of the
                              rate estimate      
            '''))
    #adding arguments 
    parser.add_argument('-j', metavar='<LEISR.json>', type=str, help='input LEISR json file')
    parser.add_argument('-r', metavar='<rates.csv>', type=str, help='output csv file with site-wise rates')

    args = parser.parse_args()

    #set up output file name if none is given
    if args.r is None:
        out_rates = "site_rates.csv"
    else:
        out_rates = args.r
        
    json_file=args.j
    rates_file=out_rates

    parse_json(json_file, rates_file)
    

if __name__ == "__main__":
    main()     
