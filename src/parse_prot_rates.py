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


	site_block =  parsed_json[ "Relative site rate estimates" ]
	final_header = "Site,Rate,lower_bound_95CI,upper_bound_95CI\n"
	with open(relprot_csv_file, "w") as f:
		f.write(final_header)
		for site in range(1, len(site_block)+1):
			rate = site_block[str(site)]
			f.write( ",".join( [ str(site), str(rate["MLE"]), str(rate["LB"]), str(rate["UB"]) ] ) + "\n")

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
	parser.add_argument('-j', metavar='<site-rates.json>', type=str, help='input FEL json file')
	parser.add_argument('-r', metavar='<rates.csv>', type=str, help='output csv file with site-wise rates')

	args = parser.parse_args()

	#set up output file name if none is given
	if args.r is None:
		out_rates = "site-rates.csv"
	else:
		out_rates = args.r
		
	json_file=args.j
	rates_file=args.r

	parse_json(json_file, rates_file)
	

if __name__ == "__main__":
	main()     
