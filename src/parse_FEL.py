#!/usr/bin/python

'''
Parse HyPhy output file =  FEL JSON into a CSV, **assuming a single partition**

Author: Stephanie J. Spielman, Dariya K. Sydykova
'''
import json, argparse

def parse_json(fel_json_file, fel_csv_file):

	with open (fel_json_file , "rU") as f:
		parsed_json = json.load(f)

	site_block =  parsed_json[ "MLE" ] 
	raw_header = site_block[ "headers" ]
	raw_content = site_block[ "content" ][ "0" ]
	
	final_header = "Site," + ",".join( [x[0].replace(" ","_") for x in raw_header] ) + "\n"
		
	with open(fel_csv_file, "w") as f:
		f.write(final_header)
		site = 1
		for row in raw_content:
			f.write( str(site) + "," + ",".join(str(x) for x in row) + "\n")
			site += 1

def main():
	'''
	Parse FEL json file to write csv file with site-wise dN/dS.
	'''
	
	#creating a parser
	parser = argparse.ArgumentParser(description='Parse FEL json file to write csv file with site-wise dN/dS.')
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
