import sys

def main(argv):

	if len(argv) != 3: # wrong number of arguments
		print("""Usage:
format_aln_id.py <fasta_file> <reformated_fasta_file>
""")
		sys.exit()

	fasta_file=argv[1]
	out_fasta_file=argv[2]
	
	f=open(fasta_file,"r")
	out=open(out_fasta_file,"w")
	for line in f:
		if line.startswith(">"):
			if "." in line:
				line=line.replace(".","_")
			if "|" in line:
				line=line.replace("|","_")
			out.write(line)
		else:
			out.write(line)
			continue

if __name__ == "__main__":
	main(sys.argv)