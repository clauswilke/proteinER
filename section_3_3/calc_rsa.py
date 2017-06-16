#!/usr/bin/python

'''
This script parses an input PDB file and returns relative solvent accessibility
values (RSA) and raw DSSP output.

Author: Benjamin R. Jack
'''

import subprocess
import argparse
import csv

# Three letter to one letter amino acid code conversion
RESDICT = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
           'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
           'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
           'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

# ASA normalization constants were taken from:
# M. Z. Tien, A. G. Meyer, D. K. Sydykova, S. J. Spielman, C. O. Wilke (2013).
# Maximum allowed solvent accessibilities of residues in proteins. PLOS ONE
# 8:e80635.
RES_MAX_ACC = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
               'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
               'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
               'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
               'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0}

def run_dssp(pdb_path, output_dssp):
    '''
    Run mkdssp.
    '''
    command = ['mkdssp', '-i', pdb_path, '-o', output_dssp]
    process = subprocess.call(command)
    return process

def parse_dssp_line(line):
    '''
    Extract values from a single line of DSSP output and calculate RSA.
    '''
    solvent_acc = int(line[35:39])  # record SA value for given AA
    amino_acid = line[13].strip()  # retrieve amino acid
    residue = line[6:10].strip()
    chain = line[11].strip()
    if amino_acid.islower():
        # if the AA is a lowercase letter, then it's a cysteine
        amino_acid = "C"
    if amino_acid in ['X', '!', '*']:
        rsa = 0
    else:
        max_acc = RES_MAX_ACC[amino_acid]  # Find max SA for residue
        rsa = solvent_acc / max_acc # Normalize SA
    return residue, amino_acid, chain, rsa

def parse_dssp(raw_dssp_output):
    '''
    Parse a DSSP output file and return a dictionary of amino acids, PDB residue
    number, chain, and RSA.
    '''
    with open(raw_dssp_output, 'r') as dssp_file:
        lines = dssp_file.readlines()
        output = {'residue': [],  # pdb residue numbers
                  'amino_acid': [],  # amino acid
                  'chain': [],  #  chain
                  'rsa': []}  #  solvent accessibiity values
    for line in lines[28:]:
        # Skip first 28 lines of header info
        residue, amino_acid, chain, rsa = parse_dssp_line(line)
        if amino_acid not in ['*', '!']:
            # Append data to lists
            output['amino_acid'].append(amino_acid)
            output['residue'].append(residue)
            output['chain'].append(chain)
            output['rsa'].append(rsa)

    return output

def main():
    '''
    Run mkdssp and parse output from input PDB.
    '''
    parser = argparse.ArgumentParser(
        description='Calculate RSA values for an input PDB.')
    parser.add_argument('pdb', metavar='<PDB path>', type=str,
                        help='input pdb file')
    parser.add_argument('prefix', metavar='<output prefix>', type=str,
                        help='prefix for output files')
    args = parser.parse_args()
    # Define output file names
    output_rsa = args.prefix + '.rsa.csv'
    asa_file = args.prefix + '.asa.txt'
    if run_dssp(args.pdb, asa_file):
        # Check for DSSP errors
        raise RuntimeError("Call to DSSP failed.")
    else:
        # DSSP succeeded, write output to CSV
        output_dict = parse_dssp(asa_file)
        with open(output_rsa, 'w') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(output_dict.keys())
            output_dict_zip = zip(*output_dict.values())
            writer.writerows(output_dict_zip)

if __name__ == "__main__":
    main()
