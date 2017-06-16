#!/usr/bin/python

'''
NEEDS DOCSTRING
'''

import subprocess
import sys


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
    command = 'mkdssp' + ' -i ' + pdb_path + ' -o ' + output_dssp
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

def parse_dssp(output_dssp):
    '''
    This calculates the RSA values for a PDB using DSSP. It returns a list of
    the amino acids and a list of their RSA values.
    '''
    with open(output_dssp, 'r') as dssp_file:
        lines = dssp_file.readlines()
        residue_number_list = []  # pdb residue numbers
        solvent_acc_list = []  # solvent accessibility (SA) values
        amino_acid_list = []  # single letter amino acids for each site
        chain_list = []
        no_rsa = 0
    for line in lines[28:]:
        solvent_acc = int(line[35:39])  # record SA value for given AA
        amino_acid = line[13]  # retrieve amino acid
        residue = line[6:10]
        chain = line[11]
        if amino_acid.islower():
            # if the AA is a lowercase letter, then it's a cysteine
            amino_acid = "C"
        if amino_acid == ('!' or '*'):
            # Ignore gaps or missing data that dssp might insert
            no_rsa = no_rsa + 1
            continue
        if amino_acid != 'X':
            max_acc = RES_MAX_ACC[amino_acid]  # Find max SA for residue
            solvent_acc_list.append(solvent_acc / max_acc)  # Normalize SA
        elif amino_acid == 'X':
            solvent_acc_list.append(0)
        # Append data to lists
        amino_acid_list.append(amino_acid)
        residue_number_list.append(residue)
        chain_list.append(chain)

    return (residue_number_list, chain_list, amino_acid_list, solvent_acc_list)


def main(argv):
    '''
    Docstring needed
    '''
    if len(argv) != 4:
        print("\n\nUsage:\n\n     "+argv[0]+" <input pdb file> <output dssp ASA> <output RSA>")
    else:
        pdb_path = argv[1]
        output_dssp = argv[2]
        output_rsa = argv[3]
        run_dssp(pdb_path, output_dssp)
        [residue_number_list,
         chain_list,
         amino_acid_list,
         solvent_acc_list] = parse_dssp(output_dssp)
        rsa_outputfile = open(output_rsa, 'w')
        rsa_outputfile.write('Residue,Amino_Acid,Chain,RSA\n')
        for i, amino_acid in enumerate(amino_acid_list):
            rsa_outputfile.write(str(residue_number_list[i]) + ',' +
                                 str(amino_acid) + ',' + chain_list[i] + ',' +
                                 str(solvent_acc_list[i]) + '\n')
        return


if __name__ == "__main__":
    main(sys.argv)
