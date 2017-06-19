#!/usr/bin/python
import sys
import os
import csv
import warnings
import argparse
from Bio.PDB import PDBParser

# Three letter to one letter amino acid code conversion
RESDICT = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
           'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
           'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
           'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

def inv_sq_distance(coord1, coord2):
    '''
    Doc string
    '''
    distance = 0.0
    for i in range(len(coord1)):
        distance += (coord1[i]-coord2[i])**2
    return 1/distance

def calculate_wcn(residues):
    '''
    Doc string
    '''
    for residue in residues:
        wcn_ca = 0
        wcn_sc = 0
        for other_residue in residues:
            if residue != other_residue:
                wcn_ca += inv_sq_distance(residue['coord_ca'],
                                          other_residue['coord_ca'])
                wcn_sc += inv_sq_distance(residue['sidechain_center'],
                                          other_residue['sidechain_center'])
        residue['wcn_ca'] = wcn_ca
        residue['wcn_sc'] = wcn_sc

    return residues

def process_residue(residue, structure):
    '''
    Doc string
    '''
    output_dict = {}
    atoms_seen = []
    output_dict['amino_acid'] = RESDICT[residue.resname]
    output_dict['residue'] = residue.get_id()[1]
    output_dict['chain'] = residue.get_full_id()[2]
    sidechain_coords = []
    for atom in structure.get_atoms():
        # atom.name is equivalent to atom.get_id()
        if atom.parent.id != residue.id:
            continue
        atoms_seen.append(atom.name)
        if atom.name == 'CA':
            output_dict['coord_ca'] = atom.get_coord()
        if atom.name not in ['C', 'CA', 'O', 'N']:
            sidechain_coords.append(atom.get_coord())

    warning_message = "Missing {} in residue (" + \
                        str(output_dict['residue']) + ", " + \
                        str(output_dict['amino_acid']) + ")"

    for mainchain_atom in ['N', 'C', 'O']:
        if mainchain_atom not in atoms_seen:
            warnings.warn(warning_message.format(mainchain_atom),
                          RuntimeWarning)
    if 'coord_ca' not in output_dict:
        raise RuntimeError(warning_message.format('CA') +
                           '. Cannot calculate C-alpha WCN.')

    if len(sidechain_coords) == 0:
        if output_dict['amino_acid'] != 'G':
            warnings.warn(warning_message.format('sidechain') +
                          '. Using CA instead.', RuntimeWarning)
        sidechain_coords.append(output_dict['coord_ca'])

    # Calculate side chain properties:
    output_dict['sidechain_size'] = len(sidechain_coords)
    output_dict['sidechain_center'] = sum(sidechain_coords)/\
                                        output_dict['sidechain_size']

    return output_dict

def collect_coordinates(structure):
    '''
    Doc string
    '''
    output_list = []
    for residue in structure.get_residues():
        het_flag = residue.get_id()[0]
        if het_flag[0:2] == 'H_' or het_flag[0] == 'W':
            # Skip hetatoms and waters in PDB
            continue
        output_list.append(process_residue(residue, structure))
    return output_list

def main():
    '''
    Doc string
    '''
    parser = argparse.ArgumentParser(
        description='Calculate WCN values for an input PDB.')
    parser.add_argument('pdb', metavar='<PDB path>', type=str,
                        help='input pdb file')
    parser.add_argument('-o', metavar='<output prefix>', type=str,
                        help='prefix for output files')
    args = parser.parse_args()
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    # Define output file names
    if args.o is None:
        # If no output prefix given, assign prefix using input filename
        args.o = pdb_name
    output_wcn = args.o + '.wcn.csv'

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_name.upper(), args.pdb)
    output_list = collect_coordinates(structure)
    output_list = calculate_wcn(output_list)

    with open(output_wcn, 'w') as csvfile:
        writer = csv.DictWriter(csvfile,
                                fieldnames=['residue', 'chain', 'amino_acid',
                                            'wcn_sc', 'wcn_ca'],
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_list)

if __name__ == "__main__":
    main()
