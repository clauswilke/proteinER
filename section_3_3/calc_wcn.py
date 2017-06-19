#!/usr/bin/python
import sys, subprocess, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.PDB.Polypeptide import *
from Bio.PDB import *
#from Bio.PDB import PDBParser
#from Bio.PDB import PDBIO
#from Bio.PDB import Dice
from collections import OrderedDict
import math
import numpy
# import get_wcn_invsq

# This code reads in a given input PDB file and outputs, in a given output file, the coordinates of backbone atoms: C, CA, O, N and the CB atom and the Center-Of-Mass (COM) of the side chains and the COM of the BackBone of each residue.
# Also on the output is the Bfactors for each of the corresponding atoms and the average Bfactor for the case Side-Chain (SC) COM and the entire Amino Acid (AA) COM (including backbone atoms).
#
# INPUT:  pdb files in ../structures/*  and the name of the output file.
# OUTPUT: A file containing all relevant residue properties described above, in chronological order of pdbs in the directory.
# Amir Shahmoradi, 12:10 PM, Tuesday Aug 12 2014, Wilke Lab, iCMB, The University of Texas at Austin.

def main(argv):
    if len( argv ) != 3:
        print("     ", argv[0], "<input PDB file>", "<output summary file>", '\n')
        sys.exit('Program stopped.\n')
    else:
        pdb_in  = argv[1]  # path for the input PDB file to be read
        sum_out = argv[2]   # summary file containing all residue coordinate & Bfactor data


    sum_out_file = open(sum_out,'w')
    sum_out_file.write('resnam'+ ',' + 'resnum,' +
                      'wcnSC' + ','
                      'wcnCA\n') 


    p = PDBParser()
    pdb_name = os.path.basename(pdb_in).split('.')[0].upper()
    structure = p.get_structure(pdb_name,pdb_in)

    output_list = []

    resnam   = []     # A list containing all Residue Names
    resnum   = []     # A list containing all Residue Numbers
    reschain = []     # A list containing the chain name in which the residues lie
    crdN    = []
    crdCA   = []
    crdC    = []
    crdO    = []
    crdCB   = []
    crdAA   = []
    crdSC   = []
    bfN     = []
    bfCA    = []
    bfC     = []
    bfO     = []
    bfCB    = []
    bfAA    = []
    bfSC    = []
    sizeSC  = []   # A list containing the total number of atoms in each residue Side Chain (SC).
    sizeAA  = []   # A list containing the total number of atoms in each Amino Acid (AA).

    for residue in structure.get_residues():
        output_dict = {}
        #print residue
        output_dict['residue_name'] = residue.resname
        output_dict['residue_num'] = residue.get_full_id()[3][1]
        output_dict['chain'] = residue.get_full_id()[2]
        no_n  = True
        no_c_a = True
        no_c  = True
        no_o  = True
        no_c_b = True
        no_sidechain = True
        coord_ca = None
        sidechain_coords = []  # A list containing the coordinates of all side chain atoms of the current residue. Will be used to calculate the COM of the side chain.
        amino_acid_coords = []  # A list containing the coordinates of all atoms of the current Amino Acid. Will be used to calculate the COM of the Amino Acid.
        for atom in structure.get_atoms():
            # atom.name is equivalent to atom.get_id()
            if atom.parent.id != residue.id:
                continue
            if atom.name == 'N':
                no_n = False
            elif atom.name == 'CA':
                no_c_a = False
                coord_ca = atom.get_coord()
            elif atom.name == 'C':
                no_c = False
            elif atom.name == 'O':
                no_o = False
            elif atom.name == 'CB':
                no_c_b = False

            if atom.name not in ['C','CA','O','N']:
                no_sidechain = False
                sidechain_coords.append(atom.get_coord())

            amino_acid_coords.append(atom.get_coord())
        
        warning_message = "Missing {} backbone atom in residue ({}, {}) in PDB {}."

        if no_n:
            warnings.warn(warning_message.format('N',
                                                 output_dict['residue_num'], 
                                                 output_dict['residue_name'],
                                                 pdb_name.upper()),
                          RuntimeWarning)
        if no_c_a:
            raise RuntimeError('Missing CA backbone atom in residue (' + 
                               str(output_dict['residue_num']) + ', ' +
                               str(output_dict['residue_name']) +
                               ') in PDB ' + pdb_name.upper() + '. Cannot' 
                               ' calculate C-alpha WCN.')
        if no_c:
            warnings.warn(warning_message.format('C',
                                                 output_dict['residue_num'], 
                                                 output_dict['residue_name'],
                                                 pdb_name.upper()),
                          RuntimeWarning)
        if no_o:
            warnings.warn(warning_message.format('O',
                                                 output_dict['residue_num'], 
                                                 output_dict['residue_name'],
                                                 pdb_name.upper()),
                          RuntimeWarning)
        if no_c_b:
            warnings.warn(warning_message.format('CB',
                                                 output_dict['residue_num'], 
                                                 output_dict['residue_name'],
                                                 pdb_name.upper()),
                          RuntimeWarning)
        if no_sidechain:
            warnings.warn(warning_message.format('sidechain',
                                                 output_dict['residue_num'], 
                                                 output_dict['residue_name'],
                                                 pdb_name.upper() + 
                                                 '. Using CA instead.'),
                          RuntimeWarning)
            sidechain_coords.append(coord_ca)
        else:
            # Calculate side chain properties:
            output_dict['sidechain_size'] = len(sidechain_coords)
            output_dict['sidechain_center'] = sum(sidechain_coords)/output_dict['sidechain_size']

    # Now calcualte the Contact numbers for differnt sets of coordinates and output the results :

    wcn_sidechain = []
    wcn_ca = []

    for residue in output_list:

        wcnCAi = 0.      # WCN for atom CA of the ith residue in the PDB file.
        for other_residue in output_list:
           if residue != other_residue:
               wcnCAi += 1./( (residue['coord_ca'][0]-other_residue['coord_ca'][0])**2 + (residue['coord_ca'][1]-other_residue['coord_ca'][1])**2 + (residue['coord_ca'][2]-other_residue['coord_ca'][2])**2 )
        wcn_ca.append(wcnCAi)


        # Now write out (or append to) the ouput file
        sum_out_file.write(residue['residue_name'] + ',' + 
                           str(wcnCA[i])+'\n')

if __name__ == "__main__":
   main(sys.argv)

