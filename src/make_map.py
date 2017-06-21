#!/usr/bin/python

'''
This script parses an input PDB file and FASTA sequence and maps the sequence 
to the PDB file using mafft. The script produces a CSV with columns for with 
PDB residue numbering, corresponding FASTA sequence numbering, and the amino 
acid. The map produced will have gaps because not all amino acids may be 
respresented in teh PDB structure.

Author: Benjamin R. Jack
'''

import os
import csv
import warnings
import argparse
import tempfile
import subprocess
from Bio.Data import SCOPData
from Bio.PDB import PDBParser
from Bio.PDB import is_aa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import StringIO

def get_aa_seq(structure):
    aa_list = []
    residue_numbers = []
    for residue in structure.get_residues():
        if is_aa(residue):
            aa_list.append(SCOPData.protein_letters_3to1[residue.resname])
            residue_numbers.append(str(residue.get_id()[1]) + \
                                   residue.get_id()[2].strip())
    aa_seq = SeqRecord(Seq(''.join(aa_list)), id='pdb_seq', description='')
    return aa_seq, residue_numbers

def run_mafft(fasta_seq, pdb_seq):
    sequences = [fasta_seq, pdb_seq]
    with tempfile.NamedTemporaryFile() as temp_fasta:
        with tempfile.NamedTemporaryFile() as temp_aln:
            SeqIO.write(sequences, temp_fasta.name, "fasta")
            subprocess.call(['linsi', temp_fasta.name], stdout = temp_aln)
            alignment = AlignIO.read(temp_aln.name, 'fasta')
    return alignment

def main():
    parser = argparse.ArgumentParser(
        description='Calculate WCN values for an input PDB.')
    parser.add_argument('fasta', metavar='<FASTA path>', type=str,
                        help='input FASTA file')
    parser.add_argument('pdb', metavar='<PDB path>', type=str,
                        help='input PDB file')
    parser.add_argument('-c', metavar='<PDB chain>', type=str,
                        help='if there are multiple chains in PDB, map FASTA '
                             'sequence to this chain')
    parser.add_argument('-o', metavar='<output prefix>', type=str,
                        help='prefix for output files')
    args = parser.parse_args()

    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    # Define output file names
    if args.o is None:
        # If no output prefix given, assign prefix using input filename
        args.o = pdb_name
    output_map = args.o + '.map.csv'

    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(pdb_name.upper(), args.pdb)

    pdb_record, residue_numbers = get_aa_seq(structure)

    fasta_seq = SeqIO.read('./structural_features/HRH1.fasta', 'fasta')
    fasta_record = SeqRecord(fasta_seq.seq, id='fasta_seq', description='')
    alignment = run_mafft(fasta_record, pdb_record)

    align_list = [list(rec) for rec in alignment]
    fasta_aln = align_list[0]
    pdb_aln = align_list[1]
    print(len(residue_numbers))
    print(len(pdb_aln))
    pdb_counter = 0
    fasta_position = 1
    out_list = []
    for fasta_aa, pdb_aa in zip(fasta_aln, pdb_aln):
        out_dict = {}
        out_dict['fasta_aa'] = fasta_aa
        out_dict['pdb_aa'] = pdb_aa
        if pdb_aa != '-':
            out_dict['pdb_position'] = residue_numbers[pdb_counter]
            pdb_counter += 1
        else:
            out_dict['pdb_position'] = ''
        if fasta_aa != '-':
            out_dict['fasta_position'] = fasta_position
            fasta_position += 1
        else:
            out_dict['fasta_position'] = ''
        out_list.append(out_dict)
    
    with open(output_map, 'w') as csvfile:
        writer = csv.DictWriter(csvfile,
                                fieldnames=['fasta_aa', 'pdb_aa', 'pdb_position', 'fasta_position'],
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(out_list)
    
if __name__ == "__main__":
    main()