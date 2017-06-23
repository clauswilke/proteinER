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
import argparse
import tempfile
import subprocess

from Bio import SeqIO, AlignIO
from Bio.Data import SCOPData
from Bio.PDB import PDBParser, is_aa
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_aa_seq(chain):
    '''
    Extract amino acid sequence from a PDB chain object and return sequence as
    Bio.SeqRecord object.
    '''
    aa_list = []
    residue_numbers = []
    for residue in chain:
        if is_aa(residue):
            aa_list.append(SCOPData.protein_letters_3to1[residue.resname])
            residue_numbers.append(str(residue.get_id()[1]) + \
                                   residue.get_id()[2].strip())
    aa_seq = SeqRecord(Seq(''.join(aa_list)), id='pdb_seq', description='')
    return aa_seq, residue_numbers

def run_mafft(fasta_seq, pdb_seq):
    '''
    Align two Bio.SeqRecord sequences and return an a biopython alignment
    object.
    '''
    sequences = [fasta_seq, pdb_seq]
    with tempfile.NamedTemporaryFile() as temp_fasta:
        # Temporary file for input pdb_seq fasta
        with tempfile.NamedTemporaryFile() as temp_aln:
            # Temporary file for mafft output
            SeqIO.write(sequences, temp_fasta.name, "fasta")
            try:
                subprocess.call(['mafft-linsi', temp_fasta.name],
                                stdout=temp_aln,
                                stderr=open(os.devnull, 'wb'))
            except:
                raise RuntimeError('Call to mafft failed. Check that mafft is '
                                   'in your PATH.')
            alignment = AlignIO.read(temp_aln.name, 'fasta')
    return alignment

def load_pdb_chain(name, pdb_file, model_name, chain_name):
    '''
    Load a specified chain from a PDB, with error checking.
    '''
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(name, pdb_file)
    try:
        model = structure[model_name]
    except KeyError:
        raise RuntimeError('PDB model could not be found. Please inspect PDB'
                           'file and specify model for map.')
    try:
        chain = model[chain_name]
    except KeyError:
        raise RuntimeError('PDB chain could not be found. Please inspect PDB'
                           'file and specify chain for map.')
    return chain

def make_map(alignment, residue_numbers):
    '''
    Make a amino acid to PDB residue number map using an alignment and a list of
    PDB residue numbers. Return a list of dictionaries designed to be converted
    to a CSV.
    '''
    # Split aligned sequences into two lists
    align_list = [list(rec) for rec in alignment]
    fasta_aln = align_list[0]
    pdb_aln = align_list[1]
    # Track *index* of where we are in the PDB amino acid sequence
    pdb_index = 0
    # Track FASTA amino acid sequence position (starts at 1, not an index!)
    fasta_position = 1
    out_list = []
    for fasta_aa, pdb_aa in zip(fasta_aln, pdb_aln):
        out_dict = {}
        if pdb_aa != '-':
            out_dict['pdb_position'] = residue_numbers[pdb_index]
            out_dict['pdb_aa'] = pdb_aa
            pdb_index += 1
        else:
            out_dict['pdb_position'] = ''
            out_dict['pdb_aa'] = ''
        if fasta_aa != '-':
            out_dict['fasta_position'] = fasta_position
            out_dict['fasta_aa'] = fasta_aa
            fasta_position += 1
        else:
            out_dict['fasta_position'] = ''
            out_dict['fasta_aa'] = ''
        if fasta_aa != '-' and pdb_aa != '-' and fasta_aa != pdb_aa:
            out_dict['mismatch'] = 1
        else:
            out_dict['mismatch'] = 0
        out_list.append(out_dict)
    return out_list

def main():
    '''
    Make a PDB amino acid to FASTA amino acid sequence map. This map is required
    to align evolutionary rates to positions in the structure.
    '''
    parser = argparse.ArgumentParser(
        description='Generate sequence-to-structure map for aligning '
                    'evolutionary rates to a PDB structure.')
    parser.add_argument('fasta', metavar='<FASTA path>', type=str,
                        help='input FASTA file')
    parser.add_argument('pdb', metavar='<PDB path>', type=str,
                        help='input PDB file')
    parser.add_argument('-c', metavar='<PDB chain>', type=str,
                        default='A',
                        help='if there are multiple chains in PDB, map FASTA '
                             'sequence to this chain')
    parser.add_argument('-m', metavar='<PDB model>', type=int,
                        default=0,
                        help='if there are multiple models in PDB, map FASTA '
                             'sequence to this model')
    parser.add_argument('-o', metavar='<output prefix>', type=str,
                        help='prefix for output files')
    args = parser.parse_args()
    # Grab PDB name from filename
    pdb_name = os.path.splitext(os.path.basename(args.pdb))[0]
    # Define output file names
    if args.o is None:
        # If no output prefix given, assign prefix using input filename
        args.o = pdb_name
    output_map = args.o + '.map.csv'
    # Load chain
    chain = load_pdb_chain(pdb_name.upper(), args.pdb, args.m, args.c)
    # Extract PDB numbering and amino acid sequence
    pdb_record, residue_numbers = get_aa_seq(chain)
    # Load FASTA sequence
    fasta_seq = SeqIO.read(args.fasta, 'fasta')
    fasta_record = SeqRecord(fasta_seq.seq, id='fasta_seq', description='')
    # Align PDB and FASTA sequence
    alignment = run_mafft(fasta_record, pdb_record)
    # Generate map
    output_list = make_map(alignment, residue_numbers)
    # Write map to CSV
    with open(output_map, 'w') as csvfile:
        writer = csv.DictWriter(csvfile,
                                fieldnames=['fasta_aa', 'pdb_aa',
                                            'pdb_position', 'fasta_position',
                                            'mismatch'],
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(output_list)

if __name__ == "__main__":
    main()
