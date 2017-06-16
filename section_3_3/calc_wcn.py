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
        print '''

Usage:'''
        print "     ", argv[0], "<input PDB file>", "<output summary file>", '\n'
        sys.exit('Program stopped.\n')
    else:
        pdb_in  = argv[1]  # path for the input PDB file to be read
        sum_out = argv[2]   # summary file containing all residue coordinate & Bfactor data


    sum_out_file = open(sum_out,'w')
    sum_out_file.write('resnam'+ ',' + 'resnum,' + # '\t' + 'sizeSC' + '\t' + 'sizeAA' + '\t'  \
                      'wcnSC' + ',' + # '\t' + 'bfSC'  + '\t' + \
                     # 'wcnAA' + '\t' + 'bfAA'  + '\t' + \
                     # 'wcnN'  + '\t' + 'bfN'   + '\t' + \
                      'wcnCA\n') # + '\t' + 'bfCA'  + '\t' + \
                     # 'wcnC'  + '\t' + 'bfC'   + '\t' + \
                     # 'wcnO'  + '\t' + 'bfO'   + '\t' + \
                     # 'wcnCB' + '\t' + 'bfCB'  + '\n' )
                      #'bfc' + '\t' + 'bfca' + '\t' + 'bfo' + '\t' + 'bfn' + '\t' + 'bfcb' + '\t' + 'bfaa' + '\t' + 'bfsc' + '\n')
                      #'N.x' + '\t' + 'N.y' + '\t' + 'N.z' + '\t' + \
                      #'CA.x' + '\t' + 'CA.y' + '\t' + 'CA.z' + '\t' + \
                      #'C.x' + '\t' + 'C.y' + '\t' + 'C.z' + '\t' + \
                      #'O.x' + '\t' + 'O.y' + '\t' + 'O.z' + '\t' + \
                      #'CB.x' + '\t' + 'CB.y' + '\t' + 'CB.z' + '\t' + \
                      #'SC.x' + '\t' + 'SC.y' + '\t' + 'SC.z' + '\t' + \
                      #'AA.x' + '\t' + 'AA.y' + '\t' + 'AA.z' + '\t' + \
                      #'bfc' + '\t' + 'bfca' + '\t' + 'bfo' + '\t' + 'bfn' + '\t' + 'bfcb' + '\t' + 'bfaa' + '\t' + 'bfsc' + '\n')

       #sum_out_file.write('pdb' + '\t' + 'chain' + '\t' + 'resnam' + '\t' + 'resnum' + '\t' + 'wcnc' + '\t' + 'wcnca' + '\t' + 'wcno' + '\t' + 'wcnn' + '\t' + 'wcncb' + '\t' + 'wcnaa' + '\t' + 'wcnsc' + '\t' + 'bfc' + '\t' + 'bfca' + '\t' + 'bfo' + '\t' + 'bfn' + '\t' + 'bfcb' + '\t' + 'bfaa' + '\t' + 'bfsc' + '\t' + 'sizeSC' + '\n')

    p = PDBParser()
    pdb_name = os.path.basename(pdb_in).split('.')[0].upper()
    # pdb_chain = pdb_in[-5:-4]
    structure = p.get_structure(pdb_name,pdb_in)

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

    #Ncounter  = 0
    #CAcounter = 0
    #Ccounter  = 0
    #Ocounter  = 0
    #CBcounter = 0

    for residue in structure.get_residues():
        #print residue
        resnam.append(residue.resname)
        resnum.append(residue.get_full_id()[3][1])
        reschain.append(residue.get_full_id()[2])
        noN  = True
        noCA = True
        noC  = True
        noO  = True
        noCB = True
        noSC = True
        rescrd_SC = []  # A list containing the coordinates of all side chain atoms of the current residue. Will be used to calculate the COM of the side chain.
        rescrd_AA = []  # A list containing the coordinates of all atoms of the current Amino Acid. Will be used to calculate the COM of the Amino Acid.
        resbf_SC  = []  # A list containing the Bfactors of all side chain atoms of the current residue. Will be used to calculate the side chain average Bfactor.
        resbf_AA  = []  # A list containing the Bfactors of all Amino Acid atoms of the current Amino Acid. Will be used to calculate the Amino Acid average Bfactor.
        for atom in structure.get_atoms():
            # atom.name is equivalent to atom.get_id()
            if atom.parent.id == residue.id and atom.name == 'N':
                noN = False
                crdN.append(atom.get_coord())
                bfN.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'CA':
                noCA = False
                crdCA.append(atom.get_coord())
                bfCA.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'C':
                noC = False
                #Ccounter += 1
                #print Ccounter
                crdC.append(atom.get_coord())
                bfC.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'O':
                noO = False
                crdO.append(atom.get_coord())
                bfO.append(atom.get_bfactor())
            elif atom.parent.id == residue.id and atom.name == 'CB':
                noCB = False
                crdCB.append(atom.get_coord())
                bfCB.append(atom.get_bfactor())

            if atom.parent.id == residue.id and atom.name not in ['C','CA','O','N']:
                noSC = False
                rescrd_SC.append(atom.get_coord())
                resbf_SC.append(atom.get_bfactor())

            if atom.parent.id == residue.id:
                rescrd_AA.append(atom.get_coord())
                resbf_AA.append(atom.get_bfactor())

        if noN:
            print 'missing N backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_name.upper()
            crdN.append(crdCA[-1])
            bfN.append(bfCA[-1])
        if noCA:
            print 'FATAL: missing CA backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_name.upper()
            crdCA.append(['NA','NA','NA'])
            bfCA.append('NA')
            sys.exit()
        if noC:
            print 'missing C backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_name.upper()
            crdC.append(crdCA[-1])
            bfC.append(bfCA[-1])
        if noO:
            print 'missing O backbone atom in residue: ', resnum[-1], resnam[-1], 'in PDB:',pdb_name.upper()
            crdO.append(crdCA[-1])
            bfO.append(bfCA[-1])
        if noCB:
            print 'missing CB backbone atom in residue: ', resnum[-1], resnam[-1], ', using CA instead, in PDB:',pdb_name.upper()
            crdCB.append(crdCA[-1])
            bfCB.append(bfCA[-1])
        if noSC:
            print 'missing side chain in residue: ', resnum[-1], resnam[-1], ', using CA instead in PDB:',pdb_name.upper()
            crdSC.append(crdCA[-1])
            bfSC.append(bfCA[-1])
            sizeSC.append(0)
        else:
            # Calculate side chain properties:
            sizeSC.append(len(rescrd_SC))
            crdSC.append(sum(rescrd_SC)/float(sizeSC[-1]))
            bfSC.append(sum(resbf_SC)/float(sizeSC[-1]))
            if sizeSC[-1] != len(resbf_SC):
                print 'something is terribly wrong with the code!: sizeSC[-1] != len(resbf_SC)', sizeSC[-1], len(resbf_SC)
                sys.exit()

        # Now calculate the Amino Acid properties:
        sizeAA.append(len(rescrd_AA))
        crdAA.append(sum(rescrd_AA)/float(sizeAA[-1]))
        bfAA.append(sum(resbf_AA)/float(sizeAA[-1]))
        if sizeAA[-1] != len(resbf_AA):
            print 'something is terribly wrong with the code!: sizeSC[-1] != len(resbf_SC)', sizeSC[-1], len(resbf_SC)
            sys.exit()

    # Now calcualte the Contact numbers for differnt sets of coordinates and output the results :

    # wcnN     = get_wcn_invsq.get_wcn_invsq(crdN)
    # wcnCA    = get_wcn_invsq.get_wcn_invsq(crdCA)
    # wcnC     = get_wcn_invsq.get_wcn_invsq(crdC)
    # wcnO     = get_wcn_invsq.get_wcn_invsq(crdO)
    # wcnCB    = get_wcn_invsq.get_wcn_invsq(crdCB)
    # wcnSC    = get_wcn_invsq.get_wcn_invsq(crdSC)
    # wcnAA    = get_wcn_invsq.get_wcn_invsq(crdAA)

    wcnSC = []
    wcnCA = []

    for i in range(len(resnam)):

        #wcnNi = 0.      # WCN for atom N of the ith residue in the PDB file.
        #for j in range(len(crdN)) :
        #    if i != j :
        #        wcnNi += 1./( (crdN[i][0]-crdN[j][0])**2 + (crdN[i][1]-crdN[j][1])**2 + (crdN[i][2]-crdN[j][2])**2 )
        #wcnN.append(wcnNi)

        wcnCAi = 0.      # WCN for atom CA of the ith residue in the PDB file.
        for j in range(len(crdCA)) :
           if i != j :
               wcnCAi += 1./( (crdCA[i][0]-crdCA[j][0])**2 + (crdCA[i][1]-crdCA[j][1])**2 + (crdCA[i][2]-crdCA[j][2])**2 )
        wcnCA.append(wcnCAi)

        #wcnCi = 0.      # WCN for atom C of the ith residue in the PDB file.
        #for j in range(len(crdC)) :
        #    if i != j :
        #        wcnCi += 1./( (crdC[i][0]-crdC[j][0])**2 + (crdC[i][1]-crdC[j][1])**2 + (crdC[i][2]-crdC[j][2])**2 )
        #wcnC.append(wcnCi)
        #
        #wcnOi = 0.      # WCN for atom O of the ith residue in the PDB file.
        #for j in range(len(crdO)) :
        #    if i != j :
        #        wcnOi += 1./( (crdO[i][0]-crdO[j][0])**2 + (crdO[i][1]-crdO[j][1])**2 + (crdO[i][2]-crdO[j][2])**2 )
        #wcnO.append(wcnOi)
        #
        #wcnCBi = 0.      # WCN for atom CB of the ith residue in the PDB file.
        #for j in range(len(crdCB)) :
        #    if i != j :
        #        wcnCBi += 1./( (crdCB[i][0]-crdCB[j][0])**2 + (crdCB[i][1]-crdCB[j][1])**2 + (crdCB[i][2]-crdCB[j][2])**2 )
        #wcnCB.append(wcnCBi)

        wcnSCi = 0.      # WCN for atom SC of the ith residue in the PDB file.
        for j in range(len(crdSC)) :
           if i != j :
               wcnSCi += 1./( (crdSC[i][0]-crdSC[j][0])**2 + (crdSC[i][1]-crdSC[j][1])**2 + (crdSC[i][2]-crdSC[j][2])**2 )
        wcnSC.append(wcnSCi)

        #wcnAAi = 0.      # WCN for atom N of the ith residue in the PDB file.
        #for j in range(len(crdAA)) :
        #    if i != j :
        #        wcnAAi += 1./( (crdAA[i][0]-crdAA[j][0])**2 + (crdAA[i][1]-crdAA[j][1])**2 + (crdAA[i][2]-crdAA[j][2])**2 )
        #wcnAA.append(wcnAAi)

        # Now write out (or append to) the ouput file
        sum_out_file.write(resnam[i] + ',' + str(resnum[i]) + ',' + # '\t' + str(sizeSC[i]) + '\t' + str(sizeAA[i]) + '\t' + \
                           str(wcnSC[i])  + ',' + # '\t' + str(bfSC[i])  + '\t' + \
                           #str(wcnAA[i])  + '\t' + str(bfAA[i])  + '\t' + \
                           #str( wcnN[i])  + '\t' + str( bfN[i])  + '\t' + \
                           str(wcnCA[i])+'\n')  # + '\t' + str(bfCA[i])  + '\t' + \
                           #str( wcnC[i])  + '\t' + str( bfC[i])  + '\t' + \
                           #str( wcnO[i])  + '\t' + str( bfO[i])  + '\t' + \
                           #str(wcnCB[i])  + '\t' + str(bfCB[i])  + '\n' )
                           #str(crdN[i])  + '\t' + str(bfN[i])  + '\t' + \
                           #str(crdCA[i]) + '\t' + str(bfCA[i]) + '\t' + \
                           #str(crdC[i])  + '\t' + str(bfC[i])  + '\t' + \
                           #str(crdO[i])  + '\t' + str(bfO[i])  + '\t' + \
                           #str(crdCB[i]) + '\t' + str(bfCB[i]) + '\t' + \
                           #str(crdSC[i]) + '\t' + str(bfSC[i]) + '\t' + \
                           #str(crdAA[i]) + '\t' + str(bfAA[i]) + '\t' + '\n')

if __name__ == "__main__":
   main(sys.argv)



#crd.append( numpy.array( [ float(record.split()[-3]) , float(record.split()[-2]) , float(record.split()[-1]) ] ) ) # stores CA atom triplet coordinates as individual elements of the list crd.
## Now calculate residue wcn :
#wcn = []   # A list containing all CA-atom WCN in the pdb file
#for i in range(len(crd)) :
#    sum_terms = 0.
#    for j in range(len(crd)) :
#        if i != j :
#            sum_terms += 1./( (crd[i][0]-crd[j][0])**2 + (crd[i][1]-crd[j][1])**2 + (crd[i][2]-crd[j][2])**2 )
#    wcn.append(sum_terms)
#input.close()
