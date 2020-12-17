#!/usr/bin/env python3:q
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import argparse
import sys

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','UNK':'X'}

def interface_extract(strs, strc, dssps, dsspc, chain):
    '''
    extracts residues on the protein surface and classify them
    in interacting or not interacting ones.
    '''
    surface = []
    for residue in strs:
        resnumb = residue.get_id()[1]
        try: rasas = dssps[(chain, resnumb)][3]
        except: rasas = None
        try: rasac = dsspc[(chain, resnumb)][3]
        except: rasac = None

        if rasas == 'NA': rasas = None
        if rasac == 'NA': rasac = None
        if rasas == None or rasac == None: surface.append([residue, None, 0.5])
        elif rasac < rasas: surface.append([residue, rasas, 1.0])
        elif rasac == rasas: surface.append([residue, rasas, 0.0])

    return surface

if __name__ == "__main__":
    p = argparse.ArgumentParser(description = 'Converts PDB complex in ISPRED4 format')
    p.add_argument('-s', required= True, help='structure')
    p.add_argument('-c', required= True, help='complex structure')
    p.add_argument('-o', required= True, help='output ISPRED4-formatted file')
    ns = p.parse_args()

    p = PDBParser(QUIET=True)
    str1 = p.get_structure('', ns.s)
    chains = [c.get_id() for c in str1[0]]
    strc = p.get_structure('', ns.c)

    try: dssp1 = DSSP(str1[0], ns.s)
    except: 
        print("# {} DSSP calculation failed!".format(ns.s))
        sys.exit()
    try: dsspc = DSSP(strc[0], ns.c)
    except: 
        print("# {} DSSP calculation failed!".format(ns.c))
        sys.exit()

    outfile = open(ns.o,'w')
    res_listc = Selection.unfold_entities(strc[0], 'R')
    for ch in chains:
        res_list1 = Selection.unfold_entities(str1[0][ch], 'R')
        surf = interface_extract(res_list1, res_listc, dssp1, dsspc, ch)
        for res, rasa, score in surf: 
            outfile.write('{}\t{}\t{}\tX\t{}\t{}\n'\
                          .format(ns.s, ch, three2one[res.get_resname()], rasa, score))
    outfile.close()

    
