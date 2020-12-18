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
             'MSE':'M','GLX':'Q'}

def interface_extract(stru, strb, strc, dsspu, dsspb, dsspc, chain):
    '''
    extracts residues on the protein surface and classify them
    in interacting or not interacting ones.
    '''
    count = 0
    surface = []
    upositions = [residue.get_id()[1] for residue in stru[chain]]
    bpositions = [residue.get_id()[1] for residue in strb[chain]]

    for resnumb in upositions:
        if resnumb in bpositions: 
            residue = strb[chain][resnumb]
            try: rasab = round(dsspb[(chain, resnumb)][3], 3)
            except: rasab = None
            try: rasac = round(dsspc[(chain, resnumb)][3], 3)
            except: rasac = None
            if rasab == 'NA': rasab = None
            if rasac == 'NA': rasac = None
            if rasab == None or rasac == None: surface.append([residue, None, 0.5])
            elif rasac < rasab: surface.append([residue, rasab, 1.0])
            elif rasac == rasab: surface.append([residue, rasab, 0.0])
        else:
            residue = stru[chain][resnumb]
            try: rasau = round(dsspu[(chain, resnumb)][3], 3)
            except: rasau = None
            surface.append([residue, rasau, 0.5])
    return surface

if __name__ == "__main__":
    p = argparse.ArgumentParser(description = 'Converts PDB complex in ISPRED4 format')
    p.add_argument('-b', required= True, help='bound structure')
    p.add_argument('-u', required= True, help='unbound structure')
    p.add_argument('-c', required= True, help='complex structure')
    p.add_argument('-o', required= True, help='output ISPRED4-formatted file')
    ns = p.parse_args()

    p = PDBParser(QUIET=True)
    strb = p.get_structure('', ns.b)
    stru = p.get_structure('', ns.u)
    strc = p.get_structure('', ns.c)
    for c in stru[0]: chain = c.get_id()
    positions = [res.get_id()[1] for res in stru[0][chain]]
    

    try: dsspb = DSSP(strb[0], ns.b)
    except: 
        print("# {} DSSP calculation failed!".format(ns.b))
        sys.exit()
    try: dsspu = DSSP(stru[0], ns.u)
    except:
        print("# {} DSSP calculation failed!".format(ns.u))
        sys.exit()
    try: dsspc = DSSP(strc[0], ns.c)
    except: 
        print("# {} DSSP calculation failed!".format(ns.c))
        sys.exit()

    outfile = open(ns.o, 'w')
    surf = interface_extract(stru[0], strb[0], strc[0], dsspu, dsspb, dsspc, chain)
    if len(surf) != len(positions): print ('In {} extracted prot len ({}) different from original unbound len ({})'\
                                           .format(ns.c, len(surf), len(positions)))
    for res, rasa, score in surf: 
        outfile.write('{}\t{}\t{}\tX\t{}\t{}\n'\
                      .format(ns.b.split('\t')[-1], chain, three2one[res.get_resname()], rasa, score))
    outfile.close()

    
