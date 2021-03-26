#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np
from Bio.PDB import *
from Bio import pairwise2
from Bio.Align import substitution_matrices
blosum62 = substitution_matrices.load("BLOSUM62")

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','GLX':'Q'}

def match_alignment(infile, aligned, mask, chain, newpos):
    
    pdblines = [line for line in open(infile) if line.startswith('ATOM')]

    count = 0
    aligned_pdblines = {}
    previous = pdblines[0][22:27]
    previousname = pdblines[0][17:20]
    for res in aligned:
        if res == '-' or len(pdblines) == 0:
            aligned_pdblines[count] = []
            count += 1
            continue
        
        while pdblines[0][22:27] == previous:
            aligned_pdblines[count] = aligned_pdblines.get(count, [])
            if pdblines[0][17:20] == previousname:
                aligned_pdblines[count].append(pdblines.pop(0))
            else: pdblines.pop(0)
            if len(pdblines) == 0: break

        if len(pdblines) > 0: 
            previous = pdblines[0][22:27]
            previousname = pdblines[0][17:20]
        count += 1

    outlist = []
    for p, res in enumerate(mask):
        if res == '-': continue
        renumb = ' '*(4-len(str(newpos)))+str(newpos)+' '
        for line in aligned_pdblines[p]:
            outline = line[:21]+chain+renumb+line[27:].rstrip()
            outlist.append(outline+'\n')
        newpos += 1
    outlist.append('TER\n')
    return newpos, outlist

def renumber(infileu, infileb, chain, firstpos):
    stru = p.get_structure('', infileu)
    strb = p.get_structure('', infileb)

    ##### align chains 
    for BIOchain in stru[0]: uchain = BIOchain
    for BIOchain in strb[0]: bchain = BIOchain
    useq = ''.join([three2one[res.get_resname()] for res in uchain])
    bseq = ''.join([three2one[res.get_resname()] for res in bchain])
    alignments = pairwise2.align.globalxx(useq, bseq)
    ualigned = alignments[0][0]
    baligned = alignments[0][1]
    
    ##### renumber
    lastposu, outlistu = match_alignment(infileu, ualigned, ualigned, chain, firstpos)
    lastposb, outlistb = match_alignment(infileb, baligned, ualigned, chain, firstpos)
    return lastposu, outlistu, outlistb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Analyze ISPRED4 predictions")
    parser.add_argument("-c", type=str, help="code list")
    parser.add_argument("-i", type=str, help="input folder path")
    parser.add_argument("-o", type=str, help="out folder path")
    ns = parser.parse_args()

    if not os.path.exists(ns.o): os.mkdir(ns.o)
    p = PDBParser(QUIET=True)

    for code in open(ns.c):
        code = code.rstrip()
        print (code)
        pathu1 = ns.i+'/'+code+'_u1.pdb'
        pathu2 = ns.i+'/'+code+'_u2.pdb'
        pathb1 = ns.i+'/'+code+'_b1.pdb'
        pathb2 = ns.i+'/'+code+'_b2.pdb'

        if os.path.exists(ns.o+'/str_b/'+code+'.pdb'): continue
        outu1 = open(ns.o+'/'+code+'_u1.pdb','w')
        outu2 = open(ns.o+'/'+code+'_u2.pdb','w')
        outu = open(ns.o+'/'+code+'_uc.pdb','w')
        outb1 = open(ns.o+'/'+code+'_b1.pdb','w')
        outb2 = open(ns.o+'/'+code+'_b2.pdb','w')
        outb = open(ns.o+'/'+code+'_bc.pdb','w')
    
        new, linesu1, linesb1 = renumber(pathu1, pathb1, 'A', 1)
        new, linesu2, linesb2 = renumber(pathu2, pathb2, 'B', 1)
   
        for line in linesu1: outu1.write(line) 
        outu1.write('END')
        outu1.close()

        for line in linesu2: outu2.write(line)
        outu2.write('END')
        outu2.close()

        for line in linesu1+linesu2: outu.write(line)
        outu.write('END')
        outu.close()
    
        for line in linesb1: outb1.write(line)
        outb1.write('END')
        outb1.close()

        for line in linesb2: outb2.write(line)
        outb2.write('END')
        outb2.close()

        for line in linesb1+linesb2: outb.write(line)
        outb.write('END')
        outb.close() 
        
