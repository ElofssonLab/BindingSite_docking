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
             'MSE':'M'}

def rewrite(infileu, infileb, chain, outlistu, outlistb, newpos):
    diff = newpos
    prev_orig = ''
    stru = p.get_structure('', infileu)
    strb = p.get_structure('', infileb)

    ##### match chains 
    for chain in stru[0]: uchain = chain
    for chain in strb[0]: bchain = chain
    useq = ''.join([three2one[res.get_resname()] for res in uchain])
    bseq = ''.join([three2one[res.get_resname()] for res in bchain])
    alignments = pairwise2.align.globalds(useq, bseq, blosum62, -10, -0.5)
    ualigned = alignments[0][0]
    baligned = alignments[0][1]

    ##### renumber unbound structure
    outlistu = []
    for line in open(infileu):
        if line.startswith('ATOM'):
            if prev_orig != '' and prev_orig != line[22:27]:
                while three2one[line[17:20]] != ualigned[newpos-diff]: newpos+=1
            renumb = ' '*(4-len(str(newpos)))+str(newpos)+' '
            outline = line[:21]+chain+renumb+line[27:].rstrip()
            outlistu.append(outline+'\n')
        prev_orig = line[22:27]
    outlistu.append('TER\n')

    ##### renumber bound structure
    newpos = diff
    outlistb = []
    for line in open(infileb):
        if line.startswith('ATOM'):
            if prev_orig != '' and prev_orig != line[22:27]: 
                while three2one[line[17:20]] != baligned[newpos-diff]: newpos+=1
            renumb = ' '*(4-len(str(newpos)))+str(newpos)+' '
            outline = line[:21]+chain+renumb+line[27:].rstrip()
            outlistb.append(outline+'\n')
        prev_orig = line[22:27]
    outlistb.append('TER\n')

    return newpos, outlistu, outlistb

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Analyze ISPRED4 predictions")
    parser.add_argument("-c", type=str, help="code list")
    parser.add_argument("-i", type=str, help="input folder path")
    parser.add_argument("-o", type=str, help="out folder path")
    ns = parser.parse_args()

    if not os.path.exists(ns.o): os.mkdir(ns.o)
    if not os.path.exists(ns.o+'/str_u'): os.mkdir(ns.o+'/str_u')
    if not os.path.exists(ns.o+'/str_b'): os.mkdir(ns.o+'/str_b')
    p = PDBParser(QUIET=True)

    for code in open(ns.c):
        code = code.rstrip()
        pathu1 = ns.i+'/'+code+'_u1.pdb'
        pathu2 = ns.i+'/'+code+'_u2.pdb'
        pathb1 = ns.i+'/'+code+'_b1.pdb'
        pathb2 = ns.i+'/'+code+'_b2.pdb'
    
        outu1 = open(ns.o+'/str_u/'+code+'_u1.pdb','w')
        outu2 = open(ns.o+'/str_u/'+code+'_u2.pdb','w')
        outu = open(ns.o+'/str_u/'+code+'.pdb','w')
        outb1 = open(ns.o+'/str_b/'+code+'_b1.pdb','w')
        outb2 = open(ns.o+'/str_b/'+code+'_b2.pdb','w')
        outb = open(ns.o+'/str_b/'+code+'.pdb','w')
    
        new, linesu1, linesb1 = rewrite(pathu1, pathb1, 'A', [], [], 0)
        new, linesu2, linesb2 = rewrite(pathu2, pathb2, 'B', [], [], new)
   
        for line in linesu1: outu1.write(line) 
        outu1.write('END')
        outu1.close()

        for line in linesu1: outu2.write(line)
        outu2.write('END')
        outu2.close()

        for line in linesu1+linesu2: outu.write(line)
        outu.write('END')
        outu.close()
    
        for line in linesb1: outb1.write(line)
        outb1.write('END')
        outb1.close()

        for line in linesb1: outb2.write(line)
        outb2.write('END')
        outb2.close()

        for line in linesb1+linesb2: outb.write(line)
        outb.write('END')
        outb.close() 
        