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
             'MSE':'M','GLX':'Z'}

def match_numbering(code, partner):
    pathu = '{}/{}_{}_u.pdb'.format(ns.i, code.upper(), partner)
    pathb = '{}/{}_{}_b.pdb'.format(ns.i, code.upper(), partner)
    if partner == 'r': chain = 'A'
    else: chain = 'B'
    if partner == 'r': partner = 1
    else: partner = 2

    ##### fetch atomic coordinates
    ucoord = [line for line in open(pathu) if line.startswith('ATOM')]
    bcoord = [line for line in open(pathb) if line.startswith('ATOM')]
    if partner == 1: 
        outbc = open('{}/{}_bc.pdb'.format(ns.o, code),'w')
        outuc = open('{}/{}_uc.pdb'.format(ns.o, code),'w')
    else: 
        outbc = open('{}/{}_bc.pdb'.format(ns.o, code),'a')
        outuc = open('{}/{}_uc.pdb'.format(ns.o, code),'a')
    outb = open('{}/{}_b{}.pdb'.format(ns.o, code, partner),'w')
    outu = open('{}/{}_u{}.pdb'.format(ns.o, code, partner),'w')

    ##### align sequences from the two structures
    sequ = ''
    previous = ''
    for line in ucoord: 
        if line[22:27] != previous: 
            sequ += three2one[line[17:20]]
            previous = line[22:27]
    seqb = ''
    for line in bcoord:
        if line[22:27] != previous:
            seqb += three2one[line[17:20]]
            previous = line[22:27]

    alignments = pairwise2.align.globalds(sequ, seqb, blosum62, -10, -0.5)
    alignedu = alignments[0][0]
    alignedb = alignments[0][1]

    print (alignedu)
    print (alignedb)

    renumber = 1
    ucount = bcount = 0
    ##### cycle over positions in the alignment
    for u, b in zip(alignedu, alignedb):
        ##### skip alignment positions with gaps in unbound structures
        ##### without increasing the new residue number
        if u == '-': continue 
        ##### write out renumbered coordinates matching the residue number
        ##### of atom at index ucount and increment ucount
        res_name = three2one[ucoord[ucount][17:20]]
        res_number = ucoord[ucount][22:27]
        while ucoord[ucount][22:27] == res_number\
        and three2one[ucoord[ucount][17:20]] == res_name\
        and ucount < len(ucoord):
            matchu = [three2one[ucoord[ucount][17:20]], renumber]
            rs = str(renumber).rjust(4)
            renumbered = ucoord[ucount][:21]+chain+rs+' '+ucoord[ucount][27:]
            outuc.write(renumbered)
            outu.write(renumbered)
            if ucount < len(ucoord)-1: ucount += 1
            else: break
        ##### skip eventual overlapping residues with same res number
        while ucoord[ucount][22:27] == res_number and ucount < len(ucoord)-1:
            ucount += 1

        ##### skip alignment positions with gaps in bound structures
        ##### increasing the new residue number
        if b == '-':
            renumber += 1
            continue
        ##### increment bcount without writing out, when residue name in 
        ##### indexed line doesnt match residue indicated in the alignment
        ##### (most likely due to gap occurring in unbound)
        while b != three2one[bcoord[bcount][17:20]]: 
            print (b, three2one[bcoord[bcount][17:20]])
            bcount += 1
        ##### write out renumbered coordinates matching the residue number
        ##### of atom at index bcount and increment bcount
        res_name = three2one[bcoord[bcount][17:20]]
        res_number = bcoord[bcount][22:27]
        while bcoord[bcount][22:27] == res_number\
        and three2one[bcoord[bcount][17:20]] == res_name\
        and bcount < len(bcoord):
            matchb = [three2one[bcoord[bcount][17:20]], renumber]
            rs = str(renumber).rjust(4)
            renumbered = bcoord[bcount][:21]+chain+rs+' '+bcoord[bcount][27:]
            outbc.write(renumbered)
            outb.write(renumbered)
            if bcount < len(bcoord)-1: bcount += 1
            else: break
        ##### skip eventual overlapping residues with same res number
        while bcoord[bcount][22:27] == res_number and bcount < len(bcoord)-1: 
            bcount += 1
        ##### increment the new residue number
        if matchu != matchb: print (matchu, matchb)
        renumber += 1
    
    outb.write('TER\nEND')
    outu.write('TER\nEND')
    if partner == 1:
        outuc.write('TER\n')
        outbc.write('TER\n')
    else:
        outbc.write('TER\nEND')
        outbc.write('TER\nEND')
    outbc.close()
    outuc.close()
    outu.close()
    outb.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Prepare Benchmark4 structures")
    parser.add_argument("-c", type=str, help="code list")
    parser.add_argument("-i", type=str, help="input folder path")
    parser.add_argument("-o", type=str, help="out folder path")
    ns = parser.parse_args()

    if not os.path.exists(ns.o): os.mkdir(ns.o)

    p = PDBParser(QUIET=True)
    for code in open(ns.c):
        code = code.rstrip()
        print (code)
        match_numbering(code, 'r')
        match_numbering(code, 'l')

