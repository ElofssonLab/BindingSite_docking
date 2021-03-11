#!/usr/bin/env python3
import os
import sys
import argparse
import numpy as np

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','GLX':'Z'}

def match_numbering(code, partner):
    pathu = '{}/{}_u{}.pdb'.format(ns.i, code, partner)
    pathb = '{}/{}_b{}.pdb'.format(ns.i, code, partner)
    if partner == 1: chain = 'A'
    else: chain = 'B'

    ##### fetch alignment from unbound file header
    previous = 0
    align = False
    useq = bseq = ''
    with open(pathu) as infile: 
        for line in infile:
            if line.startswith('REMARK   Unbound'): align = True
            if align and line.startswith('REMARK   Unbound'):
                if previous != 0 and int(line.split()[2]) != previous+1: break
                useq += line.split()[-2]
                previous = int(line.split()[-1])
            if align and line.startswith('REMARK   Bound'):
                bseq += line.split()[-2]
    
    ##### fetch atomic coordinates
    ucoord = [line for line in open(pathu) if line.startswith('ATOM')]
    bcoord = [line for line in open(pathb) if line.startswith('ATOM')]
    if partner == 1: outbc = open('{}/{}_b.pdb'.format(ns.o, code, partner),'w')
    else: outbc = open('{}/{}_b.pdb'.format(ns.o, code, partner),'a')
    outu = open('{}/{}_u{}.pdb'.format(ns.o, code, partner),'w')
    outb = open('{}/{}_b{}.pdb'.format(ns.o, code, partner),'w')

    renumber = 1
    ucount = bcount = 0
    ##### cycle over positions in the alignment
    for u, b in zip(useq, bseq):
        ##### skip alignment positions with gaps in unbound structures
        ##### without increasing the new residue number
        if u == '-': continue 
        ##### write out all modified coordinates matching the residue number
        ##### of atom at index ucount and increment ucount
        res_number = ucoord[ucount][22:27]
        while ucoord[ucount][22:27] == res_number:
            rs = str(renumber).rjust(4)
            renumbered = ucoord[ucount][:21]+chain+rs+' '+ucoord[ucount][27:]
            outu.write(renumbered)
            ucount += 1
            if ucount == len(ucoord): break
        ##### skip alignment positions with gaps in bound structures
        ##### increasing the new residue number
        if b == '-':
            renumber += 1
            continue
        ##### increment bcount without writing out, when residue name in 
        ##### indexed line doesnt match residue indicated in the alignment
        ##### (most likely due to gap occurring in unbound)
        while b != three2one[bcoord[bcount][17:20]]: bcount += 1
        ##### write out all modified coordinates matching the residue number
        ##### of atom at index bcount and increment bcount
        res_number = bcoord[bcount][22:27]
        while bcoord[bcount][22:27] == res_number:
            rs = str(renumber).rjust(4)
            renumbered = bcoord[bcount][:21]+chain+rs+' '+bcoord[bcount][27:]
            outbc.write(renumbered)
            outb.write(renumbered)
            bcount += 1
            if bcount == len(bcoord): break
        ##### increment the new residue number
        renumber += 1

    outu.write('TER')
    outb.write('TER')
    outbc.write('TER')
    if partner != 1: outbc.write('END')
    outbc.close()
    outu.close()
    outb.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Prepare Benchmark4 structures")
    parser.add_argument("-c", type=str, help="code list")
    parser.add_argument("-i", type=str, help="input folder path")
    parser.add_argument("-o", type=str, help="out folder path")
    ns = parser.parse_args()

    if not os.path.exists(ns.o): os.mkdir(ns.o)

    for code in open(ns.c):
        code = code.rstrip()
        match_numbering(code, 1)
        match_numbering(code, 2)

