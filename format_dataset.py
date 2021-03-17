#!/usr/bin/env python3:q
import argparse
import random
import sys

from Bio.PDB import *
from Bio import pairwise2
from Bio.Align import substitution_matrices
from Bio.PDB.Selection import unfold_entities
from Bio.PDB.NeighborSearch import NeighborSearch

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','GLX':'Q'}

def extract_structure(spath, apath):
    seq = ''
    str_map = {}
    ##### load pdb file
    struc = p.get_structure('', spath)
    
    ##### load naccess file
    rasas = [float(line.split()[5])/100 if line.split()[5] != 'N/A' 
             else 0.0 for line in open(apath) if line.startswith('RES') 
             and line.split()[2] == chain]
    
    ##### join structure to relative accessibility
    for n, res in enumerate(struc[0][chain]):
        resnum = res.get_id()[1]
        resname = three2one[res.get_resname()]
        str_map[resnum] = [res, rasas[n]]
        seq += resname

    assert len(rasas) == len(str_map.keys()),\
           'Structure and rasa have different lengths!'

    return str_map, seq

def interface_extract():
    ##### load structural data
    str_map = {}
    str_mapu, sequ = extract_structure(pathu, pathua)
    str_mapb, seqb = extract_structure(pathb, pathba)
    str_mapuc, sequc = extract_structure(pathuc, pathuca)
    str_mapbc, seqbc = extract_structure(pathbc, pathbca)

    ##### iterate over all aligned positions
    maxpos = max(list(str_mapu.keys()))
    for pos in range(1, maxpos+1):
        ##### if position is mapped to bound structure
        if pos in str_mapb: 
            ##### mark residue as interface if rasa monomer > rasa dimer
            if str_mapb[pos][1] > str_mapbc[pos][1]:
                str_map[pos] = str_mapu[pos]+[1.0]
            else:
                str_map[pos] = str_mapu[pos]+[0.0]
        ##### same with unbound structure if position has no 
        ##### bound structure match
        else:
            if str_mapu[pos][1] > str_mapuc[pos][1]:
                str_map[pos] = str_mapu[pos]+[1.0]
            else:
                str_map[pos] = str_mapu[pos]+[0.0]

    ##### convert residue object to one letter aminoacid code
    for key in str_map:
        resname = three2one[str_map[key][0].get_resname()]
        str_map[key][0] = resname

    return str_map, sequ

def randomize_interfaces(tpr, ppv):
    randomized = [[intmap[pos][1], intmap[pos][2]] for pos in keylist]

    ##### compute number of real binding site (BS) residues
    positives = 0
    for acc, score in randomized:
        if score == 1.0: positives += 1
    ##### calculate how many false non-BS we need for desired TPR
    FNneeded = int(positives*(1-tpr))
    ##### mark random BS residues to be turned in false non-BS
    while FNneeded > 0:
        for p, (acc, score) in enumerate(randomized):
            if score == 1.0 and random.randint(0,1) == 1:
                randomized[p][1] = 0.75
                FNneeded -= 1
            if FNneeded <= 0: break

    ##### compute number of remaining residues correctly marked as BS
    positives = 0
    for acc, score in randomized:
        if score == 1.0: positives += 1
    ##### calculate how many false BS we need for desired PPV
    FPneeded = int(((positives/(ppv*100))*100)-positives)
    ##### mark random true non-BS residues to be turned in false BS
    while FPneeded > 0:
        nonegative = True
        for p, (acc, score) in enumerate(randomized):
            if score == 0.0 and acc > 0.2 : nonegative = False
            if score == 0.0 and acc > 0.2 and random.randint(0,1) == 1:
                randomized[p][1] = 0.25
                FPneeded -= 1
            if FPneeded <= 0: break
        if nonegative: break

    ##### turn all marked residues in the relative class
    for p, (acc, score) in enumerate(randomized):
        if score == 0.25: randomized[p][1] = 1.0
        if score == 0.75: randomized[p][1] = 0.0
    
    randomized = {pos:randomized[idx] for idx, pos in enumerate(keylist)}
    return randomized

def format_ISPRED4():
    seq = ''
    predmap = {}
    for n, line in enumerate(open(isp)):
        predmap[n+1] = [line.split()[3], line.split()[-1]]
        seq += line.split()[3]
    return predmap, seq

def format_SPPIder():
    seq = ''
    count = 1
    predmap = {}
    previous = ''
    for line in open(spp):
        if not line.startswith('ATOM'): continue
        resnum = str(int(line[22:26]))
        resname = three2one[line[17:20]]
        score = round(float(line[61:66])/100, 2)
        if resnum != previous:
            predmap[count] = [resname, score]
            previous = resnum
            seq += resname
            count += 1
    return predmap, seq

def format_dynJET2():
    seq = ''
    count = 1
    predmap = {}
    for line in open(dyn):
        if not line.startswith('AA'):
            score = max(float(line.split()[12]), 
                        float(line.split()[13]), 
                        float(line.split()[14]))
            resname = three2one[line.split()[0]]
            if line.split()[15] == '1.0':
                predmap[count] = [resname, score]
            else: 
                predmap[count] = [resname, 0.0]
            seq += resname
            count += 1
    return predmap, seq

def align_predictions(pseq, pdic):
    alignments = pairwise2.align.globalxx(realseq, pseq)

    posr = posp = 1
    alignedr = alignments[0][0]
    alignedp = alignments[0][1]
    for r, p in zip(alignedr, alignedp):
        ##### if there is a match or mismatch
        if r != '-' and p != '-':
            if r != p: print ('Mismatch!', r, p)
            intmap[posr].append(pdic[posp][1])
            posr += 1
            posp += 1
        ##### if there is a gap in the prediction
        elif r != '-':
            intmap[posr].append(0.0)
            posr += 1
        ##### skip positions gapped in unbound sequence
        elif p != '-':
            posp += 1


if __name__ == "__main__":
    p = argparse.ArgumentParser(description = 'Converts PDB complex in ISPRED4 format')
    p.add_argument('-c', required= True, help='code')
    p.add_argument('-s', required= True, help='structure path')
    p.add_argument('-r', required= True, help='freesasa rsa path')
    p.add_argument('-isp', required= True, help='ISPRED4 binding site prediction')
    p.add_argument('-spp', required= True, help='SPPIder binding site prediction')
    p.add_argument('-dyn', required= True, help='dynJET2 binding site prediction')
    p.add_argument('-o', required= True, help='output formatted file')
    ns = p.parse_args()

    random.seed(42)
    p = PDBParser(QUIET=True)

    pdb = ns.c[:4]
    numb = ns.c[6]
    pathu = '{}{}_u{}.pdb'.format(ns.s, pdb, numb)
    pathua = '{}{}_u{}.rsa'.format(ns.r, pdb, numb)
    pathuc = '{}{}_uc.pdb'.format(ns.s, pdb)
    pathuca = '{}{}_uc.rsa'.format(ns.r, pdb)

    pathb = '{}{}_b{}.pdb'.format(ns.s, pdb, numb)
    pathba = '{}{}_b{}.rsa'.format(ns.r, pdb, numb)
    pathbc = '{}{}_bc.pdb'.format(ns.s, pdb)
    pathbca = '{}{}_bc.rsa'.format(ns.r, pdb)
    
    isp = '{}{}_u{}.ispred4.tsv'.format(ns.isp, pdb, numb)
    dyn = '{}{}_u{}_jet.res'.format(ns.dyn, pdb, numb)
    spp = '{}{}_u{}.pdb'.format(ns.spp, pdb, numb)

    if '_u1' in ns.c: 
        chain = 'A'
        oppchain = 'B'
    else: 
        chain = 'B'
        oppchain = 'A'

    print ('####################', ns.c[:-4])
    ##### binding sites ground truth
    intmap, realseq = interface_extract()
    keylist = intmap.keys()
    
    ##### ground truth controlled randomization
    for tpr in [0.25, 0.5]:
        for ppv in [0.5, 0.75]:
            rand = randomize_interfaces(tpr, ppv)
            for key in keylist:
                intmap[key].append(rand[key][1])

    ##### mapping predictions to unbound sequence 
    predisp, seqisp = format_ISPRED4()
    align_predictions(seqisp, predisp)

    predspp, seqspp = format_SPPIder()
    align_predictions(seqspp, predspp)
    
    preddyn, seqdyn = format_dynJET2()
    align_predictions(seqdyn, preddyn)

    ##### write out
    outfile = open(ns.o+ns.c[:-4]+'.site', 'w')
    for pos in keylist:
        for n, val in enumerate(intmap[pos][1:]):
            if n != len(intmap[pos][1:])-1: 
                outfile.write(str(round(float(val), 2))+'\t')
            else: 
                outfile.write(str(round(float(val), 2))+'\n')
    outfile.close()

    
