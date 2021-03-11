#!/usr/bin/env python3:q
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import argparse
import random
import sys

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','GLX':'Q'}

def get_struct_data(path):
    struct = p.get_structure('', path)
    dssp = DSSP(struct[0], path)
    return struct, dssp

def join_str_dssp(struc, dssp):
    joined = []
    for residue in struc[chain]:
        resnumb = residue.get_id()[1]
        try: rasa = round(dssp[(chain, resnumb)][3], 3)
        except: rasa = 0.0
        if rasa == 'NA': rasa = 0.0
        joined.append([residue, rasa])
    return joined

def interface_extract():
    posc = 0
    interface = []
    jointu = join_str_dssp(stru[0], dsspu)
    jointc = join_str_dssp(strc[0], dsspc)
    for posu in range(len(jointu)):
        resu = jointu[posu][0]
        rasau = jointu[posu][1]
        rasac = jointu[posc][1]
        resnumbu = jointu[posu][0].get_id()[1]
        resnumbc = jointc[posc][0].get_id()[1]
        if resnumbu < resnumbc: posc += 1
        if rasac < rasau*0.99: interface.append([resu, rasau, 1.0])
        else: interface.append([resu, rasau, 0.0])
        posc += 1
    return interface

def randomize_interfaces(surf, tpr, ppv):
    randomized = [[residue[1], residue[2]] for residue in surf]

    ##### compute number of real binding site (BS) residues
    positives = 0
    for acc, score in randomized:
        if score == 1.0: positives += 1
    ##### calculate how many false non-BS we need for desired TPR
    FNactual = 0
    FNneeded = int(positives*(1-tpr))
    ##### mark random BS residues to be turned in false non-BS
    while FNactual < FNneeded:
        for p, (acc, score) in enumerate(randomized):
            if score == 1.0 and random.randint(0,1) == 1:
                randomized[p][1] = 0.75
                FNactual += 1
            if FNactual == FNneeded: break

    ##### compute number of remaining residues correctly marked as BS
    positives = 0
    for acc, score in randomized:
        if score == 1.0: positives += 1
    ##### calculate how many false BS we need for desired PPV
    FPactual = 0
    FPneeded = int(((positives/(ppv*100))*100)-positives)
    ##### mark random true non-BS residues to be turned in false BS
    while FPactual < FPneeded:
        nonegative = True
        for p, (acc, score) in enumerate(randomized):
            if score == 0.0 and acc > 0.2 : nonegative = False
            if score == 0.0 and acc > 0.2 and random.randint(0,1) == 1:
                randomized[p][1] = 0.25
                FPactual += 1
            if FPactual == FPneeded: break
        if nonegative: break

    ##### turn all marked residues in the relative class
    for p, (acc, score) in enumerate(randomized):
        if score == 0.25: randomized[p][1] = 1.0
        if score == 0.75: randomized[p][1] = 0.0
        
    return randomized

def match_pred(pred, real, predseq, realseq):
    if len(pred) > len(real):
        substr = realseq[:5]
        idx = predseq.find(substr)
        pred = pred[idx:]
        pred = pred[:len(real)]
        return pred, real
    
    if len(pred) < len(real):
        substr = predseq[:5]
        idx = realseq.find(substr)
        for n in range(idx): pred = [0.0] + pred
        while len(pred) < len(real): pred.append(0.0)
        return pred, real

    return pred, real

def format_ISPRED4(interface, realseq):
    pred = [line.split()[-1] for line in open(ns.isp)]
    predseq = ''.join([line.split()[3] for line in open(ns.isp)])
    if len(pred) != len(interface):
        pred, interface = match_pred(pred, interface, predseq, realseq) 
    if len(pred) != len(interface):
        print ('RealLen:\t', len(interface))
        print (realseq)
        print ('ISPLen:\t', len(pred))
        print (predseq)
        sys.exit()
    return pred

def format_SPPIder(interface, realseq):
    count = 0
    pred = []
    predseq = ''
    previous = ''
    for line in open(ns.spp):
            if not line.startswith('ATOM'): continue
            resnum = line[22:26]
            residue = line[17:20]
            score = round(float(line[61:66])/100, 2)
            if resnum != previous:
                pred.append(score)
                predseq += three2one[residue]
                previous = resnum
                count += 1
    if len(pred) != len(interface):
        pred, interface = match_pred(pred, interface, predseq, realseq)
    if len(pred) != len(interface):
        print ('RealLen:\t', len(interface))
        print (realseq)
        print ('SPPILen:\t', len(pred))
        print (predseq)
        sys.exit()
    return pred

def format_dynJET2(interface, realseq):
    #pred = [line.split()[11] for line in open(ns.dyn) if not line.startswith('AA')]
    pred = [max(float(line.split()[12]), float(line.split()[13]), float(line.split()[14]))\
            if line.split()[15] == 1.0 else 0.0 for line in open(ns.dyn) if not line.startswith('AA')]
    predseq = ''.join([three2one[line.split()[0]] for line in open(ns.dyn) if not line.startswith('AA')])
    if len(pred) != len(interface):
        pred, interface = match_pred(pred, interface, predseq, realseq)
    if len(pred) != len(interface):
        print ('RealLen:\t', len(interface))
        print (realseq)
        print ('dynLen:\t', len(pred))
        print (predseq)
        sys.exit()
    return pred

if __name__ == "__main__":
    p = argparse.ArgumentParser(description = 'Converts PDB complex in ISPRED4 format')
    p.add_argument('-b', required= True, help='bound structure')
    p.add_argument('-u', required= True, help='unbound structure')
    p.add_argument('-c', required= True, help='complex structure')
    p.add_argument('-isp', required= True, help='ISPRED4 binding site prediction')
    p.add_argument('-spp', required= True, help='SPPIder binding site prediction')
    p.add_argument('-dyn', required= True, help='dynJET2 binding site prediction')
    p.add_argument('-o', required= True, help='output formatted file')
    ns = p.parse_args()

    random.seed(42)
    p = PDBParser(QUIET=True)
    strb, dsspb = get_struct_data(ns.b)
    stru, dsspu = get_struct_data(ns.u)
    strc, dsspc = get_struct_data(ns.c)
    for c in strb[0]: chain = c.get_id()
    positions = [res.get_id()[1] for res in strb[0][chain]]
    print ('######################', ns.u)

    ##### binding sites ground truth
    outfile = open(ns.o, 'w')
    interface = interface_extract()
    realseq = ''.join([three2one[res[0].get_resname()]\
                       for res in interface]) 

    ##### ground truth controlled randomization
    for tpr in [0.25, 0.5]:
        for ppv in [0.5, 0.75]:
            rand = randomize_interfaces(interface, tpr, ppv)
            for n, (acc, lbl) in enumerate(rand):
                interface[n].append(lbl)

    pred = format_ISPRED4(interface, realseq)
    for n, score in enumerate(pred): interface[n].append(score)
    
    pred = format_SPPIder(interface, realseq)
    for n, score in enumerate(pred): interface[n].append(score)
    
    pred = format_dynJET2(interface, realseq)
    for n, score in enumerate(pred): interface[n].append(score)

    ##### write out
    for residue in interface: 
        for n, val in enumerate(residue[1:]):
            if n != len(residue)-2: 
                outfile.write(str(round(float(val), 2))+'\t')
            else: 
                outfile.write(str(round(float(val), 2))+'\n')
    outfile.close()

    
