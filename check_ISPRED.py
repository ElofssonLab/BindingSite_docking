#!/usr/bin/env python3:q
import os
import sys
import math
import random
from sklearn.metrics import auc
import numpy as np
from Bio.PDB import *
from Bio.PDB.DSSP import DSSP
import matplotlib.pyplot as plt
import seaborn as sb
import pandas as pd

three2one = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D',
             'CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
             'HIS':'H','ILE':'I','LEU':'L','LYS':'K',
             'MET':'M','PHE':'F','PRO':'P','SER':'S',
             'THR':'T','TRP':'W','TYR':'Y','VAL':'V',
             'MSE':'M','UNK':'X'}

def interface_extract(strs, strc, dssps, dsspc, ispred, chain):

    TP = FP = 0

    count = 0
    for residue in strs:
        count += 1
        resnumb = residue.get_id()[1]
        resname = residue.get_resname()

        if count not in ispred: continue

        if (chain, resnumb) not in dssps: 
            print ('Missing residue {}{} in DSSP for {}'.format(resnumb, chain, code))
            continue

        if ispred[count][1] != dssps[(chain, resnumb)][1]: 
            print ('Error in {}: ispred and dssp does not match! '.format(code))
            return None
        if three2one[resname] != dssps[(chain, resnumb)][1]:
            print ('Error in {}: pdb and dssp does not match! '.format(code))
            return None
        if three2one[resname] != ispred[count][1]:
            print ('Error in {}: pdb and ispred does not match! '.format(code))
            return None

        rasas = dssps[(chain, resnumb)][3]
        rasac = dsspc[(chain, resnumb)][3]
        if rasas == 'NA' or rasac == 'NA': continue

        if abs(rasas-rasac)!=0: TP += 1
        else: FP += 1
    
    return [TP/(TP+FP), TP, TP+FP]

def get_auc(strs, strc, dssps, dsspc, ispred, chain, sort_list):
    count = 0
    for residue in strs:
        count += 1
        resnumb = residue.get_id()[1]
        resname = residue.get_resname()

        if count not in ispred: continue

        if (chain, resnumb) not in dssps:
            print ('Missing residue {}{} in DSSP for {}'.format(resnumb, chain, code))
            continue

        if ispred[count][1] != dssps[(chain, resnumb)][1]:
            print ('Error in {}: ispred and dssp does not match! '.format(code))
            return None, None
        if three2one[resname] != dssps[(chain, resnumb)][1]:
            print ('Error in {}: pdb and dssp does not match! '.format(code))
            return None, None
        if three2one[resname] != ispred[count][1]:
            print ('Error in {}: pdb and ispred does not match! '.format(code))
            return None, None

        rasas = dssps[(chain, resnumb)][3]
        rasac = dsspc[(chain, resnumb)][3]
        if rasas == 'NA' or rasac == 'NA': continue

        if abs(rasas-rasac)!=0: ispred[count].append(1)
        else: ispred[count].append(0)

    prev = 1
    FN = len([el for el in ispred if len(ispred[el])==3 and ispred[el][-1] == 1])
    TN = len([el for el in ispred if len(ispred[el])==3 and ispred[el][-1] == 0])
    if FN == 0 or TN == 0: return None, None
    TP = FP = 0

    PPVs = []
    TPRs = []
    FPRs = []
    for key in sort_list:
        if len(ispred[key])<3: continue
        if ispred[key][-1] == 1:
            FN -= 1
            TP += 1
        else:
            TN -= 1
            FP += 1
        if ispred[key][0] != prev:
            PPVs.append(TP/(TP+FP))
            TPRs.append(TP/(TP+FN))
            FPRs.append(FP/(FP+TN))
            prev = ispred[key][0]
    if len(PPVs) < 2 or len(TPRs) < 2 or len(FPRs) < 2: return None, None
    else: return auc(FPRs, TPRs), auc(TPRs, PPVs) 

def ispred_extract(path):
    top = 10
    pred = {}
    for line in open(path):
        pred[int(line.split('\t')[2])]=[float(line.split('\t')[6].rstrip()), line.split('\t')[3]]
    
    sortpred = []
    for key in pred: 
        if sortpred == []: sortpred.append(key)
        else:
            pos = 0
            for stored in sortpred: 
                if pred[key] > pred[stored]: sortpred = sortpred[:pos]+[key]+sortpred[pos:]
                pos += 1
                break
        if key not in sortpred: sortpred.append(key)

    maxpred = {}
    if top != None:
        if len(sortpred)>top: sortpred = sortpred[:top]
    for key in sortpred: maxpred[key] = pred[key]
    return maxpred, sortpred

if __name__ == "__main__":
    codelist = [line.rstrip() for line in open('codes')]

    codes = []
    counts = {'both':0, 'one':0, 'none':0}
    outfile = open('aucs','w')
    for code in codelist:

        p = PDBParser(QUIET=True)
        path1 = 'binding_spots/trDocking/clean_benchmark/b4/structures/'+code+'_u1.pdb'
        ispath1 = 'binding_spots/ISPRED/'+code+'_u1.ispred4.tsv'
        path2 = 'binding_spots/trDocking/clean_benchmark/b4/structures/'+code+'_u2.pdb'
        ispath2 = 'binding_spots/ISPRED/'+code+'_u2.ispred4.tsv'
        path3 = 'binding_spots/trDocking/clean_benchmark/b4/structures/'+code+'.pdb'
        str1 = p.get_structure('', path1)
        str2 = p.get_structure('', path2)
        strc = p.get_structure('', path3)

        try: dssp1 = DSSP(str1[0], path1)
        except: 
            print("# {} DSSP calculation failed!".format(path1))
            continue
        try: dssp2 = DSSP(str2[0], path2)
        except: 
            print("# {} DSSP calculation failed!".format(path2))
            continue
        try: dsspc = DSSP(strc[0], path3)
        except: 
            print("# {} DSSP calculation failed!".format(path3))
            continue

        ispred1, sortlist1 = ispred_extract(ispath1)
        ispred2, sortlist2 = ispred_extract(ispath2)
        res_list1 = Selection.unfold_entities(str1[0], 'R')
        res_list2 = Selection.unfold_entities(str2[0], 'R')
        res_listc = Selection.unfold_entities(strc[0], 'R')

#        stats1 = interface_extract(res_list1, res_listc, dssp1, dsspc, ispred1, 'A')
#        print (code, 'U1: ', stats1)
#        stats2 = interface_extract(res_list2, res_listc, dssp2, dsspc, ispred2, 'B')
#        print (code, 'U2: ', stats2)
#        if stats1 == None or stats2 == None: counts['none'] += 1
#        elif isinstance(stats1[0], float) and stats1[0] and isinstance(stats2[0], float) and stats2[0]!= 0.0: 
#            counts['both'] += 1
#            codes.append(code)
#        elif (isinstance(stats1[0], float) and stats1[0]) or (isinstance(stats2[0], float) and stats2[0]!= 0.0): counts['one'] += 1
#        else: counts['none'] += 1
#    
#    with open('top10pos','w') as f:
#        for code in codes: f.write(code+'\n')
#    print (counts)
    
        ROCauc1, PRauc1 = get_auc(res_list1, res_listc, dssp1, dsspc, ispred1, 'A', sortlist1)
        ROCauc2, PRauc2 = get_auc(res_list2, res_listc, dssp2, dsspc, ispred2, 'B', sortlist2)
        outfile.write('{} {} {} {} {}\n'.format(code, ROCauc1, ROCauc2, PRauc1, PRauc2))


